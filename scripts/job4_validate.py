"""Job 4 validation of the split feedback channels (lkspin.secular P_FB_L, P_FB_E). Writes results/job4/validation.json.

1. Bit identity: the default multipliers (1, 1) reproduce the stored step 6 case (b) pilot chunks exactly (every array).
2. Conserved quantities of each variant, LK off (spin-orbit problem alone) and LK on (quadrupole, e_b = 0):
   both channels (case b), L only, E only, neither.  Reported: dJ/J of J = S + L, dS/S, dL/L, d(S.L_hat), d(J.z), d(L_hat),
   and the orthogonality l.e (a channel that rotates e without rotating L, or the reverse, leaves the orbital plane).
Run:  uv run python scripts/job4_validate.py
"""

import json

import numpy as np

import step6_run as W
import step6_validate as V
from lkspin import secular as S

OUT = S.ROOT / "results" / "job4"
PILOT_SYSTEMS = (3, 4, 5)  # hot-Jupiter systems of the stored b pilot (0, 1, 2 are mixed / timeout and have no series)
# (fb_l, fb_e, fb_s, fb_p); the last two are the supplementary Hamiltonian-consistent split (b_P apsidal only, b_F the rest)
VARIANTS = {"both": (1.0, 1.0, 1.0, 1.0), "L_only": (1.0, 0.0, 1.0, 1.0), "E_only": (0.0, 1.0, 1.0, 1.0), "neither": (0.0, 0.0, 1.0, 1.0),
            "P_apsidal_only": (0.0, 1.0, 0.0, 1.0), "F_all_but_apsidal": (1.0, 1.0, 1.0, 0.0)}
COMPARE_KEYS = ("status", "t_end", "steps", "n_events", "theta_f", "theta_f_twin", "series", "series_members")


def bit_identity() -> dict:
    table = W.step5_table()
    items = [it for it in W.build_items("b", 20, table, 20, 0.0) if it[3] in PILOT_SYSTEMS]
    out = {}
    for item in items:
        stored = dict(np.load(W.chunk_path("b", 20, item[3])))
        new = W.integrate_system("b", *item[4:9][:4], item[7], item[8]) if False else W.integrate_system("b", item[4], item[5], item[6], item[7], item[8], item[9])
        same = {k: bool(np.array_equal(stored[k], new[k], equal_nan=True)) for k in COMPARE_KEYS if k in stored}
        out[f"system_{item[3]}"] = {"outcome": W.OUTCOME_NAMES[new["outcome"]], "all_arrays_identical": all(same.values()), "per_array": same}
    return out


def variant_par(fb_l: float, fb_e: float, fb_s: float, fb_p: float, lk: bool) -> np.ndarray:
    base = S.make_params(V.CFG4, srf=False, tides=False, lk=lk)
    par = S.extend_params(base, fb_l=fb_l, fb_e=fb_e, fb_s=fb_s, fb_p=fb_p)
    par[S.P_RTIDE], par[S.P_STOPA] = 0.0, 0.0
    return par


def drift_summary(run: dict) -> dict:
    l_vec, s_vec = run["l"], run["s"]
    j_vec = l_vec + s_vec
    lh = l_vec / np.linalg.norm(l_vec, axis=1)[:, None]
    smag = np.linalg.norm(s_vec, axis=1)
    lmag = np.linalg.norm(l_vec, axis=1)
    return {"max_dJ_over_J": float(np.abs(j_vec - j_vec[0]).max() / np.linalg.norm(j_vec[0])),
            "max_dJz_over_J": float(np.abs(j_vec[:, 2] - j_vec[0, 2]).max() / np.linalg.norm(j_vec[0])),
            "max_d_absS_rel": float(np.abs(smag / smag[0] - 1).max()), "max_d_absL_rel": float(np.abs(lmag / lmag[0] - 1).max()),
            "max_dLhat_rad": float(np.arccos(np.clip(lh @ lh[0], -1, 1)).max()),
            "max_d_S_dot_Lhat_rel": float(np.abs(np.einsum("ij,ij->i", s_vec, lh) / smag[0] - (s_vec[0] @ lh[0]) / smag[0]).max()),
            "max_abs_l_dot_e_over_l_e": float(np.abs(run["l_dot_e"]).max())}


def l_dot_e_track(par: np.ndarray, theta: float, t_total_tk: float) -> np.ndarray:
    y = V.tilt_spin(S.initial_state_ext(par, V.CFG4, theta, V.P_STAR_STRONG_DAYS), V.TILT_DEG)
    tk = S.t_k_years(par, V.CFG4["a0_au"])
    vals, t = [], 0.0
    for k in range(V.N_SEGMENTS + 1):
        t_next = k * t_total_tk * tk / V.N_SEGMENTS
        if k:
            y = S.run_ext(y, par, t_next, t_start_yr=t, ev_max=10)["y"]
            t = t_next
        vals.append(float(y[:3] @ y[3:6] / (np.linalg.norm(y[:3]) * np.linalg.norm(y[3:6]))))
    return np.array(vals)


def conservation() -> dict:
    out = {}
    span_off = V.LK_OFF_SPAN_YR
    for name, mults in VARIANTS.items():
        par_off, par_on = variant_par(*mults, False), variant_par(*mults, True)
        tk = S.t_k_years(par_off, V.CFG4["a0_au"])
        run_off = V.track(par_off, 89.0, V.P_STAR_STRONG_DAYS, span_off / tk, V.N_SEGMENTS, V.TILT_DEG)
        run_off["l_dot_e"] = l_dot_e_track(par_off, 89.0, span_off / tk)
        run_on = V.track(par_on, V.THETA_LK_DEG, V.P_STAR_STRONG_DAYS, V.N_SEGMENTS * V.SEGMENT_TK, V.N_SEGMENTS, V.TILT_DEG)
        run_on["l_dot_e"] = l_dot_e_track(par_on, V.THETA_LK_DEG, V.N_SEGMENTS * V.SEGMENT_TK)
        out[name] = {"lk_off": drift_summary(run_off), "lk_on": drift_summary(run_on),
                     "lk_on_e_max": float(run_on["e"].max()), "lk_on_e_min": float(run_on["e"].min())}
    return out


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    res = {"bit_identity_vs_step6_b_pilot": bit_identity(), "conservation": conservation()}
    (OUT / "validation.json").write_text(json.dumps(res, indent=2) + "\n")
    print(json.dumps(res, indent=2))


if __name__ == "__main__":
    main()
