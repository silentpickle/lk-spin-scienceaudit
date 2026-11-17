"""Step 6 validation of lkspin.secular's opt-in physics (feedback, octupole, spin-down). Writes results/step6/validation.json.

1. feedback multiplier = 0 reproduces the step 5 theta_f (first 10 systems that step 5 labels non-chaotic, and the first 10
   systems regardless of label, to show what chaos does);
2. angular momentum: S* + L with LK off; (S* + L).L_b with LK on (e_b = 0); the energy H_Q + H_SL of the quadrupole problem
   with feedback; with octupole, the quadrupole-level j_z and H_Q are broken at order eps_oct;
3. A16 Eq. Lvec as printed against the e<->j mirror of Eq. evec (the form implemented), and the octupole check against an
   N-body integration is scripts/step6_nbody_oct.py;
4. Skumanich spin-down against Omega0 / sqrt(1 + 2 alpha Omega0^2 t).
Run:  uv run python scripts/step6_validate.py
"""

import json
import math

import numpy as np

import step6_run as W
from lkspin import secular as S

CFG4 = S.load_config()
OUT = S.ROOT / W.CFG6["results"]
N_SYSTEMS = 10
TILT_DEG = 30.0
THETA_LK_DEG = 75.0
P_STAR_STRONG_DAYS = 1.0
N_SEGMENTS = 48
SEGMENT_TK = 0.75
N_RANDOM_STATES = 4
LK_OFF_SPAN_YR = 1.0e9  # span of the LK-off conservation run
ALPHA_MB_YR = 1.5e-14  # A16 l.1903
SKUMANICH_TIMES_YR = (1e7, 1e8, 1e9, 5e9)
OCT_E_B = 0.6
OCT_PHI = 0.7
OCT_THETA_DEG = 80.0
OCT_FACTORS = (1.0, 0.5, 0.25)
OCT_CYCLES_TK = 12.0
LITERAL_FACTOR = 75.0 / 64.0
CHAOTIC_LABEL = "III"
LABELS_PATH = S.ROOT / "results" / "step5" / "labels_v2.npz"
TIGHT_FACTOR = 1e-2  # tolerance multiplier of the tight same-algorithm comparison


def tilt_spin(y: np.ndarray, deg: float) -> np.ndarray:
    """y with S* rotated by deg away from L_hat (toward +x), magnitude kept."""
    out = y.copy()
    lhat = out[:3] / np.linalg.norm(out[:3])
    smag = np.linalg.norm(out[S.N_ORBIT:S.N_ORBIT + S.N_SPIN])
    d = math.radians(deg)
    out[S.N_ORBIT:S.N_ORBIT + S.N_SPIN] = smag * (math.cos(d) * lhat + math.sin(d) * np.array([1.0, 0.0, 0.0]))
    return out


def hamiltonians(y: np.ndarray, par: np.ndarray) -> tuple[float, float]:
    """(H_Q, H_SL) per unit time, in the normalisation dj/dt = -(j x grad_j H + e x grad_e H), j = L / L_circ.

    H_Q = (3 / 8 t_k)(5 e_z^2 - j_z^2 - 2 e^2) (quadrupole, L_b = z); H_SL = -(1/6) alpha (S/L_circ)(3 cos^2 th_sl - 1).
    The A16 feedback equations (dL/dt|SL, de/dt|SL) follow from H_SL; the sum is conserved without SRFs, tides, octupole."""
    lmag, j, a = S.orbit_elements(y, par)
    n = math.sqrt(par[S.P_G] * par[S.P_MTOT] / a**3)
    tk_inv = par[S.P_MB] / par[S.P_MTOT] * (a / par[S.P_AB]) ** 3 * n
    jz, e = j * y[2] / lmag, y[3:6]
    h_q = 0.375 * tk_inv * (5.0 * e[2] ** 2 - jz**2 - 2.0 * (e @ e))
    s = y[S.N_ORBIT:S.N_ORBIT + S.N_SPIN]
    smag = np.linalg.norm(s)
    cos_sl = (y[:3] / lmag) @ s / smag
    omega = smag * par[S.P_LSCALE] / par[S.P_ISTAR]
    alpha = par[S.P_STAR_COEF] * omega / (a**3 * j**3)
    h_sl = -alpha * (smag * j / lmag) * (3.0 * cos_sl**2 - 1.0) / 6.0
    return float(h_q), float(h_sl)


def lk_params(*, lk: bool = True, fb_mult: float = 1.0, e_b: float = 0.0, octupole: bool = False, alpha_mb: float = 0.0) -> np.ndarray:
    base = S.make_params(CFG4, srf=False, tides=False, lk=lk)
    par = S.extend_params(base, fb_mult=fb_mult, e_b=e_b, phi_b=OCT_PHI, octupole=octupole, alpha_mb=alpha_mb)
    par[S.P_RTIDE], par[S.P_STOPA] = 0.0, 0.0
    return par


def track(par: np.ndarray, theta_lb0: float, p_star: float, t_total_tk: float, n_seg: int, tilt: float) -> dict:
    """Integrate in n_seg segments; return per-segment arrays of time/t_k, J, S, L vectors, e and the two Hamiltonians."""
    y = tilt_spin(S.initial_state_ext(par, CFG4, theta_lb0, p_star), tilt)
    tk = S.t_k_years(par, CFG4["a0_au"])
    rows, t = [], 0.0
    for k in range(n_seg + 1):
        t_next = k * t_total_tk * tk / n_seg
        if k:
            y = S.run_ext(y, par, t_next, t_start_yr=t, ev_max=10)["y"]
            t = t_next
        hq, hsl = hamiltonians(y, par)
        rows.append((t / tk, *y[:3], *y[S.N_ORBIT:S.N_ORBIT + S.N_SPIN], np.linalg.norm(y[3:6]), hq, hsl))
    arr = np.array(rows)
    return {"t_tk": arr[:, 0], "l": arr[:, 1:4], "s": arr[:, 4:7], "e": arr[:, 7], "h_q": arr[:, 8], "h_sl": arr[:, 9]}


def rel_range(x: np.ndarray, scale: float) -> float:
    return float(np.ptp(x) / scale)


def precession_periods(par: np.ndarray, run: dict, span_yr: float) -> float:
    """Number of precession periods of L and S* about J in span_yr: rate = |Omega_ps| |J| / |L|, Omega_ps = alpha cos(theta_sl)."""
    l_vec, s_vec = run["l"][0], run["s"][0]
    a, j = par[S.P_A0], math.sqrt(1.0 - CFG4["e0"] ** 2)  # initial state: |l| = j, so a = a0
    omega = np.linalg.norm(s_vec) * par[S.P_LSCALE] / par[S.P_ISTAR]
    alpha = par[S.P_STAR_COEF] * omega / (a**3 * j**3)
    cos_sl = (l_vec @ s_vec) / (np.linalg.norm(l_vec) * np.linalg.norm(s_vec))
    rate = abs(alpha * cos_sl) * np.linalg.norm(l_vec + s_vec) / np.linalg.norm(l_vec)
    return float(rate * span_yr / (2.0 * math.pi))


def check_conservation() -> dict:
    par_off = lk_params(lk=False)
    run = track(par_off, 89.0, P_STAR_STRONG_DAYS, LK_OFF_SPAN_YR / S.t_k_years(par_off, CFG4["a0_au"]), N_SEGMENTS, TILT_DEG)
    j_vec = run["l"] + run["s"]
    j0 = np.linalg.norm(j_vec[0])
    periods = precession_periods(par_off, run, LK_OFF_SPAN_YR)
    out = {"lk_off_max_dJ_over_J": float(np.abs(j_vec - j_vec[0]).max() / j0),
           "lk_off_max_d_absS": float(np.abs(np.linalg.norm(run["s"], axis=1) / np.linalg.norm(run["s"][0]) - 1).max()),
           "lk_off_max_d_absL": float(np.abs(np.linalg.norm(run["l"], axis=1) / np.linalg.norm(run["l"][0]) - 1).max()),
           "lk_off_S_over_L": float(np.linalg.norm(run["s"][0]) / np.linalg.norm(run["l"][0])),
           "lk_off_precession_periods": periods}
    par = lk_params()
    run = track(par, THETA_LK_DEG, P_STAR_STRONG_DAYS, N_SEGMENTS * SEGMENT_TK, N_SEGMENTS, TILT_DEG)
    j_vec = run["l"] + run["s"]
    h = run["h_q"] + run["h_sl"]
    out.update({"lk_on_max_dJz_over_Jz": float(np.abs(j_vec[:, 2] - j_vec[0, 2]).max() / abs(j_vec[0, 2])),
                "lk_on_max_dJ_over_J": float(np.abs(j_vec - j_vec[0]).max() / np.linalg.norm(j_vec[0])),
                "lk_on_max_dLz_over_Lz_NOT_conserved": float(np.abs(run["l"][:, 2] - run["l"][0, 2]).max() / abs(run["l"][0, 2])),
                "lk_on_e_max": float(run["e"].max()),
                "energy_range_HQ": float(np.ptp(run["h_q"])), "energy_range_HSL": float(np.ptp(run["h_sl"])),
                "energy_range_HQ_plus_HSL": float(np.ptp(h)),
                "energy_conservation_ratio": float(np.ptp(h) / np.ptp(run["h_q"]))})
    return out


def check_octupole_scaling() -> dict:
    """Quadrupole-level invariants j_z and H_Q drift at order eps_oct (ratio of drifts for eps, eps/2, eps/4)."""
    cfg_par = lk_params(fb_mult=0.0, e_b=OCT_E_B, octupole=True)
    eps0 = float(cfg_par[S.P_EOCT] * CFG4["a0_au"])
    drift_j, drift_h = [], []
    for factor in OCT_FACTORS:
        par = cfg_par.copy()
        par[S.P_EOCT] *= factor
        run = track(par, OCT_THETA_DEG, P_STAR_STRONG_DAYS, OCT_CYCLES_TK, N_SEGMENTS, 0.0)
        jz = np.array([j_z(l_vec, e) for l_vec, e in zip(run["l"], run["e"], strict=True)])
        drift_j.append(float(np.ptp(jz)))
        drift_h.append(float(np.ptp(run["h_q"])))
    return {"eps_oct": [eps0 * f for f in OCT_FACTORS], "drift_jz": drift_j, "drift_HQ": drift_h,
            "drift_jz_ratios": [drift_j[i] / drift_j[i + 1] for i in range(len(drift_j) - 1)],
            "drift_HQ_ratios": [drift_h[i] / drift_h[i + 1] for i in range(len(drift_h) - 1)]}


def j_z(l_vec: np.ndarray, e_abs: float) -> float:
    """z-component of j = sqrt(1 - e^2) L_hat."""
    return math.sqrt(1.0 - e_abs**2) * l_vec[2] / np.linalg.norm(l_vec)


def literal_octupole(y: np.ndarray, par: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """A16 Eqs. Lvec and evec octupole terms transcribed AS PRINTED (Lvec's X and Y carry (e.Lb) where evec has (e.eb))."""
    lmag, jm, a = S.orbit_elements(y, par)
    e = y[3:6]
    j = jm * y[:3] / lmag
    b, q = np.array([0.0, 0.0, 1.0]), np.array([par[S.P_QX], par[S.P_QY], 0.0])
    n = math.sqrt(par[S.P_G] * par[S.P_MTOT] / a**3)
    tk_inv = par[S.P_MB] / par[S.P_MTOT] * (a / par[S.P_AB]) ** 3 * n
    eps = par[S.P_EOCT] * a
    eb, jb, eq, jq, e2 = e @ b, j @ b, e @ q, j @ q, e @ e
    w = 1.6 * e2 - 0.2 - 7.0 * eb**2 + jb**2
    d_l = -LITERAL_FACTOR * eps * lmag / jm * tk_inv * (np.cross(2 * (eb * jb + eb * jq) * j + 2 * (jq * jb - 7 * eb * eb) * e, b)
                                                       + np.cross(2 * eb * jb * j + w * e, q))
    d_e = -LITERAL_FACTOR * eps * tk_inv * (np.cross(2 * eb * jb * e + w * j, q)
                                           + np.cross(2 * (eq * jb + eb * jq) * e + 2 * (jb * jq - 7 * eb * eq) * j, b)
                                           + 3.2 * eq * np.cross(j, e))
    return d_l, d_e


def constraint_residuals(y: np.ndarray, d_l: np.ndarray, d_e: np.ndarray) -> tuple[float, float]:
    """a is constant: d ln|l|^2 = d ln(1 - e^2); and l.e = 0 stays: both relative to the size of the derivative."""
    l_vec, e = y[:3], y[3:6]
    scale = np.linalg.norm(d_e)
    a_resid = ((l_vec @ d_l) / (l_vec @ l_vec) + (e @ d_e) / (1.0 - e @ e)) / scale
    return float(a_resid), float((d_l @ e + l_vec @ d_e) / scale)


def random_state(rng: np.random.Generator) -> np.ndarray:
    l_hat = rng.normal(size=3)
    l_hat /= np.linalg.norm(l_hat)
    e_hat = np.cross(l_hat, rng.normal(size=3))
    e_hat /= np.linalg.norm(e_hat)
    e_abs = rng.uniform(0.1, 0.9)
    y = np.zeros(S.N_EXT)
    y[:3], y[3:6] = math.sqrt(1.0 - e_abs**2) * l_hat, e_abs * e_hat
    return y


def check_octupole_consistency() -> dict:
    par = lk_params(e_b=0.7, octupole=True)
    rng = np.random.default_rng(1)
    res = {"mirrored_a_resid": [], "mirrored_le_resid": [], "printed_a_resid": [], "printed_le_resid": []}
    for _ in range(N_RANDOM_STATES):
        y = random_state(rng)
        lmag, j, a = S.orbit_elements(y, par)
        n = math.sqrt(par[S.P_G] * par[S.P_MTOT] / a**3)
        tk_inv = par[S.P_MB] / par[S.P_MTOT] * (a / par[S.P_AB]) ** 3 * n
        out = np.zeros(6)
        S.oct_terms(y, par, lmag, j, a, tk_inv, out)
        res["mirrored_a_resid"].append(constraint_residuals(y, out[:3], out[3:])[0])
        res["mirrored_le_resid"].append(constraint_residuals(y, out[:3], out[3:])[1])
        printed = literal_octupole(y, par)
        res["printed_a_resid"].append(constraint_residuals(y, *printed)[0])
        res["printed_le_resid"].append(constraint_residuals(y, *printed)[1])
    return {k: float(np.abs(v).max()) for k, v in res.items()}


def check_skumanich() -> dict:
    par = lk_params(lk=False, alpha_mb=ALPHA_MB_YR)
    y0 = S.initial_state_ext(par, CFG4, 89.0, P_STAR_STRONG_DAYS)
    om0 = S.omega_star_of_period(P_STAR_STRONG_DAYS)
    rel = []
    for t in SKUMANICH_TIMES_YR:
        y = S.run_ext(y0, par, t, ev_max=10)["y"]
        om = np.linalg.norm(y[S.N_ORBIT:S.N_ORBIT + S.N_SPIN]) * par[S.P_LSCALE] / par[S.P_ISTAR]
        rel.append(float(om / (om0 / math.sqrt(1.0 + 2.0 * ALPHA_MB_YR * om0**2 * t)) - 1.0))
    return {"times_yr": list(SKUMANICH_TIMES_YR), "rel_error": rel, "omega_ratio_last": float(om / om0)}


def theta_members(table: dict, i: int) -> tuple[np.ndarray, np.ndarray]:
    """(theta_f of the 11 spins run with fb_mult = 0, step 5 theta_f of the same 11) for table row i."""
    m_p, p_star, a_b, theta = (float(table[k][i]) for k in W.INPUT_FIELDS)
    cfg = W.R.system_config(CFG4, m_p, a_b)
    par = S.extend_params(S.make_params(cfg), fb_mult=0.0)
    periods, _ = W.member_periods(p_star)
    got = np.array([S.theta_sl_ext(W.run_member(par, cfg, theta, p, twin=k == len(periods) - 1)["y"])
                    for k, p in enumerate(periods)])
    return got, np.append(table["theta_f"][i], table["theta_f_twin"][i])


def check_feedback_zero() -> dict:
    table = W.step5_table()
    labels = np.load(LABELS_PATH)["label"]
    regular = [i for i in range(labels.size) if labels[i] != CHAOTIC_LABEL][:N_SYSTEMS]
    out = {}
    for name, rows in (("first_non_chaotic", regular), ("first_any", list(range(N_SYSTEMS)))):
        diffs = []
        for i in rows:
            got, ref = theta_members(table, i)
            diffs.append(float(np.abs(got - ref).max()))
        out[name] = {"rows": rows, "labels": [str(labels[i]) for i in rows], "max_abs_diff_deg": diffs}
    return out


def single_spin_theta(cfg: dict, theta: float, p_days: float, rtol: float, atol: float) -> float:
    """theta_f from the existing integrator (run_diag, one spin, no feedback)."""
    par = S.make_params(cfg)
    y0 = S.initial_state(par, cfg, theta, 1)
    res = S.run_diag(y0, par, [p_days], cfg["stop_t_gyr"] * 1e9, rtol=rtol, atol=atol, ev_max=W.CFG5["ev_max"])
    return float(S.theta_sl_final(res["y"], 1)[0])


def ext_theta(cfg: dict, theta: float, p_days: float, rtol: float, atol: float) -> float:
    par = S.extend_params(S.make_params(cfg), fb_mult=0.0)
    y0 = S.initial_state_ext(par, cfg, theta, p_days)
    return S.theta_sl_ext(S.run_ext(y0, par, cfg["stop_t_gyr"] * 1e9, rtol=rtol, atol=atol, ev_max=W.CFG5["ev_max"])["y"])


def check_feedback_zero_same_algorithm(rows: list) -> dict:
    """Member 0 of each non-chaotic row: run_ext(fb_mult = 0) against run_diag with ONE spin, at the production tolerance and
    at TIGHT_FACTOR times tighter (both converge to the same solution, so the difference shrinks)."""
    table = W.step5_table()
    out = {"production": [], "tight": []}
    for i in rows:
        m_p, p_star, a_b, theta = (float(table[k][i]) for k in W.INPUT_FIELDS)
        cfg = W.R.system_config(CFG4, m_p, a_b)
        p_days = W.member_periods(p_star)[0][0]
        for key, factor in (("production", 1.0), ("tight", TIGHT_FACTOR)):
            rtol, atol = S.RTOL_DEFAULT * factor, S.ATOL_DEFAULT * factor
            out[key].append(abs(ext_theta(cfg, theta, p_days, rtol, atol) - single_spin_theta(cfg, theta, p_days, rtol, atol)))
    return {k: [float(v) for v in vals] for k, vals in out.items()}


def main() -> None:
    result = {"skumanich": check_skumanich(), "octupole_consistency": check_octupole_consistency(),
              "conservation": check_conservation(), "octupole_scaling": check_octupole_scaling(),
              "feedback_zero_vs_step5": (zero := check_feedback_zero()),
              "feedback_zero_same_algorithm": check_feedback_zero_same_algorithm(zero["first_non_chaotic"]["rows"])}
    OUT.mkdir(parents=True, exist_ok=True)
    (OUT / "validation.json").write_text(json.dumps(result, indent=2) + "\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
