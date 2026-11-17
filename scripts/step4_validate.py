"""Step 4 validation 1-8 for lkspin.secular. Prints numbers, writes results/step4/validation.json and orbit_89.png.

Run:  uv run python scripts/step4_validate.py [--skip-plot]
"""

import argparse
import json
import math
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import brentq

from lkspin import secular as S
from lkspin.lk_orbit import cos_theta_lb_of_x, make_orbit, x_of_t
from lkspin.spin import ATOL_DEFAULT, RTOL_DEFAULT, run_cycles
from lkspin.spin import initial_state as spin_initial_state

CFG = S.load_config()
OUT_DIR = S.ROOT / CFG["results"]
THETA_89 = 89.0
THETA_87 = 87.0
THETA_85 = 85.0
RTOL_TIGHT = 1e-11
ATOL_TIGHT = 1e-13
REC_MAX = 4_000_000
REC_STRIDE = 3
N_LK_CYCLES_V1 = 2
N_CYCLES_V3 = 6
N_MIN_V4 = 50
EPS_V4 = (0.015, 0.03)
P_STAR_SPOTS = (30.0, 29.87, 7.07, 1.67)
P_STAR_TOL = (30.0, 7.0, 1.67)
TOL_FACTOR = 10.0
DEG = 180.0 / math.pi
MYR = 1e6
STOP_T_YR = CFG["stop_t_gyr"] * 1e9
PLOT_T_RANGE_MYR = (1.0, 3000.0)


def omega_p0() -> float:
    return 2.0 * math.pi * S.HOURS_PER_YR / CFG["p_planet0_hours"]


def v1_pure_lk() -> dict:
    par = S.make_params(CFG, srf=False, tides=False)
    tk = S.t_k_years(par, CFG["a0_au"])
    orb = make_orbit(THETA_85, CFG["e0"])
    y0 = S.initial_state(par, CFG, THETA_85, 1)
    r = S.run(y0, par, [30.0], N_LK_CYCLES_V1 * orb.period * tk, RTOL_TIGHT, ATOL_TIGHT, rec_max=REC_MAX)
    rec = r["rec"]
    t = rec[:, 0] / tk
    x_an = np.array([x_of_t(ti, orb.params()) for ti in t])
    cos_an = np.array([cos_theta_lb_of_x(xi, orb.params()) for xi in x_an])
    return {"t_k_yr": tk, "max_dx": float(np.abs(1 - rec[:, 2] ** 2 - x_an).max()),
            "max_dcos_theta_lb": float(np.abs(np.cos(rec[:, 3]) - cos_an).max()), "steps": r["steps"]}


def emax_a16(par: np.ndarray, theta_deg: float, bulge_factor: float = 1.0) -> float:
    """Solve A16 Eq. emax. bulge_factor rescales eps_Tide (1.0 as printed)."""
    a0 = CFG["a0_au"]
    eg, et, er = S.srf_epsilons(par, a0, omega_p0())
    c2 = math.cos(math.radians(theta_deg)) ** 2

    def f(e):
        j = math.sqrt(1 - e * e)
        lhs = (eg * (1 / j - 1) + bulge_factor * et / 15 * ((1 + 3 * e * e + 3 / 8 * e**4) / j**9 - 1)
               + er / 3 * (1 / j**3 - 1))
        return lhs - 9 * e * e / (8 * j * j) * (j * j - 5 / 3 * c2)

    return brentq(f, 0.02, 1 - 1e-10)


def first_cycle_emax(par: np.ndarray, theta_deg: float) -> float:
    tk = S.t_k_years(par, CFG["a0_au"])
    y0 = S.initial_state(par, CFG, theta_deg, 1)
    r = S.run(y0, par, [30.0], 8 * tk, RTOL_TIGHT, ATOL_TIGHT, rec_max=REC_MAX)
    return float(r["rec"][:, 2].max())


def v2_emax() -> dict:
    out = {}
    for label, bulge in (("A16_eom_7.5", S.BULGE_COEF_A16_EOM), ("potential_15", 2 * S.BULGE_COEF_A16_EOM)):
        par = S.make_params(CFG, tides=False, bulge_coef=bulge)
        for th in (THETA_89, THETA_87):
            sim, eq = first_cycle_emax(par, th), emax_a16(par, th)
            out[f"{label}_theta{int(th)}"] = {"sim_emax": sim, "eq_emax": eq, "sim_1me": 1 - sim, "eq_1me": 1 - eq}
    return out


def conservation_srf_only() -> dict:
    par = S.make_params(CFG, tides=False)
    tk = S.t_k_years(par, CFG["a0_au"])
    y = S.initial_state(par, CFG, THETA_89, 1)
    om = np.array([S.omega_star_of_period(1.0)])
    y0 = y.copy()
    t, h = 0.0, S.H_INITIAL
    worst = {"e_dot_l": 0.0, "lz": 0.0, "a_rel": 0.0, "norm_S": 0.0}
    for k in range(1, 20 * N_CYCLES_V3 + 1):
        t, h, *_ = S.advance(y, t, k * tk / 3.0, h, par, om, RTOL_TIGHT, ATOL_TIGHT, 10**9, np.zeros((0, 7)), 1)
        lmag, j, a = S.orbit_elements(y, par)
        worst["e_dot_l"] = max(worst["e_dot_l"], abs(np.dot(y[3:6], y[0:3]) / lmag))
        worst["lz"] = max(worst["lz"], abs(y[2] - y0[2]))
        worst["a_rel"] = max(worst["a_rel"], abs(a / CFG["a0_au"] - 1))
        worst["norm_S"] = max(worst["norm_S"], abs(np.linalg.norm(y[7:10]) - 1))
    return worst


def v3_tides_J(rec: np.ndarray) -> dict:
    """A16 Eq. kozai_constant_diss: J = sqrt(a) j cos(theta_lb) over the full migration."""
    a, e, th = rec[:, 1], rec[:, 2], rec[:, 3]
    j_const = np.sqrt(a) * np.sqrt(1 - e**2) * np.cos(th)
    j0 = j_const[0]
    return {"J0": float(j0), "max_rel_drift": float(np.abs(j_const / j0 - 1).max()),
            "final_rel_drift": float(j_const[-1] / j0 - 1), "end_a": float(a[-1])}


def v4_spin_crosscheck() -> dict:
    par = S.make_params(CFG, srf=False, tides=False)
    a0 = CFG["a0_au"]
    tk = S.t_k_years(par, a0)
    orb = make_orbit(THETA_85, CFG["e0"])
    out = {}
    for eps in EPS_V4:
        om_star = orb.alpha0(eps) / tk * a0**3 / par[S.P_STAR_COEF]
        p_days = 2 * math.pi * S.DAYS_PER_YR / om_star
        y, t, h = S.initial_state(par, CFG, THETA_85, 1), 0.0, S.H_INITIAL
        sec = []
        for k in range(1, N_MIN_V4 + 1):
            t, h, *_ = S.advance(y, t, k * orb.period * tk, h, par, np.array([om_star]), RTOL_TIGHT, ATOL_TIGHT,
                                 10**9, np.zeros((0, 7)), 1)
            sec.append(S.theta_sl_final(y, 1)[0])
        ref = run_cycles(spin_initial_state(1.0, 0.0, orb.params()), N_MIN_V4, orb.period, orb.alpha0(eps),
                         orb.params(), RTOL_DEFAULT, ATOL_DEFAULT)
        ref_th = np.degrees(np.arccos(np.clip(ref[:, 0], -1, 1)))
        out[f"eps_{eps}"] = {"P_star_days": p_days, "max_diff_deg": float(np.abs(np.array(sec) - ref_th).max()),
                             "theta_sl_range_deg": [float(ref_th.min()), float(ref_th.max())]}
    return out


def first_cycle_slice(rec: np.ndarray) -> slice:
    """Indices of the first LK cycle: t = 0 to the second eccentricity minimum."""
    e = rec[:, 2]
    interior = np.where((e[1:-1] < e[:-2]) & (e[1:-1] <= e[2:]))[0] + 1
    return slice(0, interior[0] + 1)


def v5_shape_numbers(par: np.ndarray) -> dict:
    """N_max,0 = alpha_bar / n_e and A_bar,0 = alpha_bar / (2 <Omega_pl sin theta_lb>) over the first LK cycle."""
    tk = S.t_k_years(par, CFG["a0_au"])
    y0 = S.initial_state(par, CFG, THETA_89, 1)
    r = S.run(y0, par, [1.0], 20 * tk, RTOL_TIGHT, ATOL_TIGHT, rec_max=REC_MAX)
    rec = r["rec"]
    sl = first_cycle_slice(rec)
    t, a, e, th, node = (rec[sl, i] for i in range(5))
    period = t[-1] - t[0]
    n_e = 2 * math.pi / period
    j = np.sqrt(1 - e**2)
    alpha = par[S.P_STAR_COEF] * S.omega_star_of_period(1.0) / (a**3 * j**3)  # P_star = 1 d
    alpha_bar = np.trapezoid(alpha, t) / period
    alpha_bar_half = np.trapezoid(alpha[::2], t[::2]) / (t[::2][-1] - t[0])
    dnode = np.diff(np.unwrap(node))
    mean_opl_sin = abs(np.sum(0.5 * (np.sin(th[1:]) + np.sin(th[:-1])) * dnode)) / period
    return {"period_yr": period, "n_e_per_yr": n_e, "period_over_tk": period / tk,
            "Nmax0_times_Pstar": alpha_bar / n_e, "Abar0_times_Pstar": alpha_bar / (2 * mean_opl_sin),
            "Abar_over_Nmax": (alpha_bar / (2 * mean_opl_sin)) / (alpha_bar / n_e),
            "alpha_bar_decimation_change": float(alpha_bar_half / alpha_bar - 1),
            "paper": {"Abar0_times_Pstar": 25.5, "Nmax0_times_Pstar": 3.34, "ratio": 7.63}}


def v6_tolerance(par: np.ndarray) -> dict:
    out = {}
    y0 = S.initial_state(par, CFG, THETA_89, 1)
    for p in P_STAR_TOL:
        base = S.run(y0, par, [p], STOP_T_YR)
        tight = S.run(y0, par, [p], STOP_T_YR, RTOL_DEFAULT / TOL_FACTOR, ATOL_DEFAULT / TOL_FACTOR)
        th0, th1 = S.theta_sl_final(base["y"], 1)[0], S.theta_sl_final(tight["y"], 1)[0]
        out[str(p)] = {"theta_sl_f": th0, "theta_sl_f_tight": th1, "change_deg": th1 - th0,
                       "t_end_rel_change": tight["t"] / base["t"] - 1, "steps": base["steps"], "steps_tight": tight["steps"]}
    return out


def cycle_extrema(rec: np.ndarray) -> dict:
    """Per-cycle maxima of e (min 1-e) and theta_lb there."""
    e = rec[:, 2]
    idx = np.where((e[1:-1] >= e[:-2]) & (e[1:-1] > e[2:]))[0] + 1
    idx = idx[e[idx] > 0.5]
    return {"n_cycles": len(idx), "min_1me": (1 - e[idx]).tolist(), "t_myr": (rec[idx, 0] / MYR).tolist(),
            "theta_lb_at_emax_deg": (rec[idx, 3] * DEG).tolist()}


def v7_timeline(rec: np.ndarray, t_end: float) -> dict:
    ex = cycle_extrema(rec)
    m = np.array(ex["min_1me"])
    th = np.array(ex["theta_lb_at_emax_deg"])
    a = rec[:, 1]
    t = rec[:, 0]
    return {"t_end_myr": t_end / MYR, "n_cycles": ex["n_cycles"], "min_1me_overall": float(m.min()),
            "min_1me_first": float(m[0]), "min_1me_median": float(np.median(m)), "min_1me_last": float(m[-1]),
            "theta_lb_at_emax_first": float(th[0]), "theta_lb_at_emax_median": float(np.median(th)),
            "theta_lb_at_emax_last": float(th[-1]), "theta_lb_range_deg": [float(rec[:, 3].min() * DEG),
                                                                          float(rec[:, 3].max() * DEG)],
            "t_a_below_1.0_myr": float(t[np.argmax(a < 1.0)] / MYR), "t_a_below_0.5_myr": float(t[np.argmax(a < 0.5)] / MYR),
            "t_first_emax_myr": ex["t_myr"][0]}


def plot_orbit(rec: np.ndarray, p_stars: list, out_path) -> None:
    fig, ax = plt.subplots(2, 3, figsize=(13, 7))
    t = rec[:, 0] / MYR
    ax[0, 0].loglog(t[1:], 1 - rec[1:, 2], "k", lw=0.3)
    ax[0, 0].set_ylabel("1 - e")
    ax[0, 0].set_ylim(5e-3, 1.2)
    ax[0, 1].semilogx(t[1:], rec[1:, 3] * DEG, "k", lw=0.3)
    ax[0, 1].set_ylabel(r"$\theta_{lb}$ (deg)")
    ax[0, 1].set_ylim(0, 180)
    ax[0, 2].loglog(t[1:], rec[1:, 1], "k")
    ax[0, 2].set_ylabel("a (AU)")
    colors = {30.0: "r", 7.07: "b", 1.67: "purple"}
    for p in (30.0, 7.07, 1.67):
        k = p_stars.index(p)
        ax[1, (30.0, 7.07, 1.67).index(p)].semilogx(t[1:], rec[1:, S.N_REC_FIXED + k] * DEG, color=colors[p], lw=0.4)
        ax[1, (30.0, 7.07, 1.67).index(p)].set_title(f"P* = {p} d")
        ax[1, (30.0, 7.07, 1.67).index(p)].set_ylim(0, 180)
        ax[1, (30.0, 7.07, 1.67).index(p)].set_ylabel(r"$\theta_{sl}$ (deg)")
    for a_ in ax.ravel():
        a_.set_xlabel("t (Myr)")
        a_.set_xlim(*PLOT_T_RANGE_MYR)
    fig.suptitle(r"$\theta_{lb,0}=89^\circ$, 5 $M_J$, a0 = 1.5 AU (cf. S16 Fig. 9)")
    fig.tight_layout()
    fig.savefig(out_path, dpi=130)


def v8_spots(par: np.ndarray, a_bar_coef: float, nmax_coef: float) -> dict:
    y0 = S.initial_state(par, CFG, THETA_89, 1)
    out = {}
    for p in P_STAR_SPOTS:
        r = S.run(y0, par, [p], STOP_T_YR)
        a0, n0 = a_bar_coef / p, nmax_coef / p
        regime = "I" if a0 < 1 else ("II" if n0 < 0.5 else "IV")
        out[str(p)] = {"theta_sl_f": float(S.theta_sl_final(r["y"], 1)[0]), "regime_by_A0_Nmax0": regime,
                       "Abar0": a0, "Nmax0": n0, "status": S.STATUS_NAMES[r["status"]], "t_end_myr": r["t"] / MYR}
    return out


def v_frozen(par: np.ndarray) -> dict:
    """Does theta_sl change if we stop at a < 0.05 AU instead of 0.1 AU?"""
    out = {}
    y0 = S.initial_state(par, CFG, THETA_89, len(P_STAR_SPOTS))
    p2 = par.copy()
    p2[S.P_STOPA] = par[S.P_STOPA] / 2
    r1 = S.run(y0, par, list(P_STAR_SPOTS), STOP_T_YR)
    r2 = S.run(y0, p2, list(P_STAR_SPOTS), STOP_T_YR)
    for k, p in enumerate(P_STAR_SPOTS):
        out[str(p)] = {"stop_0.1": float(S.theta_sl_final(r1["y"], len(P_STAR_SPOTS))[k]),
                       "stop_0.05": float(S.theta_sl_final(r2["y"], len(P_STAR_SPOTS))[k]), "t_0.05_myr": r2["t"] / MYR}
    return out


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--skip-plot", action="store_true")
    args = ap.parse_args()
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    par = S.make_params(CFG)
    res = {}
    for name, fn in (("v1_pure_lk", v1_pure_lk), ("v2_emax", v2_emax), ("v3_conservation_srf_only", conservation_srf_only)):
        res[name] = fn()
        print(name, json.dumps(res[name], indent=1))
    p_list = list(P_STAR_SPOTS)
    y0 = S.initial_state(par, CFG, THETA_89, len(p_list))
    t0 = time.time()
    run = S.run(y0, par, p_list, STOP_T_YR, rec_max=REC_MAX, stride=REC_STRIDE)
    print("full migration, K=4:", time.time() - t0, "s", run["steps"], "steps", S.STATUS_NAMES[run["status"]])
    res["v3_tides_J"] = v3_tides_J(run["rec"])
    res["v4_spin_crosscheck"] = v4_spin_crosscheck()
    res["v5_shape"] = v5_shape_numbers(par)
    res["v6_tolerance"] = v6_tolerance(par)
    res["v7_timeline"] = v7_timeline(run["rec"], run["t"])
    res["v8_spots"] = v8_spots(par, res["v5_shape"]["Abar0_times_Pstar"], res["v5_shape"]["Nmax0_times_Pstar"])
    res["frozen_check"] = v_frozen(par)
    for k in list(res)[3:]:
        print(k, json.dumps(res[k], indent=1))
    (OUT_DIR / "validation.json").write_text(json.dumps(res, indent=1))
    if not args.skip_plot:
        plot_orbit(run["rec"], p_list, OUT_DIR / "orbit_89.png")


if __name__ == "__main__":
    main()
