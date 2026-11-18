"""Prep job 2 (SQ4): chaos maps in the non-dissipative, SRF-modified setting (config [prep_sq4]).

Setting: lkspin.secular with SRFs on, tides off, a fixed at a0, no feedback, e_b = 0 (the orbit is the same for every
P_star, so one orbit per theta_lb0 carries all spins). Spins are sampled where the eccentricity has a minimum: the orbit's
minima times come from run_diag once, then every spin set is integrated with advance(), landing exactly on those times.

Stages (each resumable; chunk files in results/prep_sq4/chunks/):
  ftle      (a) aligned + twin (1e-9 rad) at 300 P_star in [1, 40] d, theta_lb0 = 87 and 89 deg, 1000 cycles
  sections  (b) 40 initial conditions x 8 P_star at 87 deg, 800 cycles, (p, phi) at e minima
  omega     (c) FTLE against Omega_hat_L at fixed N_max,0
  validate  SRFs off at 85 deg: FTLE vs eps against the step 3 bifurcation (S15 Fig. 12), and against lkspin.spin
  plot      figures and summary.json
Run:  uv run python scripts/prep_sq4.py --stage all
"""

import argparse
import json
import math
import multiprocessing as mp
import os
import time
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import prep_common as C
from lkspin import regime as R
from lkspin import secular as S
from lkspin import spin as SP
from lkspin.lk_orbit import make_orbit

CFG4 = S.load_config()
CFG = C.load_table("prep_sq4")
OUT = S.ROOT / CFG["results"]
CHUNKS = OUT / "chunks"
DPI = 150
RTOL, ATOL = CFG["rtol"], CFG["atol"]
SECTION_DOT = 1.0
TWIN = CFG["twin_angle_rad"]
EPS_TOL_CHECK = 10.0           # the tolerance check runs this much tighter
N_TOL_CHECK = 6                # P_star values in the tolerance check
SPIN_ORBIT_ALIGNED_P = 1.0
LOG_FLOOR = 1e-12


# ------------------------------------------------------------------ file helpers
def save_atomic(path: Path, **arrays) -> None:
    tmp = path.with_suffix(f".{os.getpid()}.tmp.npz")
    np.savez(tmp, **arrays)
    tmp.replace(path)


def orbit_times(theta: float, par: np.ndarray, tag: str, n_cycles: int) -> np.ndarray:
    """Cached e-minimum times for the orbit of (par, theta)."""
    path = OUT / f"orbit_{tag}.npz"
    if path.exists():
        return np.load(path)["times"][:n_cycles]
    times = C.orbit_minima_times(par, CFG4, theta, n_cycles)
    save_atomic(path, times=times, theta=theta)
    return times


# ------------------------------------------------------------------ one FTLE integration (aligned + twin per P_star)
def ftle_integrate(par: np.ndarray, theta: float, p_list: np.ndarray, times: np.ndarray, rtol: float = RTOL,
                   atol: float = ATOL) -> dict:
    """Aligned and twin spin for every P in p_list on one orbit. Returns per-P arrays and the series."""
    n_p = p_list.size
    y0 = S.initial_state(par, CFG4, theta, 2 * n_p)
    for i in range(n_p):
        y0 = _twin_in_place(y0, 2 * i + 1)
    periods = np.repeat(p_list, 2)
    states, status = C.run_to_times(y0, par, periods, times, rtol, atol)
    spins = C.spin_block(states, 2 * n_p)
    sep = C.angle_between(spins[:, 0::2], spins[:, 1::2])          # [n_cycles, n_p]
    theta_al = C.theta_sl_deg_series(states, 2 * n_p)[:, 0::2]     # [n_cycles, n_p]
    return summarize_ftle(p_list, sep, theta_al, status)


def _twin_in_place(y0: np.ndarray, spin_index: int) -> np.ndarray:
    return R.twin_state(y0, spin_index, TWIN)


def summarize_ftle(p_list: np.ndarray, sep: np.ndarray, theta_al: np.ndarray, status: int) -> dict:
    n_p = p_list.size
    fits = [C.ftle_fit(sep[:, i], TWIN, CFG["sat_rad"], CFG["min_fit_points"]) for i in range(n_p)]
    last = theta_al[-CFG["spread_cycles"]:]
    return {"P_star": p_list, "ftle": np.array([f["ftle"] for f in fits]), "k_sat": np.array([f["k_sat"] for f in fits]),
            "n_fit": np.array([f["n_fit"] for f in fits]), "theta_std_last": np.nanstd(last, axis=0),
            "theta_range_last": np.nanmax(last, axis=0) - np.nanmin(last, axis=0),
            "sep_final": sep[-1], "status": np.full(n_p, status), "sep": sep.T.astype(np.float32),
            "theta_sl": theta_al.T.astype(np.float32)}


# ------------------------------------------------------------------ (a) FTLE scan
def p_grid() -> np.ndarray:
    lo, hi = CFG["p_star_range_days"]
    return np.geomspace(lo, hi, CFG["n_p_star"])


def ftle_chunk(args: tuple) -> str:
    theta, idx, p_list = args
    path = CHUNKS / f"ftle_th{theta:g}_c{idx:03d}.npz"
    if path.exists():
        return f"skip {path.name}"
    t0 = time.time()
    par = C.nondissipative_params(CFG4)
    times = orbit_times(theta, par, f"th{theta:g}", CFG["n_cycles"])
    res = ftle_integrate(par, theta, np.asarray(p_list), times)
    save_atomic(path, wall_s=time.time() - t0, **res)
    return f"done {path.name} {time.time() - t0:.1f}s"


def ftle_items() -> list:
    grid, n = p_grid(), CFG["chunk_p"]
    return [(th, i, grid[s:s + n].tolist()) for th in CFG["theta_lb0_deg"] for i, s in enumerate(range(0, grid.size, n))]


def merge_ftle(theta: float) -> dict:
    parts = [np.load(CHUNKS / f"ftle_th{theta:g}_c{i:03d}.npz") for i in range(math.ceil(CFG["n_p_star"] / CFG["chunk_p"]))]
    keys = ("P_star", "ftle", "k_sat", "n_fit", "theta_std_last", "theta_range_last", "sep_final", "status")
    merged = {k: np.concatenate([p[k] for p in parts]) for k in keys}
    series = {"sep": np.concatenate([p["sep"] for p in parts]), "theta_sl": np.concatenate([p["theta_sl"] for p in parts])}
    merged["wall_s_total"] = np.array(sum(float(p["wall_s"]) for p in parts))
    merged["chaotic"] = (merged["ftle"] > CFG["ftle_chaos_per_cycle"]).astype(int)
    np.savez(OUT / f"ftle_scan_{theta:g}.npz", **merged, twin_angle_rad=TWIN, sat_rad=CFG["sat_rad"])
    np.savez_compressed(OUT / f"series_ftle_{theta:g}.npz", **series)
    return merged


def runs_of_true(mask: np.ndarray, x: np.ndarray) -> list:
    """Contiguous x-intervals [lo, hi] where mask is true (x ascending)."""
    out, start = [], None
    for i, m in enumerate(mask):
        if m and start is None:
            start = i
        if not m and start is not None:
            out.append([float(x[start]), float(x[i - 1])])
            start = None
    if start is not None:
        out.append([float(x[start]), float(x[-1])])
    return out


def ftle_bands(m: dict) -> dict:
    x = m["P_star"]
    return {"ftle_gt_threshold_P_ranges_days": runs_of_true(m["ftle"] > CFG["ftle_chaos_per_cycle"], x),
            "spread_gt_threshold_P_ranges_days": runs_of_true(m["theta_std_last"] > CFG["spread_chaos_deg"], x),
            "n_ftle_chaotic": int((m["ftle"] > CFG["ftle_chaos_per_cycle"]).sum()),
            "n_spread": int((m["theta_std_last"] > CFG["spread_chaos_deg"]).sum()),
            "n_both": int(((m["ftle"] > CFG["ftle_chaos_per_cycle"]) & (m["theta_std_last"] > CFG["spread_chaos_deg"])).sum()),
            "ftle_min": float(m["ftle"].min()), "ftle_max": float(m["ftle"].max())}


def stage_ftle(pool) -> dict:
    CHUNKS.mkdir(parents=True, exist_ok=True)
    for th in CFG["theta_lb0_deg"]:
        orbit_times(th, C.nondissipative_params(CFG4), f"th{th:g}", CFG["n_cycles"])
    for msg in pool.imap_unordered(ftle_chunk, ftle_items(), chunksize=1):
        print(msg, flush=True)
    return {f"{th:g}": ftle_bands(merge_ftle(th)) for th in CFG["theta_lb0_deg"]}


# ------------------------------------------------------------------ tolerance check of the FTLE method
def stage_tolerance() -> dict:
    """Repeat a few P_star at 10x tighter tolerances and compare FTLE and final theta_sl (one orbit, 87 deg)."""
    theta = 87.0
    par = C.nondissipative_params(CFG4)
    times = orbit_times(theta, par, f"th{theta:g}", CFG["n_cycles"])
    scan = np.load(OUT / f"ftle_scan_{theta:g}.npz")
    pick = np.linspace(0, scan["P_star"].size - 1, N_TOL_CHECK).astype(int)
    tight = ftle_integrate(par, theta, scan["P_star"][pick], times, RTOL / EPS_TOL_CHECK, ATOL / EPS_TOL_CHECK)
    return {"theta": theta, "P_star": scan["P_star"][pick].tolist(), "ftle_default": scan["ftle"][pick].tolist(),
            "ftle_tight": tight["ftle"].tolist(), "theta_std_default": scan["theta_std_last"][pick].tolist(),
            "theta_std_tight": tight["theta_std_last"].tolist()}


# ------------------------------------------------------------------ (b) sections
def section_ics() -> list:
    p0s = np.linspace(-CFG["sections_p0_max"], CFG["sections_p0_max"], CFG["sections_n_p0"])
    return [(float(p0), float(phi0)) for phi0 in (0.0, math.pi) for p0 in p0s]


def section_chunk(p_star: float) -> str:
    path = CHUNKS / f"section_P{p_star:g}.npz"
    if path.exists():
        return f"skip {path.name}"
    t0 = time.time()
    theta = CFG["sections_theta_lb0_deg"]
    par = C.nondissipative_params(CFG4)
    times = orbit_times(theta, par, f"th{theta:g}", max(CFG["n_cycles"], CFG["sections_n_cycles"]))[:CFG["sections_n_cycles"]]
    ics = section_ics()
    y0 = S.initial_state(par, CFG4, theta, len(ics))
    for k, (p0, phi0) in enumerate(ics):
        o = S.N_ORBIT + S.N_SPIN * k
        y0[o:o + S.N_SPIN] = C.initial_spin_vector(theta, p0, phi0)
    states, status = C.run_to_times(y0, par, [p_star] * len(ics), times, RTOL, ATOL)
    lvec = np.broadcast_to(states[:, None, 0:3], (times.size, len(ics), S.N_SPIN))
    p, phi = C.phase_and_p(lvec, C.spin_block(states, len(ics)))
    save_atomic(path, p=p.astype(np.float32), phi=phi.astype(np.float32), ics=np.array(ics), status=status, P_star=p_star)
    return f"done {path.name} {time.time() - t0:.1f}s"


def plot_sections() -> None:
    ps = CFG["sections_p_star_days"]
    fig, axes = plt.subplots(2, 4, figsize=(18, 9))
    for ax, p_star in zip(axes.ravel(), ps, strict=True):
        d = np.load(CHUNKS / f"section_P{p_star:g}.npz")
        ax.plot(d["phi"].ravel(), d["p"].ravel(), "k.", ms=SECTION_DOT, markeredgewidth=0, rasterized=True)
        ax.set_xlim(0.0, C.TWO_PI)
        ax.set_ylim(-1.03, 1.03)
        ax.set_title(rf"$\theta_{{lb}}^0=87^\circ$, $P_\star={p_star:g}$ d")
        ax.set_xlabel(r"$\phi$")
        ax.set_ylabel(r"$\cos\theta_{\rm sl}$")
    fig.tight_layout()
    fig.savefig(OUT / "prep_sq4_sections.png", dpi=DPI)
    plt.close(fig)


def section_fill(p_star: float) -> dict:
    """Occupancy of the section: fraction of the 40 orbits whose p range exceeds 1 (rough chaos-sea indicator)."""
    d = np.load(CHUNKS / f"section_P{p_star:g}.npz")
    rng = np.nanmax(d["p"], axis=0) - np.nanmin(d["p"], axis=0)
    return {"P_star": p_star, "median_p_range": float(np.median(rng)), "n_orbits_range_gt_1": int((rng > 1.0).sum()),
            "n_status_ok": int(d["status"] == S.ST_TIMEOUT)}


# ------------------------------------------------------------------ (c) Omega_hat_L test
def omega_thetas() -> np.ndarray:
    lo, hi = CFG["omega_theta_range_deg"]
    return np.linspace(lo, hi, CFG["omega_n_theta"])


def omega_chunk(args: tuple) -> str:
    idx, theta = args
    path = CHUNKS / f"omega_{idx:03d}.npz"
    if path.exists():
        return f"skip {path.name}"
    t0 = time.time()
    par = C.nondissipative_params(CFG4)
    shape = C.first_cycle_shape(par, CFG4, theta)
    ns = np.array(CFG["omega_n_max0"])
    p_list = shape["nmax_coef"] / ns
    times = orbit_times(theta, par, f"om{idx:03d}", CFG["n_cycles"])
    res = ftle_integrate(par, theta, p_list, times)
    keep = {k: res[k] for k in ("P_star", "ftle", "k_sat", "theta_std_last", "theta_range_last", "status")}
    save_atomic(path, N_max0=ns, theta_lb0=theta, omega_hat_l=shape["omega_hat_l"], nmax_coef=shape["nmax_coef"],
                abar_coef=shape["abar_coef"], e_max=shape["e_max"], period_yr=shape["period_yr"],
                a_bar0=shape["abar_coef"] / p_list, wall_s=time.time() - t0, **keep)
    return f"done {path.name} {time.time() - t0:.1f}s"


def merge_omega() -> dict:
    parts = [np.load(CHUNKS / f"omega_{i:03d}.npz") for i in range(CFG["omega_n_theta"])]
    keys = ("N_max0", "P_star", "ftle", "k_sat", "theta_std_last", "theta_range_last", "a_bar0", "status")
    out = {k: np.stack([p[k] for p in parts]) for k in keys}          # [n_theta, 3]
    for k in ("theta_lb0", "omega_hat_l", "nmax_coef", "abar_coef", "e_max", "period_yr"):
        out[k] = np.array([float(p[k]) for p in parts])
    np.savez(OUT / "omega_test.npz", **out)
    return out


def spearman(x: np.ndarray, y: np.ndarray) -> float:
    rx, ry = np.argsort(np.argsort(x)), np.argsort(np.argsort(y))
    return float(np.corrcoef(rx, ry)[0, 1])


def omega_summary(o: dict) -> dict:
    res = {}
    for j, n in enumerate(CFG["omega_n_max0"]):
        f, x, th = o["ftle"][:, j], o["omega_hat_l"], o["theta_lb0"]
        chaotic = f > CFG["ftle_chaos_per_cycle"]
        res[f"{n:g}"] = {"spearman_ftle_vs_omega_hat_l": spearman(x, f), "spearman_ftle_vs_theta_lb0": spearman(th, f),
                         "n_chaotic": int(chaotic.sum()), "omega_hat_l_chaotic_range": [float(x[chaotic].min()), float(x[chaotic].max())] if chaotic.any() else None,
                         "omega_hat_l_regular_range": [float(x[~chaotic].min()), float(x[~chaotic].max())] if (~chaotic).any() else None,
                         "a_bar0_range": [float(o["a_bar0"][:, j].min()), float(o["a_bar0"][:, j].max())],
                         "P_star_range_days": [float(o["P_star"][:, j].min()), float(o["P_star"][:, j].max())]}
    return res


# ------------------------------------------------------------------ validation: SRFs off, 85 deg, vs step 3 and lkspin.spin
def eps_grid() -> np.ndarray:
    lo, hi = CFG["validate_eps_range"]
    return np.geomspace(lo, hi, CFG["validate_n_eps"])


def p_of_eps(par: np.ndarray, eps: np.ndarray, theta: float) -> np.ndarray:
    """P_star (d) giving S15's epsilon = Omega_pl,0 / alpha_0 on the secular orbit with SRFs off (a0 fixed)."""
    orb = make_orbit(theta, CFG4["e0"])
    tk = S.t_k_years(par, CFG4["a0_au"])
    alpha0_per_yr = orb.alpha0(1.0) / tk / eps
    omega_star = alpha0_per_yr * CFG4["a0_au"] ** 3 / par[S.P_STAR_COEF]
    return C.TWO_PI * S.DAYS_PER_YR / omega_star


def validate_chunk(args: tuple) -> str:
    idx, eps_list = args
    path = CHUNKS / f"validate_c{idx:03d}.npz"
    if path.exists():
        return f"skip {path.name}"
    theta = CFG["validate_theta_lb0_deg"]
    par = S.make_params(CFG4, srf=False, tides=False)
    times = orbit_times(theta, par, f"val{theta:g}", CFG["n_cycles"])
    eps = np.asarray(eps_list)
    res = ftle_integrate(par, theta, p_of_eps(par, eps, theta), times)
    sp = spin_py_ftle(eps)
    save_atomic(path, eps=eps, **{k: res[k] for k in ("P_star", "ftle", "theta_std_last", "k_sat")},
                spin_ftle=sp["ftle"], spin_theta_std_last=sp["theta_std_last"], spin_k_sat=sp["k_sat"])
    return f"done {path.name}"


def spin_py_ftle(eps_list: np.ndarray) -> dict:
    """Independent check: lkspin.spin on the analytic orbit, aligned + twin, same cycles and fit."""
    orbit = make_orbit(CFG["validate_theta_lb0_deg"], CFG4["e0"])
    orb = orbit.params()
    res = {"ftle": [], "theta_std_last": [], "k_sat": []}
    for eps in eps_list:
        sep, theta_al = spin_py_series(orbit, orb, float(eps))
        fit = C.ftle_fit(sep, TWIN, CFG["sat_rad"], CFG["min_fit_points"])
        res["ftle"].append(fit["ftle"])
        res["k_sat"].append(fit["k_sat"])
        res["theta_std_last"].append(float(np.std(theta_al[-CFG["spread_cycles"]:])))
    return {k: np.array(v) for k, v in res.items()}


def spin_py_series(orbit, orb: np.ndarray, eps: float) -> tuple[np.ndarray, np.ndarray]:
    alpha0 = orbit.alpha0(eps)
    ya = SP.initial_state(1.0, 0.0, orb)
    yb = ya.copy()
    xhat = np.array([1.0, 0.0, 0.0])
    yb[0:3] = math.cos(TWIN) * ya[0:3] + math.sin(TWIN) * xhat
    n = CFG["n_cycles"]
    sep, theta_al = np.empty(n), np.empty(n)
    ha = hb = SP.H_INITIAL
    for k in range(n):
        ha = SP.integrate_to(ya, k * orbit.period, (k + 1) * orbit.period, ha, alpha0, orb, RTOL, ATOL)
        hb = SP.integrate_to(yb, k * orbit.period, (k + 1) * orbit.period, hb, alpha0, orb, RTOL, ATOL)
        sep[k] = float(C.angle_between(ya[0:3], yb[0:3]))
        theta_al[k] = math.degrees(math.acos(min(1.0, max(-1.0, SP.observe(ya, (k + 1) * orbit.period, orb)[0]))))
    return sep, theta_al


def validate_items() -> list:
    eps, n = eps_grid(), CFG["chunk_p"]
    return [(i, eps[s:s + n].tolist()) for i, s in enumerate(range(0, eps.size, n))]


def merge_validate() -> dict:
    parts = [np.load(CHUNKS / f"validate_c{i:03d}.npz") for i in range(len(validate_items()))]
    out = {k: np.concatenate([p[k] for p in parts]) for k in parts[0].files}
    f12 = np.load(S.ROOT / "results/step3/bifurcation_fig12.npz")
    last = f12["theta_deg"][:, -CFG["validate_spread_cycles"]:]
    out["step3_eps"] = f12["eps"]
    out["step3_theta_std_last"] = last.std(axis=1)
    np.savez(OUT / "validation_srf_off.npz", **out)
    return out


def validation_summary(v: dict) -> dict:
    thr_f, thr_s = CFG["ftle_chaos_per_cycle"], CFG["spread_chaos_deg"]
    eps = v["eps"]
    step3_on = np.interp(np.log(eps), np.log(v["step3_eps"]), v["step3_theta_std_last"]) > thr_s
    sec_on, spin_on = v["ftle"] > thr_f, v["spin_ftle"] > thr_f
    spread_on = v["theta_std_last"] > thr_s
    return {"secular_ftle_eps_ranges": runs_of_true(sec_on, eps), "spin_py_ftle_eps_ranges": runs_of_true(spin_on, eps),
            "step3_spread_eps_ranges": runs_of_true(step3_on, eps), "secular_spread_eps_ranges": runs_of_true(spread_on, eps),
            "agree_secular_ftle_vs_step3_spread": float((sec_on == step3_on).mean()),
            "agree_spin_ftle_vs_step3_spread": float((spin_on == step3_on).mean()),
            "agree_secular_vs_spin_ftle": float((sec_on == spin_on).mean()),
            "ftle_median_abs_diff_secular_vs_spin": float(np.median(np.abs(v["ftle"] - v["spin_ftle"])))}


# ------------------------------------------------------------------ plots and summary
def plot_ftle(scans: dict) -> None:
    fig, axes = plt.subplots(2, len(scans), figsize=(7 * len(scans), 7), sharex=True)
    for col, (th, m) in enumerate(scans.items()):
        axes[0, col].semilogx(m["P_star"], m["ftle"], "k.", ms=3)
        axes[0, col].axhline(CFG["ftle_chaos_per_cycle"], color="C3", lw=0.8)
        axes[0, col].set_ylabel("FTLE per LK cycle")
        axes[0, col].set_title(rf"$\theta_{{lb}}^0={th}^\circ$ (SRFs on, tides off)")
        axes[1, col].loglog(m["P_star"], np.maximum(m["theta_std_last"], 1e-3), "k.", ms=3)
        axes[1, col].axhline(CFG["spread_chaos_deg"], color="C3", lw=0.8)
        axes[1, col].set_ylabel(r"std $\theta_{\rm sl}$ at e minima, last 500 cycles (deg)")
        axes[1, col].set_xlabel(r"$P_\star$ (d)")
    fig.tight_layout()
    fig.savefig(OUT / "prep_sq4_ftle.png", dpi=DPI)
    plt.close(fig)


def plot_omega(o: dict) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))
    for j, n in enumerate(CFG["omega_n_max0"]):
        axes[0].plot(o["omega_hat_l"], o["ftle"][:, j], "o-", ms=3, label=rf"$N_{{\max,0}}={n:g}$")
        axes[1].semilogy(o["omega_hat_l"], np.maximum(o["theta_std_last"][:, j], 1e-3), "o-", ms=3, label=rf"$N_{{\max,0}}={n:g}$")
    axes[0].axhline(CFG["ftle_chaos_per_cycle"], color="k", lw=0.6)
    axes[1].axhline(CFG["spread_chaos_deg"], color="k", lw=0.6)
    axes[0].set_xlabel(r"$\hat\Omega_L$")
    axes[1].set_xlabel(r"$\hat\Omega_L$")
    axes[0].set_ylabel("FTLE per LK cycle")
    axes[1].set_ylabel(r"std $\theta_{\rm sl}$, last 500 cycles (deg)")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(OUT / "prep_sq4_omega_test.png", dpi=DPI)
    plt.close(fig)


def plot_validation(v: dict) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(9, 8), sharex=True)
    axes[0].semilogx(v["eps"], v["ftle"], "k.", ms=4, label="secular.py (SRFs off)")
    axes[0].semilogx(v["eps"], v["spin_ftle"], "r.", ms=3, label="lkspin.spin")
    axes[0].axhline(CFG["ftle_chaos_per_cycle"], color="C0", lw=0.8)
    axes[0].set_ylabel("FTLE per LK cycle")
    axes[0].legend()
    axes[1].loglog(v["step3_eps"], np.maximum(v["step3_theta_std_last"], 1e-3), "b-", lw=1, label="step 3 (S15 Fig. 12), last 250 of 500 cycles")
    axes[1].loglog(v["eps"], np.maximum(v["theta_std_last"], 1e-3), "k.", ms=4, label="secular.py, last 500 of 1000")
    axes[1].axhline(CFG["spread_chaos_deg"], color="C0", lw=0.8)
    axes[1].set_xlabel(r"$\epsilon$")
    axes[1].set_ylabel(r"std $\theta_{\rm sl}$ at minima (deg)")
    axes[1].legend()
    fig.tight_layout()
    fig.savefig(OUT / "prep_sq4_validation_srf_off.png", dpi=DPI)
    plt.close(fig)


def shear_baseline() -> dict:
    """Pure shear (linear growth) FTLE over the same window, to justify the chaos threshold."""
    k = np.arange(1, CFG["n_cycles"] + 1)
    sep = TWIN * (1.0 + k)
    return C.ftle_fit(sep, TWIN, CFG["sat_rad"], CFG["min_fit_points"])


def load_summary() -> dict:
    path = OUT / "summary.json"
    return json.loads(path.read_text()) if path.exists() else {}


def write_summary(**sections) -> None:
    summary = load_summary()
    summary.update(sections)
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1, default=float))


# ------------------------------------------------------------------ stages
def run_stage(stage: str, pool) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    CHUNKS.mkdir(parents=True, exist_ok=True)
    t0 = time.time()
    if stage == "ftle":
        bands = stage_ftle(pool)
        plot_ftle({f"{th:g}": dict(np.load(OUT / f"ftle_scan_{th:g}.npz")) for th in CFG["theta_lb0_deg"]})
        write_summary(ftle=bands, shear_baseline_fit=shear_baseline(), tolerance_check=stage_tolerance(),
                      thresholds={"ftle_chaos_per_cycle": CFG["ftle_chaos_per_cycle"], "spread_chaos_deg": CFG["spread_chaos_deg"]})
    if stage == "sections":
        theta = CFG["sections_theta_lb0_deg"]
        orbit_times(theta, C.nondissipative_params(CFG4), f"th{theta:g}", CFG["n_cycles"])
        for msg in pool.imap_unordered(section_chunk, CFG["sections_p_star_days"], chunksize=1):
            print(msg, flush=True)
        plot_sections()
        write_summary(sections=[section_fill(p) for p in CFG["sections_p_star_days"]])
    if stage == "omega":
        for msg in pool.imap_unordered(omega_chunk, list(enumerate(omega_thetas())), chunksize=1):
            print(msg, flush=True)
        o = merge_omega()
        plot_omega(o)
        write_summary(omega_test=omega_summary(o))
    if stage == "validate":
        orbit_times(CFG["validate_theta_lb0_deg"], S.make_params(CFG4, srf=False, tides=False), f"val{CFG['validate_theta_lb0_deg']:g}", CFG["n_cycles"])
        for msg in pool.imap_unordered(validate_chunk, validate_items(), chunksize=1):
            print(msg, flush=True)
        v = merge_validate()
        plot_validation(v)
        write_summary(validation_srf_off=validation_summary(v))
    print(f"stage {stage} wall {time.time() - t0:.0f}s", flush=True)
    write_summary(**{f"wall_s_{stage}_last_invocation": time.time() - t0})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all", choices=["ftle", "sections", "omega", "validate", "all"])
    ap.add_argument("--procs", type=int, default=CFG["n_procs"])
    args = ap.parse_args()
    stages = ["validate", "ftle", "sections", "omega"] if args.stage == "all" else [args.stage]
    with mp.Pool(min(args.procs, CFG["n_procs"])) as pool:
        for st in stages:
            run_stage(st, pool)


if __name__ == "__main__":
    main()
