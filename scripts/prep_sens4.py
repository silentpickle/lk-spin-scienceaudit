"""Prep job 3: sensitivity of the step 4 reproduction to the ASSUMED parameters, one at a time (config [sens4]).

Per (variant, theta_lb0): first-cycle shape coefficients, 1 - e_max, then ONE integration carrying a stellar spin for every
P_star (grouped spins, the orbit does not depend on P_star). Resumable: results/sens4/<variant>_th<theta>.npz.
Run:  uv run python scripts/prep_sens4.py [--stage run|analyze|all]
"""

import argparse
import json
import multiprocessing as mp
import os
import time
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import prep_common as C
from lkspin import secular as S

CFG4 = S.load_config()
CFG = C.load_table("sens4")
OUT = S.ROOT / CFG["results"]
BASE = "baseline"
BASE_BULGE = S.BULGE_COEF_A16_EOM
DPI = 150
GYR = 1e9
MYR = 1e6
STEP4_DIR = S.ROOT / "results/step4"
N_DEG_MAX = 180.0


def variants() -> dict:
    out = {BASE: {}}
    out.update(CFG["variant"])
    return out


def p_grid() -> tuple[np.ndarray, int, int]:
    """(all P_star, n_main, n_extra): 200 points picked from the step 4 1000-point grid, then extras and the fine grid."""
    grid = np.geomspace(CFG["p_star_max_days"], CFG["p_star_min_days"], CFG["step4_grid_points"])
    idx = np.round(np.linspace(0, grid.size - 1, CFG["n_p_star"])).astype(int)
    lo, hi = CFG["fine_p_star_range_days"]
    extra = np.concatenate([CFG["extra_p_star_days"], np.linspace(lo, hi, CFG["fine_p_star_n"])])
    return np.concatenate([grid[idx], extra]), idx.size, extra.size


def job_path(name: str, theta: float) -> Path:
    return OUT / f"{name}_th{theta:g}.npz"


def variant_setup(name: str) -> tuple[dict, np.ndarray]:
    over = dict(variants()[name])
    bulge = over.pop("bulge_coef", BASE_BULGE)
    cfg = dict(CFG4, **over)
    return cfg, S.make_params(cfg, bulge_coef=bulge)


def run_variant(args: tuple) -> str:
    name, theta = args
    path = job_path(name, theta)
    if path.exists():
        return f"skip {path.name}"
    t0 = time.time()
    cfg, par = variant_setup(name)
    shape = C.first_cycle_shape(par, cfg, theta)
    periods, n_main, _ = p_grid()
    y0 = S.initial_state(par, cfg, theta, periods.size)
    r = S.run(y0, par, list(periods), cfg["stop_t_gyr"] * GYR)
    extended = 0
    if r["status"] == S.ST_TIMEOUT:
        r, extended = S.run(y0, par, list(periods), CFG["stop_t_gyr_extended"] * GYR), 1
    _, j, a_end = S.orbit_elements(r["y"], par)
    tmp = path.with_suffix(f".{os.getpid()}.tmp.npz")
    np.savez(tmp, P_star=periods, theta_sl_f=S.theta_sl_final(r["y"], periods.size), status=r["status"], t_end_yr=r["t"],
             a_end=a_end, extended=extended, n_main=n_main, wall_s=time.time() - t0, **shape)
    tmp.replace(path)
    return f"done {path.name} {S.STATUS_NAMES[r['status']]} t_end {r['t'] / MYR:.0f} Myr {time.time() - t0:.0f}s"


# ------------------------------------------------------------------ analysis
def load(name: str, theta: float) -> dict:
    return dict(np.load(job_path(name, theta)))


def chaos_mask(theta: float, p: np.ndarray) -> np.ndarray:
    lo, hi = CFG["chaos_band_87_days"]
    return (p >= lo) & (p <= hi) if theta == 87.0 else np.zeros(p.size, bool)


def point_value(d: dict, p_target: float) -> float:
    return float(d["theta_sl_f"][int(np.argmin(np.abs(d["P_star"] - p_target)))])


def compare_to_baseline(d: dict, b: dict, theta: float) -> dict:
    n = int(d["n_main"])
    p = d["P_star"][:n]
    diff = np.abs(d["theta_sl_f"][:n] - b["theta_sl_f"][:n])
    keep = ~chaos_mask(theta, p)
    bimodal = p >= CFG["bimodal_min_days"][0 if theta == 89.0 else 1]
    reg = keep & ~bimodal
    return {"max_dtheta_excl_chaos": float(diff[keep].max()), "max_dtheta_regular_only": float(diff[reg].max()),
            "median_dtheta_regular_only": float(np.median(diff[reg])),
            "bimodal_branch_agreement": float(((d["theta_sl_f"][:n][bimodal] > CFG["branch_split_deg"])
                                               == (b["theta_sl_f"][:n][bimodal] > CFG["branch_split_deg"])).mean()),
            "P_at_max_regular_days": float(p[reg][np.argmax(diff[reg])])}


def fine_branch_fraction(d: dict) -> float:
    n = int(d["n_main"]) + len(CFG["extra_p_star_days"])
    return float((d["theta_sl_f"][n:] > CFG["branch_split_deg"]).mean())


def row(name: str, theta: float, b: dict) -> dict:
    d = load(name, theta)
    out = {"variant": name, "theta_lb0": theta, "status": S.STATUS_NAMES[int(d["status"])], "extended_run": int(d["extended"]),
           "t_hj_Myr": float(d["t_end_yr"]) / MYR, "one_minus_e_max": float(d["one_minus_e_max"]),
           "nmax_coef": float(d["nmax_coef"]), "abar_coef": float(d["abar_coef"]),
           "abar_over_nmax": float(d["abar_coef"] / d["nmax_coef"]), "omega_hat_l": float(d["omega_hat_l"]),
           "period_Myr": float(d["period_yr"]) / MYR,
           "theta_30d": point_value(d, 30.0), "theta_29_87d": float(d["theta_sl_f"][int(d["n_main"])]),
           "high_branch_fraction_24_30d": fine_branch_fraction(d)}
    if theta == 89.0:
        out["nmax_ratio_to_S16"] = out["nmax_coef"] / CFG["s16_nmax_coef"]
        out["abar_ratio_to_S16"] = out["abar_coef"] / CFG["s16_abar_coef"]
    out.update(compare_to_baseline(d, b, theta))
    return out


def sensitivities(rows: list) -> dict:
    """Log-derivative of the two 89 deg coefficients and 1 - e_max with respect to each two-sided parameter pair."""
    by = {(r["variant"], r["theta_lb0"]): r for r in rows}
    pairs = {"chi": ("chi_5", "chi_20", 5.0, 20.0), "r_p": ("rp_0p9", "rp_1p2", 0.9, 1.2), "e0": ("e0_0p001", "e0_0p05", 0.001, 0.05),
             "p_planet0": ("pp_5h", "pp_20h", 5.0, 20.0), "r_star": ("rs_0p9", "rs_1p1", 0.9, 1.1)}
    out = {}
    for key, (lo, hi, x_lo, x_hi) in pairs.items():
        a, c = by[(lo, 89.0)], by[(hi, 89.0)]
        dlog = np.log(x_hi / x_lo)
        out[key] = {f: float(np.log(c[f] / a[f]) / dlog) for f in ("nmax_coef", "abar_coef", "one_minus_e_max", "t_hj_Myr")}
    return out


def plot_variant(name: str, base: dict[float, dict]) -> None:
    fig, axes = plt.subplots(1, 2, figsize=(14, 4.5))
    for ax, theta in zip(axes, CFG["theta_lb0_deg"], strict=True):
        d, b = load(name, theta), base[theta]
        n = int(d["n_main"])
        ax.semilogx(b["P_star"][:n], b["theta_sl_f"][:n], "k.", ms=4, label="baseline")
        ax.semilogx(d["P_star"][:n], d["theta_sl_f"][:n], "r.", ms=2.5, label=name)
        ax.set_ylim(-5, N_DEG_MAX)
        ax.set_xlabel(r"$P_\star$ (d)")
        ax.set_ylabel(r"$\theta_{\rm sl,f}$ (deg)")
        ax.set_title(rf"$\theta_{{lb}}^0={theta:g}^\circ$  ({S.STATUS_NAMES[int(d['status'])]}, t={float(d['t_end_yr']) / MYR:.0f} Myr)")
        ax.legend()
    fig.tight_layout()
    fig.savefig(OUT / f"sens4_fig8_{name}.png", dpi=DPI)
    plt.close(fig)


def step4_check(base: dict[float, dict]) -> dict:
    """Grouped baseline against the step 4 one-spin-per-run scan at the same P_star values."""
    out = {}
    grid = np.geomspace(CFG["p_star_max_days"], CFG["p_star_min_days"], CFG["step4_grid_points"])
    idx = np.round(np.linspace(0, grid.size - 1, CFG["n_p_star"])).astype(int)
    for theta, b in base.items():
        ref = np.load(STEP4_DIR / f"scan_{theta:g}.npz")
        n = int(b["n_main"])
        diff = np.abs(ref["theta_sl_f"][idx] - b["theta_sl_f"][:n])
        keep = ~chaos_mask(theta, b["P_star"][:n])
        out[f"{theta:g}"] = {"max_abs_diff_deg_all": float(diff.max()), "max_abs_diff_deg_excl_chaos": float(diff[keep].max()),
                             "n_points_diff_gt_1deg": int((diff > 1.0).sum()), "n_points": n, "P_grid_matches": bool(np.allclose(ref["P_star"][idx], b["P_star"][:n]))}
    return out


def analyze() -> None:
    base = {th: load(BASE, th) for th in CFG["theta_lb0_deg"]}
    rows = [row(v, th, base[th]) for v in variants() for th in CFG["theta_lb0_deg"]]
    for v in variants():
        plot_variant(v, base)
    summary = {"rows": rows, "log_derivatives_89deg": sensitivities(rows), "baseline_vs_step4_scan": step4_check(base),
               "s16": {"theta_30d": CFG["s16_theta_30d_deg"], "theta_29_87d": CFG["s16_theta_29_87d_deg"],
                       "nmax_coef": CFG["s16_nmax_coef"], "abar_coef": CFG["s16_abar_coef"]}}
    (OUT / "summary.json").write_text(json.dumps(summary, indent=1))
    print_table(rows)


def print_table(rows: list) -> None:
    print("variant theta status t_hj[Myr] 1-emax N/S16 A/S16 th30 th29.87 hiFrac maxdTh(excl chaos) maxdTh(regular) medDTh(reg)")
    for r in rows:
        print(f"{r['variant']:10s} {r['theta_lb0']:4.0f} {r['status'][:4]} {r['t_hj_Myr']:8.0f} {r['one_minus_e_max']:.5f} "
              f"{r.get('nmax_ratio_to_S16', float('nan')):.3f} {r.get('abar_ratio_to_S16', float('nan')):.3f} {r['theta_30d']:6.1f} "
              f"{r['theta_29_87d']:6.1f} {r['high_branch_fraction_24_30d']:.2f} {r['max_dtheta_excl_chaos']:6.1f} "
              f"{r['max_dtheta_regular_only']:6.2f} {r['median_dtheta_regular_only']:6.3f}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all", choices=["run", "analyze", "all"])
    ap.add_argument("--procs", type=int, default=CFG["n_procs"])
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage in ("run", "all"):
        items = [(v, th) for th in CFG["theta_lb0_deg"] for v in variants()]
        with mp.Pool(min(args.procs, CFG["n_procs"])) as pool:
            for msg in pool.imap_unordered(run_variant, items, chunksize=1):
                print(msg, flush=True)
    if args.stage in ("analyze", "all"):
        analyze()


if __name__ == "__main__":
    main()
