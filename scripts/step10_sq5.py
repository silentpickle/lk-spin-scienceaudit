"""Step 10 (SQ5): how large must S/L be before spin-orbit feedback changes the outcome?  Pairs case (a) (results/step5) with
case (b) (results/step6/b) per system (hot Jupiters in both), and tabulates, in log-spaced bins of S/L at the first
eccentricity maximum (results/step5 features: s_over_l):
  * the class-change fraction (v2 families, (a) vs (b)),
  * the median over systems of |Delta theta_sl,f|, itself the median over the 10 ensemble members of |theta_b - theta_a|
    (members are matched: same P_star offset),
  * the median within-ensemble spread of case (a) (standard deviation of the 10 theta_f of (a)),
each with a bootstrap interval (systems resampled within the bin).  It then finds the S/L where the class-change fraction first
exceeds 10% and 50%, and where the median |Delta theta| first exceeds the median within-ensemble spread (linear interpolation in
log S/L between bin centres; bootstrap interval from resampling all paired systems and recomputing the bins).
Writes results/step10/sq5.json and sq5.png (sq5_pilot.* with --pilot).
Run:  uv run python scripts/step10_sq5.py [--pilot]
"""

import argparse
import json
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import step5_baselines as B
import step6_analyze as A
import step6_run as W

CFG6 = W.CFG6
OUT = W.S.ROOT / CFG6["results_sq5"]
MIN_BIN_COUNT = 5  # ASSUMED: bins with fewer paired systems are not used to locate a crossing
FIG_SIZE = (7.5, 9.0)
DPI = 150
N_PANELS = 3


def paired_arrays(case: str, pilot: bool) -> dict:
    a_data, c_data = A.load_case("a", False), A.load_case(case, pilot)
    pos_a = {int(i): k for k, i in enumerate(a_data["feats"]["index"])}
    rows_c = np.array([k for k, i in enumerate(c_data["feats"]["index"]) if int(i) in pos_a])
    rows_a = np.array([pos_a[int(c_data["feats"]["index"][k])] for k in rows_c])
    theta_a, theta_c = a_data["feats"]["theta_f"][rows_a], c_data["feats"]["theta_f"][rows_c]
    return {"s_over_l": a_data["feats"]["s_over_l"][rows_a], "changed": a_data["families"][rows_a] != c_data["families"][rows_c],
            "dtheta": np.median(np.abs(theta_c - theta_a), axis=1), "spread_a": np.std(theta_a, axis=1),
            "n_a_hj": int(a_data["feats"]["index"].size), "n_case_hj": int(c_data["feats"]["index"].size)}


def log_edges(s_over_l: np.ndarray) -> np.ndarray:
    return np.geomspace(s_over_l.min(), s_over_l.max() * (1.0 + 1e-9), CFG6["sq5_n_bins"] + 1)


def bin_stats(data: dict, edges: np.ndarray, rows: np.ndarray | None = None) -> dict:
    """Per bin: count, class-change fraction, median |dtheta|, median spread of (a). rows: optional (resampled) system rows."""
    rows = np.arange(data["s_over_l"].size) if rows is None else rows
    which = np.clip(np.searchsorted(edges, data["s_over_l"][rows], side="right") - 1, 0, edges.size - 2)
    out = {k: np.full(edges.size - 1, np.nan) for k in ("frac", "dtheta", "spread")}
    out["n"] = np.zeros(edges.size - 1, dtype=int)
    for b in range(edges.size - 1):
        r = rows[which == b]
        out["n"][b] = r.size
        if r.size:
            out["frac"][b] = data["changed"][r].mean()
            out["dtheta"][b] = np.median(data["dtheta"][r])
            out["spread"][b] = np.median(data["spread_a"][r])
    return out


def first_crossing(x_log: np.ndarray, y: np.ndarray, level: float, n: np.ndarray) -> float:
    """S/L (not log) where y first rises above level, interpolated linearly in log S/L between bin centres; NaN if never."""
    ok = n >= MIN_BIN_COUNT
    xs, ys = x_log[ok], y[ok]
    above = np.where(ys > level)[0]
    if above.size == 0:
        return math.nan
    k = int(above[0])
    if k == 0:
        return float(10.0 ** xs[0])
    t = (level - ys[k - 1]) / (ys[k] - ys[k - 1])
    return float(10.0 ** (xs[k - 1] + t * (xs[k] - xs[k - 1])))


def crossings(stats: dict, centres_log: np.ndarray) -> dict:
    out = {f"frac_gt_{level:g}": first_crossing(centres_log, stats["frac"], level, stats["n"]) for level in CFG6["sq5_fraction_levels"]}
    out["dtheta_gt_spread"] = first_crossing(centres_log, stats["dtheta"] - stats["spread"], 0.0, stats["n"])
    return out


def boot_bins(data: dict, edges: np.ndarray, rng: np.random.Generator) -> dict:
    """Bootstrap of the binned statistics (systems resampled WITHIN each bin) as percentile arrays, and of the crossings (all)."""
    n_sys = data["s_over_l"].size
    which = np.clip(np.searchsorted(edges, data["s_over_l"], side="right") - 1, 0, edges.size - 2)
    members = [np.where(which == b)[0] for b in range(edges.size - 1)]
    draws = {k: np.full((B.N_BOOT, edges.size - 1), np.nan) for k in ("frac", "dtheta")}
    for it in range(B.N_BOOT):
        for b, m in enumerate(members):
            if m.size:
                r = m[rng.integers(0, m.size, m.size)]
                draws["frac"][it, b] = data["changed"][r].mean()
                draws["dtheta"][it, b] = np.median(data["dtheta"][r])
    return {k: np.nanpercentile(v, B.CI, axis=0) for k, v in draws.items()} | {"n_sys": n_sys}


def boot_crossings(data: dict, edges: np.ndarray, centres_log: np.ndarray, rng: np.random.Generator) -> dict:
    n_sys = data["s_over_l"].size
    samples = [crossings(bin_stats(data, edges, rng.integers(0, n_sys, n_sys)), centres_log) for _ in range(B.N_BOOT)]
    out = {}
    for key in samples[0]:
        vals = np.array([s[key] for s in samples])
        good = vals[np.isfinite(vals)]
        out[key] = {"finite_fraction": float(good.size / vals.size),
                    "ci": [float(v) for v in np.percentile(good, B.CI)] if good.size else [math.nan, math.nan]}
    return out


def plot(edges: np.ndarray, centres: np.ndarray, stats: dict, boots: dict, cross: dict, path) -> None:
    fig, axes = plt.subplots(N_PANELS, 1, figsize=FIG_SIZE, sharex=True)
    ok = stats["n"] > 0  # every populated bin is drawn; only bins with MIN_BIN_COUNT or more count for the crossings
    axes[0].errorbar(centres[ok], stats["frac"][ok], yerr=[stats["frac"][ok] - boots["frac"][0][ok], boots["frac"][1][ok] - stats["frac"][ok]],
                     fmt="o-")
    for level in CFG6["sq5_fraction_levels"]:
        axes[0].axhline(level, color="k", ls="--", lw=1)
        x = cross[f"frac_gt_{level:g}"]
        if np.isfinite(x):
            axes[0].axvline(x, color="gray", ls=":", lw=1)
    axes[0].set(ylabel="class-change fraction (a) -> (b)", ylim=(0, 1))
    axes[1].errorbar(centres[ok], stats["dtheta"][ok], yerr=[stats["dtheta"][ok] - boots["dtheta"][0][ok], boots["dtheta"][1][ok] - stats["dtheta"][ok]],
                     fmt="o-", label=r"median $|\Delta\theta_{\rm sl,f}|$ (a vs b)")
    axes[1].plot(centres[ok], stats["spread"][ok], "s--", label="median ensemble spread of (a) (std of 10 members)")
    axes[1].set(ylabel="deg", yscale="log")
    axes[1].legend(fontsize=7)
    axes[2].bar(centres, stats["n"], width=0.8 * np.diff(edges), align="center")
    axes[2].set(ylabel="systems per bin", xlabel=r"$S_\star/L$ at the first $e$ maximum (step 5 $P_\star$)", xscale="log")
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    data = paired_arrays("b", args.pilot)
    edges = log_edges(data["s_over_l"])
    centres = np.sqrt(edges[:-1] * edges[1:])
    centres_log = np.log10(centres)
    rng = np.random.default_rng(CFG6["boot_seed"])
    stats = bin_stats(data, edges)
    cross = crossings(stats, centres_log)
    boots = boot_bins(data, edges, rng)
    suffix = "_pilot" if args.pilot else ""
    res = {"n_paired": int(data["s_over_l"].size), "n_a_hot_jupiters": data["n_a_hj"], "n_b_hot_jupiters": data["n_case_hj"],
           "n_boot": B.N_BOOT, "min_bin_count": MIN_BIN_COUNT, "edges": edges.tolist(), "centres": centres.tolist(),
           "n": stats["n"].tolist(), "class_change_fraction": stats["frac"].tolist(), "class_change_ci": boots["frac"].tolist(),
           "median_abs_dtheta_deg": stats["dtheta"].tolist(), "median_abs_dtheta_ci": boots["dtheta"].tolist(),
           "median_spread_a_deg": stats["spread"].tolist(), "global_median_spread_a_deg": float(np.median(data["spread_a"])),
           "global_median_abs_dtheta_deg": float(np.median(data["dtheta"])), "crossings_s_over_l": cross,
           "crossings_bootstrap": boot_crossings(data, edges, centres_log, rng),
           "overall_class_change_fraction": float(data["changed"].mean())}
    (OUT / f"sq5{suffix}.json").write_text(json.dumps(A.nan_to_none(res), indent=2) + "\n")
    plot(edges, centres, stats, boots, cross, OUT / f"sq5{suffix}.png")
    print(json.dumps(A.nan_to_none({k: res[k] for k in ("n_paired", "n", "class_change_fraction", "median_abs_dtheta_deg", "median_spread_a_deg",
                                          "crossings_s_over_l", "crossings_bootstrap")}), indent=2))


if __name__ == "__main__":
    main()
