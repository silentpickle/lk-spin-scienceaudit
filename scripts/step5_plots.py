"""Step 5 outputs from results/step5/features.npz and series.npz: regime map, confusion matrix, hand-check montages.

Labels are recomputed from config [step5] on every run, so changing a threshold and re-running this script is enough.
Writes regime_map.png, confusion_first_pass.txt, handcheck_01..05.png, handcheck_cases.csv.
Run:  uv run python scripts/step5_plots.py
"""

import csv
import math
import tomllib
from collections import Counter

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lkspin import regime as R
from lkspin import secular as S

with S.CONFIG_PATH.open("rb") as _fh:
    CFG5 = tomllib.load(_fh)["step5"]
OUT_DIR = S.ROOT / CFG5["results"]
COLORS = {"I": "#3d3d9e", "II": "#b0357a", "III": "#2f9e44", "IV.1": "#b59a1f", "IV.2": "#b59a1f", "IV.3": "#b59a1f",
          "IV.4+": "#b59a1f", "V": "#e8710a", "unclassified": "#888888"}
LEGEND_ROWS = (("I", "I. bimodal"), ("II", "II. adiabatic stationary"), ("III", "III. chaotic"),
               ("IV.1", "IV. adiabatic advected"), ("V", "V. fully adiabatic"), ("unclassified", "unclassified"))
MAP_XLIM = (0.05, 100.0)
MAP_YLIM = (0.02, 40.0)
MAP_SPLITS_N = (0.5, 1.0, 2.0)
MAP_FIGSIZE = (7.5, 5.2)
MAP_DPI = 150
DOT_AREA = 9.0
MONTAGE_ROWS = 10
MONTAGE_FIGSIZE = (20.0, 24.0)
MONTAGE_DPI = 100
MONTAGE_WIDTHS = (3.2, 1.0, 2.4)
THETA_AXIS_MAX = 180.0
JITTER = 0.15
N_MONTAGES = 5
DEG = 180.0 / math.pi
MYR = 1e6
CSV_FIELDS = ("montage", "row", "index", "pick_seed", "m_p_mjup", "p_star_days", "a_b_au", "theta_lb0_deg", "a_bar0",
              "n_max0", "s_over_l", "t_hj_yr", "theta_f_min", "theta_f_max", "range", "gap", "n_low", "n_high",
              "range_low", "range_high", "twin_div_deg", "twin_div_m0_deg", "m0_n_cycles", "m0_locked_frac_last",
              "m0_locked_frac_first", "m0_n_dom", "m0_w_median_last", "label", "predicted")


def load_table() -> dict:
    table = dict(np.load(OUT_DIR / "features.npz"))
    table["label"] = np.array(R.classify_table(table, CFG5))
    table["predicted"] = np.array(R.predicted_table(table, CFG5))
    return table


# ------------------------------------------------------------------ regime map
def draw_map_guides(ax) -> None:
    ax.axvline(CFG5["map_a_bar_split"], color="k", ls="--", lw=0.9)
    for n in MAP_SPLITS_N:
        ax.axhline(n, xmin=_axfrac(CFG5["map_a_bar_split"]), color="k", ls="--", lw=0.9)
    for text, x, y in (("I", 0.12, 0.3), ("II", 5.0, 0.07), ("IV.1", 30.0, 0.65), ("IV.2", 30.0, 1.35),
                       ("IV.3", 30.0, 2.8)):
        ax.text(x, y, text, fontsize=13, fontweight="bold", ha="center")


def _axfrac(x: float) -> float:
    lo, hi = MAP_XLIM
    return (math.log(x) - math.log(lo)) / (math.log(hi) - math.log(lo))


def plot_map(table: dict) -> None:
    fig, ax = plt.subplots(figsize=MAP_FIGSIZE)
    for label in R.CLASS_ORDER:
        sel = table["label"] == label
        if sel.any():
            ax.scatter(table["a_bar0"][sel], table["n_max0"][sel], s=DOT_AREA, color=COLORS[label], lw=0, alpha=0.85)
    for key, text in LEGEND_ROWS:
        n = int((np.array([R.label_family(x) for x in table["label"]]) == R.label_family(key)).sum())
        ax.scatter([], [], s=30, color=COLORS[key], label=f"{text} ({n})")
    draw_map_guides(ax)
    ax.set(xscale="log", yscale="log", xlim=MAP_XLIM, ylim=MAP_YLIM, xlabel=r"$\bar{\mathcal{A}}_0$",
           ylabel=r"$N_{\max,0}$")
    ax.set_title(f"First-pass outcome class, {len(table['index'])} hot-Jupiter systems (cf. S16 Fig. 14)", fontsize=10)
    ax.legend(loc="lower right", fontsize=8, framealpha=0.9)
    fig.tight_layout()
    fig.savefig(OUT_DIR / "regime_map.png", dpi=MAP_DPI)
    plt.close(fig)


# ------------------------------------------------------------------ confusion matrix
def write_confusion(table: dict) -> None:
    pred_rows = (R.PRED_I, R.PRED_II, R.PRED_IV)
    cols = R.CLASS_ORDER
    pairs = Counter(zip(table["predicted"], table["label"], strict=True))
    width = 14
    lines = ["Predicted (S16 map: I if A<1; II if A>=1, N<0.5; IV if A>=1, N>=0.5)  vs  first-pass class (rows: predicted)",
             f"{len(table['index'])} hot-Jupiter systems; thresholds from config.toml [step5]", "",
             "predicted".ljust(width) + "".join(c.rjust(width) for c in cols) + "total".rjust(width)]
    for p in pred_rows:
        counts = [pairs[(p, c)] for c in cols]
        lines.append(p.ljust(width) + "".join(str(c).rjust(width) for c in counts) + str(sum(counts)).rjust(width))
    totals = [sum(pairs[(p, c)] for p in pred_rows) for c in cols]
    lines.append("total".ljust(width) + "".join(str(t).rjust(width) for t in totals) + str(sum(totals)).rjust(width))
    lines += ["", "Same, with IV.1/IV.2/IV.3/IV.4+ merged into IV:", ""]
    fam = Counter((p, R.label_family(c)) for (p, c), v in pairs.items() for _ in range(v))
    fam_cols = ("I", "II", "III", "IV", "V", "unclassified")
    lines.append("predicted".ljust(width) + "".join(c.rjust(width) for c in fam_cols))
    for p in pred_rows:
        lines.append(p.ljust(width) + "".join(str(fam[(p, c)]).rjust(width) for c in fam_cols))
    agree = sum(fam[(p, c)] for p, c in (("I", "I"), ("II", "II"), ("IV", "IV")))
    lines += ["", f"agreement (predicted family == first-pass family): {agree} / {len(table['index'])}"]
    (OUT_DIR / "confusion_first_pass.txt").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


# ------------------------------------------------------------------ hand-check montages
def pick_handcheck(table: dict) -> np.ndarray:
    """Row numbers (into the HJ table) of handcheck_n systems, drawn without replacement with seed handcheck_seed."""
    rng = np.random.default_rng(CFG5["handcheck_seed"])
    return rng.choice(len(table["index"]), size=CFG5["handcheck_n"], replace=False)


def series_theta_w(series: np.ndarray, n_cycles: int) -> dict:
    """theta_sl (deg) at minima for member 0, middle member, twin; W_k for member 0; time of the end of LK activity."""
    t = series[:, 0] / MYR
    n_active_end = t[min(n_cycles, t.size - 1)]
    psi0 = series[:, 5]
    return {"t": t, "theta0": np.degrees(np.arccos(np.clip(series[:, 4], -1, 1))),
            "theta_mid": np.degrees(np.arccos(np.clip(series[:, 6], -1, 1))),
            "theta_twin": np.degrees(np.arccos(np.clip(series[:, 8], -1, 1))), "w": np.diff(psi0) / R.TWO_PI,
            "t_active_end": n_active_end}


def panel_theta(ax, d: dict) -> None:
    ax.plot(d["t"], d["theta_mid"], color="0.6", lw=3.5, label="middle member")
    ax.plot(d["t"], d["theta_twin"], color="#d62728", lw=1.0, label="twin")
    ax.plot(d["t"], d["theta0"], color="#1f4fa8", lw=1.0, marker=".", ms=2.5, label="member 0")
    ax.axvline(d["t_active_end"], color="k", ls=":", lw=0.8)
    ax.set(ylim=(0, THETA_AXIS_MAX), xlabel="t (Myr)", ylabel=r"$\theta_{sl}$ at e minima (deg)")
    ax.legend(fontsize=7, loc="upper right", ncol=3)


def panel_strip(ax, theta_f: np.ndarray, twin: float, seed: int) -> None:
    rng = np.random.default_rng(seed)
    ax.scatter(rng.uniform(-JITTER, JITTER, theta_f.size), theta_f, s=22, color="#1f4fa8")
    ax.scatter([0.0], [twin], s=60, marker="x", color="#d62728")
    ax.set(xlim=(-0.5, 0.5), ylim=(0, THETA_AXIS_MAX), xticks=[], ylabel=r"final $\theta_{sl}$ (deg)")


def panel_winding(ax, d: dict, n_cycles: int) -> None:
    w = d["w"]
    ax.plot(np.arange(w.size), w, color="#1f4fa8", lw=0.8, marker=".", ms=2.5)
    lo, hi = math.floor(w.min()) - 1, math.ceil(w.max()) + 1
    for k in range(lo, hi + 1):
        ax.axhline(k, color="0.7", lw=0.5, zorder=0)
    ax.axvline(n_cycles, color="k", ls=":", lw=0.8)
    ax.set(ylim=(lo, hi), xlabel="cycle index k", ylabel=r"$W_k$ (member 0)")


def row_title(table: dict, r: int, pos: int) -> str:
    return (f"#{pos}  system {int(table['index'][r])}   A0={table['a_bar0'][r]:.3g}  Nmax0={table['n_max0'][r]:.3g}  "
            f"P*={table['p_star_days'][r]:.3g} d  Mp={table['m_p_mjup'][r]:.2g}  ab={table['a_b_au'][r]:.0f}  "
            f"th_lb0={table['theta_lb0_deg'][r]:.1f}   first-pass: {table['label'][r]}  (map: {table['predicted'][r]})")


def draw_montage_row(fig, gs_row, table: dict, series, r: int, pos: int) -> None:
    axes = [fig.add_subplot(gs_row[c]) for c in range(3)]
    sys_i = int(table["index"][r])
    n_cycles = int(table["m0_n_cycles"][r])
    d = series_theta_w(series[f"s{sys_i}"], n_cycles)
    panel_theta(axes[0], d)
    panel_strip(axes[1], table["theta_f"][r], float(table["theta_f_twin"][r]), sys_i)
    panel_winding(axes[2], d, n_cycles)
    fig.text(0.01, axes[0].get_position().y1 + 0.004, row_title(table, r, pos), fontsize=10, fontweight="bold")


def draw_montage(table: dict, series, rows: np.ndarray, number: int) -> None:
    fig = plt.figure(figsize=MONTAGE_FIGSIZE, dpi=MONTAGE_DPI)
    grid = fig.add_gridspec(MONTAGE_ROWS, 3, width_ratios=MONTAGE_WIDTHS, left=0.05, right=0.99, top=0.975,
                            bottom=0.03, hspace=0.75, wspace=0.22)
    for k, r in enumerate(rows):
        draw_montage_row(fig, [grid[k, c] for c in range(3)], table, series, int(r), (number - 1) * MONTAGE_ROWS + k + 1)
    fig.savefig(OUT_DIR / f"handcheck_{number:02d}.png", dpi=MONTAGE_DPI)
    plt.close(fig)


def csv_row(table: dict, r: int, montage: int, row: int) -> dict:
    out = {"montage": montage, "row": row, "pick_seed": CFG5["handcheck_seed"],
           "theta_f_min": float(table["theta_f"][r].min()), "theta_f_max": float(table["theta_f"][r].max())}
    for k in CSV_FIELDS:
        if k not in out:
            out[k] = table[k][r].item() if hasattr(table[k][r], "item") else table[k][r]
    return out


def write_handcheck(table: dict, picks: np.ndarray) -> None:
    series = np.load(OUT_DIR / "series.npz")
    rows = []
    for m in range(N_MONTAGES):
        chunk = picks[m * MONTAGE_ROWS:(m + 1) * MONTAGE_ROWS]
        draw_montage(table, series, chunk, m + 1)
        rows += [csv_row(table, int(r), m + 1, k + 1) for k, r in enumerate(chunk)]
    with (OUT_DIR / "handcheck_cases.csv").open("w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=CSV_FIELDS)
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    table = load_table()
    plot_map(table)
    write_confusion(table)
    write_handcheck(table, pick_handcheck(table))
    print("class counts:", dict(Counter(table["label"])))


if __name__ == "__main__":
    main()
