"""Plot theta_sl,f vs P_star in the style of S16 Fig. 8 from whatever results/step4/scan_*.npz exist.

One panel per theta_lb0 (top = 89, bottom = 87). Pilot files are drawn as orange diamonds over the full scans.
If results/step4/validation.json exists, the top axis shows N_max,0 = coefficient / P_star with our measured coefficient
and a dashed line marks A_bar,0 = 1.

Run:  uv run python scripts/step4_plot.py
"""

import json
import re

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lkspin import secular as S

CFG = S.load_config()
OUT_DIR = S.ROOT / CFG["results"]
Y_MAX_DEG = 180.0
NMAX_TICKS = (0.5, 1, 2, 3, 4, 5, 10)
X_PAD = 1.05
PANEL_HEIGHT = 4.0
PANEL_WIDTH = 6.0
FILE_RE = re.compile(r"scan_([0-9.]+?)(_pilot)?\.npz$")


def find_scans() -> dict:
    """{theta: {"scan": path or None, "pilot": path or None}} for files in OUT_DIR."""
    found: dict = {}
    for path in sorted(OUT_DIR.glob("scan_*.npz")):
        m = FILE_RE.search(path.name)
        if not m:
            continue
        entry = found.setdefault(float(m.group(1)), {"scan": None, "pilot": None})
        entry["pilot" if m.group(2) else "scan"] = path
    return dict(sorted(found.items(), reverse=True))


def shape_coefficients() -> tuple[float, float] | None:
    path = OUT_DIR / "validation.json"
    if not path.exists():
        return None
    v5 = json.loads(path.read_text())["v5_shape"]
    return v5["Abar0_times_Pstar"], v5["Nmax0_times_Pstar"]


def draw_panel(ax, theta: float, entry: dict, coefs) -> None:
    if entry["scan"] is not None:
        d = np.load(entry["scan"])
        ax.plot(d["P_star"], d["theta_sl_f"], "k.", ms=3, label="scan")
    if entry["pilot"] is not None:
        d = np.load(entry["pilot"])
        ax.plot(d["P_star"], d["theta_sl_f"], "D", color="darkorange", ms=4, mfc="none", label="pilot (12 values)")
    ax.set_xscale("log")
    ax.set_xlim(CFG["p_star_max_days"] * X_PAD, CFG["p_star_min_days"] / X_PAD)
    ax.set_ylim(0, Y_MAX_DEG)
    ax.set_ylabel(r"$\theta_{\rm sl,f}$ (deg)")
    ax.set_xlabel(r"$P_\star$ (days)")
    ax.set_xticks([30, 10, 5, 2, 1, 0.5, 0.2, 0.1])
    ax.set_xticklabels(["30", "10", "5", "2", "1", "0.5", "0.2", "0.1"])
    ax.text(0.03, 0.93, rf"$\theta_{{lb,0}}={theta:g}^\circ$", transform=ax.transAxes, va="top")
    if coefs is not None and theta == 89.0:
        add_top_axis(ax, coefs)
    ax.legend(loc="upper right", fontsize=8)


def add_top_axis(ax, coefs) -> None:
    a_coef, n_coef = coefs
    ax.axvline(a_coef, color="r", ls="--", lw=0.8)
    top = ax.secondary_xaxis("top", functions=(lambda p: n_coef / p, lambda n: n_coef / n))
    top.set_xticks(list(NMAX_TICKS))
    top.set_xticklabels([str(n) for n in NMAX_TICKS])
    top.set_xlabel(r"$N_{\max,0}$ (our coefficient)")


def main() -> None:
    scans = find_scans()
    if not scans:
        raise SystemExit(f"no scan_*.npz in {OUT_DIR}")
    fig, axes = plt.subplots(len(scans), 1, figsize=(PANEL_WIDTH, PANEL_HEIGHT * len(scans)), squeeze=False)
    coefs = shape_coefficients()
    for ax, (theta, entry) in zip(axes[:, 0], scans.items()):
        draw_panel(ax, theta, entry, coefs)
    fig.tight_layout()
    out = OUT_DIR / "fig8_repro.png"
    fig.savefig(out, dpi=140)
    print("saved", out)


if __name__ == "__main__":
    main()
