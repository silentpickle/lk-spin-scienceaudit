"""Reproduces S15 Fig. 12 (log eps 0.008-1), Fig. 13 (0.04-0.07) and Fig. 14 (0.025-0.0286, vs 1/eps).

Aligned start (p0 = 1), N_CYCLES LK cycles per epsilon, theta_sl (deg) recorded at every eccentricity minimum.
Run: uv run python scripts/step3_bifurcation.py [--replot]   (--replot reuses the saved .npz; data -> results/step3/bifurcation_fig*.npz, plots -> fig*.png)
"""

import os
import sys
import time
from multiprocessing import Pool

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from lkspin.lk_orbit import make_orbit
from lkspin.spin import ATOL_DEFAULT, RTOL_DEFAULT, initial_state, run_cycles

N_CYCLES = 500
OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "results", "step3")
FIGURES = {  # name: (eps_min, eps_max, n_eps, log_spaced)
    "fig12": (0.008, 1.0, 500, True),
    "fig13": (0.04, 0.07, 400, False),
    "fig14": (0.025, 0.0286, 400, False),
}
DOT_SIZE = 1.5
FIG_SIZE = (8.5, 4.5)
THETA_MAX_DEG = 180.0
YLIM_PAD = 5.0
FIG12_TICKS = (0.01, 0.02, 0.05, 0.1, 0.2, 0.5, 1.0)
DPI = 200


def theta_series(eps):
    orbit = make_orbit()
    orb = orbit.params()
    y0 = initial_state(1.0, 0.0, orb)
    res = run_cycles(y0, N_CYCLES, orbit.period, orbit.alpha0(eps), orb, RTOL_DEFAULT, ATOL_DEFAULT)
    return np.degrees(np.arccos(np.clip(res[:, 0], -1.0, 1.0))).astype(np.float32)


def eps_grid(spec):
    lo, hi, n, log = spec
    return np.geomspace(lo, hi, n) if log else np.linspace(lo, hi, n)


def run_figure(name, pool):
    eps = eps_grid(FIGURES[name])
    t = time.time()
    theta = np.array(pool.map(theta_series, eps, chunksize=1))
    np.savez_compressed(os.path.join(OUT_DIR, f"bifurcation_{name}.npz"), eps=eps, theta_deg=theta)
    print(f"{name}: {len(eps)} eps x {N_CYCLES} cycles in {time.time() - t:.0f} s")
    return eps, theta


def plot_figure(name, eps, theta):
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    x = 1.0 / eps if name == "fig14" else eps
    ax.plot(np.repeat(x, theta.shape[1]), theta.ravel(), "k.", ms=DOT_SIZE, markeredgewidth=0, rasterized=True)
    ax.set_ylim(-YLIM_PAD, THETA_MAX_DEG)
    ax.set_ylabel(r"$\theta_{\rm sl}$ (deg)")
    ax.set_xlabel(r"$1/\epsilon$" if name == "fig14" else r"$\epsilon$")
    if name == "fig12":
        ax.set_xscale("log")
        ax.set_xticks(FIG12_TICKS)
        ax.set_xticklabels([str(v) for v in FIG12_TICKS])
        ax.minorticks_off()
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, f"bifurcation_{name}.png"), dpi=DPI)
    plt.close(fig)


def replot_all():
    for name in FIGURES:
        data = np.load(os.path.join(OUT_DIR, f"bifurcation_{name}.npz"))
        plot_figure(name, data["eps"], data["theta_deg"])


def main():
    t0 = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)
    if "--replot" in sys.argv:
        replot_all()
        return
    with Pool(os.cpu_count()) as pool:
        for name in FIGURES:
            plot_figure(name, *run_figure(name, pool))
    print(f"total runtime {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
