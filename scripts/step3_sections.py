"""Reproduces S15 Fig. 9 top row: surfaces of section (cos theta_sl vs phi at eccentricity minima), eps = 0.1, 0.05, 0.02.

40 initial conditions (20 values of p0 in [-0.98, 0.98] x phi0 in {0, pi}), N_CYCLES LK cycles each.
Run: uv run python scripts/step3_sections.py   (data -> results/step3/sections.npz, plot -> sections.png)
"""

import os
import sys
import time
from multiprocessing import Pool

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from lkspin.lk_orbit import make_orbit
from lkspin.spin import ATOL_DEFAULT, COL_P, COL_PHI, RTOL_DEFAULT, initial_state, run_cycles

OUT_DIR = os.path.join(ROOT, "results", "step3")
EPSILONS = (0.1, 0.05, 0.02)
N_CYCLES = 1000
N_P0 = 20
P0_MAX = 0.98
PHI0_VALUES = (0.0, np.pi)
DOT_SIZE = 0.8
PANEL_SIZE = (5.0, 5.0)
DPI = 200
PHI_TICKS = (0.0, np.pi / 2, np.pi, 3 * np.pi / 2, 2 * np.pi)
PHI_LABELS = ("0", r"$\pi/2$", r"$\pi$", r"$3\pi/2$", r"$2\pi$")
Y_PAD = 0.03


def initial_conditions():
    return [(p0, phi0) for phi0 in PHI0_VALUES for p0 in np.linspace(-P0_MAX, P0_MAX, N_P0)]


def section_points(job):
    eps, p0, phi0 = job
    orbit = make_orbit()
    orb = orbit.params()
    res = run_cycles(initial_state(p0, phi0, orb), N_CYCLES, orbit.period, orbit.alpha0(eps), orb,
                     RTOL_DEFAULT, ATOL_DEFAULT)
    return res[:, (COL_P, COL_PHI)].astype(np.float32)


def plot_panel(ax, eps, pts):
    ax.plot(pts[..., 1].ravel(), pts[..., 0].ravel(), "k.", ms=DOT_SIZE, markeredgewidth=0, rasterized=True)
    ax.set_xlim(0.0, 2 * np.pi)
    ax.set_ylim(-1.0 - Y_PAD, 1.0 + Y_PAD)
    ax.set_xticks(PHI_TICKS)
    ax.set_xticklabels(PHI_LABELS)
    ax.set_xlabel(r"$\phi$")
    ax.set_ylabel(r"$\cos\theta_{\rm sl}$")
    ax.set_title(rf"$\epsilon={eps}$")


def main():
    t0 = time.time()
    os.makedirs(OUT_DIR, exist_ok=True)
    ics = initial_conditions()
    jobs = [(eps, p0, phi0) for eps in EPSILONS for p0, phi0 in ics]
    with Pool(os.cpu_count()) as pool:
        flat = np.array(pool.map(section_points, jobs, chunksize=1))
    data = flat.reshape(len(EPSILONS), len(ics), N_CYCLES, 2)
    np.savez_compressed(os.path.join(OUT_DIR, "sections.npz"), eps=np.array(EPSILONS), ics=np.array(ics), points=data)
    fig, axes = plt.subplots(1, len(EPSILONS), figsize=(PANEL_SIZE[0] * len(EPSILONS), PANEL_SIZE[1]))
    for ax, eps, pts in zip(axes, EPSILONS, data, strict=True):
        plot_panel(ax, eps, pts)
    fig.tight_layout()
    fig.savefig(os.path.join(OUT_DIR, "sections.png"), dpi=DPI)
    print(f"sections: {len(jobs)} runs x {N_CYCLES} cycles, NaN rows: {int(np.isnan(data).any(axis=3).sum())}, "
          f"runtime {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
