"""Exact single-resonance locations/widths (S15 Fig. 5 filled circles) and the Fig. 11 boundary rebuilt from them.

Run:  uv run python scripts/step3_exact.py
Writes results/step3/fourier_widths_exact.png, chaos_boundary_exact.npz/.png and bifurcation_with_boundary.png (the
last reads results/step3/bifurcation_fig12.npz without modifying it), and prints the Fig. 5 comparison table.
"""

import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ROOT = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(ROOT))

from lkspin.lk_orbit import make_orbit  # noqa: E402
from lkspin.shape import (  # noqa: E402
    EXACT_FIELDS,
    chaos_boundary_exact,
    cos_coefficients,
    coupling,
    exact_resonances,
    script_n,
    shape_functions,
    sin_coefficients,
    widths,
)

OUT = ROOT / "results" / "step3"
M_COEF = 500
EPS_FIG5 = (0.1, 0.05, 0.02)
COLORS = ("k", "r", "g")
N_PLOT = 50
EPS_MIN, EPS_MAX, N_EPS = 0.003, 1.0, 1500
DPI = 130
# Values read off S15 Fig. 5 (filled black = exact, eps = 0.1) by pixel measurement; N = -9..9.
PAPER_W_EPS01 = {-9: 0.082, -8: 0.150, -7: 0.210, -6: 0.269, -5: 0.328, -4: 0.387, -3: 0.451, -2: 0.524, -1: 0.629,
                 0: 0.826, 1: 0.351, 2: 0.227, 3: 0.220, 4: 0.205, 5: 0.188, 6: 0.167, 7: 0.142, 8: 0.113, 9: 0.080}
BIF_DOT = 1.0
BIF_BOUND = "any_edge"


def resonance_set(c_all, sn, psi0, eps):
    n_top = min(int(np.floor(sn / eps)), M_COEF)
    n = np.arange(-n_top, n_top + 1)
    cc = c_all[n + M_COEF]
    return n, exact_resonances(cc, n, eps, sn, psi0), widths(cc, n, eps, sn)


def save_width_figure(c_all, sn, psi0) -> None:
    fig, (ax_p, ax_w) = plt.subplots(2, 1, figsize=(8, 8), sharex=True)
    for eps, col in zip(EPS_FIG5, COLORS):
        n, ex, wa = resonance_set(c_all, sn, psi0, eps)
        ax_p.plot(n, ex.centre, "o", color=col, ms=3.5, label=rf"exact $\epsilon={eps}$")
        ax_p.plot(n, -n * eps / sn, "o", mfc="none", color="b", ms=6, mew=0.5)
        ax_w.plot(n, ex.width, "o", color=col, ms=3.5)
        ax_w.plot(n, wa, "o", mfc="none", color="b", ms=6, mew=0.5)
    ax_p.set(ylabel=r"$p_N$", ylim=(-1.05, 1.05), xlim=(-N_PLOT, N_PLOT))
    ax_p.legend()
    ax_w.set(xlabel="N", ylabel=r"$w_N$ (filled exact, open blue approximate)")
    fig.savefig(OUT / "fourier_widths_exact.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def print_table(c_all, sn, psi0) -> None:
    n, ex, wa = resonance_set(c_all, sn, psi0, 0.1)
    print("eps = 0.1:   N  paper(exact)  w_exact  w_approx  | p_exact  p_nominal")
    for k, nn in enumerate(n):
        print(f"  {nn:3d}  {PAPER_W_EPS01[nn]:8.3f}  {ex.width[k]:8.4f}  {wa[k]:8.4f}  | {ex.centre[k]:+.4f}  {-nn * 0.1 / sn:+.4f}"
              f"   exact/paper {ex.width[k] / PAPER_W_EPS01[nn]:.3f}  approx/paper {wa[k] / PAPER_W_EPS01[nn]:.3f}")


def save_boundary_figure(b, eps) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(8, 9), sharex=True)
    for ax, where, title in ((axes[0], "centre", "exact centre of outermost resonance"),
                             (axes[1], "edge", "outer separatrix edge of outermost resonance")):
        ax.fill_between(eps, b[f"bottom_any_{where}"], b[f"top_any_{where}"], color="0.75", label="largest |N| anywhere")
        ax.plot(eps, b[f"top_contig_{where}"], "b-", lw=0.7, label="contiguous chain")
        ax.plot(eps, b[f"bottom_contig_{where}"], "b-", lw=0.7)
        ax.set(xscale="log", ylabel=r"$\cos\theta_{\rm sl}$", xlim=(EPS_MIN, EPS_MAX), ylim=(-1, 1), title=title)
        ax.legend(loc="center right")
    axes[1].set_xlabel(r"$\epsilon$")
    fig.savefig(OUT / "chaos_boundary_exact.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def save_overlay(b, eps) -> None:
    data = np.load(OUT / "bifurcation_fig12.npz")
    e, theta = data["eps"], data["theta_deg"]
    top, bottom = (np.degrees(np.arccos(np.clip(b[f"{s}_{BIF_BOUND}"], -1, 1))) for s in ("top", "bottom"))
    fig, ax = plt.subplots(figsize=(9, 5))
    ax.plot(np.repeat(e, theta.shape[1]), theta.ravel(), "k.", ms=BIF_DOT, markeredgewidth=0, rasterized=True)
    ax.plot(eps, top, "r-", lw=0.8, label=r"predicted chaos bound from $N^-_{\rm out}$ (upper $\cos\theta_{sl}$)")
    ax.plot(eps, bottom, "c-", lw=0.8, label=r"predicted chaos bound from $N^+_{\rm out}$")
    ax.set(xscale="log", xlim=(e[0], e[-1]), ylim=(-5, 180), xlabel=r"$\epsilon$", ylabel=r"$\theta_{\rm sl}$ (deg)")
    ax.legend(loc="upper right")
    fig.savefig(OUT / "bifurcation_with_boundary.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def main() -> None:
    t0 = time.time()
    orb = make_orbit()
    sh = shape_functions(orb)
    sn = script_n(orb)
    psi0 = cos_coefficients(sh.psi, 1)[0]
    c_all = coupling(cos_coefficients(sh.beta, M_COEF), sin_coefficients(sh.gamma, M_COEF), M_COEF)
    print(f"script_N={sn:.6f} psi_0={psi0:.6f}")
    print_table(c_all, sn, psi0)
    save_width_figure(c_all, sn, psi0)
    eps = np.logspace(np.log10(EPS_MIN), np.log10(EPS_MAX), N_EPS)
    b = chaos_boundary_exact(c_all, M_COEF, sn, psi0, eps)
    np.savez(OUT / "chaos_boundary_exact.npz", epsilon=eps, **b)
    save_boundary_figure(b, eps)
    save_overlay(b, eps)
    print(f"runtime {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
