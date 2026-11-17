"""Work plan step 3 (analytic side): reproduce S15 Figs. 2, 3, 5 and 11.

Run:  uv run python scripts/step3_fourier.py
Writes results/step3/fourier_coefficients.npz, fourier_shape_functions.png, fourier_coefficients.png,
fourier_widths.png, chaos_boundary.npz and chaos_boundary.png, and prints the S16 A_bar cross-check and the runtime.
"""

import sys
import time
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from lkspin.lk_orbit import make_orbit  # noqa: E402
from lkspin.shape import (  # noqa: E402
    BOUNDARY_FIELDS,
    a_bar_from_beta0,
    a_bar_from_time_average,
    chaos_boundary,
    cos_coefficients,
    coupling,
    script_n,
    shape_functions,
    sin_coefficients,
    widths,
)

OUT = Path(__file__).resolve().parent.parent / "results" / "step3"
M_SAVE = 200  # Fourier orders saved
N_SAVE = 200  # |N| saved for C_N
M_BOUNDARY = 500  # orders needed for epsilon down to 0.003 (N_max = 327)
EPS_WIDTHS = 0.1
EPS_MIN, EPS_MAX, N_EPS = 0.003, 1.0, 1500
DPI = 130
FIG_W, FIG_H = 8.0, 4.5


def save_shape_figure(sh) -> None:
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.plot(sh.tau, sh.beta, label=r"$\beta$")
    ax.plot(sh.tau, sh.gamma, label=r"$\gamma$")
    ax.plot(sh.tau, sh.psi, label=r"$\psi$")
    ax.set(xlabel=r"$\tau$", ylabel="shape function", xlim=(0, 2 * np.pi))
    ax.legend()
    fig.savefig(OUT / "fourier_shape_functions.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def save_coefficient_figure(n: np.ndarray, c: np.ndarray) -> None:
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    pos = c > 0
    ax.semilogy(n[pos], np.abs(c[pos]), "o", color="k", ms=4, label="positive")
    ax.semilogy(n[~pos], np.abs(c[~pos]), "o", mfc="none", color="k", ms=4, label="negative")
    ax.set(xlabel="N", ylabel=r"$|C_N|$  ($\beta_N+\gamma_N$ for $N>0$)", xlim=(-100, 100), ylim=(1e-5, 1))
    ax.legend()
    fig.savefig(OUT / "fourier_coefficients.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def save_width_figure(n: np.ndarray, c: np.ndarray, sn: float) -> None:
    keep = np.abs(n) <= np.floor(sn / EPS_WIDTHS)
    w = widths(c, n, EPS_WIDTHS, sn)
    fig, ax = plt.subplots(figsize=(FIG_W, FIG_H))
    ax.plot(n[keep], w[keep], "o-")
    ax.set(xlabel="N", ylabel=r"$w_N$", title=rf"approximate widths, $\epsilon={EPS_WIDTHS}$")
    fig.savefig(OUT / "fourier_widths.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def save_boundary_figure(b) -> None:
    fig, axes = plt.subplots(2, 1, figsize=(FIG_W, 2 * FIG_H), sharex=True)
    for ax, suffix, title in ((axes[0], "", "p = -N eps / script_N"), (axes[1], "_floor", "p = -N / floor(N_max)")):
        ax.fill_between(b.epsilon, getattr(b, "bottom_any" + suffix), getattr(b, "top_any" + suffix), color="0.75",
                        label="largest |N| anywhere (shaded)")
        ax.plot(b.epsilon, getattr(b, "top_contig" + suffix), "b-", lw=0.7, label="contiguous chain")
        ax.plot(b.epsilon, getattr(b, "bottom_contig" + suffix), "b-", lw=0.7)
        ax.set(xscale="log", ylabel=r"$\cos\theta_{\rm sl}$", xlim=(EPS_MIN, EPS_MAX), ylim=(-1, 1), title=title)
        ax.legend(loc="center right")
    axes[1].set_xlabel(r"$\epsilon$")
    fig.savefig(OUT / "chaos_boundary.png", dpi=DPI, bbox_inches="tight")
    plt.close(fig)


def report_a_bar(orb, beta0: float, sn: float) -> None:
    a_tau = a_bar_from_beta0(beta0, EPS_WIDTHS)
    a_time = a_bar_from_time_average(orb, EPS_WIDTHS, sn)
    print(f"A_bar at eps={EPS_WIDTHS}: 1/(2 eps beta_0) = {a_tau:.8g};  N_max/(2 Omega_hat_L) = {a_time:.8g};"
          f"  ratio = {a_tau / a_time:.3e}")


def main() -> None:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    orb = make_orbit()
    sh = shape_functions(orb)
    sn = script_n(orb)
    beta_m = cos_coefficients(sh.beta, M_BOUNDARY)
    gamma_m = sin_coefficients(sh.gamma, M_BOUNDARY)
    psi_m = cos_coefficients(sh.psi, M_BOUNDARY)
    c_all = coupling(beta_m, gamma_m, M_BOUNDARY)
    n_all = np.arange(-M_BOUNDARY, M_BOUNDARY + 1)
    sel = np.abs(n_all) <= N_SAVE
    np.savez(OUT / "fourier_coefficients.npz", beta_M=beta_m[: M_SAVE + 1], gamma_M=gamma_m[: M_SAVE + 1],
             psi_M=psi_m[: M_SAVE + 1], N=n_all[sel], C_N=c_all[sel], script_N=sn)
    print(f"script_N = {sn:.6f}  beta_0 = {beta_m[0]:.6f}  C_0 = {c_all[M_BOUNDARY]:.6f}  psi_0 = {psi_m[0]:.6f}")
    save_shape_figure(sh)
    save_coefficient_figure(n_all, c_all)
    save_width_figure(n_all, c_all, sn)
    report_a_bar(orb, beta_m[0], sn)
    eps = np.logspace(np.log10(EPS_MIN), np.log10(EPS_MAX), N_EPS)
    b = chaos_boundary(c_all, M_BOUNDARY, sn, eps)
    np.savez(OUT / "chaos_boundary.npz", epsilon=eps, **{f: getattr(b, f) for f in BOUNDARY_FIELDS})
    save_boundary_figure(b)
    print(f"runtime {time.time() - t0:.1f} s")


if __name__ == "__main__":
    main()
