"""Fast checks for lkspin.shape (S15 shape functions, Fourier coefficients, widths, chaos boundary).

Run:  uv run pytest tests/test_shape.py
"""

import numpy as np
import pytest

from lkspin.lk_orbit import make_orbit
from lkspin.shape import (
    a_bar_from_beta0,
    a_bar_from_time_average,
    anywhere_outermost,
    contiguous_outermost,
    cos_coefficients,
    exact_resonances,
    widths,
    coupling,
    script_n,
    shape_functions,
    sin_coefficients,
)

N_DENSE_FAST = 2**18
N_TAU_FAST = 2**13
M_RECON = 1000
M_LOW = 40  # "low-order" coefficients compared across resolutions
SCRIPT_N_PUBLISHED = 0.98
SCRIPT_N_COMPUTED = 0.9792
SCRIPT_N_TOL = 1e-4
RECON_TOL = 1e-7
RES_TOL = 1e-7
MEAN_TOL = 1e-10
SYM_TOL = 1e-7
ABAR_RTOL = 1e-7
EPS_TEST = 0.1
SPACING = 0.2
WIDE = 0.3
EPS_SMALL = 1e-3
EPS_MEDIUM = 1e-2
N_COEF_TEST = 50
N_CHECK = (-3, 0, 3)
WIDTH_REL_TOL = 1e-3  # exact/approx - 1 measured <= 2.2e-4 at eps = 1e-3
WIDTH_SHRINK = 5.0  # discrepancy shrinks by at least this when eps drops 10x (measured about 10)
CENTRE_SLOPE_TOL = 1e-3  # |(p_c - p_N)/eps - psi_0| measured <= 5e-4 at eps = 1e-3
CENTRE_SECOND_ORDER = 1.0  # |p_c - (p_N + eps psi_0)| / eps^2 measured <= 0.46
NARROW = 0.01


@pytest.fixture(scope="module")
def orb():
    return make_orbit()


@pytest.fixture(scope="module")
def shapes(orb):
    return shape_functions(orb, N_TAU_FAST, N_DENSE_FAST)


def test_script_n_matches_paper(orb):
    value = script_n(orb, N_DENSE_FAST)
    assert abs(value - SCRIPT_N_COMPUTED) < SCRIPT_N_TOL
    assert abs(value - SCRIPT_N_PUBLISHED) < 2e-3


def test_fourier_series_reconstructs_shapes(shapes):
    tau = shapes.tau
    m = np.arange(M_RECON + 1)[:, None] * tau[None, :]
    beta = cos_coefficients(shapes.beta, M_RECON) @ np.cos(m)
    gamma = sin_coefficients(shapes.gamma, M_RECON) @ np.sin(m)
    psi = cos_coefficients(shapes.psi, M_RECON) @ np.cos(m)
    assert abs(beta - shapes.beta).max() < RECON_TOL
    assert abs(gamma - shapes.gamma).max() < RECON_TOL
    assert abs(psi - shapes.psi).max() < RECON_TOL


def test_gamma_has_zero_mean(shapes):
    assert abs(shapes.gamma.mean()) < MEAN_TOL


def test_symmetry_about_tau_pi(shapes):
    """beta, psi even and gamma odd about tau = pi (grid index n -> N - n)."""
    assert abs(shapes.beta[1:] - shapes.beta[:0:-1]).max() < SYM_TOL
    assert abs(shapes.psi[1:] - shapes.psi[:0:-1]).max() < SYM_TOL
    assert abs(shapes.gamma[1:] + shapes.gamma[:0:-1]).max() < SYM_TOL


def test_low_order_coefficients_resolution_stable(orb, shapes):
    fine = shape_functions(orb, 2 * N_TAU_FAST, 2 * N_DENSE_FAST)
    for coarse_f, fine_f, fn in ((shapes.beta, fine.beta, cos_coefficients), (shapes.gamma, fine.gamma, sin_coefficients),
                                 (shapes.psi, fine.psi, cos_coefficients)):
        assert abs(fn(coarse_f, M_LOW) - fn(fine_f, M_LOW)).max() < RES_TOL


def test_coupling_sign_bookkeeping():
    beta_m = np.array([2.0, 5.0, 7.0])
    gamma_m = np.array([0.0, 1.0, 3.0])
    c = coupling(beta_m, gamma_m, 2)  # N = -2..2
    assert list(c) == [7.0 - 3.0, 5.0 - 1.0, 2 * 2.0, 5.0 + 1.0, 7.0 + 3.0]


def test_a_bar_tau_and_time_averages_agree(orb, shapes):
    beta0 = cos_coefficients(shapes.beta, 1)[0]
    sn = script_n(orb, N_DENSE_FAST)
    tau_avg = a_bar_from_beta0(beta0, EPS_TEST)
    time_avg = a_bar_from_time_average(orb, EPS_TEST, sn, N_DENSE_FAST)
    assert tau_avg == pytest.approx(time_avg, rel=ABAR_RTOL)


def test_chain_readings_differ_on_broken_chain():
    w = np.array([1.0, WIDE, WIDE, NARROW, NARROW, WIDE, WIDE])  # chain breaks at k=3, k=6 overlaps again
    assert contiguous_outermost(w, SPACING) == 2
    assert anywhere_outermost(w, SPACING) == 6


@pytest.fixture(scope="module")
def coeffs(orb, shapes):
    beta = cos_coefficients(shapes.beta, N_COEF_TEST)
    gamma = sin_coefficients(shapes.gamma, N_COEF_TEST)
    psi0 = cos_coefficients(shapes.psi, 1)[0]
    return coupling(beta, gamma, N_COEF_TEST), psi0, script_n(orb, N_DENSE_FAST)


def width_discrepancy(coeffs, n: int, eps: float) -> float:
    c, psi0, sn = coeffs
    arr = np.array([n])
    exact = exact_resonances(c[arr + N_COEF_TEST], arr, eps, sn, psi0).width[0]
    approx = widths(c[arr + N_COEF_TEST], arr, eps, sn)[0]
    return abs(exact / approx - 1.0)


@pytest.mark.parametrize("n", N_CHECK)
def test_exact_width_tends_to_approximate(coeffs, n):
    small = width_discrepancy(coeffs, n, EPS_SMALL)
    assert small < WIDTH_REL_TOL
    assert small < width_discrepancy(coeffs, n, EPS_MEDIUM) / WIDTH_SHRINK


@pytest.mark.parametrize("n", N_CHECK)
def test_exact_centre_first_order(coeffs, n):
    c, psi0, sn = coeffs
    arr = np.array([n])
    pc = exact_resonances(c[arr + N_COEF_TEST], arr, EPS_SMALL, sn, psi0).centre[0]
    nominal = -n * EPS_SMALL / sn
    assert abs((pc - nominal) / EPS_SMALL - psi0) < CENTRE_SLOPE_TOL
    assert abs(pc - nominal - EPS_SMALL * psi0) < CENTRE_SECOND_ORDER * EPS_SMALL**2
