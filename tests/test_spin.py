"""Fast checks of lkspin.spin: norm drift, limits, scipy cross-check, convergence, (p, phi) round trip."""

import math

import numpy as np

from lkspin.lk_orbit import make_orbit
from lkspin.spin import ATOL_DEFAULT, H_INITIAL, RTOL_DEFAULT, initial_state, integrate_to, observe, run_cycles
from lkspin.spin_check import scipy_p_at_minima

N_CYCLES = 20
NORM_TOL = 1e-8
EPS_REGULAR = 0.015
EPS_NONADIABATIC = 1e3
EPS_ADIABATIC = 1e-3
ADIABATIC_MAX_DEG = 0.5
CROSS_TOL = 1e-6
CONVERGE_TOL_DEG = 1e-4
TIGHTEN = 10.0
THETA_SB_TOL_DEG = 0.2  # measured drift is about 0.07 deg over 20 cycles (physical, alpha/Omega_pl = 1e-3)
PHASE_TOL = 1e-12
EPS_NORM_CASES = (1.0, 0.1, 0.01)
ROUND_TRIP_CASES = ((0.3, 0.7), (-0.8, 4.0), (0.95, 0.0), (0.0, math.pi))

ORBIT = make_orbit()
ORB = ORBIT.params()


def run(eps, n=N_CYCLES, p0=1.0, phi0=0.0, scale=1.0):
    return run_cycles(initial_state(p0, phi0, ORB), n, ORBIT.period, ORBIT.alpha0(eps), ORB,
                      RTOL_DEFAULT * scale, ATOL_DEFAULT * scale)


def test_norm_drift_small():
    for eps in EPS_NORM_CASES:
        assert np.max(np.abs(run(eps)[:, 2] - 1.0)) < NORM_TOL


def test_nonadiabatic_keeps_theta_sb():
    y = initial_state(1.0, 0.0, ORB)
    sz0 = y[2]
    h = H_INITIAL
    for k in range(N_CYCLES):
        h = integrate_to(y, k * ORBIT.period, (k + 1) * ORBIT.period, h, ORBIT.alpha0(EPS_NONADIABATIC), ORB,
                         RTOL_DEFAULT, ATOL_DEFAULT)
        assert abs(math.degrees(math.acos(y[2]) - math.acos(sz0))) < THETA_SB_TOL_DEG


def test_adiabatic_stays_aligned():
    theta = np.degrees(np.arccos(np.clip(run(EPS_ADIABATIC, n=N_CYCLES)[:, 0], -1, 1)))
    assert theta.max() < ADIABATIC_MAX_DEG


def test_matches_scipy_regular():
    y0 = initial_state(1.0, 0.0, ORB)
    ref = scipy_p_at_minima(y0, N_CYCLES, ORBIT.period, ORBIT.alpha0(EPS_REGULAR), ORB)
    assert np.max(np.abs(run(EPS_REGULAR)[:, 0] - ref)) < CROSS_TOL


def test_tolerance_convergence():
    a = np.degrees(np.arccos(run(EPS_REGULAR)[:, 0]))
    b = np.degrees(np.arccos(run(EPS_REGULAR, scale=1.0 / TIGHTEN)[:, 0]))
    assert np.max(np.abs(a - b)) < CONVERGE_TOL_DEG


def test_p_phi_round_trip():
    for p0, phi0 in ROUND_TRIP_CASES:
        p, phi, _ = observe(initial_state(p0, phi0, ORB), 0.0, ORB)
        dphi = abs((phi - phi0 + math.pi) % (2 * math.pi) - math.pi)
        assert abs(p - p0) < PHASE_TOL
        assert dphi < PHASE_TOL or abs(p0) == 1.0
