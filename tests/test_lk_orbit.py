"""The analytic LK orbit agrees with scipy's cn and with direct integration of S15 Eqs. (ek)-(omk)."""

import math

import numpy as np
from scipy.integrate import solve_ivp
from scipy.special import ellipj

from lkspin.lk_orbit import cos_theta_lb_of_x, jacobi_cn, make_orbit, omega_pl_of_x, x_of_t

CN_TOL = 1e-12
ORBIT_TOL = 1e-7
RTOL = 1e-12
ATOL = 1e-13
SAMPLES = 400
PUBLISHED_SCRIPT_N = 0.98  # S15 l.634, checked in scripts/calibrate.py


def test_cn_matches_scipy():
    orb = make_orbit()
    for u in np.linspace(-3.0 * orb.big_k, 5.0 * orb.big_k, SAMPLES):
        assert abs(jacobi_cn(u, orb.m) - ellipj(u, orb.m)[1]) < CN_TOL


def test_x_endpoints_are_extrema():
    orb = make_orbit()
    p = orb.params()
    assert abs(x_of_t(0.0, p) - orb.x0) < CN_TOL
    assert abs(x_of_t(orb.period / 2.0, p) - orb.x1) < CN_TOL
    assert abs(x_of_t(orb.period, p) - orb.x0) < CN_TOL


def lk_rhs(_t, y):
    """S15 Eqs. (ek)-(omk) with t_k = 1. y = (e, Omega, theta_lb, omega)."""
    e, _node, inc, w = y
    j = math.sqrt(1.0 - e * e)
    de = 15.0 / 8.0 * e * j * math.sin(2 * w) * math.sin(inc) ** 2
    dnode = 0.75 * math.cos(inc) * (5 * e * e * math.cos(w) ** 2 - 4 * e * e - 1) / j
    dinc = -15.0 / 16.0 * e * e * math.sin(2 * w) * math.sin(2 * inc) / j
    dw = 3.0 * (2 * (1 - e * e) + 5 * math.sin(w) ** 2 * (e * e - math.sin(inc) ** 2)) / (4 * j)
    return [de, dnode, dinc, dw]


def test_analytic_orbit_matches_direct_integration():
    orb = make_orbit()
    p = orb.params()
    times = np.linspace(0.0, orb.period, SAMPLES)
    y0 = [orb.e0, 0.0, math.radians(orb.theta_lb0_deg), 0.0]
    sol = solve_ivp(lk_rhs, (0.0, orb.period), y0, method="DOP853", t_eval=times, rtol=RTOL, atol=ATOL)
    e, node, inc = sol.y[0], sol.y[1], sol.y[2]
    x_direct = 1.0 - e**2
    x_analytic = np.array([x_of_t(t, p) for t in times])
    assert np.max(np.abs(x_direct - x_analytic)) < ORBIT_TOL
    cos_analytic = np.array([cos_theta_lb_of_x(x, p) for x in x_analytic])
    assert np.max(np.abs(np.cos(inc) - cos_analytic)) < ORBIT_TOL
    # Node: integrate the analytic Omega_pl and compare with the direct node angle.
    node_sol = solve_ivp(lambda t, _y: [omega_pl_of_x(x_of_t(t, p), p)], (0.0, orb.period), [0.0],
                         method="DOP853", t_eval=times, rtol=RTOL, atol=ATOL)
    assert np.max(np.abs(node - node_sol.y[0])) < ORBIT_TOL
