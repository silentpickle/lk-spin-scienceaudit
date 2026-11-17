"""Fast checks of lkspin.secular: pure LK vs analytic, SRF e_max, conservation, spin vs step 3, full migration."""

import math

import numpy as np
from scipy.optimize import brentq

from lkspin import secular as S
from lkspin.lk_orbit import cos_theta_lb_of_x, make_orbit, x_of_t
from lkspin.spin import ATOL_DEFAULT, RTOL_DEFAULT, run_cycles
from lkspin.spin import initial_state as spin_initial_state

CFG = S.load_config()
A0 = CFG["a0_au"]
RTOL = 1e-11
ATOL = 1e-13
REC_MAX = 500_000
THETA_LK = 85.0
THETA_HIGH = 89.0
N_LK_PERIODS = 2
LK_TOL = 1e-8
EMAX_TOL = 1e-5
SPIN_TOL_DEG = 1e-4
N_MINIMA = 10
EPS_SPIN = 0.015
CONSERVE_TOL = 1e-9
N_CONSERVE_STEPS = 12
FIRST_CYCLE_TK = 8.0
P_STAR_SLOW = 30.0
P_STAR_FAST = 1.67
MIGRATION_MIN_MYR = 1000.0
MIGRATION_MAX_MYR = 2000.0
ORBIT_INDEPENDENCE_TOL = 1e-5
P_ORDER_TOL = 1e-4
OMEGA_P0 = 2.0 * math.pi * S.HOURS_PER_YR / CFG["p_planet0_hours"]


def test_pure_lk_matches_analytic():
    par = S.make_params(CFG, srf=False, tides=False)
    tk = S.t_k_years(par, A0)
    orb = make_orbit(THETA_LK, CFG["e0"])
    y0 = S.initial_state(par, CFG, THETA_LK, 1)
    rec = S.run(y0, par, [P_STAR_SLOW], N_LK_PERIODS * orb.period * tk, RTOL, ATOL, rec_max=REC_MAX)["rec"]
    x = np.array([x_of_t(t / tk, orb.params()) for t in rec[:, 0]])
    cos_th = np.array([cos_theta_lb_of_x(xi, orb.params()) for xi in x])
    assert np.abs(1 - rec[:, 2] ** 2 - x).max() < LK_TOL
    assert np.abs(np.cos(rec[:, 3]) - cos_th).max() < LK_TOL


def emax_from_a16(par, theta_deg):
    eg, et, er = S.srf_epsilons(par, A0, OMEGA_P0)
    c2 = math.cos(math.radians(theta_deg)) ** 2

    def f(e):
        j = math.sqrt(1 - e * e)
        lhs = eg * (1 / j - 1) + et / 15 * ((1 + 3 * e * e + 3 / 8 * e**4) / j**9 - 1) + er / 3 * (1 / j**3 - 1)
        return lhs - 9 * e * e / (8 * j * j) * (j * j - 5 / 3 * c2)

    return brentq(f, 0.02, 1 - 1e-10)


def test_srf_emax_matches_a16_potential_form():
    """With the bulge coefficient implied by A16's potential (15), e_max matches A16 Eq. emax; 7.5 (EOM) does not."""
    for theta in (THETA_HIGH, 87.0):
        tk = S.t_k_years(S.make_params(CFG, tides=False), A0)
        par15 = S.make_params(CFG, tides=False, bulge_coef=2 * S.BULGE_COEF_A16_EOM)
        y0 = S.initial_state(par15, CFG, theta, 1)
        rec = S.run(y0, par15, [P_STAR_SLOW], FIRST_CYCLE_TK * tk, RTOL, ATOL, rec_max=REC_MAX)["rec"]
        assert abs(rec[:, 2].max() - emax_from_a16(par15, theta)) < EMAX_TOL


def test_conservation_srf_only():
    par = S.make_params(CFG, tides=False)
    tk = S.t_k_years(par, A0)
    y = S.initial_state(par, CFG, THETA_HIGH, 1)
    y0 = y.copy()
    t, h = 0.0, S.H_INITIAL
    for k in range(1, N_CONSERVE_STEPS + 1):
        t, h, *_ = S.advance(y, t, k * tk / 3.0, h, par, np.array([S.omega_star_of_period(1.0)]), RTOL, ATOL,
                             10**9, np.zeros((0, S.N_REC_FIXED + 1)), 1)
        lmag, _, a = S.orbit_elements(y, par)
        assert abs(np.dot(y[3:6], y[0:3])) / lmag < CONSERVE_TOL
        assert abs(y[2] - y0[2]) < CONSERVE_TOL
        assert abs(a / A0 - 1) < CONSERVE_TOL
        assert abs(np.linalg.norm(y[7:10]) - 1) < CONSERVE_TOL


def test_spin_matches_step3_regular():
    par = S.make_params(CFG, srf=False, tides=False)
    tk = S.t_k_years(par, A0)
    orb = make_orbit(THETA_LK, CFG["e0"])
    om_star = orb.alpha0(EPS_SPIN) / tk * A0**3 / par[S.P_STAR_COEF]
    y, t, h = S.initial_state(par, CFG, THETA_LK, 1), 0.0, S.H_INITIAL
    got = []
    for k in range(1, N_MINIMA + 1):
        t, h, *_ = S.advance(y, t, k * orb.period * tk, h, par, np.array([om_star]), RTOL, ATOL, 10**9,
                             np.zeros((0, S.N_REC_FIXED + 1)), 1)
        got.append(S.theta_sl_final(y, 1)[0])
    ref = run_cycles(spin_initial_state(1.0, 0.0, orb.params()), N_MINIMA, orb.period, orb.alpha0(EPS_SPIN),
                     orb.params(), RTOL_DEFAULT, ATOL_DEFAULT)
    assert np.abs(np.array(got) - np.degrees(np.arccos(ref[:, 0]))).max() < SPIN_TOL_DEG


def test_full_migration_and_orbit_independent_of_p_star():
    par = S.make_params(CFG)
    y0 = S.initial_state(par, CFG, THETA_HIGH, 1)
    slow = S.run(y0, par, [P_STAR_SLOW], CFG["stop_t_gyr"] * 1e9)
    fast = S.run(y0, par, [P_STAR_FAST], CFG["stop_t_gyr"] * 1e9)
    assert slow["status"] == fast["status"] == S.ST_HOT_JUPITER
    assert MIGRATION_MIN_MYR < slow["t"] / 1e6 < MIGRATION_MAX_MYR
    assert abs(fast["t"] / slow["t"] - 1) < ORBIT_INDEPENDENCE_TOL
