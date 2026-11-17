"""Fast checks of the opt-in step 6 physics in lkspin.secular: feedback, octupole, spin-down (see scripts/step6_validate.py)."""

import math
import sys

import numpy as np

sys.path.insert(0, "scripts")
import step6_validate as V  # noqa: E402

from lkspin import secular as S  # noqa: E402

J_TOL = 1e-12
JZ_TOL = 1e-10
ENERGY_RATIO_TOL = 1e-6
SKUMANICH_TOL = 1e-8
MIRROR_TOL = 1e-12
PRINTED_MIN_RESIDUAL = 1e-2
LINEAR_RATIO_LOW = 1.8
LINEAR_RATIO_HIGH = 2.3
ZERO_FEEDBACK_TOL_DEG = 1e-5
ZERO_FEEDBACK_ROW = 4  # a regular (class I) step 5 system
JSNAP_MIN_CHANGE = 0.1  # |J| changes by more than this fraction when LK is on (the binary torques the orbit)
NONZERO_FEEDBACK_MIN_DEG = 1e-3
EB_TEST = 0.7
GAP_FACTOR_TOL = 1e-12


def test_octupole_equations_keep_a_and_orthogonality_printed_form_does_not():
    res = V.check_octupole_consistency()
    assert res["mirrored_a_resid"] < MIRROR_TOL
    assert res["mirrored_le_resid"] < MIRROR_TOL
    assert res["printed_a_resid"] > PRINTED_MIN_RESIDUAL
    assert res["printed_le_resid"] > PRINTED_MIN_RESIDUAL


def test_skumanich_matches_analytic_solution():
    res = V.check_skumanich()
    assert max(abs(r) for r in res["rel_error"]) < SKUMANICH_TOL


def test_total_angular_momentum_conserved_with_lk_off():
    res = V.check_conservation()
    assert res["lk_off_max_dJ_over_J"] < J_TOL
    assert res["lk_off_precession_periods"] > 1.0


def test_jz_and_energy_conserved_with_lk_quadrupole_and_feedback():
    res = V.check_conservation()
    assert res["lk_on_max_dJz_over_Jz"] < JZ_TOL
    assert res["lk_on_max_dJ_over_J"] > JSNAP_MIN_CHANGE
    assert res["energy_conservation_ratio"] < ENERGY_RATIO_TOL


def test_octupole_breaks_quadrupole_invariants_linearly_in_eps():
    res = V.check_octupole_scaling()
    for ratio in res["drift_jz_ratios"] + res["drift_HQ_ratios"]:
        assert LINEAR_RATIO_LOW < ratio < LINEAR_RATIO_HIGH


def test_zero_feedback_matches_single_spin_integrator():
    table = V.W.step5_table()
    i = ZERO_FEEDBACK_ROW
    m_p, p_star, a_b, theta = (float(table[k][i]) for k in V.W.INPUT_FIELDS)
    cfg = V.W.R.system_config(V.CFG4, m_p, a_b)
    p_days = V.W.member_periods(p_star)[0][0]
    got = V.ext_theta(cfg, theta, p_days, S.RTOL_DEFAULT, S.ATOL_DEFAULT)
    ref = V.single_spin_theta(cfg, theta, p_days, S.RTOL_DEFAULT, S.ATOL_DEFAULT)
    assert abs(got - ref) < ZERO_FEEDBACK_TOL_DEG


def test_feedback_changes_the_outcome():
    table = V.W.step5_table()
    i = ZERO_FEEDBACK_ROW
    m_p, p_star, a_b, theta = (float(table[k][i]) for k in V.W.INPUT_FIELDS)
    cfg = V.W.R.system_config(V.CFG4, m_p, a_b)
    p_days = V.W.member_periods(p_star)[0][0]
    par = S.extend_params(S.make_params(cfg), fb_mult=1.0)
    y0 = S.initial_state_ext(par, cfg, theta, p_days)
    on = S.theta_sl_ext(S.run_ext(y0, par, cfg["stop_t_gyr"] * 1e9)["y"])
    off = V.ext_theta(cfg, theta, p_days, S.RTOL_DEFAULT, S.ATOL_DEFAULT)
    assert abs(on - off) > NONZERO_FEEDBACK_MIN_DEG


def test_extend_params_keeps_quadrupole_timescale_and_sets_eps_oct():
    base = S.make_params(V.CFG4)
    par = S.extend_params(base, e_b=EB_TEST, octupole=True)
    assert np.array_equal(par[:S.NPAR], base)
    a_b_real = base[S.P_AB] / math.sqrt(1.0 - EB_TEST**2)
    mass = (base[S.P_MSTAR] - base[S.P_MP]) / (base[S.P_MSTAR] + base[S.P_MP])
    expected = mass * (V.CFG4["a0_au"] / a_b_real) * EB_TEST / (1.0 - EB_TEST**2)
    assert abs(par[S.P_EOCT] * V.CFG4["a0_au"] - expected) < GAP_FACTOR_TOL
    assert S.extend_params(base)[S.P_EOCT] == 0.0
