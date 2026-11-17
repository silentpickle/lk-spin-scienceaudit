"""Fast checks of the step 5 diagnostics (lkspin.secular.run_diag) and first-pass classifier (lkspin.regime)."""

import math
import tomllib

import numpy as np

from lkspin import regime as R
from lkspin import secular as S
from lkspin.lk_orbit import make_orbit

CFG4 = S.load_config()
with S.CONFIG_PATH.open("rb") as _fh:
    CFG5 = tomllib.load(_fh)["step5"]
THETA_LK = 85.0
N_CYCLES = 6
PERIOD_TOL = 1e-6
P_SPINS = [30.0, 5.0, 1.0]
EPS_LOCK = 0.0869
EPS_ADIABATIC = 0.008
N_LOCK = 11
LOCK_DEV = 0.1
ADIABATIC_REL = 1e-3
STEP3_CYCLES = 60
SKIP = 20
SPOT_P = 7.07
TWIN_ANGLE = 1e-8
ANGLE_TOL = 1e-12


def pure_lk_run(p_list, n_cycles):
    par = S.make_params(CFG4, srf=False, tides=False)
    tk = S.t_k_years(par, CFG4["a0_au"])
    orb = make_orbit(THETA_LK, CFG4["e0"])
    y0 = S.initial_state(par, CFG4, THETA_LK, len(p_list))
    return par, tk, orb, y0, S.run_diag(y0, par, p_list, n_cycles * orb.period * tk)


def test_diag_matches_plain_bitwise():
    par, tk, orb, y0, res = pure_lk_run(P_SPINS, N_CYCLES)
    plain = S.run(y0, par, P_SPINS, N_CYCLES * orb.period * tk)
    assert np.array_equal(plain["y"], res["y"])
    assert plain["t"] == res["t"] and plain["steps"] == res["steps"]


def test_minima_spaced_by_lk_period():
    _, tk, orb, _, res = pure_lk_run(P_SPINS, N_CYCLES)
    minima, maxima = R.split_events(res["ev"])
    assert np.abs(np.diff(minima[:, 0]) / (orb.period * tk) - 1).max() < PERIOD_TOL
    assert maxima.shape[0] >= N_CYCLES - 1


def step3_winding(eps):
    par = S.make_params(CFG4, srf=False, tides=False)
    a0 = CFG4["a0_au"]
    tk = S.t_k_years(par, a0)
    orb = make_orbit(THETA_LK, CFG4["e0"])
    p = 2 * math.pi * S.DAYS_PER_YR / (orb.alpha0(eps) / tk * a0**3 / par[S.P_STAR_COEF])
    res = S.run_diag(S.initial_state(par, CFG4, THETA_LK, 1), par, [p], STEP3_CYCLES * orb.period * tk)
    minima, _ = R.split_events(res["ev"])
    return R.winding_series(minima, 0, 1)[SKIP:], par, p, tk, orb


def test_winding_locks_to_integer_in_resonance():
    w, *_ = step3_winding(EPS_LOCK)
    assert np.abs(w - N_LOCK).max() < LOCK_DEV


def test_winding_equals_nmax_when_adiabatic():
    w, par, p, tk, orb = step3_winding(EPS_ADIABATIC)
    nmax = R.shape_coefficients(par, CFG4, THETA_LK, 4 * orb.period * tk)["nmax_coef"] / p
    assert abs(w.mean() / nmax - 1) < ADIABATIC_REL


def test_twin_is_rotated_by_requested_angle():
    par = S.make_params(CFG4)
    y0 = S.initial_state(par, CFG4, THETA_LK, 2)
    y = R.twin_state(y0, 1, TWIN_ANGLE)
    lhat = y[0:3] / np.linalg.norm(y[0:3])
    ang = np.linalg.norm(np.cross(R.spin_vector(y, 1), lhat))  # sin(angle), exact at tiny angles
    assert abs(ang - TWIN_ANGLE) < ANGLE_TOL
    assert np.array_equal(R.spin_vector(y, 0), R.spin_vector(y0, 0))


def features(**kw):
    base = {"gap": 1.0, "n_low": 5, "n_high": 5, "range_low": 1.0, "range_high": 1.0, "range": 2.0,
            "twin_div_deg": 0.0, "theta_max": 10.0, "m0_locked_frac_last": 0.0, "m0_n_dom": math.nan,
            "m0_w_median_last": 0.5}
    return {**base, **kw}


def test_classifier_each_branch():
    c = CFG5
    assert R.classify(features(gap=90.0, range=95.0, range_low=3.0, range_high=3.0), c) == "I"
    assert R.classify(features(twin_div_deg=40.0), c) == "III"
    assert R.classify(features(m0_locked_frac_last=0.9, m0_n_dom=2.0, theta_max=60.0), c) == "IV.2"
    assert R.classify(features(theta_max=1.0, m0_w_median_last=9.0), c) == "V"
    assert R.classify(features(), c) == "II"
    assert R.classify(features(theta_max=80.0), c) == "unclassified"


def test_classifier_threshold_edges():
    c = CFG5
    gap = c["bimodal_gap_deg"]
    edge = features(gap=gap, range=gap + 2.0, range_low=c["bimodal_side_range_deg"], range_high=1.0)
    assert R.classify(edge, c) == "I"
    assert R.classify({**edge, "gap": gap - 1e-9}, c) != "I"
    assert R.classify({**edge, "range_low": c["bimodal_side_range_deg"] + 1e-9}, c) != "I"
    assert R.classify(features(twin_div_deg=c["chaos_twin_deg"]), c) == "III"
    assert R.classify(features(twin_div_deg=c["chaos_twin_deg"] - 1e-9), c) != "III"
    assert R.classify(features(m0_locked_frac_last=c["advect_lock_frac"], m0_n_dom=1.0), c) == "IV.1"
    assert R.classify(features(m0_locked_frac_last=c["advect_lock_frac"] - 1e-9, m0_n_dom=1.0), c) != "IV.1"


def test_predicted_class_map():
    c = CFG5
    assert R.predicted_class(0.5, 3.0, c) == "I"
    assert R.predicted_class(1.0, 0.1, c) == "II"
    assert R.predicted_class(1.0, 0.5, c) == "IV"


def test_spot_case_labels_ii():
    f = R.integrate_system(5.0, SPOT_P, 300.0, 89.0, CFG4, CFG5)
    assert R.classify(f, CFG5) == "II"
