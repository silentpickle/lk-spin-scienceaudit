"""Fast checks of the prep job helpers (scripts/prep_*.py)."""

import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import prep_common as C  # noqa: E402
import prep_sens4 as P3  # noqa: E402
import prep_sq4 as P2  # noqa: E402
from lkspin import secular as S  # noqa: E402
from lkspin.lk_orbit import make_orbit  # noqa: E402
from lkspin.spin import initial_state as spin_initial_state  # noqa: E402

CFG4 = S.load_config()
THETA = 87.0
N_CYC = 60
P0, PHI0 = 0.3, 2.0
ANGLE_TOL = 1e-12
TINY_ANGLE = 1e-9
GROWTH = 0.4
N_FIT = 40
SAT_TEST = 1e-8
FIT_TOL = 1e-9
FTLE_CHAOS_P, FTLE_REGULAR_P = 3.0, 10.0
CHAOS_MIN = 0.02
REGULAR_MAX = 1e-3
E_MIN_TOL = 1e-5
TIME_RTOL = 1e-12
N_SENS_VARIANTS = 12
N_GRID_TOTAL = 200 + 1 + 61
GAP_BIMODAL = 90.0
N_SYNTH = 100
RNG_SEED = 3
P_HIGH_TRUE = 0.3
P_HIGH_TOL = 1e-9


def test_phase_roundtrip_matches_step3_convention():
    orb = make_orbit(THETA, CFG4["e0"]).params()
    s_spin = spin_initial_state(P0, PHI0, orb)[:3]
    s_new = C.initial_spin_vector(THETA, P0, PHI0)
    assert np.abs(s_spin - s_new).max() < ANGLE_TOL
    th = math.radians(THETA)
    lvec = np.array([0.0, -math.sin(th), math.cos(th)])
    p, phi = C.phase_and_p(lvec, s_new)
    assert abs(p - P0) < ANGLE_TOL and abs(phi - PHI0) < ANGLE_TOL


def test_angle_between_tiny_angle():
    u = np.array([0.0, 0.0, 1.0])
    v = np.array([math.sin(TINY_ANGLE), 0.0, math.cos(TINY_ANGLE)])
    assert abs(C.angle_between(u, v) - TINY_ANGLE) < ANGLE_TOL


def test_ftle_fit_recovers_exponential_growth():
    k = np.arange(1, N_FIT + 1)
    sep = TINY_ANGLE * np.exp(GROWTH * k)
    fit = C.ftle_fit(sep, TINY_ANGLE, SAT_TEST, 4)
    assert abs(fit["ftle"] - GROWTH) < FIT_TOL and fit["k_sat"] > 0


def test_branch_summary_two_clusters():
    rng = np.random.default_rng(RNG_SEED)
    high = rng.random(N_SYNTH) < P_HIGH_TRUE
    theta = np.where(high, 135.0 + rng.uniform(-3, 3, N_SYNTH), 40.0 + rng.uniform(-3, 3, N_SYNTH))
    b = C.branch_summary(theta, 25.0, 30.0)
    assert b["n_branches"] == 2 and abs(b["p_high"] - high.mean()) < P_HIGH_TOL
    assert b["peak_low_deg"] < GAP_BIMODAL < b["peak_high_deg"]


def test_run_to_times_lands_on_minima():
    par = C.nondissipative_params(CFG4)
    times = C.orbit_minima_times(par, CFG4, THETA, 3)
    y0 = S.initial_state(par, CFG4, THETA, 1)
    states, status = C.run_to_times(y0, par, [10.0], times, 1e-10, 1e-12)
    e = np.linalg.norm(states[:, 3:6], axis=1)
    assert status == S.ST_TIMEOUT and np.abs(e - CFG4["e0"]).max() < E_MIN_TOL


def test_twin_ftle_chaotic_vs_regular():
    par = C.nondissipative_params(CFG4)
    times = C.orbit_minima_times(par, CFG4, THETA, 100)
    res = P2.ftle_integrate(par, THETA, np.array([FTLE_CHAOS_P, FTLE_REGULAR_P]), times)
    assert res["ftle"][0] > CHAOS_MIN and res["ftle"][1] < REGULAR_MAX


def test_sens4_variants_and_grid():
    assert len(P3.variants()) == N_SENS_VARIANTS
    p, n_main, n_extra = P3.p_grid()
    assert p.size == N_GRID_TOTAL and p[0] == P3.CFG["p_star_max_days"] and abs(p[n_main - 1] - P3.CFG["p_star_min_days"]) < FIT_TOL
    cfg, par = P3.variant_setup("bulge_15")
    assert par[S.P_BULGE] == P3.CFG["variant"]["bulge_15"]["bulge_coef"] and cfg["chi"] == CFG4["chi"]
