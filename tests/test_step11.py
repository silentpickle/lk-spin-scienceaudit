"""Fast checks of the step 11 pieces: N-body initial state, GR force against the 1PN formula, peak finder, flags."""

import math
import sys

import numpy as np

sys.path.insert(0, "scripts")
import step11_common as C  # noqa: E402
import step11_nbody as N  # noqa: E402

GR_SPAN_YR = 2000.0
GR_E0 = 0.5
GR_RATIO_TOL = 1e-3
ELEMENT_TOL = 1e-12
FIDUCIAL = C.System("T", "fiducial", 5.0, 300.0, 89.0, 1.0, 2.0)


def test_initial_orbit_matches_secular_state():
    sim, _rebx = N.build_sim(FIDUCIAL, C.CFG4["e0"])
    o = sim.particles[N.PLANET].orbit()
    assert abs(o.a - C.CFG4["a0_au"]) < ELEMENT_TOL
    assert abs(o.e - C.CFG4["e0"]) < ELEMENT_TOL
    assert abs(math.degrees(o.inc) - FIDUCIAL.theta_lb0_deg) < ELEMENT_TOL
    assert abs(o.omega) < ELEMENT_TOL


def test_gr_precession_matches_1pn_formula():
    sim, _rebx = N.build_sim(FIDUCIAL, GR_E0, m_b=1e-12)
    om0 = sim.particles[N.PLANET].orbit().pomega
    sim.integrate(GR_SPAN_YR)
    got = sim.particles[N.PLANET].orbit().pomega - om0
    mtot = C.CFG4["m_star_msun"] + FIDUCIAL.m_p_mjup * C.S.MJUP_PER_MSUN
    a = C.CFG4["a0_au"]
    n = math.sqrt(C.S.G_AU3_PER_MSUN_YR2 * mtot / a**3)
    want = 3 * C.S.G_AU3_PER_MSUN_YR2 * mtot / (C.CFG["c_au_per_yr"] ** 2 * a * (1 - GR_E0**2)) * n * GR_SPAN_YR
    assert abs(got / want - 1) < GR_RATIO_TOL


def test_peak_finder_ignores_small_wiggles():
    t = np.linspace(0, 1, 2000)
    e = 1 - 10 ** (-2.0 * np.sin(np.pi * t) ** 2 - 0.01)
    e += 1e-6 * np.sin(200 * t)
    assert len(C.peak_indices(t, e)) == 1


def test_flags_flag_close_to_roche():
    b = {"t_near_over_pout": 100.0, "rperi_over_rtide": 0.5, "rperi_over_rtide_gr": 2.0,
         **{f"porb_alpha_P{p:g}d": 0.0 for p in C.CFG["p_star_days"]}}
    assert C.flags(b)["flag_roche"]
    assert not C.flags(b)["flag_double_avg"]
