"""Fast checks of lkspin/numshape.py: hat coefficients against the analytic shape.py, and the Hamiltonian against secular."""

import math
import sys
from pathlib import Path

import numpy as np
from scipy.integrate import solve_ivp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "scripts"))
import prep_common as C  # noqa: E402
from lkspin import numshape as N  # noqa: E402
from lkspin import secular as S  # noqa: E402
from lkspin import shape as SH  # noqa: E402
from lkspin.lk_orbit import make_orbit  # noqa: E402

CFG4 = S.load_config()
THETA_ANALYTIC = 85.0
THETA_SRF = 87.0
N_T = 2**12
N_TAU = 2**11
N_DENSE_REF = 2**17
N_TAU_REF = 2**13
M_MAX = 10
COEF_TOL = 1e-6          # measured 2e-8 at these grids
TOL_RTOL = 1e-12
TOL_ATOL = 1e-14
N_CYCLES = 5
P_STAR = 12.0
P0, PHI0 = 0.5, 1.0
HAM_RTOL = 1e-11
HAM_ATOL = 1e-13
DP_TOL = 1e-6
DPHI_TOL = 1e-5
TWO_PI = 2.0 * math.pi


def test_hat_coefficients_match_analytic_shape():
    par = S.make_params(CFG4, srf=False, tides=False)
    sh = N.first_cycle_shapes(par, CFG4, THETA_ANALYTIC, N_T, N_TAU, M_MAX, TOL_RTOL, TOL_ATOL)
    orb = make_orbit(THETA_ANALYTIC, CFG4["e0"])
    ref = SH.shape_functions(orb, N_TAU_REF, N_DENSE_REF)
    sn = SH.script_n(orb, N_DENSE_REF)
    got = np.stack([sh.beta_cos, sh.gamma_sin, sh.psi_cos])
    want = sn * np.stack([SH.cos_coefficients(ref.beta, M_MAX), SH.sin_coefficients(ref.gamma, M_MAX),
                          SH.cos_coefficients(ref.psi, M_MAX)])
    assert np.abs(got - want).max() / (sn * SH.cos_coefficients(ref.beta, M_MAX)[0]) < COEF_TOL


def test_hamiltonian_reproduces_secular_spin_over_five_cycles():
    par = C.nondissipative_params(CFG4)
    sh = N.first_cycle_shapes(par, CFG4, THETA_SRF, N_T, N_TAU, M_MAX, TOL_RTOL, TOL_ATOL)
    times = C.orbit_minima_times(par, CFG4, THETA_SRF, N_CYCLES)
    y0 = S.initial_state(par, CFG4, THETA_SRF, 1)
    y0[S.N_ORBIT:S.N_ORBIT + S.N_SPIN] = C.initial_spin_vector(THETA_SRF, P0, PHI0)
    states, _ = C.run_to_times(y0, par, [P_STAR], times, TOL_RTOL, TOL_ATOL)
    lvec = states[:, None, 0:3]
    p_sec, phi_sec = C.phase_and_p(lvec, C.spin_block(states, 1))
    t_eval = TWO_PI * np.arange(1, N_CYCLES + 1)
    sol = solve_ivp(N.hamiltonian_rhs, (0.0, t_eval[-1]), [PHI0, P0], method="DOP853", t_eval=t_eval, rtol=HAM_RTOL,
                    atol=HAM_ATOL, args=(sh.nmax_coef / P_STAR, N.shape_interp(sh)))
    dphi = (sol.y[0] - phi_sec[:, 0] + math.pi) % TWO_PI - math.pi
    assert np.abs(sol.y[1] - p_sec[:, 0]).max() < DP_TOL
    assert np.abs(dphi).max() < DPHI_TOL
