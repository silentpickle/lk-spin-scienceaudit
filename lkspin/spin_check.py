"""Independent scipy DOP853 reference for the spin equations (used by tests/test_spin.py and scripts/step3_validate.py)."""

import math

import numpy as np
from scipy.integrate import solve_ivp

from lkspin.lk_orbit import alpha_of_x, cos_theta_lb_of_x, omega_pl_of_x, x_of_t

REF_RTOL = 1e-12
REF_ATOL = 1e-14


def scipy_p_at_minima(y0, n_cycles, period, alpha0, orb):
    """p = cos theta_sl at t = kP, k = 1..n_cycles, from scipy solve_ivp (DOP853)."""
    def f(t, y):
        x = x_of_t(t, orb)
        c = cos_theta_lb_of_x(x, orb)
        s = math.sqrt(1.0 - c * c)
        lx, ly, lz = s * math.sin(y[3]), -s * math.cos(y[3]), c
        k = -alpha_of_x(x, alpha0) * (y[0] * lx + y[1] * ly + y[2] * lz)
        return [k * (ly * y[2] - lz * y[1]), k * (lz * y[0] - lx * y[2]), k * (lx * y[1] - ly * y[0]),
                omega_pl_of_x(x, orb)]

    times = period * np.arange(1, n_cycles + 1)
    sol = solve_ivp(f, (0.0, times[-1]), y0, method="DOP853", t_eval=times, rtol=REF_RTOL, atol=REF_ATOL)
    out = np.empty(n_cycles)
    for i, t in enumerate(times):
        c = cos_theta_lb_of_x(x_of_t(t, orb), orb)
        s = math.sqrt(1.0 - c * c)
        y = sol.y[:, i]
        lvec = np.array([s * math.sin(y[3]), -s * math.cos(y[3]), c])
        out[i] = float(y[:3] @ lvec / np.linalg.norm(y[:3]))
    return out
