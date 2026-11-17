"""Stellar-spin evolution in the inertial frame under an analytic Lidov-Kozai orbit (Storch & Lai 2015, S15).

State y = (Sx, Sy, Sz, Omega); L_b = z-hat; dS/dt = -alpha (S.L)(L x S); dOmega/dt = Omega_pl(x(t)).
Adaptive Dormand-Prince 5(4) (numba) integrates one LK period at a time and lands exactly on t = k P (eccentricity
minima), where p = cos(theta_sl) and the precession phase phi (S15 convention) are recorded.
Used by scripts/step3_bifurcation.py and scripts/step3_sections.py; checked by tests/test_spin.py.
Run the checks:  uv run pytest tests/test_spin.py
"""

import math

import numpy as np
from numba import njit

from lkspin.lk_orbit import alpha_of_x, cos_theta_lb_of_x, omega_pl_of_x, x_of_t

TWO_PI = 2.0 * math.pi
RTOL_DEFAULT = 1e-10
ATOL_DEFAULT = 1e-12
H_INITIAL = 1e-3
H_MIN = 1e-14
SAFETY = 0.9
FAC_MIN = 0.2
FAC_MAX = 5.0
ERR_EXPONENT = 0.2  # 1/(order+1) with order 4 for the embedded estimate
MAX_STEPS = 2_000_000
NSTATE = 4

# Output columns of run_cycles
COL_P, COL_PHI, COL_NORM = range(3)
NCOL = 3

# Dormand-Prince 5(4) tableau
C2, C3, C4, C5 = 1 / 5, 3 / 10, 4 / 5, 8 / 9
A21 = 1 / 5
A31, A32 = 3 / 40, 9 / 40
A41, A42, A43 = 44 / 45, -56 / 15, 32 / 9
A51, A52, A53, A54 = 19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729
A61, A62, A63, A64, A65 = 9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656
B1, B3, B4, B5, B6 = 35 / 384, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84
E1, E3, E4, E5, E6, E7 = (71 / 57600, -71 / 16695, 71 / 1920, -17253 / 339200, 22 / 525, -1 / 40)


@njit(cache=True)
def rhs(t, y, alpha0, orb, out):
    """Inertial-frame spin equations. Writes dy/dt into out."""
    x = x_of_t(t, orb)
    cos_lb = cos_theta_lb_of_x(x, orb)
    sin_lb = math.sqrt(max(0.0, 1.0 - cos_lb * cos_lb))
    om = y[3]
    lx, ly, lz = sin_lb * math.sin(om), -sin_lb * math.cos(om), cos_lb
    sx, sy, sz = y[0], y[1], y[2]
    sdl = sx * lx + sy * ly + sz * lz
    f = -alpha_of_x(x, alpha0) * sdl
    out[0] = f * (ly * sz - lz * sy)
    out[1] = f * (lz * sx - lx * sz)
    out[2] = f * (lx * sy - ly * sx)
    out[3] = omega_pl_of_x(x, orb)


@njit(cache=True)
def _stage(t, y, h, k, coef, alpha0, orb, tmp, out):
    """out = f(t, y + h * sum(coef[i] k[i])) for the stages given in coef."""
    for j in range(NSTATE):
        s = 0.0
        for i in range(coef.shape[0]):
            s += coef[i] * k[i, j]
        tmp[j] = y[j] + h * s
    rhs(t, tmp, alpha0, orb, out)


@njit(cache=True)
def dp_step(t, y, h, k1, alpha0, orb, ynew, k7):
    """One DP5(4) step with FSAL. k1 = f(t, y) on entry. Returns the scaled error vector in k7-side buffer (err)."""
    k = np.empty((7, NSTATE))
    tmp = np.empty(NSTATE)
    k[0] = k1
    _stage(t + C2 * h, y, h, k, np.array([A21]), alpha0, orb, tmp, k[1])
    _stage(t + C3 * h, y, h, k, np.array([A31, A32]), alpha0, orb, tmp, k[2])
    _stage(t + C4 * h, y, h, k, np.array([A41, A42, A43]), alpha0, orb, tmp, k[3])
    _stage(t + C5 * h, y, h, k, np.array([A51, A52, A53, A54]), alpha0, orb, tmp, k[4])
    _stage(t + h, y, h, k, np.array([A61, A62, A63, A64, A65]), alpha0, orb, tmp, k[5])
    for j in range(NSTATE):
        ynew[j] = y[j] + h * (B1 * k[0, j] + B3 * k[2, j] + B4 * k[3, j] + B5 * k[4, j] + B6 * k[5, j])
    rhs(t + h, ynew, alpha0, orb, k7)
    err = np.empty(NSTATE)
    for j in range(NSTATE):
        err[j] = h * (E1 * k[0, j] + E3 * k[2, j] + E4 * k[3, j] + E5 * k[4, j] + E6 * k[5, j] + E7 * k7[j])
    return err


@njit(cache=True)
def err_norm(err, y, ynew, rtol, atol):
    s = 0.0
    for j in range(NSTATE):
        sc = atol + rtol * max(abs(y[j]), abs(ynew[j]))
        s += (err[j] / sc) ** 2
    return math.sqrt(s / NSTATE)


@njit(cache=True)
def integrate_to(y, t0, t1, h, alpha0, orb, rtol, atol):
    """Advance y in place from t0 to t1 (landing exactly on t1). Returns the step size to carry on with, or -1 on failure."""
    t = t0
    k1 = np.empty(NSTATE)
    ynew = np.empty(NSTATE)
    k7 = np.empty(NSTATE)
    rhs(t, y, alpha0, orb, k1)
    for _ in range(MAX_STEPS):
        if t >= t1:
            return h
        clipped = h >= t1 - t
        h = t1 - t if clipped else h
        err = dp_step(t, y, h, k1, alpha0, orb, ynew, k7)
        en = err_norm(err, y, ynew, rtol, atol)
        if en <= 1.0:
            t = t1 if clipped else t + h
            y[:] = ynew
            k1[:] = k7
        fac = SAFETY * en ** (-ERR_EXPONENT) if en > 0.0 else FAC_MAX
        h *= min(FAC_MAX, max(FAC_MIN, fac))
        if h < H_MIN:
            return -1.0
    return -1.0


@njit(cache=True)
def phase_of(y, orb_dummy_node, cos_lb):
    """Precession phase phi = atan2(S.x', -S.y') mod 2 pi with x' = node direction, y' = z' x x'."""
    om = orb_dummy_node
    sin_lb = math.sqrt(max(0.0, 1.0 - cos_lb * cos_lb))
    nx, ny = math.cos(om), math.sin(om)
    lx, ly, lz = sin_lb * math.sin(om), -sin_lb * math.cos(om), cos_lb
    # y' = L x x' with x' = (nx, ny, 0)
    yx, yy, yz = ly * 0.0 - lz * ny, lz * nx - lx * 0.0, lx * ny - ly * nx
    sxp = y[0] * nx + y[1] * ny
    syp = y[0] * yx + y[1] * yy + y[2] * yz
    return math.atan2(sxp, -syp) % TWO_PI


@njit(cache=True)
def observe(y, t, orb):
    """Return (p, phi, |S|) at time t for state y."""
    x = x_of_t(t, orb)
    cos_lb = cos_theta_lb_of_x(x, orb)
    sin_lb = math.sqrt(max(0.0, 1.0 - cos_lb * cos_lb))
    om = y[3]
    lx, ly, lz = sin_lb * math.sin(om), -sin_lb * math.cos(om), cos_lb
    norm = math.sqrt(y[0] ** 2 + y[1] ** 2 + y[2] ** 2)
    return (y[0] * lx + y[1] * ly + y[2] * lz) / norm, phase_of(y, om, cos_lb), norm


def initial_state(p0, phi0, orb):
    """Inertial S at t = 0 (Omega = 0) from p0 = cos theta_sl and phi0, using S15's phi convention."""
    th = math.acos(math.sqrt(orb[2] / orb[0]))  # cos th_lb = sqrt(h/x0)
    lvec = np.array([0.0, -math.sin(th), math.cos(th)])
    xp = np.array([1.0, 0.0, 0.0])
    yp = np.cross(lvec, xp)
    sin_sl = math.sqrt(max(0.0, 1.0 - p0 * p0))
    s = sin_sl * (math.sin(phi0) * xp - math.cos(phi0) * yp) + p0 * lvec
    return np.array([s[0], s[1], s[2], 0.0])


@njit(cache=True)
def run_cycles(y0, n_cycles, period, alpha0, orb, rtol, atol):
    """Integrate n_cycles LK periods; row k holds (p, phi, |S|) at t = (k+1) P. Row is NaN after a failure."""
    y = y0.copy()
    out = np.full((n_cycles, NCOL), np.nan)
    h = H_INITIAL
    for k in range(n_cycles):
        h = integrate_to(y, k * period, (k + 1) * period, h, alpha0, orb, rtol, atol)
        if h < 0.0:
            return out
        o = observe(y, (k + 1) * period, orb)
        out[k, COL_P], out[k, COL_PHI], out[k, COL_NORM] = o[0], o[1], o[2]
    return out
