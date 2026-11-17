"""Analytic quadrupole Lidov-Kozai orbit, test-particle limit, no short-range forces (Storch & Lai 2015, Sec. 2.1).

Time is in units of t_k throughout. The orbit is fixed by (theta_lb0, e0); the spin problem then depends only on
epsilon = Omega_pl,0 / alpha_0, so alpha_0 = OMEGA_PL0 / epsilon (S15 Eq. epsilon).

    x(t)       = x0 + (x1 - x0) cn^2(theta, m),  theta = (K/pi)(n_e t + pi),  x = 1 - e^2    (S15 Eq. x)
    cos th_lb  = sqrt(h / x),  h = x0 cos^2 theta_lb0                                       (S15 Eq. kozai)
    Omega_pl   = Omega_pl0 [1 - 2 (x0 - h)/(x - h)] < 0,  Omega_pl0 = (3/4) sqrt(h)          (S15 Eq. pl)
    alpha      = alpha_0 / x^{3/2}                                                           (S15 Sec. 2.2)

t = 0 is an eccentricity minimum (x = x0); t = P/2 is the maximum (x = x1); P = 2 pi / n_e.
The numba functions take the orbit as a flat float64 array (ORBIT_* indices) so they can be called from jitted code.
"""

import math
from dataclasses import dataclass

import numpy as np
from numba import njit
from scipy.special import ellipk

# S15 numerical example (Sec. 2.1)
SL15_THETA_LB0_DEG = 85.0
SL15_E0 = 0.01

AGM_TOL = 1e-15
AGM_MAX_ITER = 40
QUADRATIC_A = 1.0  # leading coefficient of the x1, x2 quadratic

# Layout of the flat parameter array passed to jitted code.
ORBIT_X0, ORBIT_X1, ORBIT_H, ORBIT_M, ORBIT_K, ORBIT_NE, ORBIT_OPL0 = range(7)
ORBIT_SIZE = 7


@dataclass(frozen=True)
class LKOrbit:
    theta_lb0_deg: float
    e0: float
    x0: float
    h: float
    x1: float
    x2: float
    m: float  # elliptic parameter k^2
    big_k: float  # K(m)
    n_e: float  # eccentricity-oscillation frequency, 1 / t_k
    omega_pl0: float  # |Omega_pl| at x = x0, 1 / t_k

    @property
    def period(self) -> float:
        return 2.0 * math.pi / self.n_e

    @property
    def e_max(self) -> float:
        return math.sqrt(1.0 - self.x1)

    def alpha0(self, epsilon: float) -> float:
        return self.omega_pl0 / epsilon

    def params(self) -> np.ndarray:
        out = np.empty(ORBIT_SIZE)
        out[ORBIT_X0], out[ORBIT_X1], out[ORBIT_H] = self.x0, self.x1, self.h
        out[ORBIT_M], out[ORBIT_K], out[ORBIT_NE], out[ORBIT_OPL0] = self.m, self.big_k, self.n_e, self.omega_pl0
        return out


def make_orbit(theta_lb0_deg: float = SL15_THETA_LB0_DEG, e0: float = SL15_E0) -> LKOrbit:
    x0 = 1.0 - e0**2
    h = x0 * math.cos(math.radians(theta_lb0_deg)) ** 2
    b = (5.0 + 5.0 * h - 2.0 * x0) / 3.0
    c = 5.0 * h / 3.0
    disc = math.sqrt(b * b - 4.0 * QUADRATIC_A * c)
    x1, x2 = (b - disc) / 2.0, (b + disc) / 2.0
    m = (x0 - x1) / (x2 - x1)
    big_k = float(ellipk(m))
    n_e = 6.0 * math.pi * math.sqrt(6.0) / (8.0 * big_k) * math.sqrt(x2 - x1)
    return LKOrbit(theta_lb0_deg, e0, x0, h, x1, x2, m, big_k, n_e, 0.75 * math.sqrt(h))


@njit(cache=True)
def jacobi_cn(u: float, m: float) -> float:
    """Jacobi cn(u | m) by the arithmetic-geometric mean (Abramowitz & Stegun 16.4)."""
    a = np.empty(AGM_MAX_ITER + 1)
    c = np.empty(AGM_MAX_ITER + 1)
    a[0], b, c[0] = 1.0, math.sqrt(1.0 - m), math.sqrt(m)
    n = 0
    while abs(c[n]) > AGM_TOL and n < AGM_MAX_ITER:
        a[n + 1] = 0.5 * (a[n] + b)
        c[n + 1] = 0.5 * (a[n] - b)
        b = math.sqrt(a[n] * b)
        n += 1
    phi = (2.0**n) * a[n] * u
    for j in range(n, 0, -1):
        phi = 0.5 * (phi + math.asin(c[j] * math.sin(phi) / a[j]))
    return math.cos(phi)


@njit(cache=True)
def x_of_t(t: float, orb: np.ndarray) -> float:
    theta = orb[ORBIT_K] / math.pi * (orb[ORBIT_NE] * t + math.pi)
    cn = jacobi_cn(theta, orb[ORBIT_M])
    return orb[ORBIT_X0] + (orb[ORBIT_X1] - orb[ORBIT_X0]) * cn * cn


@njit(cache=True)
def omega_pl_of_x(x: float, orb: np.ndarray) -> float:
    h = orb[ORBIT_H]
    return orb[ORBIT_OPL0] * (1.0 - 2.0 * (orb[ORBIT_X0] - h) / (x - h))


@njit(cache=True)
def cos_theta_lb_of_x(x: float, orb: np.ndarray) -> float:
    return math.sqrt(orb[ORBIT_H] / x)


@njit(cache=True)
def alpha_of_x(x: float, alpha0: float) -> float:
    return alpha0 / x**1.5
