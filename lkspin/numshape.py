"""Numerical "hat" shape functions of one Lidov-Kozai cycle of a `secular.py` orbit (with or without short-range forces).

The hat functions equal N_max * eps * (beta, gamma, psi) of Storch & Lai 2015 / Storch et al. 2016 and do not depend on
the stellar spin period P_star; N_max = nmax_coef / P_star[days], A_bar = abar_coef / P_star[days].

    alpha_1(t) = par[P_STAR_COEF] * omega_star(1 day) / (a^3 j^3)          (alpha = alpha_1 / P_star[days])
    tau(t)     = 2 pi int_{t_a}^t alpha_1 dt' / Phi_1,   Phi_1 = int over the cycle of alpha_1
    beta_hat   = -(Omega_pl sin th_lb) dt/dtau,  gamma_hat = +(th_lb_dot) dt/dtau,  psi_hat = -(Omega_pl cos th_lb) dt/dtau

Two passes: (i) land on a uniform t grid and build tau(t) with a cumulative Simpson rule; (ii) invert tau(t) (cubic Hermite
with the known dt/dtau), land on the times of a uniform tau grid and evaluate everything there from the state and `rhs`.

Hamiltonian for the spin (phi coordinate, p momentum), `hamiltonian_rhs`:
    H = n_max (-p^2/2) + psi_hat p - sqrt(1 - p^2) [beta_hat cos phi + gamma_hat sin phi]
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_simpson
from scipy.interpolate import CubicHermiteSpline, CubicSpline

from lkspin import secular as S
from lkspin import shape as SH
from lkspin.lk_orbit import make_orbit

TWO_PI = 2.0 * np.pi
P_STAR_REF_DAYS = 1.0       # alpha_1, Phi_1 and the hat functions are for this spin period
M_MAX_DEFAULT = 10
PERIOD_GUESS_FACTOR = 1.5   # first end time for the minimum search, in analytic LK periods
GUESS_GROW = 2.0
MAX_GUESS_TRIES = 6
EV_MAX = 256                # event rows kept while searching for the first minimum
SPIN_PROBE_DAYS = 10.0      # the single spin carried by the orbit integrations (it does not affect the orbit)
RTOL_DEFAULT = 1e-11
ATOL_DEFAULT = 1e-13


@dataclass(frozen=True)
class NumShapes:
    tau: np.ndarray
    beta_hat: np.ndarray
    gamma_hat: np.ndarray
    psi_hat: np.ndarray
    nmax_coef: float
    omega_hat_l: float
    psi_hat_0: float
    abar_coef: float
    period_yr: float
    e_min: float
    e_max: float
    beta_cos: np.ndarray
    beta_sin: np.ndarray
    gamma_sin: np.ndarray
    gamma_cos: np.ndarray
    psi_cos: np.ndarray
    psi_sin: np.ndarray


@dataclass(frozen=True)
class ShapeInterp:
    """Periodic cubic spline of (beta_hat, gamma_hat, psi_hat) in tau; call spline(tau) -> (3,) or (3, n)."""

    spline: CubicSpline


# ------------------------------------------------------------------ states
def probe_state(y_a: np.ndarray) -> np.ndarray:
    """Orbit part of y_a plus one aligned spin: the length `rhs` expects for one stellar spin."""
    lvec = y_a[0:3]
    spin = lvec / np.linalg.norm(lvec)
    return np.concatenate([y_a[:S.N_ORBIT], spin])


def land(y0: np.ndarray, t_a: float, times: np.ndarray, par: np.ndarray, rtol: float, atol: float) -> np.ndarray:
    """States at each of `times` (ascending, first >= t_a), integrating `y0` from t_a and landing exactly on each time."""
    om = np.array([S.omega_star_of_period(SPIN_PROBE_DAYS)])
    y = y0.copy()
    rec = np.zeros((0, S.N_REC_FIXED + om.size))
    out = np.empty((times.size, y.size))
    t, h = t_a, S.H_INITIAL
    for k, t_next in enumerate(times):
        t, h, status, _, _ = S.advance(y, t, max(t_next, t), h, par, om, rtol, atol, S.MAX_STEPS_DEFAULT, rec, 1)
        if status != S.ST_TIMEOUT:
            raise RuntimeError(f"orbit integration stopped with status {S.STATUS_NAMES[status]}")
        out[k] = y
    return out


def alpha_unit(states: np.ndarray, par: np.ndarray) -> np.ndarray:
    """alpha_1 = par[P_STAR_COEF] * omega_star(1 day) / (a^3 j^3) for each state."""
    lmag2 = np.sum(states[:, 0:3] ** 2, axis=1)
    j2 = 1.0 - np.sum(states[:, 3:6] ** 2, axis=1)
    a = par[S.P_A0] * lmag2 / j2
    return par[S.P_STAR_COEF] * S.omega_star_of_period(P_STAR_REF_DAYS) / (a**3 * j2**1.5)


def eccentricity(states: np.ndarray) -> np.ndarray:
    return np.linalg.norm(states[:, 3:6], axis=1)


# ------------------------------------------------------------------ first minimum
def first_minimum_time(par: np.ndarray, cfg: dict, theta_lb0_deg: float, rtol: float, atol: float) -> float:
    """Time (yr) of the first eccentricity minimum after t = 0 (event search of secular.run_diag)."""
    analytic = make_orbit(theta_lb0_deg, cfg["e0"]).period * S.t_k_years(par, cfg["a0_au"])
    t_stop = PERIOD_GUESS_FACTOR * analytic
    y0 = S.initial_state(par, cfg, theta_lb0_deg, 1)
    for _ in range(MAX_GUESS_TRIES):
        ev = S.run_diag(y0, par, [SPIN_PROBE_DAYS], t_stop, rtol, atol, ev_max=EV_MAX)["ev"]
        minima = ev[ev[:, 3] == S.KIND_MIN]
        if minima.shape[0] > 0:
            return float(minima[0, 0])
        t_stop *= GUESS_GROW
    raise RuntimeError("no eccentricity minimum found")


# ------------------------------------------------------------------ pass 1: tau(t)
@dataclass(frozen=True)
class TauMap:
    t: np.ndarray
    tau: np.ndarray
    dt_dtau: np.ndarray
    phi_1: float
    e_min: float
    e_max: float


def build_tau_map(y_a: np.ndarray, t_a: float, t_b: float, par: np.ndarray, n_t: int, rtol: float, atol: float) -> TauMap:
    """Uniform t grid of n_t intervals over [t_a, t_b]; tau(t) from the cumulative Simpson integral of alpha_1."""
    t = np.linspace(t_a, t_b, n_t + 1)
    states = land(probe_state(y_a), t_a, t, par, rtol, atol)
    alpha = alpha_unit(states, par)
    cum = cumulative_simpson(alpha, x=t, initial=0.0)
    phi_1 = float(cum[-1])
    ecc = eccentricity(states)
    return TauMap(t, TWO_PI * cum / phi_1, phi_1 / (TWO_PI * alpha), phi_1, float(ecc.min()), float(ecc.max()))


def times_of_tau(tmap: TauMap, tau_grid: np.ndarray) -> np.ndarray:
    """Times at the tau grid by cubic Hermite inversion (tau nodes, t values, dt/dtau slopes)."""
    inv = CubicHermiteSpline(tmap.tau, tmap.t, tmap.dt_dtau)
    return np.maximum.accumulate(inv(tau_grid))


# ------------------------------------------------------------------ pass 2: hat functions from states
def state_rates(states: np.ndarray, par: np.ndarray) -> np.ndarray:
    """dl/dt (rows of out[0:3] from secular.rhs) at each state."""
    om = np.array([S.omega_star_of_period(SPIN_PROBE_DAYS)])
    out = np.zeros(states.shape[1])
    rates = np.empty((states.shape[0], 3))
    for k, y in enumerate(states):
        S.rhs(0.0, y, par, om, out)
        rates[k] = out[0:3]
    return rates


def hat_functions(states: np.ndarray, rates: np.ndarray, dt_dtau: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """(beta_hat, gamma_hat, psi_hat) from states, dl/dt and dt/dtau (definitions 4 to 7 of the module docstring)."""
    lmag = np.linalg.norm(states[:, 0:3], axis=1)
    lhat = states[:, 0:3] / lmag[:, None]
    dlhat = (rates - lhat * np.sum(lhat * rates, axis=1)[:, None]) / lmag[:, None]
    node = np.arctan2(lhat[:, 0], -lhat[:, 1])
    xp = np.stack([np.cos(node), np.sin(node), np.zeros_like(node)], axis=1)
    yp = np.cross(lhat, xp)
    sin_lb = np.hypot(lhat[:, 0], lhat[:, 1])
    opl_sin = np.sum(dlhat * xp, axis=1)
    opl_cos = opl_sin * lhat[:, 2] / sin_lb
    thdot = -np.sum(dlhat * yp, axis=1)
    return -opl_sin * dt_dtau, thdot * dt_dtau, -opl_cos * dt_dtau


def cycle_shapes(y_a: np.ndarray, t_a: float, t_b: float, par: np.ndarray, n_t: int, n_tau: int,
                 m_max: int = M_MAX_DEFAULT, rtol: float = RTOL_DEFAULT, atol: float = ATOL_DEFAULT) -> NumShapes:
    """Hat functions of the cycle from the minimum at (t_a, y_a) to the next minimum at t_b."""
    tmap = build_tau_map(y_a, t_a, t_b, par, n_t, rtol, atol)
    tau = TWO_PI * np.arange(n_tau) / n_tau
    states = land(probe_state(y_a), t_a, times_of_tau(tmap, tau), par, rtol, atol)
    dt_dtau = tmap.phi_1 / (TWO_PI * alpha_unit(states, par))
    beta, gamma, psi = hat_functions(states, state_rates(states, par), dt_dtau)
    ecc = eccentricity(states)
    return assemble(tau, beta, gamma, psi, tmap, t_b - t_a, (min(tmap.e_min, ecc.min()), max(tmap.e_max, ecc.max())), m_max)


def assemble(tau: np.ndarray, beta: np.ndarray, gamma: np.ndarray, psi: np.ndarray, tmap: TauMap, period: float,
             e_range: tuple[float, float], m_max: int) -> NumShapes:
    beta_cos = SH.cos_coefficients(beta, m_max)
    psi_cos = SH.cos_coefficients(psi, m_max)
    omega_hat = float(np.mean(beta))
    nmax_coef = tmap.phi_1 / TWO_PI
    return NumShapes(tau, beta, gamma, psi, nmax_coef, omega_hat, float(psi_cos[0]), nmax_coef / (2.0 * omega_hat), period,
                     e_range[0], e_range[1], beta_cos, SH.sin_coefficients(beta, m_max), SH.sin_coefficients(gamma, m_max),
                     SH.cos_coefficients(gamma, m_max), psi_cos, SH.sin_coefficients(psi, m_max))


def first_cycle_shapes(par: np.ndarray, cfg: dict, theta_lb0_deg: float, n_t: int, n_tau: int,
                       m_max: int = M_MAX_DEFAULT, rtol: float = RTOL_DEFAULT, atol: float = ATOL_DEFAULT) -> NumShapes:
    """Shapes of the first cycle: t_a = 0 (initial state of `cfg`), t_b = the first eccentricity minimum."""
    y_a = S.initial_state(par, cfg, theta_lb0_deg, 1)
    t_b = first_minimum_time(par, cfg, theta_lb0_deg, rtol, atol)
    return cycle_shapes(y_a, 0.0, t_b, par, n_t, n_tau, m_max, rtol, atol)


# ------------------------------------------------------------------ Hamiltonian
def shape_interp(shapes: NumShapes) -> ShapeInterp:
    """Periodic cubic spline through the hat functions (the tau grid excludes 2 pi; the first point closes the loop)."""
    tau = np.append(shapes.tau, TWO_PI)
    values = np.stack([shapes.beta_hat, shapes.gamma_hat, shapes.psi_hat], axis=1)
    closed = np.vstack([values, values[:1]])
    return ShapeInterp(CubicSpline(tau, closed, bc_type="periodic"))


def hamiltonian_rhs(tau: float, state: np.ndarray, n_max: float, shapes_interp: ShapeInterp) -> list[float]:
    """d(phi, p)/dtau for H = n_max (-p^2/2) + psi_hat p - sqrt(1-p^2) (beta_hat cos phi + gamma_hat sin phi)."""
    phi, p = state
    beta, gamma, psi = shapes_interp.spline(tau % TWO_PI)
    root = np.sqrt(1.0 - p * p)
    drive = beta * np.cos(phi) + gamma * np.sin(phi)
    dphi = -n_max * p + psi + p / root * drive
    dp = root * (-beta * np.sin(phi) + gamma * np.cos(phi))
    return [dphi, dp]
