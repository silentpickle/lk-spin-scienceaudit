"""Shared helpers for the prep jobs (prep_sq2.py, prep_sq4.py, prep_sens4.py).

Nothing here changes lkspin. Pieces:
  * config tables ([step4] plus the prep_* / sens4 tables),
  * the S15 precession frame built from the current L vector (phase of the spin about L_hat),
  * `run_to_times`: integrate a state with lkspin.secular.advance, landing exactly on a list of times (e minima),
  * first-cycle orbit shape (e_min, e_max, period, N_max,0 x P, A_bar,0 x P, Omega_hat_L),
  * branch statistics of a theta_sl,f ensemble (largest-gap clusters as regime.cluster_count, bootstrap intervals).
"""

import math
import tomllib

import numpy as np

from lkspin import regime as R
from lkspin import secular as S

TWO_PI = 2.0 * math.pi
DEG = 180.0 / math.pi
MAX_STEPS = S.MAX_STEPS_DEFAULT
PERIOD_GUESS_FACTOR = 1.5      # first guess of the integration time, in units of (n_cycles + 1) analytic LK periods
GUESS_GROW = 2.0
MAX_GUESS_TRIES = 6
ORBIT_EV_MAX_PAD = 8           # event buffer rows per expected cycle (min + max, with slack)
ORBIT_PROBE_PERIOD_DAYS = 10.0 # the single spin carried while finding the minima times (does not affect the orbit)
SHAPE_T_FACTOR = 4.0           # [step5] shape_ncycle_tk: first-cycle integration length in analytic periods
MODE_MIN_COUNT = 1


def load_table(name: str) -> dict:
    with S.CONFIG_PATH.open("rb") as fh:
        return tomllib.load(fh)[name]


# ------------------------------------------------------------------ parameters
def nondissipative_params(cfg: dict, *, srf: bool = True) -> np.ndarray:
    """Tides off, tidal-disruption stop off (a is fixed, so a long run must not be stopped by the pericentre event)."""
    par = S.make_params(cfg, srf=srf, tides=False)
    par[S.P_RTIDE] = 0.0
    return par


# ------------------------------------------------------------------ S15 frame from the L vector
def phase_and_p(lvec: np.ndarray, svec: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(p, phi) of spin vectors svec[..., 3] about the orbital normal lvec[..., 3] (same leading shape), S15 convention
    as lkspin.spin.phase_of: x' = node direction (cos O, sin O, 0), O = atan2(L_x, -L_y); y' = L_hat x x';
    phi = atan2(S.x', -S.y') mod 2 pi, p = cos theta_sl."""
    lhat = lvec / np.linalg.norm(lvec, axis=-1, keepdims=True)
    shat = svec / np.linalg.norm(svec, axis=-1, keepdims=True)
    om = np.arctan2(lhat[..., 0], -lhat[..., 1])
    xp = np.stack([np.cos(om), np.sin(om), np.zeros_like(om)], axis=-1)
    yp = np.cross(lhat, xp)
    phi = np.arctan2(np.sum(shat * xp, axis=-1), -np.sum(shat * yp, axis=-1)) % TWO_PI
    return np.sum(shat * lhat, axis=-1), phi


def initial_spin_vector(theta_lb0_deg: float, p0: float, phi0: float) -> np.ndarray:
    """Unit spin at t = 0 for given p0 = cos theta_sl and S15 phase phi0 (as lkspin.spin.initial_state): L_hat0 =
    (0, -sin th, cos th), x' = (1, 0, 0), y' = L_hat0 x x'."""
    th = math.radians(theta_lb0_deg)
    lhat = np.array([0.0, -math.sin(th), math.cos(th)])
    xp = np.array([1.0, 0.0, 0.0])
    yp = np.cross(lhat, xp)
    sin_sl = math.sqrt(max(0.0, 1.0 - p0 * p0))
    return sin_sl * (math.sin(phi0) * xp - math.cos(phi0) * yp) + p0 * lhat


def tilted_spin_vector(lhat: np.ndarray, theta_sl_deg: float, azimuth_rad: float) -> np.ndarray:
    """Unit spin at angle theta_sl from lhat with azimuth about lhat measured from the x axis (x is perpendicular to
    L_hat0 = (0, -sin, cos))."""
    xp = np.array([1.0, 0.0, 0.0])
    yp = np.cross(lhat, xp)
    th = math.radians(theta_sl_deg)
    return math.cos(th) * lhat + math.sin(th) * (math.cos(azimuth_rad) * xp + math.sin(azimuth_rad) * yp)


def angle_between(u: np.ndarray, v: np.ndarray) -> np.ndarray:
    """Angle (rad) between vectors along the last axis, atan2 form (accurate for tiny angles)."""
    cross = np.linalg.norm(np.cross(u, v), axis=-1)
    return np.arctan2(cross, np.sum(u * v, axis=-1))


# ------------------------------------------------------------------ orbit minima and landing integration
def analytic_period_yr(par: np.ndarray, cfg: dict, theta_deg: float) -> float:
    return R.orbit_period_yr(par, cfg, theta_deg)


def orbit_minima_times(par: np.ndarray, cfg: dict, theta_deg: float, n_cycles: int) -> np.ndarray:
    """Times (yr) of the first n_cycles eccentricity minima after t = 0 (the orbit does not depend on the spins)."""
    t_stop = PERIOD_GUESS_FACTOR * (n_cycles + 1) * analytic_period_yr(par, cfg, theta_deg)
    y0 = S.initial_state(par, cfg, theta_deg, 1)
    for _ in range(MAX_GUESS_TRIES):
        ev_max = ORBIT_EV_MAX_PAD * (n_cycles + 1)
        res = S.run_diag(y0, par, [ORBIT_PROBE_PERIOD_DAYS], t_stop, ev_max=ev_max)
        minima, _ = R.split_events(res["ev"])
        if minima.shape[0] >= n_cycles:
            return minima[:n_cycles, 0].copy()
        t_stop *= GUESS_GROW
    raise RuntimeError("could not find enough eccentricity minima")


def run_to_times(y0: np.ndarray, par: np.ndarray, periods_days, times: np.ndarray, rtol: float, atol: float) -> tuple[np.ndarray, int]:
    """Integrate y0 (K spins) landing exactly on each time in `times` (ascending). Returns (states [n, len(y0)], status);
    rows after a stop event are NaN. status = ST_TIMEOUT means every time was reached."""
    om = np.array([S.omega_star_of_period(p) for p in np.atleast_1d(periods_days)])
    y = y0.copy()
    rec = np.zeros((0, S.N_REC_FIXED + om.size))
    out = np.full((times.size, y.size), np.nan)
    t, h, status = 0.0, S.H_INITIAL, S.ST_TIMEOUT
    for k, t_next in enumerate(times):
        t, h, status, _, _ = S.advance(y, t, t_next, h, par, om, rtol, atol, MAX_STEPS, rec, 1)
        if status != S.ST_TIMEOUT:
            break
        out[k] = y
    return out, status


def spin_block(states: np.ndarray, n_spin: int) -> np.ndarray:
    """Spin vectors [n_times, n_spin, 3] from a state array [n_times, len(y)]."""
    sl = states[:, S.N_ORBIT:S.N_ORBIT + S.N_SPIN * n_spin]
    return sl.reshape(states.shape[0], n_spin, S.N_SPIN)


def theta_sl_deg_series(states: np.ndarray, n_spin: int) -> np.ndarray:
    """theta_sl (deg) of every spin at every stored time: [n_times, n_spin]."""
    lvec = states[:, 0:3][:, None, :]
    return np.degrees(angle_between(np.broadcast_to(lvec, (states.shape[0], n_spin, S.N_SPIN)), spin_block(states, n_spin)))


# ------------------------------------------------------------------ first-cycle orbit shape
def first_cycle_shape(par: np.ndarray, cfg: dict, theta_deg: float, p_star_days: float = R.P_REFERENCE_DAYS) -> dict:
    """First-cycle (t = 0 to the first interior e minimum) orbit shape and S16 coefficients, as step4_validate check 5.

    Returns e_min, e_max, 1 - e_max, period_yr, nmax_coef (= N_max,0 x P_star[d]), abar_coef, omega_hat_l (= N/(2A) =
    <Omega_pl sin theta_lb> / n_e, independent of P_star), n_max0 and a_bar0 at p_star_days."""
    t_first = SHAPE_T_FACTOR * analytic_period_yr(par, cfg, theta_deg)
    rec = R.first_cycle_record(par, cfg, theta_deg, t_first)
    sl = R.first_cycle_slice(rec)
    if sl is None:
        return {"nan": True}
    shape = R.shape_coefficients(par, cfg, theta_deg, t_first)
    t, e = rec[sl, 0], rec[sl, 2]
    period = float(t[-1] - t[0])
    nmax, abar = shape["nmax_coef"], shape["abar_coef"]
    scale = R.P_REFERENCE_DAYS / p_star_days
    return {"e_min": float(e.min()), "e_max": float(e.max()), "one_minus_e_max": float(1.0 - e.max()),
            "period_yr": period, "nmax_coef": nmax, "abar_coef": abar, "omega_hat_l": nmax / (2.0 * abar),
            "n_max0": nmax * scale, "a_bar0": abar * scale}


# ------------------------------------------------------------------ branch statistics of a theta_f ensemble
def split_value(theta: np.ndarray) -> tuple[float, float]:
    """(largest gap, midpoint of it) of the sorted values."""
    srt = np.sort(theta)
    gaps = np.diff(srt)
    g = int(np.argmax(gaps))
    return float(gaps[g]), float(0.5 * (srt[g] + srt[g + 1]))


def branch_summary(theta: np.ndarray, tight_deg: float, gap_deg: float) -> dict:
    """Two-branch description by the largest gap (regime.cluster_count rules). n_branches = 1, 2, or 0 (not two tight
    branches). cut = midpoint of the largest gap; p_high = fraction above it (0 or 1 if single); peaks = medians."""
    n_clusters = R.cluster_count(theta, tight_deg, gap_deg)
    gap, cut = split_value(theta)
    high = theta > cut
    low_med = float(np.median(theta[~high])) if (~high).any() else math.nan
    high_med = float(np.median(theta[high])) if high.any() else math.nan
    return {"n_branches": int(n_clusters), "gap_deg": gap, "cut_deg": cut, "p_high": float(high.mean()),
            "peak_low_deg": low_med, "peak_high_deg": high_med, "median_all_deg": float(np.median(theta)),
            "range_deg": float(theta.max() - theta.min())}


def modes(theta: np.ndarray, mode_gap: float, min_frac: float) -> list[dict]:
    """Split sorted theta at every gap >= mode_gap; keep groups holding >= min_frac of the members.
    Each: n, fraction, median, lo, hi."""
    srt = np.sort(theta)
    cuts = np.where(np.diff(srt) >= mode_gap)[0] + 1
    groups = np.split(srt, cuts)
    return [{"n": int(g.size), "fraction": float(g.size / srt.size), "median": float(np.median(g)), "lo": float(g[0]),
             "hi": float(g[-1])} for g in groups if g.size >= max(MODE_MIN_COUNT, min_frac * srt.size)]


def bootstrap_branches(theta: np.ndarray, cut: float, n_boot: int, rng: np.random.Generator, ci: list) -> dict:
    """Bootstrap intervals of p_high and of the two peak medians with the cut fixed at the full-sample value."""
    idx = rng.integers(0, theta.size, size=(n_boot, theta.size))
    samp = theta[idx]
    high = samp > cut
    p_high = high.mean(axis=1)
    out = {"p_high_ci": np.percentile(p_high, ci).tolist()}
    for name, mask in (("peak_low", ~high), ("peak_high", high)):
        med = np.array([np.median(s[m]) if m.any() else np.nan for s, m in zip(samp, mask, strict=True)])
        out[f"{name}_ci"] = np.nanpercentile(med, ci).tolist() if np.isfinite(med).any() else [math.nan, math.nan]
    return out


# ------------------------------------------------------------------ FTLE from separations at e minima
def ftle_fit(sep: np.ndarray, delta0: float, sat: float, min_points: int) -> dict:
    """Slope of ln(separation) against cycle number k = 0, 1, ... from the start up to where the separation first reaches
    `sat` (exclusive; at least min_points points). sep[k-1] is the separation after k cycles, delta0 the one at k = 0.
    Returns ftle (per cycle), k_sat (first cycle with sep >= sat, or -1) and n_fit."""
    ln_all = np.log(np.concatenate([[delta0], np.maximum(sep, np.finfo(float).tiny)]))
    over = np.where(np.exp(ln_all) >= sat)[0]
    k_sat = int(over[0]) if over.size else -1
    n_fit = ln_all.size if k_sat < 0 else max(k_sat, min_points)
    k = np.arange(n_fit)
    slope = float(np.polyfit(k, ln_all[:n_fit], 1)[0])
    return {"ftle": slope, "k_sat": k_sat, "n_fit": int(n_fit)}
