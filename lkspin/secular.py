"""Secular vector integrator for S16 Fig. 8: LK quadrupole + SRFs + weak-friction planet tides + K stellar spins.

Equations: A16 appendix (arXiv:1510.08918, ms_arxiv.tex l.1758-1930), vector form, e_b = 0, no stellar feedback.
Units: yr, AU, Msun (G = 4 pi^2). State y = [l(3), e(3), s, S_1(3), ..., S_K(3)] where l = L / L_SCALE,
s = S_p / L_SCALE (planet spin angular momentum), S_k = unit stellar spin for the k-th P_star.
a follows from |l| and |e|:  a = a0 |l|^2 / (1 - e^2).  The orbit does not depend on P_star, so K spins share one orbit.

Adaptive Dormand-Prince 5(4) (numba) with events: a < stop_a, pericentre < r_Tide, t >= t_stop.
Checked by tests/test_secular.py and scripts/step4_validate.py.
"""

import math
import tomllib
from pathlib import Path

import numpy as np
from numba import njit

from lkspin.spin import (A21, A31, A32, A41, A42, A43, A51, A52, A53, A54, A61, A62, A63, A64, A65, B1, B3, B4, B5,
                         B6, C2, C3, C4, C5, E1, E3, E4, E5, E6, E7)

ROOT = Path(__file__).resolve().parent.parent
CONFIG_PATH = ROOT / "config.toml"

# Unit conversions (yr, AU, Msun).
G_AU3_PER_MSUN_YR2 = 4.0 * math.pi**2
C_AU_PER_YR = 63241.077
SEC_PER_YR = 3.15576e7
MJUP_PER_MSUN = 9.5458e-4
RJUP_AU = 4.7791e-4
RSUN_AU = 4.6505e-3
HOURS_PER_YR = 8766.0
DAYS_PER_YR = 365.25
R_TIDE_COEF = 2.7  # A16 Eq. Rtide with f = 1
DELTA_T_L_SEC_PER_CHI = 0.1  # A16: Delta t_L = 0.1 chi s

# Layout of the parameter array.
(P_G, P_MTOT, P_MSTAR, P_MP, P_MB, P_AB, P_RP, P_RSTAR, P_K2P, P_KQP, P_KP, P_KQS, P_KS, P_DTL, P_LSCALE, P_A0,
 P_C, P_IP, P_OMBRK, P_SRF, P_TIDE, P_LK, P_STOPA, P_RTIDE, P_STAR_COEF, P_BULGE) = range(26)
NPAR = 26
BULGE_COEF_A16_EOM = 7.5  # A16 Eq. omega_tide_p; the potential (Eq. epsitide) implies 15, see findings

N_ORBIT = 7  # l(3), e(3), s
IDX_S = 6
N_SPIN = 3

# Status codes.
ST_RUNNING, ST_HOT_JUPITER, ST_TIMEOUT, ST_DISRUPTED, ST_FAILED, ST_STEPS = range(6)
STATUS_NAMES = ("running", "hot_jupiter", "timeout", "disrupted", "failed", "step_limit")

# Integrator constants.
H_INITIAL = 1e-3
H_MIN = 1e-10
SAFETY = 0.9
FAC_MIN = 0.2
FAC_MAX = 5.0
ERR_EXPONENT = 0.2
RTOL_DEFAULT = 1e-10
ATOL_DEFAULT = 1e-12
MAX_STEPS_DEFAULT = 2_000_000_000
NSTAGE = 7
N_REC_FIXED = 6  # leading columns of the record array

# Planet-tide exponents and coefficients of A16 f_1..f_5.
J_POW_TIDE = 13
HALF = 0.5
COEF_TIDE_E = 5.5
COEF_TIDE_E_F3 = 9.0


def load_config() -> dict:
    with CONFIG_PATH.open("rb") as fh:
        return tomllib.load(fh)["step4"]


def make_params(cfg: dict, *, srf: bool = True, tides: bool = True, lk: bool = True, a0_au: float | None = None,
                bulge_coef: float = BULGE_COEF_A16_EOM) -> np.ndarray:
    """Parameter array (see P_* indices) from config [step4]."""
    par = np.zeros(NPAR)
    ms, mp = cfg["m_star_msun"], cfg["m_p_mjup"] * MJUP_PER_MSUN
    rp = cfg["r_p_rjup"] * RJUP_AU
    a0 = cfg["a0_au"] if a0_au is None else a0_au
    mtot = ms + mp
    mu = ms * mp / mtot
    par[P_G], par[P_MTOT], par[P_MSTAR], par[P_MP], par[P_MB] = G_AU3_PER_MSUN_YR2, mtot, ms, mp, cfg["m_b_msun"]
    par[P_AB], par[P_RP], par[P_RSTAR] = cfg["a_b_au"], rp, cfg["r_star_rsun"] * RSUN_AU
    par[P_K2P], par[P_KQP], par[P_KP] = cfg["k_2p"], cfg["k_qp"], cfg["k_p"]
    par[P_KQS], par[P_KS] = cfg["k_q_star"], cfg["k_star"]
    par[P_DTL] = DELTA_T_L_SEC_PER_CHI * cfg["chi"] / SEC_PER_YR
    par[P_LSCALE] = mu * math.sqrt(G_AU3_PER_MSUN_YR2 * mtot * a0)
    par[P_A0], par[P_C] = a0, C_AU_PER_YR
    par[P_IP] = cfg["k_p"] * mp * rp**2
    par[P_OMBRK] = math.sqrt(G_AU3_PER_MSUN_YR2 * mp / rp**3)
    par[P_SRF], par[P_TIDE], par[P_LK] = float(srf), float(tides), float(lk)
    par[P_STOPA] = cfg["stop_a_au"]
    par[P_RTIDE] = R_TIDE_COEF * rp * (ms / mp) ** (1.0 / 3.0)
    par[P_BULGE] = bulge_coef
    par[P_STAR_COEF] = 1.5 * cfg["k_q_star"] / cfg["k_star"] * mp / ms * par[P_RSTAR] ** 3
    return par


def omega_star_of_period(p_star_days: float) -> float:
    """Stellar spin frequency in rad/yr."""
    return 2.0 * math.pi * DAYS_PER_YR / p_star_days


def initial_state(par: np.ndarray, cfg: dict, theta_lb0_deg: float, n_spin: int) -> np.ndarray:
    """L_hat0 = (0, -sin th, cos th), e0 along +x, S_hat = L_hat0 (aligned), planet spin from config."""
    th = math.radians(theta_lb0_deg)
    e0 = cfg["e0"]
    lhat = np.array([0.0, -math.sin(th), math.cos(th)])
    y = np.zeros(N_ORBIT + N_SPIN * n_spin)
    y[0:3] = math.sqrt(1.0 - e0**2) * lhat
    y[3] = e0
    omega_p = 2.0 * math.pi * HOURS_PER_YR / cfg["p_planet0_hours"]
    y[IDX_S] = par[P_IP] * omega_p / par[P_LSCALE]
    for k in range(n_spin):
        y[N_ORBIT + N_SPIN * k:N_ORBIT + N_SPIN * (k + 1)] = lhat
    return y


# ---------------------------------------------------------------- eccentricity functions (A16 Eqs. F2-F5)
@njit(cache=True)
def f2(e2):
    return 1.0 + 7.5 * e2 + 45.0 / 8.0 * e2**2 + 5.0 / 16.0 * e2**3


@njit(cache=True)
def f3(e2):
    return 1.0 + 15.0 / 4.0 * e2 + 15.0 / 8.0 * e2**2 + 5.0 / 64.0 * e2**3


@njit(cache=True)
def f4(e2):
    return 1.0 + 1.5 * e2 + e2**2 / 8.0


@njit(cache=True)
def f5(e2):
    return 1.0 + 3.0 * e2 + 3.0 / 8.0 * e2**2


@njit(cache=True)
def cross3(a0, a1, a2, b0, b1, b2):
    return a1 * b2 - a2 * b1, a2 * b0 - a0 * b2, a0 * b1 - a1 * b0


@njit(cache=True)
def orbit_elements(y, par):
    """Return (|l|, j, a) from the state."""
    lmag = math.sqrt(y[0] ** 2 + y[1] ** 2 + y[2] ** 2)
    e2 = y[3] ** 2 + y[4] ** 2 + y[5] ** 2
    j2 = 1.0 - e2
    return lmag, math.sqrt(j2), par[P_A0] * lmag * lmag / j2


@njit(cache=True)
def lk_terms(y, par, lmag, j, tk_inv, out):
    """LK quadrupole (A16 Eqs. Lvec, evec), L_b = z-hat; adds to out[0:6]."""
    lh0, lh1, lh2 = y[0] / lmag, y[1] / lmag, y[2] / lmag
    j0, j1, j2 = j * lh0, j * lh1, j * lh2
    e0, e1, e2 = y[3], y[4], y[5]
    jz, ez = j2, e2
    # j x Lb, e x Lb, e x j with Lb = z
    jxb = (j1, -j0, 0.0)
    exb = (e1, -e0, 0.0)
    jxe = cross3(j0, j1, j2, e0, e1, e2)
    fl = 0.75 * lmag / j * tk_inv
    out[0] += fl * (jz * jxb[0] - 5.0 * ez * exb[0])
    out[1] += fl * (jz * jxb[1] - 5.0 * ez * exb[1])
    out[2] += fl * (jz * jxb[2] - 5.0 * ez * exb[2])
    fe = 0.75 * tk_inv
    out[3] += fe * (jz * exb[0] + 2.0 * jxe[0] - 5.0 * ez * jxb[0])
    out[4] += fe * (jz * exb[1] + 2.0 * jxe[1] - 5.0 * ez * jxb[1])
    out[5] += fe * (jz * exb[2] + 2.0 * jxe[2] - 5.0 * ez * jxb[2])


@njit(cache=True)
def srf_rate(par, a, e2, j, n, omega_p):
    """Total pericentre precession rate (GR + planet static tide + planet rotation), A16 Sec. A.3."""
    rpa = par[P_RP] / a
    w_gr = 3.0 * par[P_G] * par[P_MTOT] / (par[P_C] ** 2 * a * j * j) * n
    w_tide = par[P_BULGE] * par[P_K2P] * par[P_MSTAR] / par[P_MP] * rpa**5 * f4(e2) / j**10 * n
    om_hat = omega_p / par[P_OMBRK]
    w_rot = 1.5 * par[P_KQP] * rpa**2 * om_hat**2 / (j * j) ** 2 * n
    return w_gr + w_tide + w_rot


@njit(cache=True)
def tide_terms(y, par, lmag, j, a, n, omega_p, out):
    """Weak-friction planet tide (A16 Sec. A.4): ds/dt, dl/dt = -ds/dt l_hat, de/dt (planet spin along L)."""
    e2 = y[3] ** 2 + y[4] ** 2 + y[5] ** 2
    ta_inv = 6.0 * par[P_K2P] * par[P_DTL] * par[P_MSTAR] / par[P_MP] * (par[P_RP] / a) ** 5 * n * n
    pref = HALF * ta_inv / j**J_POW_TIDE
    ratio = omega_p / n
    j3 = j**3
    ds = -pref * lmag * (j3 * f5(e2) * ratio - f2(e2))  # in units of L (|l| is L / L_SCALE)
    out[IDX_S] += ds
    for i in range(3):
        out[i] += -ds * y[i] / lmag
        out[3 + i] += pref * (COEF_TIDE_E * j3 * f4(e2) * ratio - COEF_TIDE_E_F3 * f3(e2)) * y[3 + i]


@njit(cache=True)
def rhs(t, y, par, omega_star, out):
    """dy/dt for orbit + planet spin + K stellar spins."""
    out[:] = 0.0
    lmag, j, a = orbit_elements(y, par)
    e2 = 1.0 - j * j
    n = math.sqrt(par[P_G] * par[P_MTOT] / a**3)
    omega_p = y[IDX_S] * par[P_LSCALE] / par[P_IP]
    if par[P_LK] > 0.0:
        tk_inv = par[P_MB] / par[P_MTOT] * (a / par[P_AB]) ** 3 * n
        lk_terms(y, par, lmag, j, tk_inv, out)
    if par[P_SRF] > 0.0:
        w = srf_rate(par, a, e2, j, n, omega_p)
        cx, cy, cz = cross3(y[0] / lmag, y[1] / lmag, y[2] / lmag, y[3], y[4], y[5])
        out[3] += w * cx
        out[4] += w * cy
        out[5] += w * cz
    if par[P_TIDE] > 0.0:
        tide_terms(y, par, lmag, j, a, n, omega_p, out)
    spin_terms(y, par, omega_star, lmag, j, a, out)


@njit(cache=True)
def spin_terms(y, par, omega_star, lmag, j, a, out):
    """dS_k/dt = Omega_ps L_hat x S_k with Omega_ps = -(3/2)(kq/k*)(Mp/M*)(R*/a)^3 Omega* cos(th_sl) / j^3."""
    lh0, lh1, lh2 = y[0] / lmag, y[1] / lmag, y[2] / lmag
    base = -par[P_STAR_COEF] / (a**3 * j**3)
    for k in range(omega_star.shape[0]):
        o = N_ORBIT + N_SPIN * k
        cs = lh0 * y[o] + lh1 * y[o + 1] + lh2 * y[o + 2]
        w = base * omega_star[k] * cs
        c0, c1, c2 = cross3(lh0, lh1, lh2, y[o], y[o + 1], y[o + 2])
        out[o] = w * c0
        out[o + 1] = w * c1
        out[o + 2] = w * c2


@njit(cache=True)
def _stage_sum(kk, coef, i):
    s = 0.0
    for m in range(coef.shape[0]):
        s += coef[m] * kk[m, i]
    return s


@njit(cache=True)
def rk_step(t, y, h, kk, par, omega_star, tmp, ynew, err):
    """One DP5(4) step with FSAL. kk[0] = f(t, y) on entry; fills ynew, err and kk[6] = f(t + h, ynew)."""
    n = y.shape[0]
    a2 = np.array([A21])
    a3 = np.array([A31, A32])
    a4 = np.array([A41, A42, A43])
    a5 = np.array([A51, A52, A53, A54])
    a6 = np.array([A61, A62, A63, A64, A65])
    stages = ((a2, C2), (a3, C3), (a4, C4), (a5, C5), (a6, 1.0))
    for st in range(5):
        coef = stages[st][0]
        for i in range(n):
            tmp[i] = y[i] + h * _stage_sum(kk, coef, i)
        rhs(t + stages[st][1] * h, tmp, par, omega_star, kk[st + 1])
    for i in range(n):
        ynew[i] = y[i] + h * (B1 * kk[0, i] + B3 * kk[2, i] + B4 * kk[3, i] + B5 * kk[4, i] + B6 * kk[5, i])
    rhs(t + h, ynew, par, omega_star, kk[6])
    for i in range(n):
        err[i] = h * (E1 * kk[0, i] + E3 * kk[2, i] + E4 * kk[3, i] + E5 * kk[4, i] + E6 * kk[5, i] + E7 * kk[6, i])


@njit(cache=True)
def err_norm(err, y, ynew, rtol, atol):
    s = 0.0
    for i in range(y.shape[0]):
        sc = atol + rtol * max(abs(y[i]), abs(ynew[i]))
        s += (err[i] / sc) ** 2
    return math.sqrt(s / y.shape[0])


@njit(cache=True)
def event_status(y, par):
    """Stopping events: hot Jupiter (a < stop_a), tidal disruption (pericentre < r_Tide), else running."""
    lmag, j, a = orbit_elements(y, par)
    if a < par[P_STOPA]:
        return ST_HOT_JUPITER
    if a * j * j / (1.0 + j) < par[P_RTIDE] and par[P_RTIDE] > 0.0:  # a(1-e) = a j^2 / (1 + j)
        return ST_DISRUPTED
    return ST_RUNNING


@njit(cache=True)
def record_row(rec, idx, t, y, par, omega_star):
    """rec columns: t, a, e, theta_lb, node angle, planet spin (rad/yr), then theta_sl for each spin."""
    lmag, j, a = orbit_elements(y, par)
    rec[idx, 0] = t
    rec[idx, 1] = a
    rec[idx, 2] = math.sqrt(max(0.0, 1.0 - j * j))
    rec[idx, 3] = math.acos(min(1.0, max(-1.0, y[2] / lmag)))
    rec[idx, 4] = math.atan2(y[0], -y[1])
    rec[idx, 5] = y[IDX_S] * par[P_LSCALE] / par[P_IP]
    for k in range(omega_star.shape[0]):
        o = N_ORBIT + N_SPIN * k
        cs = (y[0] * y[o] + y[1] * y[o + 1] + y[2] * y[o + 2]) / lmag
        rec[idx, N_REC_FIXED + k] = math.acos(min(1.0, max(-1.0, cs)))


@njit(cache=True)
def advance(y, t0, t_stop, h, par, omega_star, rtol, atol, max_steps, rec, stride):
    """Advance y in place from t0 until an event or t_stop (landing exactly on t_stop).

    rec: (nrec_max, N_REC_FIXED + K) buffer; one row stored every `stride` accepted steps (row 0 = start).
    Returns (t, h, status, steps, nrec).
    """
    t = t0
    n = y.shape[0]
    kk = np.empty((NSTAGE, n))
    tmp = np.empty(n)
    ynew = np.empty(n)
    err = np.empty(n)
    rhs(t, y, par, omega_star, kk[0])
    nrec = 0
    if rec.shape[0] > 0:
        record_row(rec, 0, t, y, par, omega_star)
        nrec = 1
    steps = 0
    status = event_status(y, par)
    while status == ST_RUNNING:
        if t >= t_stop:
            return t, h, ST_TIMEOUT, steps, nrec
        if steps >= max_steps:
            return t, h, ST_STEPS, steps, nrec
        clipped = h >= t_stop - t
        h = t_stop - t if clipped else h
        rk_step(t, y, h, kk, par, omega_star, tmp, ynew, err)
        en = err_norm(err, y, ynew, rtol, atol)
        if en <= 1.0:
            t = t_stop if clipped else t + h
            y[:] = ynew
            kk[0][:] = kk[6]
            steps += 1
            if nrec < rec.shape[0] and steps % stride == 0:
                record_row(rec, nrec, t, y, par, omega_star)
                nrec += 1
            status = event_status(y, par)
        fac = SAFETY * en ** (-ERR_EXPONENT) if en > 0.0 else FAC_MAX
        h *= min(FAC_MAX, max(FAC_MIN, fac))
        if h < H_MIN:
            return t, h, ST_FAILED, steps, nrec
    return t, h, status, steps, nrec


def run(y0, par, p_star_days, t_stop_yr, rtol=RTOL_DEFAULT, atol=ATOL_DEFAULT, rec_max=0, stride=1,
        max_steps=MAX_STEPS_DEFAULT):
    """Convenience wrapper. Returns dict with final state, status, steps and (optional) record array."""
    om = np.array([omega_star_of_period(p) for p in np.atleast_1d(p_star_days)])
    y = y0.copy()
    rec = np.zeros((rec_max, N_REC_FIXED + om.size))
    t, h, status, steps, nrec = advance(y, 0.0, t_stop_yr, H_INITIAL, par, om, rtol, atol, max_steps, rec, stride)
    return {"t": t, "y": y, "status": status, "steps": steps, "rec": rec[:nrec], "h": h}


def theta_sl_final(y, n_spin):
    """theta_sl (deg) for each stellar spin at the end state."""
    lh = y[0:3] / np.linalg.norm(y[0:3])
    out = np.empty(n_spin)
    for k in range(n_spin):
        s = y[N_ORBIT + N_SPIN * k:N_ORBIT + N_SPIN * (k + 1)]
        out[k] = math.degrees(math.acos(np.clip(np.dot(s, lh) / np.linalg.norm(s), -1.0, 1.0)))
    return out


def t_k_years(par: np.ndarray, a: float) -> float:
    """A16 Eq. tk with n = sqrt(G Mtot / a^3), a_b,eff = a_b (e_b = 0)."""
    n = math.sqrt(par[P_G] * par[P_MTOT] / a**3)
    return par[P_MTOT] / par[P_MB] * (par[P_AB] / a) ** 3 / n


def srf_epsilons(par: np.ndarray, a: float, omega_p_rad_per_yr: float) -> tuple[float, float, float]:
    """(eps_GR, eps_Tide, eps_Rot) of A16 Eqs. epsigr, epsitide and the rot analogue."""
    g, mt, mb, ab, c = par[P_G], par[P_MTOT], par[P_MB], par[P_AB], par[P_C]
    eps_gr = 3.0 * g * mt**2 * ab**3 / (mb * a**4 * c**2)
    eps_tide = 15.0 * par[P_K2P] * par[P_MSTAR] * mt * ab**3 * par[P_RP] ** 5 / (mb * par[P_MP] * a**8)
    om_hat = omega_p_rad_per_yr / par[P_OMBRK]
    eps_rot = 1.5 * par[P_KQP] * om_hat**2 * mt / mb * (par[P_RP] / a) ** 2 * (ab / a) ** 3
    return eps_gr, eps_tide, eps_rot
