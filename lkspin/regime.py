"""Step 5: per-system features of a stellar-spin phase ensemble and a first-pass outcome classifier (S16 Fig. 14).

Winding sign convention. S16 Eq. (rescond): Omega_ps_bar = -alpha_bar cos(theta_sl) = N n_e. The integrator's winding
variable is dPsi/dt = +alpha cos(theta_sl) = -Omega_ps, so the winding per LK cycle W_k = Delta Psi / (2 pi) equals -N
(W > 0 for theta_sl < 90 deg, where S precesses retrograde about L). Everything here uses W; "IV.|N|" labels use |N| = |W|.

Cycle bookkeeping. The integrator logs every extremum of e. A cycle runs from one e minimum to the next; its amplitude
is e_max - e_min. Cycles count as LK cycles while the amplitude is at least `lk_amp_min` (first drop below it ends them).
"""

import math

import numpy as np

from lkspin import secular as S
from lkspin.lk_orbit import make_orbit

TWO_PI = 2.0 * math.pi
HALF_LAST = 0.5
MIN_CYCLES_FOR_STATS = 3
SHAPE_EXTEND_FACTOR = 2.0
SHAPE_MAX_EXTENSIONS = 6
SHAPE_REC_MAX = 400_000
P_REFERENCE_DAYS = 1.0  # shape coefficients are computed at this P_star and scaled by 1 / P_star
LABEL_I, LABEL_II, LABEL_III, LABEL_V, LABEL_UNCLASSIFIED = "I", "II", "III", "V", "unclassified"
PRED_I, PRED_II, PRED_IV = "I", "II", "IV"
CLASS_ORDER = ("I", "II", "III", "IV.1", "IV.2", "IV.3", "IV.4+", "V", "unclassified")
MAX_NAMED_ADVECTION = 3


# ------------------------------------------------------------------ the sample
def draw_sample(cfg5: dict) -> dict:
    """n_sample systems: M_p uniform, P_star log-uniform, a_b uniform, theta_lb0 uniform (one default_rng stream)."""
    rng = np.random.default_rng(cfg5["seed"])
    n = cfg5["n_sample"]
    lo, hi = cfg5["p_star_days_range"]
    return {"m_p_mjup": rng.uniform(*cfg5["m_p_mjup_range"], n),
            "p_star_days": np.exp(rng.uniform(math.log(lo), math.log(hi), n)),
            "a_b_au": rng.uniform(*cfg5["a_b_au_range"], n),
            "theta_lb0_deg": rng.uniform(*cfg5["theta_lb0_deg_range"], n)}


def period_offsets(cfg5: dict) -> np.ndarray:
    return np.linspace(-cfg5["period_spread"], cfg5["period_spread"], cfg5["n_members"])


def middle_index(offsets: np.ndarray) -> int:
    """Member whose offset d is closest to 0 (the first of a tie)."""
    return int(np.argmin(np.abs(offsets)))


def system_config(cfg4: dict, m_p_mjup: float, a_b_au: float) -> dict:
    return dict(cfg4, m_p_mjup=float(m_p_mjup), a_b_au=float(a_b_au))


def twin_state(y0: np.ndarray, spin_index: int, angle_rad: float) -> np.ndarray:
    """Copy of y0 whose spin `spin_index` is rotated by angle_rad away from L_hat (toward the x axis, which is
    perpendicular to L_hat0 = (0, -sin, cos))."""
    y = y0.copy()
    o = S.N_ORBIT + S.N_SPIN * spin_index
    lhat = y[0:3] / np.linalg.norm(y[0:3])
    perp = np.array([1.0, 0.0, 0.0])
    y[o:o + S.N_SPIN] = math.cos(angle_rad) * lhat + math.sin(angle_rad) * perp
    return y


# ------------------------------------------------------------------ initial shape parameters (step4_validate check 5)
def first_cycle_slice(rec: np.ndarray) -> slice | None:
    """t = 0 up to the first interior eccentricity minimum, or None if the record has none yet."""
    e = rec[:, 2]
    interior = np.where((e[1:-1] < e[:-2]) & (e[1:-1] <= e[2:]))[0] + 1
    return slice(0, interior[0] + 1) if interior.size else None


def first_cycle_record(par: np.ndarray, cfg: dict, theta_deg: float, t_first_yr: float) -> np.ndarray:
    """Dense record covering the first LK cycle; the integration time is doubled until a minimum is inside it."""
    y0 = S.initial_state(par, cfg, theta_deg, 1)
    t_stop = t_first_yr
    rec = S.run(y0, par, [P_REFERENCE_DAYS], t_stop, rec_max=SHAPE_REC_MAX)["rec"]
    for _ in range(SHAPE_MAX_EXTENSIONS):
        if first_cycle_slice(rec) is not None:
            return rec
        t_stop *= SHAPE_EXTEND_FACTOR
        rec = S.run(y0, par, [P_REFERENCE_DAYS], t_stop, rec_max=SHAPE_REC_MAX)["rec"]
    return rec


def shape_coefficients(par: np.ndarray, cfg: dict, theta_deg: float, t_first_yr: float) -> dict:
    """N_max,0 * P_star, A_bar,0 * P_star and (S/L at the first e maximum) * P_star, P_star in days.

    N_max,0 = alpha_bar / n_e and A_bar,0 = alpha_bar / (2 <Omega_pl sin theta_lb>_t) over the first LK cycle, as in
    scripts/step4_validate.py check 5 (trapezoid over the recorded steps).
    """
    rec = first_cycle_record(par, cfg, theta_deg, t_first_yr)
    sl = first_cycle_slice(rec)
    if sl is None:
        return {"nmax_coef": math.nan, "abar_coef": math.nan, "s_over_l_coef": math.nan}
    t, a, e, th, node = (rec[sl, i] for i in range(5))
    period = t[-1] - t[0]
    j = np.sqrt(1.0 - e**2)
    alpha = par[S.P_STAR_COEF] * S.omega_star_of_period(P_REFERENCE_DAYS) / (a**3 * j**3)
    alpha_bar = np.trapezoid(alpha, t) / period
    dnode = np.diff(np.unwrap(node))
    mean_opl_sin = abs(np.sum(0.5 * (np.sin(th[1:]) + np.sin(th[:-1])) * dnode)) / period
    imax = int(np.argmax(e))
    return {"nmax_coef": alpha_bar / (TWO_PI / period), "abar_coef": alpha_bar / (2.0 * mean_opl_sin),
            "s_over_l_coef": s_over_l(par, a[imax], e[imax], P_REFERENCE_DAYS)}


def s_over_l(par: np.ndarray, a: float, e: float, p_star_days: float) -> float:
    """S_star / L with S_star = k* M* R*^2 Omega*, L = mu sqrt(G M_tot a (1 - e^2)) (Msun, AU, yr)."""
    s_star = par[S.P_KS] * par[S.P_MSTAR] * par[S.P_RSTAR] ** 2 * S.omega_star_of_period(p_star_days)
    mu = par[S.P_MSTAR] * par[S.P_MP] / par[S.P_MTOT]
    ang_mom = mu * math.sqrt(par[S.P_G] * par[S.P_MTOT] * a * (1.0 - e * e))
    return s_star / ang_mom


# ------------------------------------------------------------------ cycles and winding
def split_events(ev: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """(minima rows, maxima rows) of the event array."""
    return ev[ev[:, 3] == S.KIND_MIN], ev[ev[:, 3] == S.KIND_MAX]


def cycle_amplitudes(minima: np.ndarray, maxima: np.ndarray) -> np.ndarray:
    """e_max - e_min for each cycle (between minimum j and minimum j+1); 0 if no maximum is logged inside."""
    n_cyc = max(0, minima.shape[0] - 1)
    amp = np.zeros(n_cyc)
    which = np.searchsorted(minima[:, 0], maxima[:, 0]) - 1  # cycle index of each maximum
    for w, e_max in zip(which, maxima[:, 2], strict=True):
        if 0 <= w < n_cyc:
            amp[w] = max(amp[w], e_max - minima[w, 2])
    return amp


def n_active_cycles(amp: np.ndarray, amp_min: float) -> int:
    """Number of leading cycles with amplitude >= amp_min (the LK-active phase)."""
    below = np.where(amp < amp_min)[0]
    return int(below[0]) if below.size else amp.size


def winding_series(minima: np.ndarray, spin: int, n_spin: int) -> np.ndarray:
    """W_k = (Psi(t_{k+1}) - Psi(t_k)) / 2 pi between successive e minima for the given spin."""
    psi = minima[:, S.EV_FIXED + n_spin + spin]
    return np.diff(psi) / TWO_PI


def locked_mask(w: np.ndarray, tol: float) -> np.ndarray:
    return np.abs(w - np.round(w)) < tol


def window_stats(w: np.ndarray, tol: float) -> dict:
    """Locked fraction, dominant integer (round of the median W over locked cycles) and median W of a window."""
    if w.size == 0:
        return {"locked_frac": math.nan, "n_dom": math.nan, "w_median": math.nan}
    lock = locked_mask(w, tol)
    n_dom = float(np.round(np.median(w[lock]))) if lock.any() else math.nan
    return {"locked_frac": float(lock.mean()), "n_dom": n_dom, "w_median": float(np.median(w))}


def winding_features(w_active: np.ndarray, cfg5: dict) -> dict:
    """Winding statistics for one spin over the LK-active cycles (last `last_frac` and first `first_frac` windows)."""
    n = w_active.size
    n_last = max(1, math.ceil(cfg5["last_frac"] * n)) if n else 0
    n_first = max(1, math.ceil(cfg5["first_frac"] * n)) if n else 0
    last = window_stats(w_active[n - n_last:], cfg5["lock_tol"])
    first = window_stats(w_active[:n_first], cfg5["lock_tol"])
    return {"n_cycles": n, "locked_frac_last": last["locked_frac"], "n_dom": last["n_dom"],
            "w_median_last": last["w_median"], "locked_frac_first": first["locked_frac"],
            "n_dom_first": first["n_dom"], "w_median_first": first["w_median"]}


# ------------------------------------------------------------------ ensemble statistics
def angle_between_deg(u: np.ndarray, v: np.ndarray) -> float:
    cs = np.dot(u, v) / (np.linalg.norm(u) * np.linalg.norm(v))
    return math.degrees(math.acos(min(1.0, max(-1.0, cs))))


def spin_vector(y: np.ndarray, k: int) -> np.ndarray:
    return y[S.N_ORBIT + S.N_SPIN * k:S.N_ORBIT + S.N_SPIN * (k + 1)]


def ensemble_stats(theta_f: np.ndarray) -> dict:
    """Sorted angles, range, largest gap, counts on each side of it and the range inside each side."""
    srt = np.sort(theta_f)
    gaps = np.diff(srt)
    g = int(np.argmax(gaps))
    return {"theta_sorted": srt, "range": float(srt[-1] - srt[0]), "gap": float(gaps[g]), "n_low": g + 1,
            "n_high": srt.size - g - 1, "range_low": float(srt[g] - srt[0]), "range_high": float(srt[-1] - srt[g + 1])}


# ------------------------------------------------------------------ first-pass classifier
def is_bimodal(f: dict, cfg5: dict) -> bool:
    """Class I: a wide gap with enough members and tight clusters on both sides."""
    return (f["gap"] >= cfg5["bimodal_gap_deg"] and min(f["n_low"], f["n_high"]) >= cfg5["bimodal_min_side"]
            and max(f["range_low"], f["range_high"]) <= cfg5["bimodal_side_range_deg"])


def is_chaotic(f: dict, cfg5: dict) -> bool:
    """Class III: tiny twin offset blows up, or the ensemble spans a wide range."""
    return f["twin_div_deg"] >= cfg5["chaos_twin_deg"] or f["range"] >= cfg5["chaos_range_deg"]


def is_advected(f: dict, cfg5: dict) -> bool:
    """Class IV: most late LK cycles are locked to a non-zero integer winding."""
    locked, n_dom = f["m0_locked_frac_last"], f["m0_n_dom"]
    return bool(locked >= cfg5["advect_lock_frac"] and n_dom != 0 and not math.isnan(n_dom))


def is_fully_adiabatic(f: dict, cfg5: dict) -> bool:
    """Class V: spin follows L (all theta_f small) and winds fast without locking."""
    return bool(f["theta_max"] < cfg5["adiabatic_theta_deg"] and abs(f["m0_w_median_last"]) > cfg5["adiabatic_min_w"])


def is_stationary(f: dict, cfg5: dict) -> bool:
    """Class II: regular and bounded (max theta_f below the stationary limit)."""
    return bool(f["theta_max"] < cfg5["stationary_theta_deg"])


def advection_label(n_dom: float) -> str:
    n = int(abs(n_dom))
    return f"IV.{n}" if n <= MAX_NAMED_ADVECTION else "IV.4+"


def classify(f: dict, cfg5: dict) -> str:
    """First-pass class in the brief's order: I, III, IV.N, V, II, else unclassified."""
    if is_bimodal(f, cfg5):
        return LABEL_I
    if is_chaotic(f, cfg5):
        return LABEL_III
    if is_advected(f, cfg5):
        return advection_label(f["m0_n_dom"])
    if is_fully_adiabatic(f, cfg5):
        return LABEL_V
    if is_stationary(f, cfg5):
        return LABEL_II
    return LABEL_UNCLASSIFIED


def predicted_class(a_bar0: float, n_max0: float, cfg5: dict) -> str:
    """S16 Fig. 14 map (findings 01): I if A < 1; else II if N_max < 0.5; else IV."""
    if a_bar0 < cfg5["map_a_bar_split"]:
        return PRED_I
    return PRED_II if n_max0 < cfg5["map_n_max_split"] else PRED_IV


def label_family(label: str) -> str:
    """IV.N -> IV, otherwise unchanged (for comparison with predicted_class)."""
    return "IV" if label.startswith("IV") else label


# ------------------------------------------------------------------ one system
def member_features(minima: np.ndarray, maxima: np.ndarray, spin: int, n_spin: int, cfg5: dict) -> dict:
    amp = cycle_amplitudes(minima, maxima)
    n_act = n_active_cycles(amp, cfg5["lk_amp_min"])
    w_all = winding_series(minima, spin, n_spin)
    feats = winding_features(w_all[:n_act], cfg5)
    feats["n_cycles_total"] = int(w_all.size)
    return feats


def ensemble_features(res: dict, n_members: int, mid: int) -> dict:
    """theta_f, twin divergence and ensemble statistics from a run_diag result with n_members + 1 spins."""
    theta_all = S.theta_sl_final(res["y"], n_members + 1)
    twin = n_members
    stats = ensemble_stats(theta_all[:n_members])
    return {"theta_f": theta_all[:n_members], "theta_f_twin": float(theta_all[twin]),
            "twin_div_deg": angle_between_deg(spin_vector(res["y"], twin), spin_vector(res["y"], mid)),
            "twin_div_m0_deg": angle_between_deg(spin_vector(res["y"], twin), spin_vector(res["y"], 0)),
            "theta_sorted": stats["theta_sorted"], "range": stats["range"], "gap": stats["gap"],
            "n_low": stats["n_low"], "n_high": stats["n_high"], "range_low": stats["range_low"],
            "range_high": stats["range_high"], "theta_max": float(theta_all[:n_members].max())}


def series_columns(minima: np.ndarray, maxima: np.ndarray, spins: tuple, n_spin: int) -> np.ndarray:
    """Raw per-cycle series: t, a, e_min, amplitude, then (p, Psi) at each e minimum for the given spins."""
    amp = np.append(cycle_amplitudes(minima, maxima), 0.0)
    cols = [minima[:, 0], minima[:, 1], minima[:, 2], amp]
    for k in spins:
        cols += [minima[:, S.EV_FIXED + k], minima[:, S.EV_FIXED + n_spin + k]]
    return np.column_stack(cols)


def integrate_system(m_p: float, p_star: float, a_b: float, theta_lb0: float, cfg4: dict, cfg5: dict) -> dict:
    """Run the K + 1 spin ensemble for one system. Returns status, timing, and (if a hot Jupiter) features + series."""
    cfg = system_config(cfg4, m_p, a_b)
    par = S.make_params(cfg)
    offsets = period_offsets(cfg5)
    mid = middle_index(offsets)
    periods = list(p_star * (1.0 + offsets)) + [p_star * (1.0 + offsets[mid])]
    n_spin = len(periods)
    y0 = S.initial_state(par, cfg, theta_lb0, n_spin)
    y_init = twin_state(y0, n_spin - 1, cfg5["twin_angle_rad"])
    res = S.run_diag(y_init, par, periods, cfg["stop_t_gyr"] * 1e9, ev_max=cfg5["ev_max"])
    out = {"status": res["status"], "t_end": res["t"], "steps": res["steps"], "n_events": res["ev"].shape[0]}
    if res["status"] != S.ST_HOT_JUPITER:
        return out
    out.update(hot_jupiter_features(res, par, cfg, theta_lb0, periods, mid, cfg5))
    return out


def hot_jupiter_features(res: dict, par: np.ndarray, cfg: dict, theta_lb0: float, periods: list, mid: int,
                         cfg5: dict) -> dict:
    n_members = cfg5["n_members"]
    n_spin = n_members + 1
    minima, maxima = split_events(res["ev"])
    feats = ensemble_features(res, n_members, mid)
    wind = member_features(minima, maxima, 0, n_spin, cfg5)
    feats.update({f"m0_{k}": v for k, v in wind.items()})
    t_first = cfg5["shape_ncycle_tk"] * orbit_period_yr(par, cfg, theta_lb0)
    shape = shape_coefficients(par, cfg, theta_lb0, t_first)
    p_mid = periods[mid]
    feats.update({"a_bar0": shape["abar_coef"] / p_mid, "n_max0": shape["nmax_coef"] / p_mid,
                  "s_over_l": shape["s_over_l_coef"] / p_mid, "p_mid_days": p_mid})
    feats["series"] = series_columns(minima, maxima, (0, mid, n_members), n_spin)
    return feats


def orbit_period_yr(par: np.ndarray, cfg: dict, theta_lb0: float) -> float:
    """Analytic (no-SRF) LK period in years, used only to size the first integration for the shape parameters."""
    return make_orbit(theta_lb0, cfg["e0"]).period * S.t_k_years(par, cfg["a0_au"])


CLASSIFIER_FIELDS = ("gap", "n_low", "n_high", "range_low", "range_high", "range", "twin_div_deg", "theta_max",
                     "m0_locked_frac_last", "m0_n_dom", "m0_w_median_last")


def classify_table(table: dict, cfg5: dict) -> list[str]:
    """First-pass label for every row of a features table (dict of equal-length arrays, e.g. features.npz)."""
    n = len(table[CLASSIFIER_FIELDS[0]])
    return [classify({k: table[k][i] for k in CLASSIFIER_FIELDS}, cfg5) for i in range(n)]


def predicted_table(table: dict, cfg5: dict) -> list[str]:
    return [predicted_class(a, n, cfg5) for a, n in zip(table["a_bar0"], table["n_max0"], strict=True)]


# ------------------------------------------------------------------ revised classifier (v2, after hand validation)
# Built on hand-labelled cases 1-30 (results/step5/handlabels.csv); cases 31-50 held out. See findings 04.
SER_P_M0, SER_PSI_M0, SER_AMP = 4, 5, 3


def active_rows(series: np.ndarray, amp_min: float) -> int:
    """Leading LK-active cycles of a per-system series (amplitude column)."""
    return n_active_cycles(series[:-1, SER_AMP], amp_min)


def max_swing_deg(series: np.ndarray, n_act: int) -> float:
    """Largest theta_sl (deg) of member 0 at e minima during the LK-active phase."""
    p = np.clip(series[: n_act + 1, SER_P_M0], -1.0, 1.0)
    return float(np.degrees(np.arccos(p)).max())


def longest_lock(w: np.ndarray, tol: float) -> tuple[int, int]:
    """(longest run of consecutive cycles with W within tol of the same non-zero integer, that integer)."""
    best, best_n, run, prev = 0, 0, 0, None
    for wk in w:
        n = int(np.round(wk))
        ok = n != 0 and abs(wk - n) < tol
        run = run + 1 if ok and n == prev else (1 if ok else 0)
        prev = n if ok else None
        if run > best:
            best, best_n = run, n
    return best, best_n


def best_lock_fraction(w: np.ndarray, tol: float) -> tuple[float, int]:
    """Largest fraction of LK-active cycles with W within tol of one non-zero integer N, and that N."""
    if w.size == 0:
        return 0.0, 0
    n_max = int(np.ceil(np.abs(w).max())) + 1
    fracs = [(float((np.abs(w - n) < tol).mean()), n) for n in range(-n_max, n_max + 1) if n != 0]
    return max(fracs)


def cluster_count(values: np.ndarray, tight_deg: float, gap_deg: float) -> int:
    """1 if all values lie within tight_deg; 2 if a gap >= gap_deg splits them into two tight groups; else 0 (scatter)."""
    srt = np.sort(values)
    if srt[-1] - srt[0] <= tight_deg:
        return 1
    g = int(np.argmax(np.diff(srt)))
    two = srt[g + 1] - srt[g] >= gap_deg and srt[g] - srt[0] <= tight_deg and srt[-1] - srt[g + 1] <= tight_deg
    return 2 if two else 0


def v2_features(theta_all: np.ndarray, series: np.ndarray, cfg5: dict) -> dict:
    n_act = active_rows(series, cfg5["lk_amp_min"])
    w = np.diff(series[: n_act + 1, SER_PSI_M0]) / TWO_PI
    run, n_lock = longest_lock(w, cfg5["v2_lock_tol"])
    frac, n_frac = best_lock_fraction(w, cfg5["v2_lock_tol"])
    return {"clusters": cluster_count(theta_all, cfg5["v2_tight_deg"], cfg5["v2_gap_deg"]),
            "swing": max_swing_deg(series, series.shape[0] - 1), "lock_run": run, "lock_n": n_lock,
            "lock_frac": frac, "lock_frac_n": n_frac,
            "w_median_active": float(np.median(np.abs(w))) if w.size else math.nan,
            "theta_f_max": float(np.max(theta_all))}


def classify_v2(g: dict, cfg5: dict) -> str:
    """III (scatter) -> I (<=2 tight clusters after a large swing) -> IV.N (locked fraction) -> V -> II."""
    if g["clusters"] == 0:
        return LABEL_III
    if g["swing"] >= cfg5["v2_swing_deg"]:
        return LABEL_I
    if g["lock_frac"] >= cfg5["v2_lock_frac"]:
        return advection_label(g["lock_frac_n"])
    if g["theta_f_max"] < cfg5["adiabatic_theta_deg"] and g["w_median_active"] >= cfg5["adiabatic_min_w"]:
        return LABEL_V
    return LABEL_II
