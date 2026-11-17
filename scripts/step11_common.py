"""Step 11 shared pieces: config, the 20 systems, the secular runs (GR-only and full SRFs) and the breakdown indicators.

Secular side: lkspin.secular with tides off. GR-only is made by zeroing k_2p and k_qp in the parameter array (the
bulge and rotation precession terms are proportional to them), and by disabling the disruption event so the run
continues past r_Tide. Nothing in lkspin is modified.
"""

import json
import math
import tomllib
from dataclasses import asdict, dataclass

import numpy as np
from scipy.signal import find_peaks

from lkspin import secular as S

with (S.ROOT / "config.toml").open("rb") as _fh:
    _ALL = tomllib.load(_fh)
CFG4 = _ALL["step4"]
CFG = _ALL["step11"]
OUT_DIR = S.ROOT / CFG["results"]
CHUNK_DIR = OUT_DIR / "chunks"

TWO_PI = 2.0 * math.pi
FIRST_SPAN_TK = 8.0          # first secular integration span, in t_k; doubled until two e_max are found
MAX_DOUBLINGS = 4
N_PEAKS_WANTED = 2
POOL_EXTRA_SEED_OFFSET = 0   # extras continue the same rng stream


@dataclass(frozen=True)
class System:
    name: str
    kind: str                 # fiducial, random or extreme
    m_p_mjup: float
    a_b_au: float
    theta_lb0_deg: float
    mean_anom_planet: float   # rad
    mean_anom_comp: float     # rad


def _draw(rng: np.random.Generator, n: int) -> np.ndarray:
    """n rows of (M_p, a_b, theta, M_planet, M_companion), in that column order."""
    cols = [rng.uniform(*CFG["m_p_mjup_range"], n), rng.uniform(*CFG["a_b_au_range"], n),
            rng.uniform(*CFG["theta_lb0_deg_range"], n), rng.uniform(0.0, TWO_PI, n), rng.uniform(0.0, TWO_PI, n)]
    return np.column_stack(cols)


def _system(name: str, kind: str, row: np.ndarray) -> System:
    return System(name, kind, float(row[0]), float(row[1]), float(row[2]), float(row[3]), float(row[4]))


def fiducial_systems(rng: np.random.Generator) -> list[System]:
    angles = CFG["fiducial_theta_lb0_deg"]
    phases = rng.uniform(0.0, TWO_PI, (len(angles), 2))
    return [System(f"F{int(th)}", "fiducial", CFG4["m_p_mjup"], CFG4["a_b_au"], th, float(ph[0]), float(ph[1]))
            for th, ph in zip(angles, phases, strict=True)]


def params_for(sys_: System, *, mode: str) -> np.ndarray:
    """Secular parameter array. mode 'gr': GR only; 'full': GR + planet bulge + rotation. Tides always off."""
    cfg = dict(CFG4, m_p_mjup=sys_.m_p_mjup, a_b_au=sys_.a_b_au)
    par = S.make_params(cfg, srf=True, tides=False)
    par[S.P_RTIDE] = 0.0  # no disruption event: compare the equations, flag the Roche limit separately
    if mode == "gr":
        par[S.P_K2P] = 0.0
        par[S.P_KQP] = 0.0
    return par


def peak_indices(t: np.ndarray, e: np.ndarray) -> np.ndarray:
    """Indices of e_max peaks: peaks of -log10(1-e) at least peak_prominence_dex above their surroundings."""
    y = -np.log10(np.maximum(1.0 - e, 1e-12))
    return find_peaks(y, prominence=CFG["peak_prominence_dex"])[0]


def secular_series(sys_: System, mode: str) -> dict:
    """Integrate until two e_max are found. Returns t (yr), e, theta_lb (rad), peak indices, t_k."""
    par = params_for(sys_, mode=mode)
    cfg = dict(CFG4, m_p_mjup=sys_.m_p_mjup, a_b_au=sys_.a_b_au)
    tk = S.t_k_years(par, cfg["a0_au"])
    y0 = S.initial_state(par, cfg, sys_.theta_lb0_deg, 1)
    span = FIRST_SPAN_TK * tk
    for _ in range(MAX_DOUBLINGS):
        r = S.run(y0, par, [CFG["p_star_days"][-1]], span, CFG["secular_rtol"], CFG["secular_atol"],
                  rec_max=CFG["secular_rec_max"])
        rec = r["rec"]
        peaks = peak_indices(rec[:, 0], rec[:, 2])
        if len(peaks) >= N_PEAKS_WANTED:
            break
        span *= 2.0
    return {"t": rec[:, 0], "e": rec[:, 2], "theta": rec[:, 3], "peaks": peaks, "tk": tk, "par": par,
            "a0": cfg["a0_au"]}


def secular_summary(sys_: System, mode: str) -> dict:
    """First/second e_max, their times and theta_lb at the first, for one secular run."""
    s = secular_series(sys_, mode)
    pk = s["peaks"]
    out = {"tk": s["tk"], "n_peaks": len(pk)}
    for i, key in enumerate(("1", "2")):
        ok = i < len(pk)
        out[f"emax{key}"] = s["e"][pk[i]] if ok else math.nan
        out[f"tmax{key}"] = s["t"][pk[i]] if ok else math.nan
    out["theta_at_emax1"] = math.degrees(s["theta"][pk[0]]) if len(pk) else math.nan
    out["period"] = out["tmax2"] - out["tmax1"]
    return out


def p_out_years(sys_: System) -> float:
    mtot = CFG4["m_star_msun"] + sys_.m_p_mjup * S.MJUP_PER_MSUN + CFG4["m_b_msun"]
    return math.sqrt(sys_.a_b_au**3 / mtot)


def p_orb_years(sys_: System) -> float:
    mtot = CFG4["m_star_msun"] + sys_.m_p_mjup * S.MJUP_PER_MSUN
    return math.sqrt(CFG4["a0_au"] ** 3 / mtot)


def r_tide_au(sys_: System) -> float:
    mp = sys_.m_p_mjup * S.MJUP_PER_MSUN
    return CFG["r_tide_coef"] * CFG4["r_p_rjup"] * S.RJUP_AU * (CFG4["m_star_msun"] / mp) ** (1.0 / 3.0)


def alpha_per_year(sys_: System, e: float, p_star_days: float) -> float:
    """Stellar spin precession rate alpha at eccentricity e (A16 Omega_ps without cos(theta_sl)), rad/yr."""
    mp = sys_.m_p_mjup * S.MJUP_PER_MSUN
    rstar = CFG4["r_star_rsun"] * S.RSUN_AU
    j = math.sqrt(1.0 - e * e)
    coef = 1.5 * CFG4["k_q_star"] / CFG4["k_star"] * mp / CFG4["m_star_msun"]
    return coef * (rstar / CFG4["a0_au"]) ** 3 * S.omega_star_of_period(p_star_days) / j**3


def breakdown(sys_: System, emax: float, emax_gr: float) -> dict:
    """Analytic breakdown indicators at the secular e_max (full SRFs); Roche ratio also for the GR-only e_max."""
    tk = S.t_k_years(S.make_params(dict(CFG4, m_p_mjup=sys_.m_p_mjup, a_b_au=sys_.a_b_au)), CFG4["a0_au"])
    t_near = tk * math.sqrt(1.0 - emax**2)
    out = {"t_k": tk, "p_out": p_out_years(sys_), "p_orb": p_orb_years(sys_)}
    out["t_near_over_pout"] = t_near / out["p_out"]
    for p in CFG["p_star_days"]:
        out[f"porb_alpha_P{p:g}d"] = out["p_orb"] * alpha_per_year(sys_, emax, p)
    out["rperi_over_rtide"] = CFG4["a0_au"] * (1.0 - emax) / r_tide_au(sys_)
    out["rperi_over_rtide_gr"] = CFG4["a0_au"] * (1.0 - emax_gr) / r_tide_au(sys_)
    return out


def flags(b: dict) -> dict:
    frac = CFG["flag_fraction"]
    out = {"flag_double_avg": b["t_near_over_pout"] < 1.0 / frac,
           "flag_roche": min(b["rperi_over_rtide"], b["rperi_over_rtide_gr"]) < 1.0}
    for p in CFG["p_star_days"]:
        out[f"flag_porb_alpha_P{p:g}d"] = b[f"porb_alpha_P{p:g}d"] > frac
    return out


def build_systems() -> list[System]:
    """2 fiducials + n_random plain draws + n_extreme systems with the highest secular e_max out of n_pool extra draws."""
    rng = np.random.default_rng(CFG["seed"])
    systems = fiducial_systems(rng)
    plain = _draw(rng, CFG["n_random"])
    systems += [_system(f"R{i:02d}", "random", row) for i, row in enumerate(plain)]
    pool = _draw(rng, CFG["n_pool"])
    cands = [_system("pool", "extreme", row) for row in pool]
    emax = np.array([secular_summary_emax(c) for c in cands])
    top = np.argsort(emax)[::-1][:CFG["n_extreme"]]
    for k, idx in enumerate(top):
        c = cands[idx]
        systems.append(System(f"X{k}", "extreme", c.m_p_mjup, c.a_b_au, c.theta_lb0_deg, c.mean_anom_planet,
                              c.mean_anom_comp))
    return systems


def secular_summary_emax(sys_: System) -> float:
    """Cheap first-cycle e_max with full SRFs: integrate to the first peak only (8 t_k covers it)."""
    par = params_for(sys_, mode="full")
    cfg = dict(CFG4, m_p_mjup=sys_.m_p_mjup, a_b_au=sys_.a_b_au)
    tk = S.t_k_years(par, cfg["a0_au"])
    y0 = S.initial_state(par, cfg, sys_.theta_lb0_deg, 1)
    r = S.run(y0, par, [CFG["p_star_days"][-1]], FIRST_SPAN_TK * tk, CFG["secular_rtol"], CFG["secular_atol"],
              rec_max=CFG["secular_rec_max"])
    return float(r["rec"][:, 2].max())


SYSTEMS_PATH = OUT_DIR / "systems.json"


def save_systems(systems: list[System]) -> None:
    OUT_DIR.mkdir(parents=True, exist_ok=True)
    SYSTEMS_PATH.write_text(json.dumps([asdict(s) for s in systems], indent=1))


def load_systems() -> list[System]:
    return [System(**d) for d in json.loads(SYSTEMS_PATH.read_text())]
