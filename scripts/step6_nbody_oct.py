"""Step 6 validation 3: the octupole secular equations (lkspin.secular.oct_terms, A16 Eqs. Lvec/evec) against a direct
REBOUND N-body integration (IAS15, point masses, no GR) of star + planet + eccentric companion.

SRFs, tides and feedback are off. Fiducial: M* = 1, Mp = 1 MJ, a = 1 AU, a_b = 100 AU, e_b = 0.7 (so eps_oct = 0.0137,
a_b,eff = 71.4 AU), theta_lb0 = 85 deg, e0 = 0.01, e_b-hat at azimuth phi from the planet's initial node line (+x).
Both sides give a block-averaged e(t) (window = one outer period for the N-body; the secular curve is sampled at the block
midpoints) and the e_max of every LK cycle. The N-body and the secular run use the same initial elements (the N-body's
osculating e0, so a tiny mismatch in the mean e0 is expected).
Writes results/step6/validation/octupole_nbody.{json,npz,png}.
Run:  uv run python scripts/step6_nbody_oct.py [--t-end-tk 22] [--phi 0.9] [--time-only]
"""

import argparse
import json
import math
import time

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import rebound
from scipy.signal import find_peaks

from lkspin import secular as S

CFG4 = S.load_config()
OUT = S.ROOT / "results" / "step6" / "validation"
M_P_MJUP = 1.0
A_AU = 1.0
A_B_AU = 100.0
E_B = 0.7
THETA_DEG = 85.0
E0 = 0.01
P_STAR_DAYS = 10.0  # irrelevant (feedback off); the extended state needs a spin
SAMPLES_PER_BLOCK = 100
PROMINENCE_DEX = 0.3  # a peak of -log10(1 - e) must stand this far above its surroundings to count as an e_max
PLANET = 1
FIG_SIZE = (9.0, 6.0)
DPI = 140


def secular_setup(phi: float) -> tuple[np.ndarray, dict, np.ndarray]:
    cfg = dict(CFG4, m_p_mjup=M_P_MJUP, a0_au=A_AU, e0=E0, a_b_au=A_B_AU * math.sqrt(1.0 - E_B**2))
    base = S.make_params(cfg, srf=False, tides=False, a0_au=A_AU)
    par = S.extend_params(base, fb_mult=0.0, e_b=E_B, phi_b=phi, octupole=True)
    par[S.P_RTIDE], par[S.P_STOPA] = 0.0, 0.0  # no disruption or hot-Jupiter stop: point masses, as in the N-body
    return par, cfg, S.initial_state_ext(par, cfg, THETA_DEG, P_STAR_DAYS)


def outer_period_yr(cfg: dict) -> float:
    m_in = cfg["m_star_msun"] + cfg["m_p_mjup"] * S.MJUP_PER_MSUN
    return 2.0 * math.pi * math.sqrt(A_B_AU**3 / (S.G_AU3_PER_MSUN_YR2 * (m_in + cfg["m_b_msun"])))


def secular_curve(par: np.ndarray, y0: np.ndarray, t_mid: np.ndarray) -> np.ndarray:
    """e at the given times (chunked integration; the restart of the step size at each chunk is harmless here)."""
    y, t, out = y0.copy(), 0.0, np.empty(t_mid.size)
    for i, t_next in enumerate(t_mid):
        res = S.run_ext(y, par, t_next, t_start_yr=t, ev_max=10)
        y, t = res["y"], t_next
        out[i] = np.linalg.norm(y[3:6])
    return out


def build_sim(cfg: dict, phi: float, mean_anom: float) -> rebound.Simulation:
    sim = rebound.Simulation()
    sim.G = S.G_AU3_PER_MSUN_YR2
    sim.integrator = "ias15"
    sim.add(m=cfg["m_star_msun"])
    sim.add(m=cfg["m_p_mjup"] * S.MJUP_PER_MSUN, a=A_AU, e=E0, inc=math.radians(THETA_DEG), Omega=0.0, omega=0.0, M=mean_anom)
    sim.add(m=cfg["m_b_msun"], a=A_B_AU, e=E_B, inc=0.0, Omega=0.0, omega=phi, M=mean_anom)
    sim.move_to_com()
    return sim


def nbody_blocks(sim: rebound.Simulation, p_out: float, n_blocks: int) -> np.ndarray:
    """(n_blocks, 4): block midpoint, mean e, min e, mean theta_lb (deg, planet orbit vs the initial companion plane)."""
    dt = p_out / SAMPLES_PER_BLOCK
    out = np.empty((n_blocks, 4))
    t = 0.0
    for b in range(n_blocks):
        samples = np.empty((SAMPLES_PER_BLOCK, 2))
        for i in range(SAMPLES_PER_BLOCK):
            sim.integrate(t + (i + 1) * dt)
            o = sim.particles[PLANET].orbit()
            samples[i] = (o.e, o.inc)
        out[b] = (t + 0.5 * p_out, samples[:, 0].mean(), samples[:, 0].min(), math.degrees(samples[:, 1].mean()))
        t += p_out
    return out


def find_maxima(t: np.ndarray, e: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Indices of e maxima as peaks of -log10(1 - e)."""
    peaks, _ = find_peaks(-np.log10(1.0 - np.minimum(e, 1.0 - 1e-12)), prominence=PROMINENCE_DEX)
    return t[peaks], e[peaks]


def compare(t: np.ndarray, e_nb: np.ndarray, e_sec: np.ndarray, tk: float) -> dict:
    tn, en = find_maxima(t, e_nb)
    ts, es = find_maxima(t, e_sec)
    n = min(tn.size, ts.size)
    return {"n_max_nbody": int(tn.size), "n_max_secular": int(ts.size),
            "t_max_nbody_tk": (tn / tk).tolist(), "t_max_secular_tk": (ts / tk).tolist(),
            "e_max_nbody": en.tolist(), "e_max_secular": es.tolist(),
            "one_minus_emax_ratio": ((1.0 - en[:n]) / (1.0 - es[:n])).tolist(),
            "dt_max_tk": ((tn[:n] - ts[:n]) / tk).tolist()}


def plot(t: np.ndarray, e_nb: np.ndarray, e_sec: np.ndarray, tk: float, path) -> None:
    fig, axes = plt.subplots(2, 1, figsize=FIG_SIZE, sharex=True)
    axes[0].plot(t / tk, e_nb, label="N-body (block mean)")
    axes[0].plot(t / tk, e_sec, "--", label="secular + octupole")
    axes[0].set(ylabel="e", ylim=(0, 1))
    axes[1].plot(t / tk, 1.0 - e_nb)
    axes[1].plot(t / tk, 1.0 - e_sec, "--")
    axes[1].set(yscale="log", ylabel="1 - e", xlabel=r"$t / t_k$")
    axes[0].legend()
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--t-end-tk", type=float, default=22.0)
    ap.add_argument("--phi", type=float, default=0.9, help="azimuth of e_b-hat from the initial node line (rad)")
    ap.add_argument("--mean-anom", type=float, default=0.3)
    ap.add_argument("--time-only", action="store_true", help="run 1 t_k to time the N-body, write nothing")
    args = ap.parse_args()
    par, cfg, y0 = secular_setup(args.phi)
    tk = S.t_k_years(par, A_AU)
    p_out = outer_period_yr(cfg)
    n_blocks = max(1, int((1.0 if args.time_only else args.t_end_tk) * tk / p_out))
    sim = build_sim(cfg, args.phi, args.mean_anom)
    t0 = time.time()
    blocks = nbody_blocks(sim, p_out, n_blocks)
    wall = time.time() - t0
    print(f"N-body {n_blocks} blocks of {p_out:.1f} yr ({n_blocks * p_out / tk:.2f} t_k = {n_blocks * p_out:.3g} yr): "
          f"{wall:.0f} s wall, {sim.steps_done} steps; eps_oct = {par[S.P_EOCT] * A_AU:.5f}, t_k = {tk:.0f} yr")
    if args.time_only:
        return
    t = blocks[:, 0]
    e_sec = secular_curve(par, y0, t)
    cmp = compare(t, blocks[:, 1], e_sec, tk)
    OUT.mkdir(parents=True, exist_ok=True)
    np.savez(OUT / "octupole_nbody.npz", t_yr=t, e_nbody=blocks[:, 1], e_min_nbody=blocks[:, 2], theta_nbody=blocks[:, 3],
             e_secular=e_sec, t_k=tk)
    res = {"eps_oct": float(par[S.P_EOCT] * A_AU), "t_k_yr": tk, "wall_s": wall, "phi": args.phi, "t_end_tk": n_blocks * p_out / tk,
           "steps": int(sim.steps_done), **cmp}
    (OUT / "octupole_nbody.json").write_text(json.dumps(res, indent=2) + "\n")
    plot(t, blocks[:, 1], e_sec, tk, OUT / "octupole_nbody.png")
    for k, v in cmp.items():
        print(k, np.round(v, 5).tolist() if isinstance(v, list) else v)


if __name__ == "__main__":
    main()
