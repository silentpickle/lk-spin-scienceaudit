"""Step 11 N-body runs: star + planet + companion, IAS15 + REBOUNDx "gr", block-averaged output.

Output per system: results/step11/chunks/<name>.npz with one row per averaging block (one outer binary period):
t (block midpoint, yr), e_mean, e_max (osculating, within the block), e_min, theta_mean (deg), a_mean (AU).
The planet's orbit is Jacobi relative to the star (REBOUND's Particle.orbit() default, Jacobi COM of the lower indices).
Run:  uv run python scripts/step11_nbody.py [--only NAME ...] [--validate-gr]
"""

import argparse
import math
import multiprocessing as mp
import os
import time
from pathlib import Path

import numpy as np
import rebound
import reboundx

import step11_common as C
from lkspin import secular as S

SECONDS_PER_HOUR = 3600.0
STATE_DIR = C.OUT_DIR / "state"
PLANET = 1
COLUMNS = ("t", "e_mean", "e_max", "e_min", "theta_mean", "a_mean")


def build_sim(sys_: C.System, e0: float, m_b: float | None = None) -> tuple[rebound.Simulation, reboundx.Extras]:
    """Star, planet (omega = 0, inclination theta_lb0 to the companion's plane, node on +x), companion on a circular orbit."""
    sim = rebound.Simulation()
    sim.G = S.G_AU3_PER_MSUN_YR2
    sim.integrator = "ias15"
    sim.add(m=C.CFG4["m_star_msun"])
    sim.add(m=sys_.m_p_mjup * S.MJUP_PER_MSUN, a=C.CFG4["a0_au"], e=e0, inc=math.radians(sys_.theta_lb0_deg),
            Omega=0.0, omega=math.radians(C.CFG["omega_deg"]), M=sys_.mean_anom_planet)
    sim.add(m=C.CFG4["m_b_msun"] if m_b is None else m_b, a=sys_.a_b_au, e=0.0, inc=0.0, M=sys_.mean_anom_comp)
    sim.move_to_com()
    rebx = reboundx.Extras(sim)
    gr = rebx.load_force(C.CFG["gr_force"])
    rebx.add_force(gr)
    gr.params["c"] = C.CFG["c_au_per_yr"]
    return sim, rebx


def sample_block(sim: rebound.Simulation, t_start: float, n_samples: int, dt: float) -> np.ndarray:
    """(n_samples, 3) of e, inclination to the initial companion plane (rad), a for the planet's Jacobi orbit."""
    out = np.empty((n_samples, 3))
    for i in range(n_samples):
        sim.integrate(t_start + (i + 1) * dt)
        o = sim.particles[PLANET].orbit()
        out[i] = (o.e, o.inc, o.a)
    return out


def block_row(t_mid: float, s: np.ndarray) -> tuple[float, ...]:
    return (t_mid, s[:, 0].mean(), s[:, 0].max(), s[:, 0].min(), math.degrees(s[:, 1].mean()), s[:, 2].mean())


def state_paths(name: str) -> tuple[Path, Path]:
    return STATE_DIR / f"{name}.bin", STATE_DIR / f"{name}.rebx"


def save_checkpoint(name: str, sim: rebound.Simulation, rebx: reboundx.Extras, rows: list, extra: dict) -> None:
    """Rows npz first, then the REBOUND/REBOUNDx state; each written to .tmp and renamed. On resume the rows are cut back
    to the state's time, so a crash between the two renames loses nothing."""
    C.CHUNK_DIR.mkdir(parents=True, exist_ok=True)
    STATE_DIR.mkdir(parents=True, exist_ok=True)
    arr = np.array(rows).reshape(-1, len(COLUMNS))
    rows_tmp = C.CHUNK_DIR / f"{name}.tmp.npz"
    np.savez(rows_tmp, **dict(zip(COLUMNS, arr.T, strict=True)), **{k: np.array(v) for k, v in extra.items()})
    bin_path, rebx_path = state_paths(name)
    bin_tmp, rebx_tmp = bin_path.with_suffix(".bin.tmp"), rebx_path.with_suffix(".rebx.tmp")
    bin_tmp.unlink(missing_ok=True)
    rebx_tmp.unlink(missing_ok=True)
    sim.save_to_file(str(bin_tmp))
    rebx.save(str(rebx_tmp))
    os.replace(rows_tmp, C.CHUNK_DIR / f"{name}.npz")
    os.replace(bin_tmp, bin_path)
    os.replace(rebx_tmp, rebx_path)


def start_or_resume(sys_: C.System, block_len: float) -> tuple[rebound.Simulation, reboundx.Extras, list, float]:
    """(sim, rebx, rows, wall seconds used by earlier sessions). Resumes from results/step11/state if present."""
    bin_path, rebx_path = state_paths(sys_.name)
    if not (bin_path.exists() and rebx_path.exists()):
        sim, rebx = build_sim(sys_, C.CFG4["e0"])
        return sim, rebx, [], 0.0
    sim = rebound.Simulation(str(bin_path))
    rebx = reboundx.Extras(sim, str(rebx_path))
    saved = np.load(C.CHUNK_DIR / f"{sys_.name}.npz")
    n_blocks = round(sim.t / block_len)
    rows = [tuple(saved[c][i] for c in COLUMNS) for i in range(n_blocks)]
    return sim, rebx, rows, float(saved["wall_s"])


def block_time(n_blocks: int, block_len: float) -> float:
    """Start time of block n_blocks, by the same accumulation an uninterrupted run uses."""
    t = 0.0
    for _ in range(n_blocks):
        t += block_len
    return t


def run_session(sys_: C.System, t_end: float, budget_s: float, max_blocks: int | None = None) -> dict:
    """Run until t_end, the session wall budget, or max_blocks; checkpoint; return the status."""
    p_orb = C.p_orb_years(sys_)
    dt = p_orb / C.CFG["n_per_orbit"]
    n_samples = round(C.CFG["mean_window_p_out"] * C.p_out_years(sys_) / p_orb) * C.CFG["n_per_orbit"]
    block_len = n_samples * dt
    sim, rebx, rows, wall_prev = start_or_resume(sys_, block_len)
    t = block_time(len(rows), block_len)
    t0 = time.time()
    n_run = 0
    while t < t_end and time.time() - t0 < budget_s and (max_blocks is None or n_run < max_blocks):
        rows.append(block_row(t + 0.5 * block_len, sample_block(sim, t, n_samples, dt)))
        t += block_len
        n_run += 1
        if len(rows) % C.CFG["chunk_every_blocks"] == 0:
            save_checkpoint(sys_.name, sim, rebx, rows,
                            {"complete": False, "wall_s": wall_prev + time.time() - t0, "t_end_target": t_end})
    done = t >= t_end
    wall = wall_prev + time.time() - t0
    comp_inc = math.degrees(sim.particles[2].orbit().inc)
    save_checkpoint(sys_.name, sim, rebx, rows, {"complete": done, "wall_s": wall, "t_end_target": t_end,
                                                 "steps": sim.steps_done, "comp_inc_deg_final": comp_inc})
    return {"name": sys_.name, "complete": done, "t": t, "t_end": t_end, "wall": wall, "steps": sim.steps_done}


def run_system(args: tuple) -> str:
    sys_, t_end = args
    r = run_session(sys_, t_end, C.CFG["max_wall_hours"] * SECONDS_PER_HOUR)
    state = "complete" if r["complete"] else "paused"
    return f"{r['name']} complete={r['complete']} ({state}) t={r['t']:.3e}/{r['t_end']:.3e} wall={r['wall']:.0f}s steps={r['steps']}"


def target_end(sys_: C.System) -> float:
    """End time: the secular second e_max plus a margin, or the first e_max plus a margin if the run is too long."""
    g = C.secular_summary(sys_, "gr")
    if g["tmax2"] > C.CFG["two_cycle_max_yr"]:
        return g["tmax1"] + C.CFG["stop_margin_period"] * g["period"]
    return g["tmax2"] + C.CFG["stop_margin_period"] * g["period"]


def is_done(name: str) -> bool:
    """Complete, or one of the legacy capped runs that have no saved state (kept as partial results, never rerun)."""
    if name in C.CFG["legacy_capped_no_state"]:
        return True
    path = C.CHUNK_DIR / f"{name}.npz"
    return path.exists() and bool(np.load(path)["complete"])


PRIORITY = {"fiducial": 0, "extreme": 1, "random": 2}


def validate_gr() -> None:
    """Check the REBOUNDx GR precession of an isolated planet-star pair against 3 G (M*+Mp) n / (c^2 a j^2)."""
    sys_ = C.System("V", "check", C.CFG4["m_p_mjup"], C.CFG4["a_b_au"], 0.0, 0.0, 0.0)
    e0, span = 0.5, 20000.0
    sim, _rebx = build_sim(sys_, e0, m_b=1e-12)
    om0 = sim.particles[PLANET].orbit().pomega
    sim.integrate(span)
    d_om = (sim.particles[PLANET].orbit().pomega - om0 + math.pi) % (2 * math.pi) - math.pi
    mtot = C.CFG4["m_star_msun"] + sys_.m_p_mjup * S.MJUP_PER_MSUN
    a, j2 = C.CFG4["a0_au"], 1.0 - e0**2
    n = math.sqrt(S.G_AU3_PER_MSUN_YR2 * mtot / a**3)
    expect = 3.0 * S.G_AU3_PER_MSUN_YR2 * mtot / (C.CFG["c_au_per_yr"] ** 2 * a * j2) * n * span
    print(f"GR precession over {span:g} yr: N-body {d_om:.6e} rad, formula {expect:.6e} rad, ratio {d_om / expect:.5f}")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only", nargs="*")
    ap.add_argument("--validate-gr", action="store_true")
    args = ap.parse_args()
    if args.validate_gr:
        validate_gr()
        return
    if not C.SYSTEMS_PATH.exists():
        C.save_systems(C.build_systems())
    systems = [s for s in C.load_systems() if not args.only or s.name in args.only]
    jobs = [(s, target_end(s)) for s in systems if not is_done(s.name)]
    jobs.sort(key=lambda j: (PRIORITY[j[0].kind], -j[1]))
    while jobs:  # a session budget pauses a run; go round again until every run is complete
        with mp.Pool(C.CFG["n_workers"]) as pool:
            for msg in pool.imap_unordered(run_system, jobs):
                print(msg, flush=True)
        jobs = [j for j in jobs if not is_done(j[0].name)]


if __name__ == "__main__":
    main()
