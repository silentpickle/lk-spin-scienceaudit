"""Step 11 N-body runs: star + planet + companion, IAS15 + REBOUNDx "gr", block-averaged output.

Output per system: results/step11/chunks/<name>.npz with one row per averaging block (one outer binary period):
t (block midpoint, yr), e_mean, e_max (osculating, within the block), e_min, theta_mean (deg), a_mean (AU).
The planet's orbit is Jacobi relative to the star (REBOUND's Particle.orbit() default, Jacobi COM of the lower indices).
Run:  uv run python scripts/step11_nbody.py [--only NAME ...] [--validate-gr]
"""

import argparse
import math
import multiprocessing as mp
import time

import numpy as np
import rebound
import reboundx

import step11_common as C
from lkspin import secular as S

SECONDS_PER_HOUR = 3600.0
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


def save_chunk(name: str, rows: list, extra: dict) -> None:
    C.CHUNK_DIR.mkdir(parents=True, exist_ok=True)
    arr = np.array(rows).reshape(-1, len(COLUMNS))
    tmp = C.CHUNK_DIR / f"{name}.tmp.npz"
    np.savez(tmp, **dict(zip(COLUMNS, arr.T, strict=True)), **{k: np.array(v) for k, v in extra.items()})
    tmp.replace(C.CHUNK_DIR / f"{name}.npz")


def run_system(args: tuple) -> str:
    sys_, t_end = args
    sim, _rebx = build_sim(sys_, C.CFG4["e0"])
    p_orb = C.p_orb_years(sys_)
    dt = p_orb / C.CFG["n_per_orbit"]
    orbits_per_block = round(C.CFG["mean_window_p_out"] * C.p_out_years(sys_) / p_orb)
    n_samples = orbits_per_block * C.CFG["n_per_orbit"]
    block_len = n_samples * dt
    wall_cap = C.CFG["max_wall_hours"] * SECONDS_PER_HOUR
    t0 = time.time()
    rows: list = []
    t = 0.0
    while t < t_end and time.time() - t0 < wall_cap:
        rows.append(block_row(t + 0.5 * block_len, sample_block(sim, t, n_samples, dt)))
        t += block_len
        if len(rows) % C.CFG["chunk_every_blocks"] == 0:
            save_chunk(sys_.name, rows, {"complete": False, "wall_s": time.time() - t0, "t_end_target": t_end})
    done = t >= t_end
    comp_inc = math.degrees(sim.particles[2].orbit().inc)
    save_chunk(sys_.name, rows, {"complete": done, "wall_s": time.time() - t0, "t_end_target": t_end,
                                 "steps": sim.steps_done, "comp_inc_deg_final": comp_inc})
    return f"{sys_.name} complete={done} t={t:.3e}/{t_end:.3e} wall={time.time() - t0:.0f}s steps={sim.steps_done}"


def target_end(sys_: C.System) -> float:
    g = C.secular_summary(sys_, "gr")
    return g["tmax2"] + C.CFG["stop_margin_period"] * g["period"]


def is_done(name: str) -> bool:
    path = C.CHUNK_DIR / f"{name}.npz"
    return path.exists() and bool(np.load(path)["complete"])


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
    jobs.sort(key=lambda j: -j[1])
    with mp.Pool(C.CFG["n_workers"]) as pool:
        for msg in pool.imap_unordered(run_system, jobs):
            print(msg, flush=True)


if __name__ == "__main__":
    main()
