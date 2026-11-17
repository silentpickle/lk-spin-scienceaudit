"""Step 4 scan: theta_sl,f vs P_star for one theta_lb0. Chunked, resumable, multiprocessing.

P_star values (log-spaced over [p_star_min, p_star_max]) are sorted and cut into groups of --group; each group is one
integration carrying a stellar spin per P_star (the orbit does not depend on P_star). Each group's result is saved to
results/step4/chunks/ as it finishes; re-running skips finished chunks. Merged into results/step4/scan_<theta>.npz
(pilot: scan_<theta>_pilot.npz).

Run:  uv run python scripts/step4_scan.py --theta 89 --n 300 --group 4
      uv run python scripts/step4_scan.py --pilot
"""

import argparse
import multiprocessing as mp
import time
from pathlib import Path

import numpy as np

from lkspin import secular as S

CFG = S.load_config()
OUT_DIR = S.ROOT / CFG["results"]
CHUNK_DIR = OUT_DIR / "chunks"
PILOT_THETA = 89.0
PILOT_N = 12
PILOT_GROUP = 1
DEFAULT_GROUP = 4
FIELDS = ("P_star", "theta_sl_f", "outcome", "t_end", "a_end", "e_end", "steps", "wall_s")
CHUNK_TAG_FMT = "{kind}_th{theta:g}_n{n}_g{group}_c{idx:04d}.npz"


def p_star_grid(n: int) -> np.ndarray:
    """Log-spaced P_star in days, descending (30 d first, as on the S16 axis)."""
    return np.geomspace(CFG["p_star_max_days"], CFG["p_star_min_days"], n)


def chunk_path(kind: str, theta: float, n: int, group: int, idx: int) -> Path:
    return CHUNK_DIR / CHUNK_TAG_FMT.format(kind=kind, theta=theta, n=n, group=group, idx=idx)


def run_group(theta: float, p_stars: np.ndarray) -> dict:
    """Integrate one group of P_star values together. Returns arrays keyed by FIELDS."""
    par = S.make_params(CFG)
    y0 = S.initial_state(par, CFG, theta, p_stars.size)
    t0 = time.time()
    r = S.run(y0, par, list(p_stars), CFG["stop_t_gyr"] * 1e9)
    wall = time.time() - t0
    _, j, a = S.orbit_elements(r["y"], par)
    k = p_stars.size
    return {"P_star": p_stars, "theta_sl_f": S.theta_sl_final(r["y"], k), "outcome": np.full(k, r["status"]),
            "t_end": np.full(k, r["t"]), "a_end": np.full(k, a), "e_end": np.full(k, np.sqrt(1 - j * j)),
            "steps": np.full(k, r["steps"]), "wall_s": np.full(k, wall)}


def work_item(args: tuple) -> str:
    kind, theta, n, group, idx, p_stars = args
    path = chunk_path(kind, theta, n, group, idx)
    if path.exists():
        return f"skip {path.name}"
    res = run_group(theta, np.asarray(p_stars))
    tmp = path.with_suffix(".tmp.npz")
    np.savez(tmp, **res)
    tmp.replace(path)
    return f"done {path.name} wall {res['wall_s'][0]:.1f}s"


def build_items(kind: str, theta: float, n: int, group: int) -> list:
    grid = p_star_grid(n)
    return [(kind, theta, n, group, i, grid[s:s + group].tolist())
            for i, s in enumerate(range(0, n, group))]


def merge(kind: str, theta: float, n: int, group: int, out_path: Path) -> dict:
    parts = [np.load(chunk_path(kind, theta, n, group, i)) for i in range(len(build_items(kind, theta, n, group)))]
    merged = {f: np.concatenate([p[f] for p in parts]) for f in FIELDS}
    np.savez(out_path, **merged)
    return merged


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--theta", type=float, default=PILOT_THETA)
    ap.add_argument("--n", type=int, default=CFG["n_p_star"])
    ap.add_argument("--group", type=int, default=DEFAULT_GROUP)
    ap.add_argument("--procs", type=int, default=mp.cpu_count())
    ap.add_argument("--pilot", action="store_true", help="theta 89, 12 P_star, one integration each; separate file")
    args = ap.parse_args()
    if args.pilot:
        args.theta, args.n, args.group = PILOT_THETA, PILOT_N, PILOT_GROUP
    kind = "pilot" if args.pilot else "scan"
    CHUNK_DIR.mkdir(parents=True, exist_ok=True)
    items = build_items(kind, args.theta, args.n, args.group)
    t0 = time.time()
    with mp.Pool(args.procs) as pool:
        for msg in pool.imap_unordered(work_item, items):
            print(msg, flush=True)
    suffix = "_pilot" if args.pilot else ""
    out = OUT_DIR / f"scan_{args.theta:g}{suffix}.npz"
    merged = merge(kind, args.theta, args.n, args.group, out)
    print(f"total wall {time.time() - t0:.1f}s; saved {out}")
    for i in range(args.n):
        print(f"{merged['P_star'][i]:9.4f} d  theta_sl_f {merged['theta_sl_f'][i]:8.3f}  "
              f"{S.STATUS_NAMES[int(merged['outcome'][i])]:12s} t_end {merged['t_end'][i] / 1e6:9.3f} Myr  "
              f"steps {int(merged['steps'][i]):9d}  wall {merged['wall_s'][i]:6.2f}s")


if __name__ == "__main__":
    main()
