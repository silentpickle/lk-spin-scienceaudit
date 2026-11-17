"""Step 5: integrate the sample, one phase ensemble per system, and merge features.

Draws n_sample systems (config [step5], seed recorded), integrates each (K + 1 spins on one orbit, with winding and
eccentricity-extremum diagnostics), saves one chunk per system in results/step5/chunks/ (resumable), and merges:
  results/step5/systems_all.npz   every drawn system: inputs, status, t_end, wall time
  results/step5/features.npz      one row per hot-Jupiter system: all features, first-pass and predicted class
  results/step5/series.npz        compressed raw per-cycle series (member 0, middle member, twin) per HJ system
Run:  uv run python scripts/step5_run.py [--n N] [--procs P]      (--n limits the sample to its first N draws)
"""

import argparse
import multiprocessing as mp
import time
import tomllib
from pathlib import Path

import numpy as np

from lkspin import regime as R
from lkspin import secular as S

CFG4 = S.load_config()
with S.CONFIG_PATH.open("rb") as _fh:
    CFG5 = tomllib.load(_fh)["step5"]
OUT_DIR = S.ROOT / CFG5["results"]
CHUNK_DIR = OUT_DIR / "chunks"
INPUT_FIELDS = ("m_p_mjup", "p_star_days", "a_b_au", "theta_lb0_deg")
SERIES_NOTE = "columns: t_yr, a_au, e_min, amplitude, then (cos theta_sl, Psi) for member 0, middle member, twin"
STRING_FIELDS = ("label", "predicted")


def chunk_path(i: int) -> Path:
    return CHUNK_DIR / f"sys_{i:05d}.npz"


def work_item(args: tuple) -> str:
    i, m_p, p_star, a_b, theta = args
    path = chunk_path(i)
    if path.exists():
        return f"skip {i}"
    t0 = time.time()
    out = R.integrate_system(m_p, p_star, a_b, theta, CFG4, CFG5)
    out["wall_s"] = time.time() - t0
    tmp = path.with_suffix(".tmp.npz")
    np.savez(tmp, **out)
    tmp.replace(path)
    return f"done {i} {S.STATUS_NAMES[int(out['status'])]} {out['wall_s']:.1f}s"


def build_items(sample: dict, n: int) -> list:
    return [(i, *(float(sample[k][i]) for k in INPUT_FIELDS)) for i in range(n)]


def load_chunks(n: int) -> list:
    return [dict(np.load(chunk_path(i))) for i in range(n)]


def stack_scalars(rows: list, keys: list) -> dict:
    return {k: np.array([r[k] for r in rows]) for k in keys}


def hj_table(sample: dict, chunks: list) -> tuple[dict, dict]:
    """Feature table (one row per HJ) and the series dict keyed s<index>."""
    hj = [i for i, c in enumerate(chunks) if int(c["status"]) == S.ST_HOT_JUPITER]
    rows = [chunks[i] for i in hj]
    skip = {"series", "status", "steps", "n_events", "wall_s", "t_end"}
    keys = [k for k in rows[0] if k not in skip]
    table = stack_scalars(rows, keys)
    table.update({"index": np.array(hj), "t_hj_yr": np.array([r["t_end"] for r in rows]),
                  "wall_s": np.array([r["wall_s"] for r in rows]), "n_events": np.array([r["n_events"] for r in rows])})
    table.update({k: sample[k][hj] for k in INPUT_FIELDS})
    series = {f"s{i}": chunks[i]["series"] for i in hj}
    return table, series


def add_labels(table: dict) -> None:
    table["label"] = np.array(R.classify_table(table, CFG5))
    table["predicted"] = np.array(R.predicted_table(table, CFG5))


def merge(sample: dict, n: int, wall_total: float) -> dict:
    chunks = load_chunks(n)
    status = np.array([int(c["status"]) for c in chunks])
    allsys = {k: sample[k][:n] for k in INPUT_FIELDS}
    allsys.update({"status": status, "t_end_yr": np.array([c["t_end"] for c in chunks]),
                   "wall_s": np.array([c["wall_s"] for c in chunks]), "seed": CFG5["seed"]})
    np.savez(OUT_DIR / "systems_all.npz", **allsys)
    table, series = hj_table(sample, chunks)
    add_labels(table)
    table["seed"] = CFG5["seed"]
    table["period_offsets"] = R.period_offsets(CFG5)
    table["middle_index"] = R.middle_index(R.period_offsets(CFG5))
    np.savez(OUT_DIR / "features.npz", **table)
    series["note"] = np.array(SERIES_NOTE)
    np.savez_compressed(OUT_DIR / "series.npz", **series)
    return {"status": status, "wall_each": allsys["wall_s"], "wall_total": wall_total, "n_hj": len(table["index"])}


def report(info: dict, n: int) -> None:
    status = info["status"]
    print(f"systems drawn {n}: hot Jupiter {info['n_hj']}, "
          + ", ".join(f"{S.STATUS_NAMES[s]} {int((status == s).sum())}" for s in range(S.ST_STEPS + 1) if s != S.ST_HOT_JUPITER))
    print(f"wall per system: mean {info['wall_each'].mean():.2f} s, max {info['wall_each'].max():.1f} s; "
          f"sum {info['wall_each'].sum():.0f} s CPU; this run {info['wall_total']:.0f} s wall")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--n", type=int, default=CFG5["n_sample"])
    ap.add_argument("--procs", type=int, default=CFG5["n_procs"])
    args = ap.parse_args()
    procs = min(args.procs, CFG5["n_procs"])
    CHUNK_DIR.mkdir(parents=True, exist_ok=True)
    sample = R.draw_sample(CFG5)
    items = build_items(sample, args.n)
    t0 = time.time()
    with mp.Pool(procs) as pool:
        for k, msg in enumerate(pool.imap_unordered(work_item, items, chunksize=1)):
            if k % 100 == 0:
                print(k, msg, flush=True)
    report(merge(sample, args.n, time.time() - t0), args.n)


if __name__ == "__main__":
    main()
