"""Step 6 (SQ1 cases b-d): re-run the step 5 hot-Jupiter systems with spin-orbit feedback, octupole and spin-down.

Case b: feedback.  Case c: feedback + octupole (e_b and the azimuth of e_b-hat drawn per system, seed in [step6]).
Case d: feedback + octupole + magnetic braking.  Each system = 10 phase-ensemble members + the twin of the middle member,
every spin on its OWN orbit (feedback makes the orbit depend on the spin): 11 integrations of lkspin.secular.run_ext.

Per case writes (results/step6/<case>/, or <case>_pilot/ with --pilot N):
  chunks/sys_<i>.npz      one per system (resumable)
  systems_all.npz         every system run: inputs, e_b, azimuth, per-member status, outcome code, wall time
  features.npz            one row per system whose 11 integrations all formed a hot Jupiter; the step 5 columns, in the
                          same format. a_bar0, n_max0, s_over_l and p_mid_days are COPIED from results/step5 (S16's predicted
                          class is unchanged by construction); twin_div_deg is the angle between the final inertial spin
                          vectors of the twin and the middle member (the two live on different orbits now)
  series.npz              s<index> = member-0 series in the step 5 format (t, a, e_min, amplitude, p, Psi of member 0, then
                          p, Psi of the middle member and of the twin at THEIR OWN k-th minima, NaN padded)
  series_members.npz      s<index> = (11, rows, 6) per-member (t, a, e_min, amplitude, p, Psi), each on its own orbit
Run:  uv run python scripts/step6_run.py --case b [--pilot 20] [--procs 10]
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
    _ALL = tomllib.load(_fh)
CFG5 = _ALL["step5"]
CFG6 = _ALL["step6"]
STEP5_DIR = S.ROOT / CFG5["results"]
OUT_DIR = S.ROOT / CFG6["results"]
CASES = {"b": {"octupole": False, "spindown": False}, "c": {"octupole": True, "spindown": False},
         "d": {"octupole": True, "spindown": True}}
OUT_MIXED = S.ST_STEPS + 1  # members of one system ended differently
OUTCOME_NAMES = (*S.STATUS_NAMES, "mixed")
INPUT_FIELDS = ("m_p_mjup", "p_star_days", "a_b_au", "theta_lb0_deg")
COPIED_FIELDS = ("a_bar0", "n_max0", "s_over_l", "p_mid_days")
SERIES_COLS = 6
SERIES_NOTE = ("member-0 columns t_yr, a_au, e_min, amplitude, p0, Psi0 on member 0's own orbit, then (cos theta_sl, Psi) of "
               "the middle member and the twin at their own k-th e minimum (NaN padded); see series_members.npz")
BOOK_KEYS = ("status", "t_end", "steps", "n_events")


def step5_table() -> dict:
    return dict(np.load(STEP5_DIR / "features.npz"))


def draw_binary(n_sample: int) -> tuple[np.ndarray, np.ndarray]:
    """(e_b, azimuth) for every step 5 sample index: one default_rng stream, e_b first, then the azimuth."""
    rng = np.random.default_rng(CFG6["seed"])
    return rng.uniform(*CFG6["e_b_range"], n_sample), rng.uniform(*CFG6["phi_b_range_rad"], n_sample)


def case_params(case: str, par: np.ndarray, e_b: float, phi_b: float) -> np.ndarray:
    flags = CASES[case]
    return S.extend_params(par, fb_mult=CFG6["fb_mult"], e_b=e_b, phi_b=phi_b, octupole=flags["octupole"],
                           alpha_mb=CFG6["alpha_mb_yr"] if flags["spindown"] else 0.0)


def twin_ext(y0: np.ndarray, angle_rad: float) -> np.ndarray:
    """y0 with S_star rotated by angle_rad away from L_hat (toward +x, perpendicular to L_hat0), magnitude kept."""
    y = y0.copy()
    lhat = y[0:3] / np.linalg.norm(y[0:3])
    smag = np.linalg.norm(y[S.N_ORBIT:S.N_ORBIT + S.N_SPIN])
    perp = np.array([1.0, 0.0, 0.0])
    y[S.N_ORBIT:S.N_ORBIT + S.N_SPIN] = smag * (np.cos(angle_rad) * lhat + np.sin(angle_rad) * perp)
    return y


def member_periods(p_star: float) -> tuple[list[float], int]:
    """P_star of the 10 members and the twin (the twin repeats the middle member's period), and the middle index."""
    offsets = R.period_offsets(CFG5)
    mid = R.middle_index(offsets)
    return list(p_star * (1.0 + offsets)) + [p_star * (1.0 + offsets[mid])], mid


def run_member(par: np.ndarray, cfg: dict, theta_lb0: float, p_days: float, twin: bool) -> dict:
    y0 = S.initial_state_ext(par, cfg, theta_lb0, p_days)
    y_init = twin_ext(y0, CFG5["twin_angle_rad"]) if twin else y0
    res = S.run_ext(y_init, par, cfg["stop_t_gyr"] * 1e9, ev_max=CFG5["ev_max"])
    return res


def outcome_of(status: np.ndarray) -> int:
    return int(status[0]) if np.all(status == status[0]) else OUT_MIXED


def pad_series(arrs: list[np.ndarray]) -> np.ndarray:
    rows = max(a.shape[0] for a in arrs)
    out = np.full((len(arrs), rows, SERIES_COLS), np.nan)
    for k, a in enumerate(arrs):
        out[k, :a.shape[0]] = a
    return out


def aligned_pairs(series: np.ndarray, n_rows: int) -> np.ndarray:
    """(p, Psi) columns of a member series cut or NaN-padded to n_rows."""
    out = np.full((n_rows, 2), np.nan)
    m = min(n_rows, series.shape[0])
    out[:m] = series[:m, 4:6]
    return out


def compat_series(per_member: list[np.ndarray], mid: int, twin: int) -> np.ndarray:
    s0 = per_member[0]
    return np.hstack([s0, aligned_pairs(per_member[mid], s0.shape[0]), aligned_pairs(per_member[twin], s0.shape[0])])


def hot_jupiter_features(results: list[dict]) -> dict:
    n_members = CFG5["n_members"]
    mid = R.middle_index(R.period_offsets(CFG5))
    twin = n_members
    theta_all = np.array([S.theta_sl_ext(r["y"]) for r in results])
    spins = [r["y"][S.N_ORBIT:S.N_ORBIT + S.N_SPIN] for r in results]
    feats = {"theta_f": theta_all[:n_members], "theta_f_twin": float(theta_all[twin]),
             "twin_div_deg": R.angle_between_deg(spins[twin], spins[mid]),
             "twin_div_m0_deg": R.angle_between_deg(spins[twin], spins[0]), "theta_max": float(theta_all[:n_members].max())}
    stats = R.ensemble_stats(theta_all[:n_members])
    feats.update({k: stats[k] for k in ("range", "gap", "n_low", "n_high", "range_low", "range_high")})
    feats["theta_sorted"] = stats["theta_sorted"]
    minima, maxima = R.split_events(results[0]["ev"])
    feats.update({f"m0_{k}": v for k, v in R.member_features(minima, maxima, 0, 1, CFG5).items()})
    return feats


def member_series(res: dict) -> np.ndarray:
    minima, maxima = R.split_events(res["ev"])
    return R.series_columns(minima, maxima, (0,), 1)


def integrate_system(case: str, m_p: float, p_star: float, a_b: float, theta_lb0: float, e_b: float, phi_b: float) -> dict:
    cfg = R.system_config(CFG4, m_p, a_b)
    par = case_params(case, S.make_params(cfg), e_b, phi_b)
    periods, mid = member_periods(p_star)
    results, walls, cpus = [], [], []
    for k, p_days in enumerate(periods):
        t0, c0 = time.time(), time.process_time()
        results.append(run_member(par, cfg, theta_lb0, p_days, twin=k == len(periods) - 1))
        walls.append(time.time() - t0)
        cpus.append(time.process_time() - c0)
    status = np.array([r["status"] for r in results])
    out = {"status": status, "outcome": outcome_of(status), "t_end": np.array([r["t"] for r in results]),
           "steps": np.array([r["steps"] for r in results]), "n_events": np.array([r["ev"].shape[0] for r in results]),
           "wall_members": np.array(walls), "wall_s": float(sum(walls)), "cpu_s": float(sum(cpus)), "e_b": e_b, "phi_b": phi_b,
           "eps_oct0": float(par[S.P_EOCT] * cfg["a0_au"])}
    if out["outcome"] != S.ST_HOT_JUPITER:
        return out
    out.update(hot_jupiter_features(results))
    per_member = [member_series(r) for r in results]
    out["series"] = compat_series(per_member, mid, len(periods) - 1)
    out["series_members"] = pad_series(per_member)
    return out


def chunk_dir(case: str, pilot: int) -> Path:
    return OUT_DIR / (f"{case}_pilot" if pilot else case) / "chunks"


def chunk_path(case: str, pilot: int, i: int) -> Path:
    return chunk_dir(case, pilot) / f"sys_{i:05d}.npz"


def work_item(args: tuple) -> str:
    case, pilot, i, m_p, p_star, a_b, theta, e_b, phi_b = args
    path = chunk_path(case, pilot, i)
    if path.exists():
        return f"skip {i}"
    out = integrate_system(case, m_p, p_star, a_b, theta, e_b, phi_b)
    tmp = path.with_suffix(".tmp.npz")
    np.savez(tmp, **out)
    tmp.replace(path)
    return f"done {i} {OUTCOME_NAMES[out['outcome']]} {out['wall_s']:.1f}s"


def build_items(case: str, pilot: int, table: dict, n: int) -> list:
    e_b_all, phi_all = draw_binary(CFG5["n_sample"])
    flags = CASES[case]
    items = []
    for i in range(n):
        idx = int(table["index"][i])
        e_b = float(e_b_all[idx]) if flags["octupole"] else 0.0
        phi_b = float(phi_all[idx]) if flags["octupole"] else 0.0
        items.append((case, pilot, i, *(float(table[k][i]) for k in INPUT_FIELDS), e_b, phi_b))
    return items


def merge(case: str, pilot: int, table: dict, n: int, wall_total: float) -> dict:
    chunks = [dict(np.load(chunk_path(case, pilot, i))) for i in range(n)]
    out_dir = OUT_DIR / (f"{case}_pilot" if pilot else case)
    outcome = np.array([int(c["outcome"]) for c in chunks])
    allsys = {k: table[k][:n] for k in INPUT_FIELDS}
    allsys.update({"index": table["index"][:n], "e_b": np.array([c["e_b"] for c in chunks]),
                   "phi_b": np.array([c["phi_b"] for c in chunks]), "eps_oct0": np.array([c["eps_oct0"] for c in chunks]),
                   "status_members": np.array([c["status"] for c in chunks]), "outcome": outcome,
                   "t_end_members": np.array([c["t_end"] for c in chunks]),
                   "wall_s": np.array([c["wall_s"] for c in chunks]), "cpu_s": np.array([c["cpu_s"] for c in chunks]),
                   "seed": CFG6["seed"]})
    np.savez(out_dir / "systems_all.npz", **allsys)
    hj = [i for i, c in enumerate(chunks) if int(c["outcome"]) == S.ST_HOT_JUPITER]
    write_features(out_dir, table, chunks, hj)
    return {"outcome": outcome, "wall_each": allsys["wall_s"], "cpu_each": allsys["cpu_s"], "wall_total": wall_total}


def write_features(out_dir: Path, table: dict, chunks: list, hj: list) -> None:
    skip = {"series", "series_members", "status", "outcome", "t_end", "steps", "n_events", "wall_members", "wall_s", "cpu_s"}
    rows = [chunks[i] for i in hj]
    keys = [k for k in rows[0] if k not in skip]
    feats = {k: np.array([r[k] for r in rows]) for k in keys}
    feats.update({"index": table["index"][hj], "t_hj_yr": np.array([r["t_end"].mean() for r in rows]),
                  "wall_s": np.array([r["wall_s"] for r in rows]), "cpu_s": np.array([r["cpu_s"] for r in rows]), "n_events": np.array([r["n_events"][0] for r in rows])})
    feats.update({k: table[k][hj] for k in (*INPUT_FIELDS, *COPIED_FIELDS)})
    feats["label"] = np.array(R.classify_table(feats, CFG5))
    feats["predicted"] = np.array(R.predicted_table(feats, CFG5))
    feats["seed"] = CFG6["seed"]
    feats["period_offsets"] = R.period_offsets(CFG5)
    feats["middle_index"] = R.middle_index(R.period_offsets(CFG5))
    np.savez(out_dir / "features.npz", **feats)
    series = {f"s{int(table['index'][i])}": chunks[i]["series"] for i in hj}
    series["note"] = np.array(SERIES_NOTE)
    np.savez_compressed(out_dir / "series.npz", **series)
    np.savez_compressed(out_dir / "series_members.npz",
                        **{f"s{int(table['index'][i])}": chunks[i]["series_members"] for i in hj})


def report(case: str, info: dict, n: int) -> None:
    out = info["outcome"]
    print(f"case {case}: {n} systems: " + ", ".join(f"{OUTCOME_NAMES[k]} {int((out == k).sum())}"
                                                    for k in range(len(OUTCOME_NAMES)) if (out == k).any()))
    print(f"per system: mean wall {info['wall_each'].mean():.2f} s, mean CPU {info['cpu_each'].mean():.2f} s, "
          f"max CPU {info['cpu_each'].max():.1f} s; sum CPU {info['cpu_each'].sum():.0f} s; this run {info['wall_total']:.0f} s wall")


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--case", choices=sorted(CASES), required=True)
    ap.add_argument("--pilot", type=int, default=0, help="run only the first N systems, into results/step6/<case>_pilot")
    ap.add_argument("--procs", type=int, default=CFG6["n_procs"])
    args = ap.parse_args()
    table = step5_table()
    n = args.pilot if args.pilot else table["index"].size
    chunk_dir(args.case, args.pilot).mkdir(parents=True, exist_ok=True)
    items = build_items(args.case, args.pilot, table, n)
    t0 = time.time()
    with mp.Pool(min(args.procs, CFG6["n_procs"])) as pool:
        for k, msg in enumerate(pool.imap_unordered(work_item, items, chunksize=1)):
            if k % 50 == 0:
                print(k, msg, flush=True)
    report(args.case, merge(args.case, args.pilot, table, n, time.time() - t0), n)


if __name__ == "__main__":
    main()
