"""Prep job 1 (SQ2): phase ensembles of class I systems (config [prep_sq2]).

Systems: 4 S16 fiducials + 16 step 5 class I (labels_v2) systems stratified in A_bar_0 (seed in config).
Ensembles per system (full [step4] physics: SRFs + weak-friction tides, no feedback, e_b = 0):
  A  P_star jitter: 200 members, P_star (1 + d), d uniform +-5%, one shared orbit, grouped spins
  B  orbital phase: 200 members, omega_0 uniform in [0, 2 pi), each its own orbit, P_star fixed
  W  supplementary: as A with a wide jitter (wide_jitter_max), for systems where +-5% does not populate both branches
  C  misaligned (fiducials only): theta_sl,0 in {10, 20, 30} deg (+ a 1 deg control), azimuth uniform, shared orbit
Chunks: results/prep_sq2/chunks/ (resumable). Analysis writes summary.json, theta_f_<system>.npz and prep_sq2_hist_<system>.png.
Run:  uv run python scripts/prep_sq2.py --stage run   then   --stage analyze
"""

import argparse
import json
import math
import multiprocessing as mp
import os
import time
from pathlib import Path

import matplotlib
import numpy as np

matplotlib.use("Agg")
import matplotlib.pyplot as plt

import prep_common as C
from lkspin import regime as R
from lkspin import secular as S

CFG4 = S.load_config()
CFG = C.load_table("prep_sq2")
OUT = S.ROOT / CFG["results"]
CHUNKS = OUT / "chunks"
STEP5 = S.ROOT / "results/step5"
DPI = 130
GYR = 1e9
CODE_A, CODE_B, CODE_C, CODE_W = 1, 2, 3, 4
BIN_DEG = 3.0
THETA_MAX = 180.0


# ------------------------------------------------------------------ systems
def pick_v2(rng: np.random.Generator) -> list[int]:
    """One class I system per log-spaced A_bar_0 bin (feature row indices); an empty bin takes the nearest system."""
    f, lab = np.load(STEP5 / "features.npz"), np.load(STEP5 / "labels_v2.npz")
    rows = np.where(lab["label"] == "I")[0]
    a = f["a_bar0"][rows]
    lo, hi = CFG["a_bar0_range"]
    edges = np.geomspace(lo, hi, CFG["n_v2_systems"] + 1)
    picks = []
    for k in range(CFG["n_v2_systems"]):
        inbin = rows[(a >= edges[k]) & (a < edges[k + 1])]
        if inbin.size:
            picks.append(int(rng.choice(inbin)))
            continue
        centre = math.sqrt(edges[k] * edges[k + 1])
        picks.append(int(rows[np.argmin(np.abs(np.log(a / centre)))]))
    return picks


def build_systems() -> list[dict]:
    systems = []
    for th, p in zip(CFG["fiducial_theta_lb0_deg"], CFG["fiducial_p_star_days"], strict=True):
        systems.append({"name": f"fid_{th:g}_P{p:g}", "m_p": CFG4["m_p_mjup"], "a_b": CFG4["a_b_au"], "theta": th, "P": p, "fiducial": True})
    f = np.load(STEP5 / "features.npz")
    for k, row in enumerate(pick_v2(np.random.default_rng(CFG["seed"]))):
        systems.append({"name": f"v2_{k:02d}", "m_p": float(f["m_p_mjup"][row]), "a_b": float(f["a_b_au"][row]),
                        "theta": float(f["theta_lb0_deg"][row]), "P": float(f["p_star_days"][row]), "fiducial": False,
                        "step5_index": int(f["index"][row]), "step5_a_bar0": float(f["a_bar0"][row])})
    return systems


def system_setup(sysd: dict) -> tuple[dict, np.ndarray]:
    cfg = R.system_config(CFG4, sysd["m_p"], sysd["a_b"])
    return cfg, S.make_params(cfg)


def rng_for(si: int, code: int, extra: int = 0) -> np.random.Generator:
    return np.random.default_rng([CFG["seed"], si, code, extra])


def save_atomic(path: Path, **arrays) -> None:
    tmp = path.with_suffix(f".{os.getpid()}.tmp.npz")
    np.savez(tmp, **arrays)
    tmp.replace(path)


# ------------------------------------------------------------------ tasks
def task_A(si: int, sysd: dict, jitter: float = 0.0, code: int = CODE_A) -> dict:
    """Ensemble A (jitter = config jitter_max) or the supplementary wide-jitter ensemble W (code CODE_W)."""
    cfg, par = system_setup(sysd)
    jitter = jitter or CFG["jitter_max"]
    d = rng_for(si, code).uniform(-jitter, jitter, CFG["n_members"])
    periods = sysd["P"] * (1.0 + d)
    y0 = S.initial_state(par, cfg, sysd["theta"], periods.size)
    r = S.run(y0, par, list(periods), cfg["stop_t_gyr"] * GYR)
    return {"theta_f": S.theta_sl_final(r["y"], periods.size), "P_star": periods, "status": r["status"], "t_end_yr": r["t"]}


def task_C(si: int, sysd: dict, ang_idx: int, angle: float) -> dict:
    cfg, par = system_setup(sysd)
    n = CFG["n_members"]
    az = rng_for(si, CODE_C, ang_idx).uniform(0.0, C.TWO_PI, n)
    y0 = S.initial_state(par, cfg, sysd["theta"], n)
    lhat = y0[0:3] / np.linalg.norm(y0[0:3])
    for k in range(n):
        o = S.N_ORBIT + S.N_SPIN * k
        y0[o:o + S.N_SPIN] = C.tilted_spin_vector(lhat, angle, az[k])
    r = S.run(y0, par, [sysd["P"]] * n, cfg["stop_t_gyr"] * GYR)
    return {"theta_f": S.theta_sl_final(r["y"], n), "azimuth": az, "status": r["status"], "t_end_yr": r["t"], "theta_sl0": angle}


def member_B(par: np.ndarray, cfg: dict, theta: float, p_star: float, omega0: float) -> tuple[float, int, float]:
    """One member: e0 vector (length e0) rotated by omega0 about L_hat in the orbit plane, own orbit."""
    y0 = S.initial_state(par, cfg, theta, 1)
    lhat = y0[0:3] / np.linalg.norm(y0[0:3])
    xp = np.array([1.0, 0.0, 0.0])
    yp = np.cross(lhat, xp)
    y0[3:6] = CFG["e0_vector_amp"] * (math.cos(omega0) * xp + math.sin(omega0) * yp)
    r = S.run(y0, par, [p_star], cfg["stop_t_gyr"] * GYR)
    return float(S.theta_sl_final(r["y"], 1)[0]), int(r["status"]), float(r["t"])


def task_B(si: int, sysd: dict, chunk: int) -> dict:
    cfg, par = system_setup(sysd)
    lo, hi = chunk * CFG["chunk_members"], (chunk + 1) * CFG["chunk_members"]
    omegas = rng_for(si, CODE_B).uniform(0.0, C.TWO_PI, CFG["n_members"])[lo:hi]
    res = [member_B(par, cfg, sysd["theta"], sysd["P"], w) for w in omegas]
    return {"theta_f": np.array([r[0] for r in res]), "status": np.array([r[1] for r in res]),
            "t_end_yr": np.array([r[2] for r in res]), "omega0": omegas}


def tasks(systems: list[dict]) -> list[tuple]:
    out = []
    n_chunks = CFG["n_members"] // CFG["chunk_members"]
    angles = CFG["misalign_deg"] + [CFG["control_misalign_deg"]]
    for si, s in enumerate(systems):
        out.append(("A", si, 0))
        out.append(("W", si, 0))
        out += [("B", si, c) for c in range(n_chunks)]
        if s["fiducial"]:
            out += [("C", si, k) for k in range(len(angles))]
    return out


def chunk_file(systems: list[dict], kind: str, si: int, idx: int) -> Path:
    tag = {"A": "A", "W": "W", "B": f"B_c{idx:02d}", "C": f"C_a{idx}"}[kind]
    return CHUNKS / f"{systems[si]['name']}_{tag}.npz"


def run_task(args: tuple) -> str:
    kind, si, idx = args
    systems = build_systems()
    path = chunk_file(systems, kind, si, idx)
    if path.exists():
        return f"skip {path.name}"
    t0 = time.time()
    s = systems[si]
    if kind == "A":
        res = task_A(si, s)
    elif kind == "W":
        res = task_A(si, s, CFG["wide_jitter_max"], CODE_W)
    elif kind == "B":
        res = task_B(si, s, idx)
    else:
        res = task_C(si, s, idx, (CFG["misalign_deg"] + [CFG["control_misalign_deg"]])[idx])
    save_atomic(path, wall_s=time.time() - t0, **res)
    return f"done {path.name} {time.time() - t0:.0f}s"


def run_stage(procs: int) -> None:
    CHUNKS.mkdir(parents=True, exist_ok=True)
    systems = build_systems()
    items = tasks(systems)
    print(f"{len(systems)} systems, {len(items)} tasks", flush=True)
    with mp.Pool(procs) as pool:
        for msg in pool.imap_unordered(run_task, items, chunksize=1):
            print(msg, flush=True)


# ------------------------------------------------------------------ analysis
def load_ensemble(systems: list[dict], si: int, kind: str, idx: int = 0) -> dict:
    if kind == "B":
        n_chunks = CFG["n_members"] // CFG["chunk_members"]
        parts = [np.load(chunk_file(systems, "B", si, c)) for c in range(n_chunks)]
        return {k: np.concatenate([p[k] for p in parts]) for k in ("theta_f", "status", "t_end_yr", "omega0")}
    d = np.load(chunk_file(systems, kind, si, idx))
    out = {k: d[k] for k in d.files}
    out["status"] = np.full(out["theta_f"].size, int(d["status"]))
    out["t_end_yr"] = np.full(out["theta_f"].size, float(d["t_end_yr"]))
    return out


def stats_for(theta: np.ndarray, rng: np.random.Generator) -> dict:
    """Branch statistics + bootstrap CIs of one set of theta_f values (hot-Jupiter members only)."""
    bs = C.branch_summary(theta, CFG["tight_deg"], CFG["gap_deg"])
    two = bs["gap_deg"] >= CFG["gap_deg"]
    boot = C.bootstrap_branches(theta, bs["cut_deg"], CFG["n_boot"], rng, CFG["ci_percent"])
    out = {"n": int(theta.size), "two_branches": bool(two), "n_clusters_rule": bs["n_branches"], "largest_gap_deg": bs["gap_deg"],
           "cut_deg": bs["cut_deg"], "range_deg": bs["range_deg"], "median_all_deg": bs["median_all_deg"], "mean_deg": float(theta.mean()),
           "std_deg": float(theta.std())}
    if two:
        out.update({"p_high": bs["p_high"], "peak_low_deg": bs["peak_low_deg"], "peak_high_deg": bs["peak_high_deg"], **boot})
    out["modes"] = C.modes(theta, CFG["mode_gap_deg"], CFG["mode_min_frac"])
    out["third_mode"] = bool(len(out["modes"]) > 2)
    return out


def ensemble_report(e: dict, rng: np.random.Generator) -> dict:
    ok = e["status"] == S.ST_HOT_JUPITER
    theta = e["theta_f"][ok]
    rep = {"n_members": int(e["theta_f"].size), "n_hot_jupiter": int(ok.sum()),
           "status_counts": {S.STATUS_NAMES[s]: int((e["status"] == s).sum()) for s in np.unique(e["status"])}}
    rep["t_hj_Myr_median"] = float(np.median(e["t_end_yr"][ok]) / 1e6) if ok.any() else math.nan
    rep["t_hj_Myr_range"] = [float(e["t_end_yr"][ok].min() / 1e6), float(e["t_end_yr"][ok].max() / 1e6)] if ok.any() else None
    if ok.sum() < CFG["convergence_sizes"][0]:
        return rep
    rep["full"] = stats_for(theta, rng)
    rep["convergence"] = {str(n): stats_for(theta[:n], rng) for n in CFG["convergence_sizes"] if n <= theta.size}
    return rep


def system_info(sysd: dict) -> dict:
    cfg, par = system_setup(sysd)
    shape = C.first_cycle_shape(par, cfg, sysd["theta"], sysd["P"])
    shape["omega_hat_l_check_n_over_2a"] = shape["n_max0"] / (2.0 * shape["a_bar0"])
    shape["lk_period_Myr"] = shape.pop("period_yr") / 1e6
    return shape


def plot_hist(name: str, ens: dict[str, np.ndarray]) -> None:
    n = len(ens)
    cols = min(n, 3)
    rows = math.ceil(n / cols)
    fig, axes = plt.subplots(rows, cols, figsize=(5 * cols, 3.2 * rows), squeeze=False)
    bins = np.arange(0.0, THETA_MAX + BIN_DEG, BIN_DEG)
    for ax, (label, th) in zip(axes.ravel(), ens.items(), strict=False):
        ax.hist(th[np.isfinite(th)], bins=bins, color="k")
        ax.set_title(label, fontsize=9)
        ax.set_xlabel(r"$\theta_{\rm sl,f}$ (deg)")
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    fig.suptitle(name)
    fig.tight_layout()
    fig.savefig(OUT / f"prep_sq2_hist_{name}.png", dpi=DPI)
    plt.close(fig)


def analyze_system(systems: list[dict], si: int) -> dict:
    sysd = systems[si]
    rng = np.random.default_rng(CFG["boot_seed"] + si)
    angles = CFG["misalign_deg"] + [CFG["control_misalign_deg"]]
    raw = {"A": load_ensemble(systems, si, "A"), "W": load_ensemble(systems, si, "W"), "B": load_ensemble(systems, si, "B")}
    if sysd["fiducial"]:
        raw.update({f"C_{angles[k]:g}": load_ensemble(systems, si, "C", k) for k in range(len(angles))})
    rep = {key: ensemble_report(e, rng) for key, e in raw.items()}
    np.savez(OUT / f"theta_f_{sysd['name']}.npz", **{f"{k}_theta_f": e["theta_f"] for k, e in raw.items()},
             **{f"{k}_status": e["status"] for k, e in raw.items()})
    labels = {"A: P_star jitter (+-5%)": raw["A"]["theta_f"], "B: omega_0 phase": raw["B"]["theta_f"], "W: P_star jitter (wide)": raw["W"]["theta_f"]}
    labels.update({f"C: theta_sl,0 = {angles[k]:g} deg": raw[f"C_{angles[k]:g}"]["theta_f"] for k in range(len(angles)) if sysd["fiducial"]})
    plot_hist(sysd["name"], {k: v[np.isfinite(v)] for k, v in labels.items()})
    return {"inputs": sysd, "orbit_shape": system_info(sysd), "ensembles": rep, "shifts_vs_A": shifts(rep)}


def shifts(rep: dict) -> dict:
    """Peak shifts of each C ensemble (and B) relative to ensemble A, where both are two-branch."""
    out = {}
    base = rep["A"].get("full")
    for key, e in rep.items():
        full = e.get("full")
        if key == "A" or full is None or base is None:
            continue
        out[key] = {f: (full[f] - base[f]) for f in ("peak_low_deg", "peak_high_deg", "p_high") if f in full and f in base}
        out[key]["median_all_shift_deg"] = full["median_all_deg"] - base["median_all_deg"]
    return out


def analyze() -> None:
    systems = build_systems()
    result = {"seed": CFG["seed"], "boot_seed": CFG["boot_seed"], "config": CFG,
              "systems": {s["name"]: analyze_system(systems, si) for si, s in enumerate(systems)}}
    (OUT / "summary.json").write_text(json.dumps(result, indent=1, default=float))
    for name, r in result["systems"].items():
        a = r["ensembles"]["A"].get("full", {})
        b = r["ensembles"]["B"].get("full", {})
        print(name, "A:", {k: round(a[k], 3) for k in ("p_high", "peak_low_deg", "peak_high_deg") if k in a}, "B:",
              {k: round(b[k], 3) for k in ("p_high", "peak_low_deg", "peak_high_deg") if k in b})


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--stage", default="all", choices=["run", "analyze", "all"])
    ap.add_argument("--procs", type=int, default=CFG["n_procs"])
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    if args.stage in ("run", "all"):
        run_stage(min(args.procs, CFG["n_procs"]))
    if args.stage in ("analyze", "all"):
        analyze()


if __name__ == "__main__":
    main()
