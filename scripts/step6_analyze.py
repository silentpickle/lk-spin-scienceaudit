"""Step 6 analysis (SQ1): for each of the cases (b, c, d) found in results/step6, and for case (a) from results/step5:
the v2 observed classes, the confusion matrix of S16's predicted class (unchanged, from the step 5 A_bar_0 and N_max,0) against
the observed class, the accuracy with bootstrap intervals, the theta_lb0-only and "map + Omega_hat_L cut" baselines (as
scripts/step5_baselines.py: fitted on the same random half of the 1028 step 5 systems, scored on the other half), the class
transition matrix from (a) to the case (paired by system), the misclassified systems in the (A_bar_0, N_max,0) plane
(results/step6/<case>_map.png) and a summary JSON (results/step6/summary.json).
Systems whose 11 integrations do not all form a hot Jupiter are excluded from classification and counted per case.
Run:  uv run python scripts/step6_analyze.py [--pilot]      (--pilot reads results/step6/<case>_pilot and writes *_pilot files)
"""

import argparse
import json
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import step5_baselines as B
import step5_classify as C5
import step6_run as W

CFG6 = W.CFG6
OUT = W.OUT_DIR
FAMILIES = C5.FAMILIES
MIN_TRAIN_PER_BIN = 3  # theta-only baseline: bins are cut to keep this many training systems per bin (pilot sizes only)
CHANGED_MARKER = "x"
CASE_LABELS = {"a": "(a) no feedback", "b": "(b) feedback", "c": "(c) feedback + octupole", "d": "(d) feedback + octupole + spin-down"}


def nan_to_none(obj):
    """Recursively replace NaN floats by None so the JSON is strict."""
    if isinstance(obj, dict):
        return {k: nan_to_none(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [nan_to_none(v) for v in obj]
    return None if isinstance(obj, float) and obj != obj else obj


def case_dir(case: str, pilot: bool) -> Path:
    return W.STEP5_DIR if case == "a" else OUT / (f"{case}_pilot" if pilot else case)


def load_case(case: str, pilot: bool) -> dict:
    d = case_dir(case, pilot)
    feats = np.load(d / "features.npz")
    series = np.load(d / "series.npz")
    labels, _ = C5.v2_labels(feats, series)
    n_total = 1028 if case == "a" else int(np.load(d / "systems_all.npz")["outcome"].size)
    outcome = None if case == "a" else np.load(d / "systems_all.npz")["outcome"]
    return {"feats": feats, "labels": np.array(labels), "families": B.families(np.array(labels)), "n_total": n_total, "outcome": outcome}


def outcome_counts(outcome: np.ndarray | None, n_hj: int) -> dict:
    if outcome is None:
        return {"hot_jupiter": n_hj}
    return {W.OUTCOME_NAMES[k]: int((outcome == k).sum()) for k in range(len(W.OUTCOME_NAMES)) if (outcome == k).any()}


def train_mask(index_all: np.ndarray) -> dict:
    """index -> in the training half, with the step 5 baselines' draw (default_rng(seed).random(1028) < 0.5)."""
    rng = np.random.default_rng(CFG6["baseline_split_seed"])
    flags = rng.random(index_all.size) < B.HALF
    return {int(i): bool(f) for i, f in zip(index_all, flags, strict=True)}


def baselines(feats, obs: np.ndarray, pred: np.ndarray, flags: dict, rng: np.random.Generator) -> dict:
    train = np.array([flags[int(i)] for i in feats["index"]])
    test = ~train
    if train.sum() == 0 or test.sum() == 0:
        return {}
    B.N_THETA_BINS = max(1, min(10, int(train.sum()) // MIN_TRAIN_PER_BIN))
    omega = feats["n_max0"] / (2.0 * feats["a_bar0"])
    th_pred = B.theta_only(feats["theta_lb0_deg"], obs, train)
    om_pred, cut = B.map_plus_omega(pred, omega, obs, train)
    return {"n_train": int(train.sum()), "n_test": int(test.sum()), "s16_map_test_half": B.boot_ci((pred == obs)[test], rng),
            "theta_only_test_half": B.boot_ci((th_pred == obs)[test], rng),
            "map_plus_omega_test_half": B.boot_ci((om_pred == obs)[test], rng), "omega_cut": cut,
            "theta_only_bins": B.N_THETA_BINS}


def confusion_text(m: np.ndarray) -> str:
    lines = ["predicted \\ observed (v2)  " + "  ".join(f"{f:>5}" for f in FAMILIES) + "  total"]
    for i, p in enumerate(C5.PREDICTED):
        lines.append(f"{p:>26}  " + "  ".join(f"{v:5d}" for v in m[i]) + f"  {m[i].sum():5d}")
    return "\n".join(lines)


def transition_matrix(fam_a: np.ndarray, fam_case: np.ndarray) -> np.ndarray:
    m = np.zeros((len(FAMILIES), len(FAMILIES)), dtype=int)
    for a, c in zip(fam_a, fam_case, strict=True):
        m[FAMILIES.index(a), FAMILIES.index(c)] += 1
    return m


def transition_text(m: np.ndarray) -> str:
    lines = ["(a) \\ case   " + "  ".join(f"{f:>5}" for f in FAMILIES) + "  total"]
    for i, f in enumerate(FAMILIES):
        lines.append(f"{f:>11}  " + "  ".join(f"{v:5d}" for v in m[i]) + f"  {m[i].sum():5d}")
    return "\n".join(lines)


def paired(case_data: dict, a_data: dict) -> tuple[np.ndarray, np.ndarray]:
    """(fam_a, fam_case) of the systems that are hot Jupiters in both, matched by sample index."""
    pos_a = {int(i): k for k, i in enumerate(a_data["feats"]["index"])}
    rows_c = [k for k, i in enumerate(case_data["feats"]["index"]) if int(i) in pos_a]
    rows_a = [pos_a[int(case_data["feats"]["index"][k])] for k in rows_c]
    return a_data["families"][rows_a], case_data["families"][rows_c]


def plot_misclassified(feats, obs: np.ndarray, pred: np.ndarray, title: str, path) -> None:
    fig, ax = plt.subplots(figsize=C5.FIG_SIZE)
    wrong = obs != pred
    for f in FAMILIES:
        sel = (obs == f) & ~wrong
        ax.scatter(feats["a_bar0"][sel], feats["n_max0"][sel], s=C5.DOT, c=C5.COLOURS[f], label=f"{f} observed, predicted right ({sel.sum()})")
    for f in FAMILIES:
        sel = (obs == f) & wrong
        ax.scatter(feats["a_bar0"][sel], feats["n_max0"][sel], s=4 * C5.DOT, c=C5.COLOURS[f], marker=CHANGED_MARKER, linewidths=1.5,
                   label=f"{f} observed, MISclassified ({sel.sum()})")
    ax.axvline(1.0, color="k", ls="--", lw=1)
    for n in C5.N_MAX_LINES:
        ax.plot([1.0, C5.MAP_LIMITS[0][1]], [n, n], "k--", lw=1)
    ax.set(xscale="log", yscale="log", xlim=C5.MAP_LIMITS[0], ylim=C5.MAP_LIMITS[1], xlabel=r"$\bar{\mathcal{A}}_0$",
           ylabel=r"$N_{\rm max,0}$", title=title)
    ax.legend(loc="lower right", fontsize=6)
    fig.tight_layout()
    fig.savefig(path, dpi=C5.DPI)
    plt.close(fig)


def analyse_case(case: str, data: dict, a_data: dict, flags: dict, rng: np.random.Generator, pilot: bool) -> dict:
    feats, obs = data["feats"], data["families"]
    pred = np.array(feats["predicted"]).astype(str)
    suffix = "_pilot" if pilot else ""
    m = C5.confusion(list(pred), list(data["labels"]))
    hits = pred == obs
    res = {"n_total": data["n_total"], "outcomes": outcome_counts(data["outcome"], obs.size), "n_classified": int(obs.size),
           "observed_family_counts": {f: int((obs == f).sum()) for f in FAMILIES}, "confusion": m.tolist(),
           "accuracy_map": B.boot_ci(hits, rng), "baselines": baselines(feats, obs, pred, flags, rng)}
    text = confusion_text(m) + f"\nfamily agreement: {hits.sum()} / {hits.size} = {hits.mean():.3f}"
    if case != "a":
        fam_a, fam_c = paired(data, a_data)
        tm = transition_matrix(fam_a, fam_c)
        changed = fam_a != fam_c
        res["transition_a_to_case"] = tm.tolist()
        res["n_paired"] = int(fam_a.size)
        res["class_change_fraction"] = B.boot_ci(changed, rng) if fam_a.size else []
        text += "\n\nclass transition (a) -> case, paired systems:\n" + transition_text(tm)
        text += f"\nchanged: {changed.sum()} / {changed.size}"
    (OUT / f"confusion_{case}{suffix}.txt").write_text(text + "\n")
    print(f"== case {CASE_LABELS[case]}: {res['outcomes']}\n{text}")
    plot_misclassified(feats, obs, pred, f"{CASE_LABELS[case]}: misclassified (x) by S16's map, {obs.size} hot Jupiters",
                       OUT / f"{case}{suffix}_map.png")
    return res


def main() -> None:
    ap = argparse.ArgumentParser()
    ap.add_argument("--pilot", action="store_true")
    args = ap.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    a_data = load_case("a", False)
    flags = train_mask(a_data["feats"]["index"])
    rng = np.random.default_rng(CFG6["boot_seed"])
    summary = {"seed_bootstrap": CFG6["boot_seed"], "n_boot": B.N_BOOT, "pilot": args.pilot}
    for case in ("a", "b", "c", "d"):
        present = case == "a" or (case_dir(case, args.pilot) / "features.npz").exists()
        if not present:
            continue
        data = a_data if case == "a" else load_case(case, args.pilot)
        summary[case] = analyse_case(case, data, a_data, flags, rng, args.pilot)
    (OUT / ("summary_pilot.json" if args.pilot else "summary.json")).write_text(json.dumps(nan_to_none(summary), indent=2) + "\n")


if __name__ == "__main__":
    main()
