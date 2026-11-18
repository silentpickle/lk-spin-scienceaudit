"""Job 4 analysis: which feedback channel drives the SQ1 class changes?  Cases b (both channels), b_L (dL/dt|SL only) and
b_E (de/dt|SL only), each against case (a) (results/step5).  Reuses step6_analyze (v2 classes, confusion, accuracy, baselines,
transitions) and step10_sq5 (class-change fraction against S/L).  Writes results/job4/summary.json and class_change_vs_sl.png.
Run:  uv run python scripts/job4_analyze.py
"""

import json

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import step5_baselines as B
import step6_analyze as A
import step6_run as W
import step10_sq5 as Q
from lkspin import secular as S

OUT = S.ROOT / "results" / "job4"
CASES = ("b", "b_L", "b_E", "b_P", "b_F")
PAIRS = (("b_L", "b_E"), ("b_P", "b_F"))  # (L, E) as specified; (P, F) the Hamiltonian-consistent split, supplementary
FAMILIES = A.FAMILIES
ADVECTION = "IV"
FIG_SIZE = (8.0, 5.0)
DPI = 150
COLOURS = {"b": "k", "b_L": "tab:blue", "b_E": "tab:red", "b_P": "tab:orange", "b_F": "tab:green"}
LABELS = {"b": "(b) both channels", "b_L": "(b_L) L only: dL/dt|SL", "b_E": "(b_E) E only: de/dt|SL",
          "b_P": "(b_P) apsidal precession only", "b_F": "(b_F) all but apsidal precession"}


def by_index(data: dict) -> dict:
    """sample index -> v2 family, for the hot-Jupiter systems of a loaded case."""
    return {int(i): f for i, f in zip(data["feats"]["index"], data["families"], strict=True)}


def outcome_by_index(case: str) -> dict:
    d = np.load(A.OUT / case / "systems_all.npz")
    return {int(i): int(o) for i, o in zip(d["index"], d["outcome"], strict=True)}


def transition_from(src: dict, dst: dict) -> dict:
    """Transition matrix src -> dst over systems classified in both, plus the system count."""
    keys = sorted(set(src) & set(dst))
    m = A.transition_matrix(np.array([src[k] for k in keys]), np.array([dst[k] for k in keys]))
    return {"matrix": m.tolist(), "n_paired": len(keys), "changed": int(sum(src[k] != dst[k] for k in keys))}


def reproduction(fam: dict, ch_a: str, ch_b: str) -> dict:
    """Systems classified in a, b, b_L and b_E and whose class changes between a and b: how often each channel alone gives the b class."""
    keys = sorted(set.intersection(*(set(fam[c]) for c in ("a", "b", ch_a, ch_b))))
    changed = [k for k in keys if fam["a"][k] != fam["b"][k]]
    out = {"n_common_all_four": len(keys), "n_changed_a_to_b": len(changed), "n_unchanged_a_to_b": len(keys) - len(changed)}
    for ch in (ch_a, ch_b):
        out[f"{ch}_gives_b_class"] = int(sum(fam[ch][k] == fam["b"][k] for k in changed))
        out[f"{ch}_leaves_a_class"] = int(sum(fam[ch][k] == fam["a"][k] for k in changed))
        out[f"{ch}_other_class"] = int(sum(fam[ch][k] not in (fam["a"][k], fam["b"][k]) for k in changed))
        out[f"{ch}_changes_unchanged_systems"] = int(sum(fam[ch][k] != fam["a"][k] for k in keys if k not in set(changed)))
    fa, fb_ = fam[ch_a], fam[ch_b]
    out["both_give_b_class"] = int(sum(fa[k] == fam["b"][k] == fb_[k] for k in changed))
    out[f"only_{ch_a}_gives_b_class"] = int(sum(fa[k] == fam["b"][k] != fb_[k] for k in changed))
    out[f"only_{ch_b}_gives_b_class"] = int(sum(fb_[k] == fam["b"][k] != fa[k] for k in changed))
    out["neither_gives_b_class"] = int(sum(fam["b"][k] not in (fa[k], fb_[k]) for k in changed))
    out[f"fraction_{ch_a}"] = out[f"{ch_a}_gives_b_class"] / max(1, len(changed))
    out[f"fraction_{ch_b}"] = out[f"{ch_b}_gives_b_class"] / max(1, len(changed))
    return out


def migration_loss(outs: dict, ch_a: str, ch_b: str) -> dict:
    """Systems that are hot Jupiters in (a) but not in (b): do they fail (timeout or mixed) under b_L, under b_E?"""
    hj = W.S.ST_HOT_JUPITER
    lost = [k for k, o in outs["b"].items() if o != hj]  # (a) forms a hot Jupiter for all 1028 systems
    fail = {c: {k: outs[c][k] != hj for k in lost} for c in (ch_a, ch_b)}
    res = {"n_lost_a_to_b": len(lost), "lost_by_b_outcome": {W.OUTCOME_NAMES[o]: int(sum(outs["b"][k] == o for k in lost)) for o in set(outs["b"].values()) if o != hj}}
    for c in (ch_a, ch_b):
        res[f"fail_under_{c}"] = int(sum(fail[c].values()))
        res[f"{c}_outcomes_of_lost"] = {W.OUTCOME_NAMES[o]: int(sum(outs[c][k] == o for k in lost)) for o in sorted(set(outs[c][k] for k in lost))}
    res["fail_under_both"] = int(sum(fail[ch_a][k] and fail[ch_b][k] for k in lost))
    res[f"fail_only_{ch_a}"] = int(sum(fail[ch_a][k] and not fail[ch_b][k] for k in lost))
    res[f"fail_only_{ch_b}"] = int(sum(fail[ch_b][k] and not fail[ch_a][k] for k in lost))
    res["fail_neither"] = int(sum(not fail[ch_a][k] and not fail[ch_b][k] for k in lost))
    # also the reverse: systems that are hot Jupiters under b but fail under a single channel
    for c in (ch_a, ch_b):
        res[f"hj_under_b_but_fail_{c}"] = int(sum(o == hj and outs[c][k] != hj for k, o in outs["b"].items()))
    res["n_hot_jupiter"] = {c: int(sum(o == hj for o in outs[c].values())) for c in CASES}
    return res


def sq5_case(case: str, edges: np.ndarray, rng: np.random.Generator) -> dict:
    data = Q.paired_arrays(case, False)
    centres_log = np.log10(np.sqrt(edges[:-1] * edges[1:]))
    stats = Q.bin_stats(data, edges)
    return {"data": data, "stats": stats, "boots": Q.boot_bins(data, edges, rng), "cross": Q.crossings(stats, centres_log),
            "cross_boot": Q.boot_crossings(data, edges, centres_log, rng)}


def sq5_json(r: dict) -> dict:
    return {"n_paired": int(r["data"]["s_over_l"].size), "overall_class_change_fraction": float(r["data"]["changed"].mean()),
            "n": r["stats"]["n"].tolist(), "class_change_fraction": r["stats"]["frac"].tolist(), "class_change_ci": r["boots"]["frac"].tolist(),
            "median_abs_dtheta_deg": r["stats"]["dtheta"].tolist(), "crossings_s_over_l": r["cross"], "crossings_bootstrap": r["cross_boot"]}


def plot_sq5(res: dict, edges: np.ndarray) -> None:
    centres = np.sqrt(edges[:-1] * edges[1:])
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    offsets = {"b": 1.0, "b_L": 1.03, "b_E": 0.97, "b_P": 1.06, "b_F": 0.94}  # small horizontal nudge so the error bars do not overlap
    for case in CASES:
        st, bt = res[case]["stats"], res[case]["boots"]
        ok = st["n"] > 0
        ax.errorbar(centres[ok] * offsets[case], st["frac"][ok], yerr=[st["frac"][ok] - bt["frac"][0][ok], bt["frac"][1][ok] - st["frac"][ok]],
                    fmt="o-" if case in ("b", "b_L", "b_E") else "s:", color=COLOURS[case], ms=4, capsize=2, label=f"a -> {LABELS[case]} (n={res[case]['data']['s_over_l'].size})")
    for level in W.CFG6["sq5_fraction_levels"]:
        ax.axhline(level, color="gray", ls="--", lw=1)
    ax.set(xscale="log", ylim=(0, 1), xlabel=r"$S_\star/L$ at the first $e$ maximum (step 5 $P_\star$)", ylabel="class-change fraction")
    ax.legend(fontsize=7, loc="upper left")
    fig.tight_layout()
    fig.savefig(OUT / "class_change_vs_sl.png", dpi=DPI)
    plt.close(fig)


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    a_data = A.load_case("a", False)
    flags = A.train_mask(a_data["feats"]["index"])
    rng = np.random.default_rng(W.CFG6["boot_seed"])
    data = {"a": a_data} | {c: A.load_case(c, False) for c in CASES}
    summary = {"n_boot": B.N_BOOT, "boot_seed": W.CFG6["boot_seed"], "cases": {}}
    for c in CASES:
        summary["cases"][c] = A.analyse_case(c, data[c], a_data, flags, rng, False)
    fam = {c: by_index(data[c]) for c in ("a", *CASES)}
    summary["transition_b_to"] = {c: transition_from(fam["b"], fam[c]) for c in CASES[1:]}
    summary["transition_a_to"] = {c: transition_from(fam["a"], fam[c]) for c in CASES}
    outs = {c: outcome_by_index(c) for c in CASES}
    summary["reproduction"] = {f"{x}_{y}": reproduction(fam, x, y) for x, y in PAIRS}
    summary["migration_loss"] = {f"{x}_{y}": migration_loss(outs, x, y) for x, y in PAIRS}
    summary["advection_class_IV"] = {c: int(sum(f == ADVECTION for f in fam[c].values())) for c in ("a", *CASES)}
    summary["a_IV_systems_destinations"] = {c: {f: int(sum(fam["a"][k] == ADVECTION and fam[c][k] == f for k in fam["a"] if k in fam[c])) for f in FAMILIES}
                                            for c in CASES}
    summary["a_IV_systems_still_hj"] = {c: int(sum(k in fam[c] for k in fam["a"] if fam["a"][k] == ADVECTION)) for c in CASES}
    edges = Q.log_edges(Q.paired_arrays("b", False)["s_over_l"])
    sq5 = {c: sq5_case(c, edges, rng) for c in CASES}
    summary["sq5"] = {c: sq5_json(sq5[c]) for c in CASES} | {"edges": edges.tolist()}
    plot_sq5(sq5, edges)
    summary["l_dot_e_diagnostic"] = {c: ld_diag(c) for c in CASES[1:]}
    (OUT / "summary.json").write_text(json.dumps(A.nan_to_none(summary), indent=2) + "\n")
    print_report(summary)


def ld_diag(case: str) -> dict:
    """|cos(l, e)| at the end of each run, over all chunks (a consistent orbit has 0)."""
    vals = np.concatenate([np.abs(np.load(p)["cos_le_end"]) for p in sorted((A.OUT / case / "chunks").glob("sys_*.npz"))])
    return {"n": int(vals.size), "median": float(np.median(vals)), "p95": float(np.percentile(vals, 95)), "max": float(vals.max()),
            "fraction_above_0.1": float((vals > 0.1).mean())}


def print_report(s: dict) -> None:
    for key in ("advection_class_IV", "migration_loss", "reproduction", "a_IV_systems_still_hj", "l_dot_e_diagnostic"):
        print(key, json.dumps(s[key]))
    for c in CASES:
        r = s["cases"][c]
        print(c, r["outcomes"], r["observed_family_counts"], "acc", r["accuracy_map"], "theta-only", r["baselines"].get("theta_only_test_half"),
              "cross", s["sq5"][c]["crossings_s_over_l"], s["sq5"][c]["crossings_bootstrap"])


if __name__ == "__main__":
    main()
