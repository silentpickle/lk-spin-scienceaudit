"""Step 5, revised classifier (v2): apply it to every hot-Jupiter system, score it against the hand labels, redraw
the regime map (cf. S16 Fig. 14) and write the predicted-vs-observed confusion matrix.

Hand labels (results/step5/handlabels.csv) for cases 1-30 were used to choose the v2 rules; cases 31-50 are held out.
Run:  uv run python scripts/step5_classify.py
"""

import csv
import json
import tomllib
from collections import Counter
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from lkspin import regime as R

ROOT = Path(__file__).resolve().parent.parent
CFG5 = tomllib.loads((ROOT / "config.toml").read_text())["step5"]
OUT = ROOT / CFG5["results"]
FAMILIES = ("I", "II", "III", "IV", "V")
PREDICTED = ("I", "II", "IV")
COLOURS = {"I": "#3b3fa0", "II": "#b03a7a", "III": "#2e9a46", "IV": "#b39a2a", "V": "#e8731a"}
MAP_LIMITS = ((0.05, 100.0), (0.02, 40.0))
N_MAX_LINES = (0.5, 1.0, 2.0)
FIG_SIZE = (8.5, 5.5)
DPI = 150
DOT = 9
AMBIGUOUS = "ambiguous"


def load():
    feats = np.load(OUT / "features.npz")
    series = np.load(OUT / "series.npz")
    return feats, series


def v2_labels(feats, series) -> tuple[list[str], list[dict]]:
    labels, rows = [], []
    for i, idx in enumerate(feats["index"]):
        theta_all = np.append(feats["theta_f"][i], feats["theta_f_twin"][i])
        g = R.v2_features(theta_all, series[f"s{int(idx)}"], CFG5)
        labels.append(R.classify_v2(g, CFG5))
        rows.append(g)
    return labels, rows


def score(hand_path: Path, idx_to_label: dict) -> dict:
    """Agreement with hand labels by split, at family level; ambiguous cases reported separately."""
    out = {}
    with hand_path.open() as fh:
        for r in csv.DictReader(fh):
            got = R.label_family(idx_to_label[int(r["index"])])
            rec = out.setdefault(r["split"], {"n": 0, "agree": 0, "ambiguous": 0, "misses": []})
            if r["hand_label"] == AMBIGUOUS:
                rec["ambiguous"] += 1
                rec["misses"].append((r["row"], r["index"], "ambiguous", got))
                continue
            rec["n"] += 1
            if got == r["hand_label"]:
                rec["agree"] += 1
            else:
                rec["misses"].append((r["row"], r["index"], r["hand_label"], got))
    return out


def confusion(pred: list[str], obs: list[str]) -> np.ndarray:
    m = np.zeros((len(PREDICTED), len(FAMILIES)), dtype=int)
    for p, o in zip(pred, obs, strict=True):
        m[PREDICTED.index(p), FAMILIES.index(R.label_family(o))] += 1
    return m


def write_confusion(m: np.ndarray, path: Path) -> None:
    lines = ["predicted \\ observed (v2)  " + "  ".join(f"{f:>5}" for f in FAMILIES) + "  total"]
    for i, p in enumerate(PREDICTED):
        lines.append(f"{p:>26}  " + "  ".join(f"{v:5d}" for v in m[i]) + f"  {m[i].sum():5d}")
    agree = sum(m[PREDICTED.index(f), FAMILIES.index(f)] for f in PREDICTED)
    lines.append(f"family agreement: {agree} / {m.sum()} = {agree / m.sum():.3f}")
    path.write_text("\n".join(lines) + "\n")
    print("\n".join(lines))


def plot_map(feats, labels: list[str], path: Path) -> None:
    fam = np.array([R.label_family(lab) for lab in labels])
    fig, ax = plt.subplots(figsize=FIG_SIZE)
    for f in FAMILIES:
        sel = fam == f
        ax.scatter(feats["a_bar0"][sel], feats["n_max0"][sel], s=DOT, c=COLOURS[f], label=f"{f} ({sel.sum()})")
    ax.axvline(1.0, color="k", ls="--", lw=1)
    for n in N_MAX_LINES:
        ax.plot([1.0, MAP_LIMITS[0][1]], [n, n], "k--", lw=1)
    ax.set(xscale="log", yscale="log", xlim=MAP_LIMITS[0], ylim=MAP_LIMITS[1], xlabel=r"$\bar{\mathcal{A}}_0$",
           ylabel=r"$N_{\rm max,0}$", title=f"Outcome class (v2 classifier), {fam.size} hot Jupiters (cf. S16 Fig. 14)")
    ax.legend(loc="lower right", fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=DPI)
    plt.close(fig)


def main() -> None:
    feats, series = load()
    labels, rows = v2_labels(feats, series)
    idx_to_label = {int(i): lab for i, lab in zip(feats["index"], labels, strict=True)}
    sc = score(OUT / "handlabels.csv", idx_to_label)
    for split, rec in sc.items():
        print(f"{split}: {rec['agree']}/{rec['n']} agree (+{rec['ambiguous']} ambiguous); misses {rec['misses']}")
    m = confusion(list(feats["predicted"]), labels)
    write_confusion(m, OUT / "confusion_v2.txt")
    print("v2 classes:", dict(Counter(labels)))
    plot_map(feats, labels, OUT / "regime_map_v2.png")
    np.savez_compressed(OUT / "labels_v2.npz", index=feats["index"], label=np.array(labels),
                        clusters=[r["clusters"] for r in rows], swing=[r["swing"] for r in rows],
                        lock_frac=[r["lock_frac"] for r in rows], lock_n=[r["lock_frac_n"] for r in rows])
    (OUT / "handcheck_score_v2.json").write_text(json.dumps(sc, indent=2) + "\n")


if __name__ == "__main__":
    main()
