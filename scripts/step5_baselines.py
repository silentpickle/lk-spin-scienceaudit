"""Step 5 / SQ1(a) baselines: accuracy of S16's (A_bar_0, N_max,0) map against the v2 observed classes, with bootstrap
errors, compared with (i) a classifier using theta_lb0 alone (the brief's negative-result test) and (ii) the map plus
Omega_hat_L = N_max,0 / (2 A_bar_0) splitting the IV region into III / IV. The baselines are fitted on a random half
and scored on the other half. Writes results/step5/baselines.json.

Run:  uv run python scripts/step5_baselines.py
"""

import json
import tomllib
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
CFG5 = tomllib.loads((ROOT / "config.toml").read_text())["step5"]
OUT = ROOT / CFG5["results"]
SEED = 5
N_BOOT = 2000
N_THETA_BINS = 10
OMEGA_GRID = np.linspace(0.02, 0.4, 77)
HALF = 0.5
CI = (2.5, 97.5)


def families(labels: np.ndarray) -> np.ndarray:
    return np.array([lab.split(".")[0] for lab in labels])


def boot_ci(hits: np.ndarray, rng: np.random.Generator) -> list[float]:
    means = [hits[rng.integers(0, hits.size, hits.size)].mean() for _ in range(N_BOOT)]
    return [float(np.mean(hits)), *[float(v) for v in np.percentile(means, CI)]]


def theta_only(theta: np.ndarray, obs: np.ndarray, train: np.ndarray) -> np.ndarray:
    """Majority observed class per theta_lb0 bin (bins from training quantiles); returns predictions for all."""
    edges = np.quantile(theta[train], np.linspace(0, 1, N_THETA_BINS + 1))
    bins = np.clip(np.searchsorted(edges, theta, side="right") - 1, 0, N_THETA_BINS - 1)
    major = {}
    for b in range(N_THETA_BINS):
        vals, counts = np.unique(obs[train & (bins == b)], return_counts=True)
        major[b] = vals[np.argmax(counts)]
    return np.array([major[b] for b in bins])


def map_plus_omega(pred: np.ndarray, omega: np.ndarray, obs: np.ndarray, train: np.ndarray) -> tuple[np.ndarray, float]:
    """S16 map, but predicted-IV systems with Omega_hat_L above a threshold (fitted on train) are called III."""
    def apply(cut: float) -> np.ndarray:
        return np.where((pred == "IV") & (omega >= cut), "III", pred)
    acc = [(apply(c)[train] == obs[train]).mean() for c in OMEGA_GRID]
    cut = float(OMEGA_GRID[int(np.argmax(acc))])
    return apply(cut), cut


def main() -> None:
    f = np.load(OUT / "features.npz")
    obs = families(np.load(OUT / "labels_v2.npz")["label"])
    pred = f["predicted"].astype(str)
    omega = f["n_max0"] / (2.0 * f["a_bar0"])
    rng = np.random.default_rng(SEED)
    train = rng.random(obs.size) < HALF
    test = ~train
    th_pred = theta_only(f["theta_lb0_deg"], obs, train)
    om_pred, cut = map_plus_omega(pred, omega, obs, train)
    res = {
        "n_systems": int(obs.size), "seed": SEED, "n_boot": N_BOOT,
        "s16_map_all": boot_ci(pred == obs, rng),
        "s16_map_test_half": boot_ci((pred == obs)[test], rng),
        "theta_only_test_half": boot_ci((th_pred == obs)[test], rng),
        "map_plus_omega_test_half": boot_ci((om_pred == obs)[test], rng),
        "omega_cut": cut,
    }
    (OUT / "baselines.json").write_text(json.dumps(res, indent=2) + "\n")
    for k, v in res.items():
        print(k, v)


if __name__ == "__main__":
    main()
