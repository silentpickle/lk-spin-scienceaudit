"""Step 11 analysis: secular vs N-body e_max table (results/step11/comparison.csv) and the 4-panel figure.

Method: the N-body e(t) is block-averaged over one outer binary period (see step11_nbody.py). Its e_max is the peak of
that running mean (peaks of -log10(1-e) with prominence >= peak_prominence_dex). The envelope e_max is the largest
osculating e in the peak block and its two neighbours. The secular e_max is the peak of the (smooth) secular e(t).
Run:  uv run python scripts/step11_analyze.py
"""

import csv
import math

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

import step11_common as C

PLOT_MAX_POINTS = 20000
N_PANELS_EXTREME = 2
ENVELOPE_HALF_WIDTH = 1


def load_nbody(name: str) -> dict | None:
    path = C.CHUNK_DIR / f"{name}.npz"
    if not path.exists():
        return None
    return {k: v for k, v in np.load(path).items()}


def rel(a: float, b: float) -> float:
    return a / b - 1.0 if b and not math.isnan(a) else math.nan


def nbody_summary(nb: dict | None) -> dict:
    """First/second e_max (running mean), envelope e_max, times, period, theta at the first peak."""
    keys = ("emax1", "emax2", "env1", "tmax1", "tmax2", "theta_at_emax1", "period", "reached_t", "complete", "wall_s",
            "steps", "comp_inc_deg")
    out = dict.fromkeys(keys, math.nan)
    if nb is None or len(nb["t"]) == 0:
        return out
    pk = C.peak_indices(nb["t"], nb["e_mean"])
    out.update(reached_t=nb["t"][-1], complete=bool(nb["complete"]), wall_s=float(nb["wall_s"]))
    out["steps"] = float(nb["steps"]) if "steps" in nb else math.nan
    out["comp_inc_deg"] = float(nb["comp_inc_deg_final"]) if "comp_inc_deg_final" in nb else math.nan
    for i, key in enumerate(("1", "2")):
        if i < len(pk):
            out[f"emax{key}"] = nb["e_mean"][pk[i]]
            out[f"tmax{key}"] = nb["t"][pk[i]]
    if len(pk):
        lo, hi = max(pk[0] - ENVELOPE_HALF_WIDTH, 0), pk[0] + ENVELOPE_HALF_WIDTH + 1
        out["env1"] = nb["e_max"][lo:hi].max()
        out["theta_at_emax1"] = nb["theta_mean"][pk[0]]
    out["period"] = out["tmax2"] - out["tmax1"]
    return out


def system_row(sys_: C.System) -> dict:
    gr, full = C.secular_summary(sys_, "gr"), C.secular_summary(sys_, "full")
    nb = nbody_summary(load_nbody(sys_.name))
    b = C.breakdown(sys_, full["emax1"], gr["emax1"])
    row = {"name": sys_.name, "kind": sys_.kind, "m_p_mjup": sys_.m_p_mjup, "a_b_au": sys_.a_b_au,
           "theta_lb0_deg": sys_.theta_lb0_deg, "mean_anom_planet_rad": sys_.mean_anom_planet,
           "mean_anom_comp_rad": sys_.mean_anom_comp, "t_k_yr": b["t_k"]}
    row |= {f"sec_gr_{k}": v for k, v in gr.items() if k not in ("tk", "n_peaks")}
    row |= {f"sec_full_{k}": v for k, v in full.items() if k not in ("tk", "n_peaks")}
    row |= {f"nb_{k}": v for k, v in nb.items()}
    row["rel_1me_full_vs_gr"] = rel(1 - full["emax1"], 1 - gr["emax1"])
    row["rel_1me_nb_vs_gr"] = rel(1 - nb["emax1"], 1 - gr["emax1"])
    row["rel_1me_nbenv_vs_gr"] = rel(1 - nb["env1"], 1 - gr["emax1"])
    row["abs_emax_nb_minus_gr"] = nb["emax1"] - gr["emax1"]
    row["rel_tmax1_nb_vs_gr"] = rel(nb["tmax1"], gr["tmax1"])
    row["rel_period_nb_vs_gr"] = rel(nb["period"], gr["period"])
    row["rel_tmax2_nb_vs_gr"] = rel(nb["tmax2"], gr["tmax2"])
    row["rel_1me2_nb_vs_gr"] = rel(1 - nb["emax2"], 1 - gr["emax2"])
    row["dtheta_at_emax1_deg_nb_minus_gr"] = nb["theta_at_emax1"] - gr["theta_at_emax1"]
    row |= b | C.flags(b)
    return row


def write_csv(rows: list[dict]) -> None:
    C.OUT_DIR.mkdir(parents=True, exist_ok=True)
    with (C.OUT_DIR / "comparison.csv").open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0]))
        w.writeheader()
        w.writerows(rows)


def decimate(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    step = max(1, len(x) // PLOT_MAX_POINTS)
    return x[::step], y[::step]


def plot_panel(ax, sys_: C.System) -> None:
    nb = load_nbody(sys_.name)
    for mode, style, label in (("gr", "-", "secular, GR only"), ("full", "--", "secular, full SRFs")):
        s = C.secular_series(sys_, mode)
        t, e = decimate(s["t"], s["e"])
        ax.plot(t / 1e6, 1 - e, style, lw=1.0, label=label)
    if nb is not None and len(nb["t"]):
        t = nb["t"] / 1e6
        ax.fill_between(t, 1 - nb["e_max"], 1 - nb["e_min"], color="gray", alpha=0.35, lw=0, label="N-body min..max in block")
        ax.plot(t, 1 - nb["e_mean"], "k", lw=0.8, label="N-body block mean")
    ax.set_yscale("log")
    ax.set_title(f"{sys_.name}: Mp={sys_.m_p_mjup:.2f} MJ, ab={sys_.a_b_au:.0f} AU, th0={sys_.theta_lb0_deg:.2f}", fontsize=9)
    ax.set_xlabel("t (Myr)")
    ax.set_ylabel("1 - e")


def make_figure(systems: list[C.System], rows: list[dict]) -> None:
    fid = [s for s in systems if s.kind == "fiducial"]
    ranked = sorted(zip(rows, systems, strict=True), key=lambda p: -p[0]["sec_gr_emax1"])
    extreme = [s for _, s in ranked[:N_PANELS_EXTREME] if s.kind != "fiducial"]
    chosen = (fid + extreme)[:4]
    fig, axes = plt.subplots(2, 2, figsize=(11, 7.5))
    for ax, s in zip(axes.flat, chosen, strict=False):
        plot_panel(ax, s)
    axes.flat[0].legend(fontsize=7)
    fig.tight_layout()
    fig.savefig(C.OUT_DIR / "nbody_vs_secular.png", dpi=130)


def main() -> None:
    systems = C.load_systems()
    rows = [system_row(s) for s in systems]
    write_csv(rows)
    make_figure(systems, rows)


if __name__ == "__main__":
    main()
