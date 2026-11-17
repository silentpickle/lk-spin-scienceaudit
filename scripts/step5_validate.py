"""Step 5 validation of the winding / event diagnostics and the first-pass labels on the S16 Fig. 8 spot cases.

Writes results/step5/validation.json.  Run:  uv run python scripts/step5_validate.py
"""

import json
import math
import tomllib

import numpy as np

from lkspin import regime as R
from lkspin import secular as S
from lkspin.lk_orbit import make_orbit

CFG4 = S.load_config()
with S.CONFIG_PATH.open("rb") as _fh:
    CFG5 = tomllib.load(_fh)["step5"]
OUT = S.ROOT / CFG5["results"] / "validation.json"
THETA_STEP3 = 85.0
N_CYC = 100
N_SKIP = 30
EPS_RESONANCE = 0.0869  # inside the N = 11 lock found by scanning eps (theta_sl up to ~25 deg; S15 Fig. 12 band 0.07-0.095)
EPS_ADIABATIC = 0.008   # lowest eps of S15 Fig. 12
SPOTS = (29.87, 7.07, 1.67)
SPOT_EXPECT = {29.87: "I", 7.07: "II", 1.67: "IV.2"}
SPOT_THETA, SPOT_MP, SPOT_AB = 89.0, 5.0, 300.0


def step3_series(eps: float) -> dict:
    par = S.make_params(CFG4, srf=False, tides=False)
    a0 = CFG4["a0_au"]
    tk = S.t_k_years(par, a0)
    orb = make_orbit(THETA_STEP3, CFG4["e0"])
    p = 2 * math.pi * S.DAYS_PER_YR / (orb.alpha0(eps) / tk * a0**3 / par[S.P_STAR_COEF])
    y0 = S.initial_state(par, CFG4, THETA_STEP3, 1)
    res = S.run_diag(y0, par, [p], N_CYC * orb.period * tk)
    minima, _ = R.split_events(res["ev"])
    w = R.winding_series(minima, 0, 1)[N_SKIP:]
    nmax = R.shape_coefficients(par, CFG4, THETA_STEP3, 4 * orb.period * tk)["nmax_coef"] / p
    return {"eps": eps, "P_star_days": p, "W_mean": float(w.mean()), "W_min": float(w.min()), "W_max": float(w.max()),
            "locked_frac_tol0.1": float(R.locked_mask(w, CFG5["lock_tol"]).mean()), "N_max": nmax,
            "theta_range_deg": [float(np.degrees(np.arccos(minima[N_SKIP:, 4])).min()),
                                float(np.degrees(np.arccos(minima[N_SKIP:, 4])).max())]}


def spot(p: float) -> dict:
    f = R.integrate_system(SPOT_MP, p, SPOT_AB, SPOT_THETA, CFG4, CFG5)
    return {"P_star": p, "label": R.classify(f, CFG5), "expected": SPOT_EXPECT[p], "a_bar0": f["a_bar0"],
            "n_max0": f["n_max0"], "theta_f_range": [float(f["theta_f"].min()), float(f["theta_f"].max())],
            "gap": f["gap"], "twin_div_deg": f["twin_div_deg"], "locked_frac_last": f["m0_locked_frac_last"],
            "n_dom": f["m0_n_dom"], "w_median_last": f["m0_w_median_last"]}


def main() -> None:
    res = {"resonance_N11": step3_series(EPS_RESONANCE), "adiabatic_limit": step3_series(EPS_ADIABATIC),
           "spots_S16_fig8": [spot(p) for p in SPOTS]}
    OUT.write_text(json.dumps(res, indent=1))
    print(json.dumps(res, indent=1))


if __name__ == "__main__":
    main()
