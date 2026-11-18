"""Step 9 T0 validation of lkspin/numshape.py; writes results/step9/t0_validation.json (config [step9]).

Gate 1  hat coefficients vs the analytic shape.py (SRFs off, 85 deg), and nmax_coef * eps_1 vs script-N.
Gate 2  SRFs on, tides off, 87 and 89 deg: vs prep_common.first_cycle_shape, grid convergence, symmetry diagnostics.
Gate 3  the Hamiltonian reproduces the secular spin integration over 20 cycles.
Run:  uv run python scripts/t0_validate.py
"""

import json
import math
import time

import numpy as np
import prep_common as C
from scipy.integrate import solve_ivp

from lkspin import numshape as N
from lkspin import regime as R
from lkspin import secular as S
from lkspin import shape as SH
from lkspin.lk_orbit import make_orbit

CFG4 = S.load_config()
CFG = C.load_table("step9")
OUT = S.ROOT / CFG["results"]
N_COEF = CFG["m_max"] + 1
TWO_PI = 2.0 * math.pi


def grids() -> dict:
    return {"n_t": CFG["n_t"], "n_tau": CFG["n_tau"], "m_max": CFG["m_max"], "rtol": CFG["rtol"], "atol": CFG["atol"]}


def shapes_for(par: np.ndarray, theta: float, scale: int = 1) -> N.NumShapes:
    g = grids()
    g["n_t"], g["n_tau"] = g["n_t"] * scale, g["n_tau"] * scale
    return N.first_cycle_shapes(par, CFG4, theta, **g)


def coefficient_stack(sh: N.NumShapes) -> np.ndarray:
    """(beta cos, gamma sin, psi cos) coefficients, shape (3, m_max + 1)."""
    return np.stack([sh.beta_cos, sh.gamma_sin, sh.psi_cos])


def summary(sh: N.NumShapes) -> dict:
    return {"nmax_coef": sh.nmax_coef, "abar_coef": sh.abar_coef, "omega_hat_l": sh.omega_hat_l, "psi_hat_0": sh.psi_hat_0,
            "psi_hat_0_over_omega_hat_l": sh.psi_hat_0 / sh.omega_hat_l, "e_min": sh.e_min, "e_max": sh.e_max,
            "period_yr": sh.period_yr, "beta_cos_0_10": sh.beta_cos.tolist(), "gamma_sin_0_10": sh.gamma_sin.tolist(),
            "psi_cos_0_10": sh.psi_cos.tolist()}


# ------------------------------------------------------------------ gate 1
def analytic_coefficients(orb) -> tuple[float, np.ndarray]:
    sn = SH.script_n(orb)
    ref = SH.shape_functions(orb)
    m_max = CFG["m_max"]
    stack = np.stack([SH.cos_coefficients(ref.beta, m_max), SH.sin_coefficients(ref.gamma, m_max),
                      SH.cos_coefficients(ref.psi, m_max)])
    return sn, sn * stack


def epsilon_one(par: np.ndarray, orb) -> float:
    """eps_1 = Omega_pl0 / alpha_0 at P_star = 1 day, Omega_pl0 = (3/4) sqrt(h) / t_k, alpha_0 = coef omega_star / a0^3."""
    t_k = S.t_k_years(par, par[S.P_A0])
    opl0 = 0.75 * math.sqrt(orb.h) / t_k
    alpha0 = par[S.P_STAR_COEF] * S.omega_star_of_period(N.P_STAR_REF_DAYS) / par[S.P_A0] ** 3
    return opl0 / alpha0


def gate1() -> dict:
    theta = CFG["gate1_theta_lb0_deg"]
    par = S.make_params(CFG4, srf=False, tides=False)
    sh = shapes_for(par, theta)
    orb = make_orbit(theta, CFG4["e0"])
    sn, expect = analytic_coefficients(orb)
    diff = float(np.abs(coefficient_stack(sh) - expect).max() / expect[0, 0])
    rel_n = abs(sh.nmax_coef * epsilon_one(par, orb) - sn) / sn
    return {"theta_lb0_deg": theta, "script_n": sn, "max_coef_diff_over_sn_beta0": diff,
            "pass_coef": bool(diff < CFG["gate1_coef_tol"]), "nmax_coef_times_eps_1": sh.nmax_coef * epsilon_one(par, orb),
            "eps_1": epsilon_one(par, orb), "nmax_rel_diff": rel_n, "pass_nmax": bool(rel_n < CFG["gate1_nmax_tol"]),
            "shapes": summary(sh)}


# ------------------------------------------------------------------ gate 2
def symmetry(sh: N.NumShapes) -> dict:
    b0 = sh.omega_hat_l
    return {"beta_sin_max_rel": float(np.abs(sh.beta_sin).max() / b0), "psi_sin_max_rel": float(np.abs(sh.psi_sin).max() / b0),
            "gamma_cos_max_rel": float(np.abs(sh.gamma_cos).max() / b0)}


def gate2_case(theta: float) -> dict:
    par = C.nondissipative_params(CFG4)
    sh = shapes_for(par, theta)
    fine = shapes_for(par, theta, 2)
    ref = C.first_cycle_shape(par, CFG4, theta)
    conv = float(np.abs(coefficient_stack(sh) - coefficient_stack(fine)).max() / sh.omega_hat_l)
    return {"theta_lb0_deg": theta, "ref_nmax_coef": ref["nmax_coef"], "ref_abar_coef": ref["abar_coef"],
            "nmax_rel_diff_to_ref": abs(sh.nmax_coef - ref["nmax_coef"]) / ref["nmax_coef"],
            "abar_rel_diff_to_ref": abs(sh.abar_coef - ref["abar_coef"]) / ref["abar_coef"],
            "convergence_max_rel": conv, "pass_convergence": bool(conv < CFG["gate2_conv_tol"]), "symmetry": symmetry(sh),
            "shapes": summary(sh)}


# ------------------------------------------------------------------ gate 3
def wrap(delta: np.ndarray) -> np.ndarray:
    return (delta + math.pi) % TWO_PI - math.pi


def secular_spins(par: np.ndarray, theta: float, p_star: float, times: np.ndarray, rtol: float,
                  atol: float) -> tuple[np.ndarray, np.ndarray]:
    """(p, phi) at the given times for every initial condition, as prep_sq4.section_chunk."""
    ics = CFG["gate3_ics"]
    y0 = S.initial_state(par, CFG4, theta, len(ics))
    for k, (p0, phi0) in enumerate(ics):
        o = S.N_ORBIT + S.N_SPIN * k
        y0[o:o + S.N_SPIN] = C.initial_spin_vector(theta, p0, phi0)
    states, status = C.run_to_times(y0, par, [p_star] * len(ics), times, rtol, atol)
    lvec = np.broadcast_to(states[:, None, 0:3], (times.size, len(ics), S.N_SPIN))
    return C.phase_and_p(lvec, C.spin_block(states, len(ics)))


def hamiltonian_spins(sh: N.NumShapes, n_max: float, ic: list, n_cycles: int) -> np.ndarray:
    """(phi, p) at tau = 2 pi k, k = 1..n_cycles: array (n_cycles, 2)."""
    interp = N.shape_interp(sh)
    p0, phi0 = ic
    t_eval = TWO_PI * np.arange(1, n_cycles + 1)
    sol = solve_ivp(N.hamiltonian_rhs, (0.0, t_eval[-1]), [phi0, p0], method="DOP853", t_eval=t_eval,
                    rtol=CFG["gate3_rtol"], atol=CFG["gate3_atol"], args=(n_max, interp))
    return sol.y.T


def gate3_case(par: np.ndarray, sh: N.NumShapes, times: np.ndarray, p_star: float, rtol: float, atol: float) -> list[dict]:
    theta = CFG["gate3_theta_lb0_deg"]
    p_sec, phi_sec = secular_spins(par, theta, p_star, times, rtol, atol)
    out = []
    for k, ic in enumerate(CFG["gate3_ics"]):
        ham = hamiltonian_spins(sh, sh.nmax_coef / p_star, ic, times.size)
        dp = float(np.abs(ham[:, 1] - p_sec[:, k]).max())
        dphi = float(np.abs(wrap(ham[:, 0] - phi_sec[:, k])).max())
        ok = bool(dp < CFG["gate3_dp_tol"] and dphi < CFG["gate3_dphi_tol"])
        out.append({"P_star_days": p_star, "p0": ic[0], "phi0": ic[1], "max_abs_dp": dp, "max_abs_dphi": dphi, "pass": ok})
    return out


def tight_minima_times(par: np.ndarray, theta: float, n_cycles: int) -> np.ndarray:
    """First n_cycles eccentricity minima from run_diag at the tight tolerance."""
    t_stop = CFG["gate3_tight_t_stop_periods"] * (n_cycles + 1) * C.analytic_period_yr(par, CFG4, theta)
    y0 = S.initial_state(par, CFG4, theta, 1)
    res = S.run_diag(y0, par, [N.SPIN_PROBE_DAYS], t_stop, CFG["gate3_tight_rtol"], CFG["gate3_tight_atol"],
                     ev_max=CFG["gate3_tight_ev_max"])
    minima, _ = R.split_events(res["ev"])
    return minima[:n_cycles, 0].copy()


def gate3(tight: bool) -> list[dict]:
    """tight False: minima times from prep_common.orbit_minima_times and spins at [step4] tolerance (as prep_sq4)."""
    theta = CFG["gate3_theta_lb0_deg"]
    par = C.nondissipative_params(CFG4)
    sh = shapes_for(par, theta)
    n_cycles = CFG["gate3_n_cycles"]
    times = tight_minima_times(par, theta, n_cycles) if tight else C.orbit_minima_times(par, CFG4, theta, n_cycles)
    rtol, atol = (CFG["gate3_tight_rtol"], CFG["gate3_tight_atol"]) if tight else (CFG["spin_rtol"], CFG["spin_atol"])
    rows = []
    for p_star in CFG["gate3_p_star_days"]:
        rows += gate3_case(par, sh, times, p_star, rtol, atol)
    return rows


# ------------------------------------------------------------------ main
def main() -> None:
    t0 = time.time()
    OUT.mkdir(parents=True, exist_ok=True)
    g1 = gate1()
    g2 = [gate2_case(th) for th in CFG["gate2_theta_lb0_deg"]]
    g3, g3_tight = gate3(tight=False), gate3(tight=True)
    result = {"config": CFG, "gate1": g1, "gate2": g2, "gate3": g3, "gate3_tight_tolerance": g3_tight,
              "wall_s": time.time() - t0}
    (OUT / "t0_validation.json").write_text(json.dumps(result, indent=2))
    brief = {"gate1": {k: v for k, v in g1.items() if k != "shapes"},
             "gate2": [{k: v for k, v in c.items() if k != "shapes"} for c in g2], "gate3": g3,
             "gate3_tight_tolerance": g3_tight, "wall_s": result["wall_s"]}
    print(json.dumps(brief, indent=2))


if __name__ == "__main__":
    main()
