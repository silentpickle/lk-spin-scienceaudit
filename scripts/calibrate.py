"""Work plan step 2: recompute each published prefactor from its definition and compare.

Run:  uv run python scripts/calibrate.py
Writes the table to stdout and to the results path in config.toml. Each check is one small function returning the
computed value; CHECKS pairs it with the published number and the equation it comes from.
"""

import hashlib
import json
import math
import tomllib
from dataclasses import dataclass
from pathlib import Path
from typing import Callable

from scipy.integrate import quad
from scipy.special import ellipj, ellipk

ROOT = Path(__file__).resolve().parent.parent

# Physical constants, SI. IAU 2015 nominal values (Resolution B3) and CODATA 2018.
GM_SUN = 1.3271244e20  # m^3 s^-2
GM_JUP = 1.2668653e17  # m^3 s^-2
R_SUN = 6.957e8  # m
AU = 1.495978707e11  # m
C_LIGHT = 299_792_458.0  # m s^-1
DAY = 86_400.0  # s
YEAR = 365.25 * DAY  # Julian year, s
MEGAYEAR = 1e6 * YEAR
TWO_PI = 2.0 * math.pi
M_JUP_OVER_M_SUN = GM_JUP / GM_SUN

# Stellar structure constants as each paper normalises them.
K_Q_STORCH = 0.05  # S15/S16/A16: k_q (Claret & Gimenez 1992)
K_STAR_STORCH = 0.1  # S15/S16/A16: k_star
K_Q_A17 = 0.01  # A17: k_q,0
K_STAR_A17 = 0.06  # A17: k_star

# Reference values the published scalings are normalised to.
PLANET_TO_STAR = 1e-3  # S15 Eq. alpha0: 10^3 M_p / M_star = 1
COMPANION_TO_PLANET = 1e3  # S15 Eq. epsilon: M_b / 10^3 M_p = 1
SPIN_HAT_ALPHA0 = 0.05  # S15 Eq. alpha0: Omega_hat / 0.05
SPIN_HAT_EPSILON = 0.1  # S15 Eq. epsilon: Omega_hat / 0.1
A_B_EPSILON = 300.0 * AU  # S15 Eq. epsilon
A_IN_A17 = 0.1 * AU  # A17 Eqs. Aparameter, maxaout_zero_ecc
A_OUT_A17 = 10.0 * AU  # A17 Eq. Aparameter
A_OUT_GR = 100.0 * AU  # A17 Eq. epsilon_SRFs: a_out,eff / 100 AU
A_B_A16 = 100.0 * AU  # A16 Table 1: a_b,eff / 100 AU
P_STAR_A17 = 5.0 * DAY  # A17 Eq. Aparameter
P_STAR_A16 = 1.0 * DAY  # A16 Table 1: P_star / day

# A17: LK oscillations stop when eps_GR > (9/4)(1 - 5/3 cos^2 I0); I0 = 90 deg gives the largest distance.
GR_SUPPRESSION = 9.0 / 4.0
EPS_GR_PUBLISHED_PREFACTOR = 3e-2

# SL15 numerical example (Section 2.1): e0 = 0.01, theta_lb0 = 85 deg.
SL15_E0 = 0.01
SL15_THETA_LB0_DEG = 85.0
QUAD_LIMIT = 500  # subintervals for the sharply peaked x^{-3/2} integrand


def mean_motion(gm: float, a: float) -> float:
    return math.sqrt(gm / a**3)


def spin_precession_s15() -> float:
    """alpha_0 in units of 2 pi / yr (S15 Eq. alpha0; also S16)."""
    k_ratio = K_Q_STORCH / K_STAR_STORCH
    alpha0 = 1.5 * k_ratio * PLANET_TO_STAR * SPIN_HAT_ALPHA0 * math.sqrt(GM_SUN) * R_SUN**1.5 / AU**3
    return alpha0 / (TWO_PI / YEAR)


def adiabaticity_epsilon_s15() -> float:
    """epsilon = Omega_pl,0 / alpha_0 at the reference values of S15 Eq. epsilon (k_star / 2 k_q = 1, cos = 1)."""
    gm_star = GM_SUN
    omega_pl0 = 0.75 * mean_motion(gm_star, AU) * (COMPANION_TO_PLANET * PLANET_TO_STAR) * (AU / A_B_EPSILON) ** 3
    k_ratio = K_Q_STORCH / K_STAR_STORCH
    alpha0 = 1.5 * k_ratio * PLANET_TO_STAR * SPIN_HAT_EPSILON * math.sqrt(gm_star) * R_SUN**1.5 / AU**3
    return omega_pl0 / alpha0


def triple_adiabaticity_a17() -> float:
    """A17 script-A = |Omega_ps t_k j^3 / cos theta_sl| with all masses 1 M_sun (A17 Eq. Aparameter)."""
    omega_star = TWO_PI / P_STAR_A17
    omega_ps0 = 1.5 * (K_Q_A17 / K_STAR_A17) * omega_star * (R_SUN / A_IN_A17) ** 3
    t_k = (A_OUT_A17 / A_IN_A17) ** 3 / mean_motion(GM_SUN, A_IN_A17)
    return omega_ps0 * t_k


def gr_strength_a17(a_in: float = AU, a_out: float = A_OUT_GR) -> float:
    """eps_GR = 3 G m01^2 a_out,eff^3 / (c^2 a_in^4 m2), masses 1 M_sun. Definition from Liu, Munoz & Lai 2015;
    A17 Eq. epsilon_SRFs prints only the 3e-2 prefactor."""
    return 3.0 * GM_SUN * a_out**3 / (C_LIGHT**2 * a_in**4)


def perturber_distance_a17() -> float:
    """Largest a_out,eff (AU) allowing LK at I0 = 90 deg: eps_GR = 9/4 (A17 Eq. maxaout_zero_ecc)."""
    unit = gr_strength_a17(A_IN_A17, A_OUT_GR)
    return A_OUT_GR / AU * (GR_SUPPRESSION / unit) ** (1.0 / 3.0)


def perturber_distance_rounded_a17() -> float:
    """Same, but using A17's rounded 3e-2 prefactor for eps_GR instead of the constants."""
    unit = EPS_GR_PUBLISHED_PREFACTOR * (AU / A_IN_A17) ** 4
    return A_OUT_GR / AU * (GR_SUPPRESSION / unit) ** (1.0 / 3.0)


def kozai_window_a17(eta: float = 0.0, eps_srf: float = 0.0) -> float:
    """Upper edge (cos I0)_+ of the LK window (A17 Eq. Icrit); eta = eps_SRF = 0 should give sqrt(3/5)."""
    return (-eta + math.sqrt(eta**2 + 60.0 - 20.0 * eps_srf)) / 10.0


def spin_precession_a16() -> float:
    """|Omega_ps| prefactor in yr^-1 at M_p = 1 M_J, P_star = 1 d, a = 1 AU (A16 Eq. OmegaPS)."""
    omega_star = TWO_PI / P_STAR_A16
    rate = 1.5 * (K_Q_STORCH / K_STAR_STORCH) * M_JUP_OVER_M_SUN * (R_SUN / AU) ** 3 * omega_star
    return rate * YEAR


def initial_adiabaticity_a16() -> float:
    """A16 A_0 = |Omega_ps / Omega_L| at e = 0, with |cos theta_sl0 / sin 2 theta_lb0| = 1 (A16 Eq. adiabatic).
    Omega_L(e=0) = Omega_pl sin theta_lb = (3 / 8 t_k) |sin 2 theta_lb|."""
    t_k = (A_B_A16 / AU) ** 3 / mean_motion(GM_SUN, AU)
    omega_ps0 = spin_precession_a16() / YEAR
    return omega_ps0 * t_k * 8.0 / 3.0


def lk_roots(x0: float, h: float) -> tuple[float, float]:
    """x1 < x2 solving x^2 - (5 + 5h - 2 x0) x / 3 + 5h/3 = 0 (S15 Section 2.1)."""
    b = (5.0 + 5.0 * h - 2.0 * x0) / 3.0
    disc = math.sqrt(b * b - 4.0 * (5.0 / 3.0) * h)
    return (b - disc) / 2.0, (b + disc) / 2.0


def shape_factor_s15(theta_deg: float = SL15_THETA_LB0_DEG, e0: float = SL15_E0) -> float:
    """S15's function script-N = epsilon * alpha_bar / n_e, from the analytic LK solution (S15 Eqs. x, nmax).

    alpha_bar / n_e = (alpha_0 / 2 pi) * integral of x^{-3/2} dt over one LK period, and
    epsilon = Omega_pl,0 / alpha_0 with Omega_pl,0 = 3 sqrt(h) / (4 t_k). Time is in units of t_k.
    """
    x0 = 1.0 - e0**2
    h = x0 * math.cos(math.radians(theta_deg)) ** 2
    x1, x2 = lk_roots(x0, h)
    m = (x0 - x1) / (x2 - x1)
    big_k = float(ellipk(m))
    n_e = 6.0 * math.pi * math.sqrt(6.0) / (8.0 * big_k) * math.sqrt(x2 - x1)

    def inv_x32(theta: float) -> float:
        cn = float(ellipj(theta, m)[1])
        return (x0 + (x1 - x0) * cn * cn) ** -1.5

    # theta = (K / pi)(n_e t + pi) runs from K to 3K over one period; dt = pi / (K n_e) dtheta.
    integral, _ = quad(inv_x32, big_k, 3.0 * big_k, points=[2.0 * big_k], limit=QUAD_LIMIT)
    period_integral = math.pi / (big_k * n_e) * integral
    return 3.0 * math.sqrt(h) / (8.0 * math.pi) * period_integral


@dataclass(frozen=True)
class Check:
    check_id: str
    quantity: str
    source: str
    published: float
    compute: Callable[[], float]


CHECKS = [
    Check("C1_spin_precession", "alpha_0 prefactor [2pi/yr]", "S15 Eq. alpha0 (=S16 Eq. alpha0)", 1.19e-8,
          spin_precession_s15),
    Check("C2_epsilon_prefactor", "epsilon prefactor", "S15 Eq. epsilon", 1.17, adiabaticity_epsilon_s15),
    Check("C3_triple_adiabaticity", "A17 script-A prefactor", "A17 Eq. Aparameter", 58.0, triple_adiabaticity_a17),
    Check("C4_gr_strength", "eps_GR prefactor", "A17 Eq. epsilon_SRFs", 3e-2, gr_strength_a17),
    Check("C5_perturber_distance", "max a_out,eff [AU] (constants)", "A17 Eq. maxaout_zero_ecc", 19.6,
          perturber_distance_a17),
    Check("C6_perturber_rounded", "max a_out,eff [AU] (3e-2 prefactor)", "A17 Eq. maxaout_zero_ecc", 19.6,
          perturber_distance_rounded_a17),
    Check("C7_kozai_window", "(cos I0)_+ at eta = eps = 0", "A17 Eq. Icrit vs sqrt(3/5)", math.sqrt(3.0 / 5.0),
          kozai_window_a17),
    Check("C8_shape_factor", "script-N(85 deg, e0 = 0.01)", "S15 Eq. nmax and text", 0.98, shape_factor_s15),
    Check("C9_anderson_precession", "A16 Omega_ps prefactor [1/yr]", "A16 Eq. OmegaPS", 1.64e-7,
          spin_precession_a16),
    Check("C10_anderson_adiabaticity", "A16 A_0 prefactor", "A16 Eq. adiabatic", 0.07, initial_adiabaticity_a16),
]


def round_sig(value: float, figs: int) -> float:
    return float(f"{value:.{figs - 1}e}")


def run_check(check: Check, figs: int) -> dict[str, object]:
    computed = check.compute()
    return {
        "id": check.check_id,
        "quantity": check.quantity,
        "source": check.source,
        "published": check.published,
        "computed": computed,
        "rel_diff": computed / check.published - 1.0,
        "agrees": round_sig(computed, figs) == round_sig(check.published, figs),
    }


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def print_table(rows: list[dict[str, object]], figs: int) -> None:
    print(f"{'check':28} {'published':>11} {'computed':>13} {'rel diff':>9}  {figs} s.f.")
    for r in rows:
        mark = "ok" if r["agrees"] else "DIFFERS"
        print(f"{r['id']:28} {r['published']:>11.4g} {r['computed']:>13.6g} {r['rel_diff']:>+9.2%}  {mark}")


def main() -> None:
    config = tomllib.loads((ROOT / "config.toml").read_text())
    figs = config["calibration"]["sig_figs"]
    brief = ROOT / config["brief"]["path"]
    rows = [run_check(c, figs) for c in CHECKS]
    print_table(rows, figs)
    out = ROOT / config["calibration"]["results"]
    out.parent.mkdir(parents=True, exist_ok=True)
    record = {"brief": brief.name, "brief_sha256": sha256(brief), "sig_figs": figs, "checks": rows}
    out.write_text(json.dumps(record, indent=2) + "\n")
    print(f"\nbrief sha256 {record['brief_sha256']}\nwrote {out.relative_to(ROOT)}")


if __name__ == "__main__":
    main()
