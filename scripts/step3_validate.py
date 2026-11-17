"""Validation numbers for lkspin.spin (work-plan step 3). Run: uv run python scripts/step3_validate.py"""

import math
import time

import numpy as np

from lkspin.lk_orbit import make_orbit
from lkspin.spin import ATOL_DEFAULT, H_INITIAL, RTOL_DEFAULT, initial_state, integrate_to, run_cycles
from lkspin.spin_check import scipy_p_at_minima

N_FULL = 500
N_SHORT = 20
N_CONV = 100
EPS_NORM = (1.0, 0.1, 0.01)
EPS_NONADIABATIC = 1e3
EPS_ADIABATIC = 1e-4
N_ADIABATIC = 20
EPS_CHAOTIC = 0.1
EPS_REGULAR = 0.015
TIGHTEN = 10.0

ORBIT = make_orbit()
ORB = ORBIT.params()


def run(eps, n, scale=1.0, p0=1.0):
    y0 = initial_state(p0, 0.0, ORB)
    return run_cycles(y0, n, ORBIT.period, ORBIT.alpha0(eps), ORB, RTOL_DEFAULT * scale, ATOL_DEFAULT * scale)


def theta_deg(res):
    return np.degrees(np.arccos(np.clip(res[:, 0], -1, 1)))


def check_norm():
    for eps in EPS_NORM:
        print(f"1. eps={eps}: max ||S|-1| over {N_FULL} cycles = {np.max(np.abs(run(eps, N_FULL)[:, 2] - 1)):.3e}")


def check_nonadiabatic():
    y = initial_state(1.0, 0.0, ORB)
    sz0, h, worst = y[2], H_INITIAL, 0.0
    for k in range(N_FULL):
        h = integrate_to(y, k * ORBIT.period, (k + 1) * ORBIT.period, h, ORBIT.alpha0(EPS_NONADIABATIC), ORB,
                         RTOL_DEFAULT, ATOL_DEFAULT)
        worst = max(worst, abs(math.degrees(math.acos(y[2]) - math.acos(sz0))))
    print(f"2. eps=1e3: max |theta_sb - theta_sb0| over {N_FULL} cycles = {worst:.4f} deg "
          f"(theta_sb0 = {math.degrees(math.acos(sz0)):.3f})")


def check_adiabatic():
    t = time.time()
    th = theta_deg(run(EPS_ADIABATIC, N_ADIABATIC))
    print(f"3. eps=1e-4, aligned start, {N_ADIABATIC} cycles: max theta_sl = {th.max():.5f} deg ({time.time() - t:.1f} s)")


def check_scipy(eps):
    y0 = initial_state(1.0, 0.0, ORB)
    t = time.time()
    ref = scipy_p_at_minima(y0, N_SHORT, ORBIT.period, ORBIT.alpha0(eps), ORB)
    mine = run(eps, N_SHORT)[:, 0]
    d = np.abs(mine - ref)
    print(f"4. eps={eps}: max |p_mine - p_scipy| over {N_SHORT} cycles = {d.max():.3e} (first cycle {d[0]:.1e}, "
          f"last {d[-1]:.1e}; scipy {time.time() - t:.0f} s)")


def check_convergence():
    a, b = theta_deg(run(EPS_REGULAR, N_CONV)), theta_deg(run(EPS_REGULAR, N_CONV, scale=1.0 / TIGHTEN))
    print(f"5. eps=0.015, {N_CONV} cycles, tol/10: max change in theta_sl = {np.max(np.abs(a - b)):.3e} deg "
          f"(range of theta_sl {a.min():.3f}..{a.max():.3f})")
    a, b = theta_deg(run(EPS_CHAOTIC, N_CONV)), theta_deg(run(EPS_CHAOTIC, N_CONV, scale=1.0 / TIGHTEN))
    print(f"5b. eps=0.1, {N_CONV} cycles, tol/10: max change in theta_sl = {np.max(np.abs(a - b)):.3e} deg")


def main():
    t0 = time.time()
    check_norm()
    check_nonadiabatic()
    check_adiabatic()
    check_scipy(EPS_REGULAR)
    check_scipy(EPS_CHAOTIC)
    check_convergence()
    print(f"total {time.time() - t0:.0f} s")


if __name__ == "__main__":
    main()
