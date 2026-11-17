"""S15 shape functions beta, gamma, psi, their Fourier coefficients, resonance widths and chaos boundary.

Reproduces the analytic side of Storch & Lai 2015 (arXiv:1411.0739): Fig. 2 (shape functions), Fig. 3 (beta_N+gamma_N),
Fig. 5 (widths) and Fig. 11 (predicted chaotic zone vs epsilon). Everything is epsilon-independent except the widths
and the boundary. Time is in units of t_k; tau is the rescaled time dtau ~ x^{-3/2} dt running 0..2 pi over one LK period.

Run the full reproduction with:  uv run python scripts/step3_fourier.py
"""

from dataclasses import dataclass

import numpy as np
from scipy.integrate import cumulative_simpson
from scipy.special import ellipj

from lkspin.lk_orbit import LKOrbit

TWO_PI = 2.0 * np.pi
N_DENSE = 2**20  # theta samples used to build tau(theta); the x^{-3/2} integrand spikes ~700x at e_max
N_TAU = 2**16  # uniform tau grid points
HALF = 0.5
TWO = 2.0
X_POWER = 1.5


@dataclass(frozen=True)
class Shapes:
    tau: np.ndarray
    theta: np.ndarray  # elliptic argument at each tau
    beta: np.ndarray
    gamma: np.ndarray
    psi: np.ndarray


def x_and_xdot(theta: np.ndarray, orb: LKOrbit) -> tuple[np.ndarray, np.ndarray]:
    """x = 1 - e^2 and dx/dt at elliptic argument theta (d cn/d theta = -sn dn)."""
    sn, cn, dn, _ = ellipj(theta, orb.m)
    x = orb.x0 + (orb.x1 - orb.x0) * cn**2
    dtheta_dt = orb.big_k * orb.n_e / np.pi
    xdot = (orb.x1 - orb.x0) * TWO * cn * (-sn * dn) * dtheta_dt
    return x, xdot


def theta_span(orb: LKOrbit) -> tuple[float, float]:
    """One LK period t in [0, P] maps to theta in [K, 3K] (x = x0 at both ends)."""
    return orb.big_k, 3.0 * orb.big_k


def dense_theta(orb: LKOrbit, n: int = N_DENSE) -> np.ndarray:
    lo, hi = theta_span(orb)
    return np.linspace(lo, hi, n + 1)


def omega_pl_ratio(x: np.ndarray, orb: LKOrbit) -> np.ndarray:
    """Omega_pl / Omega_pl0 (negative)."""
    return 1.0 - TWO * (orb.x0 - orb.h) / (x - orb.h)


def tau_of_theta_dense(orb: LKOrbit, n: int = N_DENSE) -> tuple[np.ndarray, np.ndarray]:
    """Dense theta grid and tau(theta) = 2 pi * int x^{-3/2} dt / int over a period (dt is proportional to dtheta)."""
    theta = dense_theta(orb, n)
    x, _ = x_and_xdot(theta, orb)
    cum = cumulative_simpson(x**-X_POWER, x=theta, initial=0.0)
    return theta, TWO_PI * cum / cum[-1]


def uniform_tau_thetas(orb: LKOrbit, n_tau: int = N_TAU, n_dense: int = N_DENSE) -> tuple[np.ndarray, np.ndarray]:
    """Uniform tau grid on [0, 2 pi) and the theta at each point (inverse of the monotone tau(theta))."""
    theta_d, tau_d = tau_of_theta_dense(orb, n_dense)
    tau = TWO_PI * np.arange(n_tau) / n_tau
    return tau, np.interp(tau, tau_d, theta_d)


def shape_functions(orb: LKOrbit, n_tau: int = N_TAU, n_dense: int = N_DENSE) -> Shapes:
    """beta, gamma, psi on a uniform tau grid (S15 Eqs. beta, gamma, psi divided by epsilon)."""
    tau, theta = uniform_tau_thetas(orb, n_tau, n_dense)
    x, xdot = x_and_xdot(theta, orb)
    cos_lb = np.sqrt(orb.h / x)
    sin_lb = np.sqrt(1.0 - orb.h / x)
    ratio = omega_pl_ratio(x, orb)
    thdot = xdot * cos_lb / (TWO * x * sin_lb)
    scale = x**X_POWER
    return Shapes(tau, theta, -ratio * scale * sin_lb, thdot * scale / orb.omega_pl0, -ratio * scale * cos_lb)


def cos_coefficients(f: np.ndarray, m_max: int) -> np.ndarray:
    """a_M, M = 0..m_max, for f = sum a_M cos(M tau) on a uniform periodic grid (a_0 is the mean)."""
    spec = np.fft.rfft(f) / f.size
    out = TWO * spec.real[: m_max + 1]
    out[0] = spec.real[0]
    return out


def sin_coefficients(f: np.ndarray, m_max: int) -> np.ndarray:
    """b_M, M = 0..m_max, for f = sum b_M sin(M tau); b_0 is set to 0 here (see coupling for gamma_0 = beta_0)."""
    spec = np.fft.rfft(f) / f.size
    out = -TWO * spec.imag[: m_max + 1]
    out[0] = 0.0
    return out


def coupling(beta_m: np.ndarray, gamma_m: np.ndarray, n_max: int) -> np.ndarray:
    """C_N for N = -n_max..n_max (index N + n_max): beta_N + gamma_N (N>0), 2 beta_0 (N=0), beta_|N| - gamma_|N| (N<0)."""
    n = np.arange(-n_max, n_max + 1)
    m = np.abs(n)
    sign = np.sign(n)
    c = beta_m[m] + sign * gamma_m[m]
    c[n_max] = TWO * beta_m[0]
    return c


def script_n(orb: LKOrbit, n: int = N_DENSE) -> float:
    """script-N = eps * alpha_bar / n_e = (3 sqrt(h) / 8 pi) int x^{-3/2} dt over one period."""
    theta = dense_theta(orb, n)
    x, _ = x_and_xdot(theta, orb)
    integral = orb.period * np.trapezoid(x**-X_POWER, theta) / (theta[-1] - theta[0])
    return 3.0 * np.sqrt(orb.h) / (8.0 * np.pi) * integral


def a_bar_from_beta0(beta0: float, epsilon: float) -> float:
    """S16 A_bar = 1 / (2 eps beta_0) (tau average)."""
    return 1.0 / (TWO * epsilon * beta0)


def a_bar_from_time_average(orb: LKOrbit, epsilon: float, script_n_value: float, n: int = N_DENSE) -> float:
    """S16 A_bar = N_max / (2 Omega_hat_L), Omega_hat_L = <|Omega_pl| sin th_lb>_t / n_e (time average); N_max = script_N/eps."""
    theta = dense_theta(orb, n)
    x, _ = x_and_xdot(theta, orb)
    integrand = np.abs(omega_pl_ratio(x, orb)) * orb.omega_pl0 * np.sqrt(1.0 - orb.h / x)
    omega_hat = np.trapezoid(integrand, theta) / (theta[-1] - theta[0]) / orb.n_e
    return (script_n_value / epsilon) / (TWO * omega_hat)


def resonance_momenta(n: np.ndarray, epsilon: float, script_n_value: float) -> np.ndarray:
    return -n * epsilon / script_n_value


def widths(c: np.ndarray, n: np.ndarray, epsilon: float, script_n_value: float) -> np.ndarray:
    """S15 Eq. appwidth: w_N = 2 [2 eps |C_N| sqrt(1 - p_N^2)]^{1/2}."""
    p = resonance_momenta(n, epsilon, script_n_value)
    return TWO * np.sqrt(TWO * epsilon * np.abs(c) * np.sqrt(np.clip(1.0 - p**2, 0.0, None)))


def overlaps_previous(w: np.ndarray, spacing: float) -> np.ndarray:
    """o[k] is True when resonance k+1 overlaps resonance k (w indexed by distance from N = 0 along one side)."""
    return HALF * (w[1:] + w[:-1]) > spacing


def one_side_widths(c_all: np.ndarray, n_coef: int, sign: int, n_top: int, epsilon: float, sn: float) -> np.ndarray:
    """Widths of N = 0, sign*1, ..., sign*n_top."""
    n = sign * np.arange(n_top + 1)
    return widths(c_all[n + n_coef], n, epsilon, sn)


def contiguous_outermost(w: np.ndarray, spacing: float) -> int:
    """Largest k such that resonances 1..k each overlap their inner neighbour (chain unbroken from N = 0)."""
    o = overlaps_previous(w, spacing)
    broken = np.flatnonzero(~o)
    return int(broken[0]) if broken.size else o.size


def anywhere_outermost(w: np.ndarray, spacing: float) -> int:
    """Largest k with k overlapping k-1 and not overlapping k+1 (a nonexistent k+1 counts as not overlapping)."""
    o = overlaps_previous(w, spacing)
    nxt = np.append(o[1:], False)
    hits = np.flatnonzero(o & ~nxt)
    return int(hits[-1]) + 1 if hits.size else 0


BOUNDARY_FIELDS = ("top_contig", "bottom_contig", "top_any", "bottom_any",
                   "top_contig_floor", "bottom_contig_floor", "top_any_floor", "bottom_any_floor")


@dataclass(frozen=True)
class Boundary:
    """cos th_sl bounds vs epsilon. top = from N-_out, bottom = from N+_out. contig = unbroken chain from N = 0,
    any = largest |N| satisfying both overlap conditions anywhere. *_floor uses p = -N / floor(N_max) (S15 Eq. pn
    with the floored N_max), so the last existing resonance sits exactly at p = +-1; the others use p = -N eps / script_N.
    """

    epsilon: np.ndarray
    top_contig: np.ndarray
    bottom_contig: np.ndarray
    top_any: np.ndarray
    bottom_any: np.ndarray
    top_contig_floor: np.ndarray
    bottom_contig_floor: np.ndarray
    top_any_floor: np.ndarray
    bottom_any_floor: np.ndarray


def outermost_pair(c_all: np.ndarray, n_coef: int, epsilon: float, sn: float, finder) -> tuple[int, int]:
    """(N+_out, N-_out) for one epsilon using the given chain-finding function."""
    n_top = min(int(np.floor(sn / epsilon)), n_coef)
    spacing = epsilon / sn
    wp = one_side_widths(c_all, n_coef, 1, n_top, epsilon, sn)
    wn = one_side_widths(c_all, n_coef, -1, n_top, epsilon, sn)
    return finder(wp, spacing), -finder(wn, spacing)


def boundary_point(c_all: np.ndarray, n_coef: int, epsilon: float, sn: float) -> list[float]:
    """Values in BOUNDARY_FIELDS order for one epsilon."""
    n_floor = max(np.floor(sn / epsilon), 1.0)  # N_max < 1: only N = 0 exists; p stays 0 because N_out = 0
    out = []
    for per_n in (epsilon / sn, 1.0 / n_floor):
        for finder in (contiguous_outermost, anywhere_outermost):
            n_pos, n_neg = outermost_pair(c_all, n_coef, epsilon, sn, finder)
            out += [-n_neg * per_n, -n_pos * per_n]
    return out


def chaos_boundary(c_all: np.ndarray, n_coef: int, sn: float, eps: np.ndarray) -> Boundary:
    pts = np.array([boundary_point(c_all, n_coef, e, sn) for e in eps])
    return Boundary(eps, *(pts[:, k] for k in range(len(BOUNDARY_FIELDS))))


# ---------------------------------------------------------------------------------------------------------------
# Exact single-resonance fixed points and separatrix widths (S15 Eq. oneharmonic), vectorised over N.
#   H_N(p, phi) = -p^2/2 + a p - (eps/2) sqrt(1-p^2) C_N cos(phi),   a = -N eps / script_N + eps psi_0
# With p = sin(u), dH/dp = 0 at fixed phi (cos phi = s) reads g(u) = sin u - a - k tan u = 0, k = (eps/2) C_N s.
# Stable point: s = -sign(C_N) so k < 0 and g is strictly increasing in u: exactly one root, found by bisection.
# Saddle: s = +sign(C_N), k > 0. g then has three roots in (-pi/2, pi/2): the physical one between the turning points
# cos^3(u) = k (where g is increasing), and two spurious ones near p = +-1 (the S15 footnote). Bisecting only on the
# interval between the turning points excludes the spurious ones by construction.
# ---------------------------------------------------------------------------------------------------------------
BISECT_ITER = 80
HALF_PI = np.pi / TWO


@dataclass(frozen=True)
class ExactResonances:
    centre: np.ndarray  # p of the stable fixed point
    lo: np.ndarray  # lower separatrix edge at the stable phi
    hi: np.ndarray  # upper separatrix edge
    has_saddle: np.ndarray  # False where no interior saddle exists (separatrix taken through the nearer pole)

    @property
    def width(self) -> np.ndarray:
        return self.hi - self.lo


def bisect_vec(f, lo: np.ndarray, hi: np.ndarray) -> np.ndarray:
    """Vectorised bisection for f increasing through zero on [lo, hi]."""
    lo, hi = lo.copy(), hi.copy()
    for _ in range(BISECT_ITER):
        mid = HALF * (lo + hi)
        up = f(mid) > 0.0
        hi = np.where(up, mid, hi)
        lo = np.where(up, lo, mid)
    return HALF * (lo + hi)


def g_of_u(u: np.ndarray, a: np.ndarray, k: np.ndarray) -> np.ndarray:
    return np.sin(u) - a - k * np.tan(u)


def stable_point(a: np.ndarray, kmag: np.ndarray) -> np.ndarray:
    """Stable fixed point p (unique root, k = -kmag <= 0)."""
    edge = HALF_PI * (1.0 - 1e-12)
    ones = np.ones_like(a)
    u = bisect_vec(lambda v: g_of_u(v, a, -kmag), -edge * ones, edge * ones)
    return np.sin(u)


def saddle_point(a: np.ndarray, kmag: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Physical saddle p (k = +kmag) and a mask for where it exists (kmag < 1 and g changes sign between turning points)."""
    kk = np.clip(kmag, 1e-300, 1.0 - 1e-12)
    third = 1.0 / 3.0
    u_hi = np.arccos(kk**third)
    u_lo = -u_hi
    exists = (kmag < 1.0) & (g_of_u(u_lo, a, kk) < 0.0) & (g_of_u(u_hi, a, kk) > 0.0)
    u = bisect_vec(lambda v: g_of_u(v, a, kk), u_lo, u_hi)
    return np.sin(u), exists


def h_stable_phi(p: np.ndarray, a: np.ndarray, kmag: np.ndarray) -> np.ndarray:
    """H at the stable phi: -p^2/2 + a p + kmag sqrt(1-p^2)."""
    return -HALF * p**2 + a * p + kmag * np.sqrt(np.clip(1.0 - p**2, 0.0, None))


def separatrix_edges(pc: np.ndarray, a: np.ndarray, kmag: np.ndarray, level: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    """Roots of H_stable(p) = level on each side of pc (H_stable is concave); clipped to +-1 if the level is not reached."""
    def f_lower(p):  # increasing on [-1, pc]
        return h_stable_phi(p, a, kmag) - level

    def f_upper(p):  # decreasing on [pc, 1]; negate to make it increasing for the bisection
        return -(h_stable_phi(p, a, kmag) - level)

    ones = np.ones_like(pc)
    lower = bisect_vec(f_lower, -ones, pc)
    upper = bisect_vec(f_upper, pc, ones)
    lower = np.where(f_lower(-ones) >= 0.0, -ones, lower)
    upper = np.where(h_stable_phi(ones, a, kmag) >= level, ones, upper)
    return lower, upper


def exact_resonances(c: np.ndarray, n: np.ndarray, epsilon: float, sn: float, psi0: float) -> ExactResonances:
    """Exact stable-point location and separatrix extent for resonances n with couplings c (arrays of equal shape).

    Resonances with C_N = 0 get zero width. Where no interior saddle exists (islands touching a pole, has_saddle False)
    the separatrix is taken to be the H level of the pole nearer the centre, so the island runs from that pole to the
    mirrored root.
    """
    a = -n * epsilon / sn + epsilon * psi0
    kmag = HALF * epsilon * np.abs(c)
    pc = stable_point(a, kmag)
    ps, exists = saddle_point(a, kmag)
    level = np.where(exists, -HALF * ps**2 + a * ps - kmag * np.sqrt(1.0 - ps**2), -np.inf)
    pole = np.where(pc >= 0.0, 1.0, -1.0)
    level = np.where(exists, level, h_stable_phi(pole, a, kmag))  # no saddle: the separatrix runs through the nearer pole
    lo, hi = separatrix_edges(pc, a, kmag, level)
    zero = kmag == 0.0
    lo, hi = np.where(zero, pc, lo), np.where(zero, pc, hi)
    return ExactResonances(pc, lo, hi, exists | zero)


# ---------------------------------------------------------------------------------------------------------------
# Chaos boundary from exact centres and widths.
# ---------------------------------------------------------------------------------------------------------------
EXACT_FIELDS = tuple(f"{side}_{reading}_{where}" for where in ("centre", "edge") for reading in ("contig", "any")
                     for side in ("top", "bottom"))


def exact_overlap_flags(centre: np.ndarray, width: np.ndarray) -> np.ndarray:
    """o[k-1] True when resonance k overlaps k-1: |p_k - p_{k-1}| < (w_k + w_{k-1}) / 2 (arrays ordered by |N|)."""
    return np.abs(np.diff(centre)) < HALF * (width[1:] + width[:-1])


def contiguous_from_flags(o: np.ndarray) -> int:
    broken = np.flatnonzero(~o)
    return int(broken[0]) if broken.size else o.size


def anywhere_from_flags(o: np.ndarray) -> int:
    nxt = np.append(o[1:], False)
    hits = np.flatnonzero(o & ~nxt)
    return int(hits[-1]) + 1 if hits.size else 0


def exact_side(c_all: np.ndarray, n_coef: int, sign: int, n_top: int, eps: float, sn: float, psi0: float):
    n = sign * np.arange(n_top + 1)
    return exact_resonances(c_all[n + n_coef], n, eps, sn, psi0)


def exact_boundary_point(c_all: np.ndarray, n_coef: int, eps: float, sn: float, psi0: float) -> list[float]:
    """Values in EXACT_FIELDS order for one epsilon. top comes from N-_out, bottom from N+_out.

    "centre": the stable fixed point of the outermost resonance. "edge": the outer separatrix edge of that resonance.
    """
    n_top = min(int(np.floor(sn / eps)), n_coef)
    pos = exact_side(c_all, n_coef, 1, n_top, eps, sn, psi0)
    neg = exact_side(c_all, n_coef, -1, n_top, eps, sn, psi0)
    flags_pos, flags_neg = (exact_overlap_flags(r.centre, r.width) for r in (pos, neg))
    out = []
    for where in ("centre", "edge"):
        for finder in (contiguous_from_flags, anywhere_from_flags):
            kp, kn = finder(flags_pos), finder(flags_neg)
            top = neg.centre[kn] if where == "centre" else neg.hi[kn]
            bottom = pos.centre[kp] if where == "centre" else pos.lo[kp]
            out += [top, bottom]
    return out


def chaos_boundary_exact(c_all: np.ndarray, n_coef: int, sn: float, psi0: float, eps: np.ndarray) -> dict[str, np.ndarray]:
    pts = np.array([exact_boundary_point(c_all, n_coef, e, sn, psi0) for e in eps])
    return {f: pts[:, k] for k, f in enumerate(EXACT_FIELDS)}
