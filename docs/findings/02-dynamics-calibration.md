# Step 3: calibrating the dynamics against Storch & Lai 2015

| | |
|---|---|
| Answers brief | `docs/brief/planet-formation-and-spin-orbit-dynamics--how-does-stellar-spin-evolve-during-lidov-kozai-migration.md` |
| Brief SHA-256 | `11966760d9b6fc5e008207e46cd21e6b5dfb3c2160a9453c8ac546143dd1c3c3` |
| Work plan step | 3 (reproduce the S15 bifurcation diagram and a surface of section; overlay the analytic chaotic zone) |
| Builds on | `docs/findings/01-definitions-and-calibration.md` (definitions, short names S15/S16/A16/A17) |
| Date | 2026-11-17 |
| Done by | Claude: Opus 5.5 (orbit module, review, comparisons, this write-up) with two Sonnet 5.5 subagents (spin integrator and runs; Fourier and resonance analysis). Not checked by a physicist. |
| Verify | `shasum -a 256 -c docs/SHA256SUMS` |

## Result

The numerical dynamics reproduce S15. The bifurcation diagram (S15 Fig. 12) and its zooms (Figs. 13 and 14) match
feature by feature. The surfaces of section (Fig. 9) match qualitatively. The shape functions, Fourier
coefficients and resonance widths (Figs. 2, 3 and 5) match, the widths to about 1%.

The one thing that does **not** reproduce is S15's analytic chaos-boundary curve (Fig. 11). The paper does not fully
specify its overlap rule, and no reading we tried matches the curve everywhere. Even so, the overlap theory does
predict where the aligned trajectory goes chaotic in our own numerical diagram (see the overlay).

The brief's step 3 check ("the chaos onset and the periodic islands visually match the published figure") passes.
The integrator is ready for step 4.

## Method

- **Orbit** (`lkspin/lk_orbit.py`): the analytic quadrupole LK solution of S15 §2.1 for $\theta_{lb}^0=85^\circ$ and
  $e_0=0.01$, with the Jacobi cn function computed by the AGM method. Time is in units of $t_k$. The LK period is
  $6.2748\,t_k$ and $e_{\max}=0.99365$. The spin problem then depends only on $\epsilon$, via $\alpha_0=\Omega_{pl,0}/\epsilon$.
- **Spin** (`lkspin/spin.py`): $d\hat S/dt=-\alpha(\hat S\cdot\hat L)\,\hat L\times\hat S$, integrated in the inertial frame
  (no coordinate singularity at the aligned start), with $\hat L(t)$ built from $\theta_{lb}(t)$ and the node $\Omega(t)$.
  It uses an adaptive Dormand–Prince 5(4) scheme compiled with numba, at rtol $10^{-10}$ and atol $10^{-12}$, and lands
  exactly on each eccentricity minimum $t=kP$. There it records $p=\cos\theta_{sl}$ and S15's phase
  $\phi=\operatorname{atan2}(\hat S\cdot\hat x',-\hat S\cdot\hat y')$.
- **Shape functions and resonances** (`lkspin/shape.py`): $\beta,\gamma,\psi$ are sampled on a uniform $\tau$ grid of
  $2^{16}$ points, inverted from a $2^{20}$-point grid in $t$. The resonance coefficient is
  $C_N=\beta_N+\gamma_N$ for $N>0$, $C_0=2\beta_0$, and $\beta_{|N|}-\gamma_{|N|}$ for $N<0$. There are approximate
  widths (S15 Eq. appwidth) and exact widths, the latter taken from the fixed points and separatrix of S15 Eq. oneharmonic.

| Run | Script | Grid | Runtime (16 cores) |
|---|---|---|---|
| Fig. 12 | `scripts/step3_bifurcation.py` | 500 ε, log-spaced from 0.008 to 1, 500 LK cycles, aligned start | 130 s |
| Fig. 13 | same | 400 ε in [0.04, 0.07] | 89 s |
| Fig. 14 | same | 400 ε in [0.025, 0.0286] | 168 s |
| Fig. 9 (top) | `scripts/step3_sections.py` | ε = 0.1, 0.05, 0.02; 40 initial conditions ($p_0$ from −0.98 to 0.98, $\phi_0\in\{0,\pi\}$); 1000 cycles | 78 s |
| Figs. 2, 3, 5, 11 | `scripts/step3_fourier.py`, `scripts/step3_exact.py` | Fourier coefficients up to $M=200$; boundary from ε = 0.003 to 1 | 8 s |

All runs are deterministic (no random numbers). None of the 1,420 runs hit the step limit. Data and plots are in `results/step3/`.

## Validation

| Check | Result |
|---|---|
| Analytic orbit vs direct integration of S15 Eqs. (ek)–(omk) over one period | $x$, $\cos\theta_{lb}$ and node agree to $<10^{-7}$ |
| $\lvert\hat S\rvert$ drift over 500 cycles | $6.5\times10^{-9}$ (ε = 1), $1.2\times10^{-7}$ (0.1), $4.3\times10^{-8}$ (0.01) |
| Non-adiabatic limit, ε = $10^3$: $\theta_{sb}$ conserved | drifts at most 0.093° over 500 cycles; the residual is physical at $\alpha/\Omega_{pl}\sim10^{-3}$ |
| Adiabatic limit, ε = $10^{-4}$, aligned start: $\theta_{sl}$ conserved | max 0.011° over 20 cycles |
| Same equations in scipy DOP853 (rtol $10^{-12}$), 20 cycles | max $\lvert\Delta p\rvert=2.3\times10^{-10}$ at ε = 0.015 (regular, the meaningful case); $2.0\times10^{-8}$ at ε = 0.1 |
| Tolerance 10× tighter, 100 cycles | $\theta_{sl}$ changes by $1.1\times10^{-5}$° (ε = 0.015) and $2.2\times10^{-5}$° (ε = 0.1) |
| $\mathcal N$ from the $\tau$ normalisation | 0.979220; S15 quotes 0.98 |
| S16's two forms of $\bar{\mathcal A}$: $1/(2\epsilon\beta_0)$ vs $N_{\max}/(2\hat\Omega_L)$ at ε = 0.1 | 11.277647 both, equal to $10^{-9}$; they are the same quantity algebraically, since $d\tau\propto\alpha\,dt$ |
| Fourier reconstruction of $\beta,\gamma,\psi$; doubling the grid | $4\times10^{-9}$; low-order coefficients change by $4\times10^{-9}$ |

`uv run pytest` runs 23 fast tests. `uv run python scripts/step3_validate.py` prints the slow checks above (36 s).

## Comparison with the published figures

The published figures were rendered from the arXiv EPS files and compared by eye. Their PNGs are not committed because
of publisher copyright; `scripts/fetch_sources.sh` and Ghostscript recreate them.

**Fig. 12, the bifurcation diagram: reproduced.**

| Feature | S15 | Ours |
|---|---|---|
| Spread stays near 0 | ε ≲ 0.02 | max θ_sl 6.1° for ε < 0.02 (median 0.9°) |
| Isolated spike near ε ≈ 0.027 | about 38° | 38.8° at ε = 0.02681. Only Fig. 14's finer grid resolves it; the 500-point Fig. 12 grid misses it. |
| Spikes, ε 0.04–0.07 | capped near 125–135° | capped at 122–132° (136.6° in the Fig. 13 grid) |
| Band, ε 0.08–0.1 | 20–30° | median of the per-ε maximum 24.6° |
| Spike at ε ≈ 0.1 | about 150° | about 156° |
| Widespread chaos begins | about 0.1–0.12 | 0.12: 96% of ε values in 0.12–0.2 exceed 90°, against 33% in 0.10–0.12 |

**Fig. 13 (0.04–0.07): reproduced.** The layout of spikes and gaps matches, including the first spikes near
0.0417 and 0.0436, the isolated spike near 0.0468, the gaps near 0.0585 and 0.0625, and the cap rising from about 117°
to about 135°.

**Fig. 14 (0.025–0.0286, vs 1/ε): reproduced in spacing, partly in height.** The chaotic zones are evenly spaced, about
1.0 apart in $1/\epsilon$ (at 35.25, 36.28 and 37.30, plus 38.3 and 39.3), and the sawtooth envelope matches. The heights
of the last two spikes differ: about 12° and 14° here, against roughly 17° and 10° by eye in S15.

**Fig. 9, top row, surfaces of section: qualitative match.** At ε = 0.1 there is a chaotic sea with islands near
$\cos\theta_{sl}=\pm0.97$ and −0.87 at $\phi\approx0$. At ε = 0.02 there are regular bands near 0.97 and 0.87 and a gap
near −0.6. At ε = 0.05 the empty band sits at about −0.60 to −0.65, against about −0.52 to −0.62 in S15, a shift of
about 0.05. This could be because different initial conditions were sampled; it was not investigated further.

**Fig. 2, shape functions; Fig. 3, Fourier coefficients: reproduced.** $C_N$ changes sign between the same pairs of $N$
as in S15 (−11/−10, −28/−27, −46/−45, −64/−63, −83/−82, 11/12, 29/30, 48/49, 66/67, 85/86). The filled and open circles
mark the same signs. $\beta_0=0.4434$ and $\psi_0=0.2951$.

**Fig. 5, resonance widths: reproduced to about 1%.** Values were read off the figure by pixel measurement at ε = 0.1.
Exact widths over paper values run from 0.97 to 1.01 for $N=-9$ to 9; the approximate widths are within 1–2% for
$\lvert N\rvert\le8$. Exact resonance centres sit about $\epsilon\psi_0\approx+0.03$ above the nominal $-N/N_{\max}$.

**Fig. 11, analytic chaos boundary: not reproduced.** S15 defines the outermost resonance $N^\pm_{\rm out}$ as the one
that "overlaps with the previous but not the next". It does not say whether the overlapping chain has to run
continuously from $N=0$. We built both readings, each with approximate and exact widths, and with the zone edge placed
at either the resonance centre or its outer separatrix edge.

| Reading | Where it matches S15 Fig. 11 | Where it fails |
|---|---|---|
| Largest $\lvert N\rvert$ anywhere, separatrix edge, exact widths (best) | ±0.28 at ε = 0.003; upper edge about 0.6 at 0.01, with jumps near 0.0165, 0.03 and 0.07; upper edge near 1 from ε ≈ 0.04; lower edge for ε ≲ 0.025; both edges saturate at ±1 for ε ≳ 0.14 | lower edge for ε ≈ 0.03–0.1: stays near −0.9 and only jumps near 0.06, where S15 shows a sloping edge with a jump to about −0.4 near 0.042; much finer sawtooth in ε than S15's points |
| Contiguous chain from $N=0$ | lower edge near 0.045 (−0.45) | far too narrow at small ε (±0.28 at 0.01, against S15's ±0.6); breaks at the first zero of $C_N$ |

Exact and approximate widths give nearly the same boundaries, so the width method is not what's missing. The
remaining gap probably comes from an unstated detail of S15's procedure, such as the ε grid, how islands touching the
poles are treated, or which edge is plotted.

**Overlay** (`results/step3/bifurcation_with_boundary.png`; best reading, shown as $\theta_{sl}=\arccos p$ on our Fig. 12):
the predicted upper edge of the chaotic zone tracks where the aligned trajectory goes chaotic.

- For ε < 0.02 the edge sits at 25–60° and the trajectory stays below 6°.
- The edge falls to about 0° across 0.04–0.06, where the isolated spikes appear.
- The edge rises to 30–40° over 0.07–0.1, where the trajectory stays in a 20–30° band beneath it.
- The edge reaches 0 at ε ≈ 0.12–0.14, where widespread chaos starts.

## Corrections and notes

1. **S15 Fig. 3 at $N=0$.** The point plotted is $\beta_0\approx0.44$, while the widths (Fig. 5, $w_0=0.83$) use $C_0=2\beta_0=0.89$. This is a labelling inconsistency in the paper and changes no result.
2. **Cycle count.** The brief asks for 2000 LK cycles per ε. S15's caption says about 500, and we used 500 to match.
3. **Which S15 figure.** The brief's target is S15 Fig. 12 (`\label{aligned}`), which has no short-range forces. S15 Fig. 1 is a different diagram, with short-range forces and sampled at eccentricity maxima.
4. **Overlap rule.** S15's resonance-overlap boundary is not fully specified. That matters for SQ4, which proposes to explain the remaining chaos with this theory. Any SQ4 claim built on the overlap criterion should state which reading it uses.

## Not done

- S15's surfaces of section with forcing truncated to $\lvert N\rvert\le N_{\max}$ harmonics (Fig. 9, bottom row).
- The ε = 0.05 section offset, and the heights of the last two spikes in Fig. 14, were not chased.
- No quantitative image-difference metric; comparisons are by eye, plus the Fig. 5 pixel measurement.
- Work plan steps 4 onward.
