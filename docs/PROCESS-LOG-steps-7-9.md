# Process log: steps 7–9

Brief SHA-256: `11966760d9b6fc5e008207e46cd21e6b5dfb3c2160a9453c8ac546143dd1c3c3`

One entry per unit of work, newest last. Each entry says what was asked, what was run, what came back, and what it
decided. Status words used throughout:

- **checked**: confirmed by a computation in this repo, with the script named.
- **lead**: plausible, not yet tested. Do not build on it without the named test.
- **failed**: tested and wrong. Kept here so it is not tried again.

The plan these entries follow is the "Fable plan" table in `HANDOFF-steps-7-9.md`.

## 2026-11-18, Stage 0: orient

- **Done.** `uv sync`, `scripts/fetch_sources.sh`, `uv run pytest` (57 pass), `shasum -a 256 -c docs/SHA256SUMS` (all OK).
  Fetched and checksum-verified `prep_sq4`, `prep_sq2` and `step5` archives.
- **Read.** README, both handoff files, the brief, findings 01–04 and 07, `docs/DATA.md`, `lkspin/shape.py`,
  `lkspin/spin.py`, `scripts/prep_sq4.py`, `scripts/prep_common.py`, and S16 `ms.tex` l.312–481 and l.625–1305.
- **Written.** The "Fable plan" section of `HANDOFF-steps-7-9.md`.
- **Not written.** A longer reasoning notebook was started and abandoned; this log and the findings files replace it.
- **Status of the plan's two starting ideas.** Both are **leads**. The only evidence so far is hand arithmetic on two
  of the four fiducials (plan table, row 7b), which is why those fiducials are in the development set and the
  odd-numbered `v2_*` systems are held out.
- **Deviation from the kickoff prompt.** The numerical shape functions (9d's construction) are built first as unit T0,
  because the FTLE data that 9a compares against has short-range forces on and `shape.py` has them off.
- **Observation from existing data, no new run.** In `results/prep_sq4/omega_test.npz`, the $N_{\max,0}=0.3$ column is
  chaotic (FTLE 0.05–1.15 per LK cycle) for the 14 orbits with $\bar{\mathcal A}_0$ = 0.53–0.72 and regular for all 26
  with $\bar{\mathcal A}_0\ge0.74$. Whether a separatrix of the $N=0$ Hamiltonian even exists at those values is what
  unit 9b.1 computes.

## 2026-11-18, T0: numerical shape functions

- **Question.** Can $\beta,\gamma,\psi(\tau)$ and the coefficients $N_{\max}$, $\bar{\mathcal A}$, $\bar{\mathcal B}$,
  $C_N$ be built from any one LK cycle of a `secular.py` orbit, accurately enough to stand in for `shape.py`?
- **Outcomes and what each means.**
  - All three gates pass: the construction and its sign conventions are right; 9a can run.
  - Gate 1 (match `shape.py`, SRFs off) fails: a sampling or normalisation error in the construction; fix before going on.
  - Gate 3 (integrating S16 Eq. kozaiham with the built shapes reproduces `secular.py`) fails while gate 1 passes: a
    sign or phase-origin convention is wrong ($\gamma$'s sign, $\phi$'s zero, or where $\tau=0$ sits); fix before going on.
- **Who.** Code and runs by a Sonnet 5.5 subagent from a written spec; Fable wrote the definitions and Hamilton's
  equations, then read `lkspin/numshape.py` line by line against them (signs of $\hat\beta,\hat\gamma,\hat\psi$, the
  node frame, $\dot\phi=+\partial H/\partial p$).
- **Result: checked.** `uv run python scripts/t0_validate.py` (3 s) writes `results/step9/t0_validation.json`.
  - Gate 1 (85°, SRFs off, against `shape.py`): largest coefficient difference $5\times10^{-9}$ of $\mathcal N\beta_0$
    (limit $10^{-4}$); $N_{\max}\epsilon$ against $\mathcal N$: $1.3\times10^{-9}$.
  - Gate 2 (87° and 89°, SRFs on): coefficients agree with `prep_common.first_cycle_shape` to $2\times10^{-4}$ or
    better; doubling both grids changes nothing above $3\times10^{-12}$; odd/even symmetry holds to $10^{-8}$.
  - Gate 3 (Hamiltonian against `secular.py`, 6 cases, 20 cycles): $\lvert\Delta p\rvert\le2\times10^{-7}$,
    $\lvert\Delta\phi\rvert\le10^{-6}$ **when the reference lands on minima found at rtol $10^{-12}$**. With minima
    found at the production rtol $10^{-10}$ the differences are $1$–$3\times10^{-5}$ in $p$ and up to $1.4\times10^{-4}$
    in $\phi$, over the limit I had set ($10^{-5}$, $10^{-4}$) in 5 of 6 cases.
  - **Judgement call on gate 3.** I count it as passed, on the tight-tolerance variant. Reason: a wrong sign or phase
    origin would give differences of order 1, not $10^{-5}$, and the residual shrinks 100-fold when only the reference's
    landing times are made more accurate. The limit was not loosened; both variants are in the JSON.
- **By-product.** Eccentricity-minimum times found at rtol $10^{-10}$ jitter by about 30 yr in $1.41\times10^{7}$ yr
  per cycle ($2\times10^{-6}$), because the minimum is very flat at $e_0=0.01$. This is harmless for the existing FTLE
  scans (both twins land on the same times). Anything that needs an exact one-cycle map (unit 9b.2) should integrate the
  Hamiltonian in $\tau$, which is exactly periodic, rather than land `secular.py` on minima.
- **Numbers that later units use** (first cycle, $P_\star$-independent):

  | Orbit | $N_{\max}P_\star$ [d] | $\bar{\mathcal A}P_\star$ [d] | $\hat\Omega_L$ | $\hat\psi_0/\hat\Omega_L$ |
  |---|---|---|---|---|
  | 85°, SRFs off | 8.5993 | 9.9038 | 0.4341 | 0.6655 |
  | 87°, SRFs on | 3.0026 | 7.9018 | 0.1900 | 0.1894 |
  | 89°, SRFs on | 3.2719 | 24.3448 | 0.0672 | 0.0660 |

- **Decided.** T0 is done; 9a runs next on these shapes.
