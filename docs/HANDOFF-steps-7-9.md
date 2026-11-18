# Handoff: work-plan steps 7–9 (SQ2, SQ3, SQ4)

For the agent that supervises the theory steps; Claude Fable 5.2 is suggested. Written 2026-11-17 by the Opus 5.5
session that did steps 1–6, 10 and 11.

Brief SHA-256: `11966760d9b6fc5e008207e46cd21e6b5dfb3c2160a9453c8ac546143dd1c3c3`

## Read first, in this order

1. `README.md`, especially "For agents": the conventions, signing, and what never to commit.
2. The brief: `docs/brief/…md`, sections "What is being asked", "What a finished answer looks like" (SQ2–SQ4), and work-plan steps 7–9.
3. `docs/findings/01` through `04`, plus 05 and 06 once written for steps 6, 10 and 11. Findings 01 lists where the brief is wrong. Don't inherit those errors.

## What these three steps are

These steps are mostly theory, and the compute is cheap. Each one asks for an analytic or semi-analytic prediction, checked
against ensembles of integrations. A wrong derivation will look just as confident as a right one, so the supervisor should
own the derivations and the comparisons, and hand the code and runs to cheaper subagents (Sonnet 5.5 has worked well
here), reviewing their physics term by term. Every number in this repo so far came from an integrator that was validated
against an independent method first. Keep that standard.

| Step | Sub-question | Deliverable (from the brief) |
|---|---|---|
| 7 | SQ2: bimodal peaks | A formula $\theta_{sl,f}^\pm=F(\text{separatrix area at }t=0)$, overlaid on $\theta_{sl,f}$ histograms from phase ensembles (about 200 per system). It must also predict how the peaks shift for $\theta_{sl,0}\neq0$; S16 says only that they move toward 90°. |
| 8 | SQ3: advection capture | Capture probability into resonance $N$ against the drift rate of the Hamiltonian coefficients, from separatrix-crossing theory (Neishtadt-type; check the textbook result), compared with Monte Carlo fractions. Plus $\theta_{sl,f}$ against $P_\star$ with the striations labelled by predicted $N$. |
| 9 | SQ4: unexplained chaos | Surfaces of section at fixed $\bar{\mathcal A}$ with $N_{\max}<1$ in the chaotic period range. Either identify the overlapping resonances, or show there is no overlap and find the mechanism (for example repeated separatrix crossing of the $N=0$ islands, or a Mathieu-type parametric resonance near $N_{\max}\sim\tfrac12$). Finite-time Lyapunov exponents against $P_\star$ should match the spread in $\theta_{sl,f}$. Otherwise, measure the boundary empirically. |

## What exists and is validated

| Tool | What it does | Validated against |
|---|---|---|
| `lkspin/lk_orbit.py` | analytic quadrupole LK orbit (S15 §2.1), no SRFs, time in $t_k$ | direct integration of S15's orbit equations, to 1e-7 |
| `lkspin/spin.py` | non-dissipative spin on the analytic orbit; samples $(p,\phi)$ at eccentricity minima, so it gives surfaces of section directly | scipy to 2e-10; limiting cases; reproduces S15 Figs. 9 and 12–14 (findings 02) |
| `lkspin/shape.py` | shape functions $\beta,\gamma,\psi(\tau)$, Fourier coefficients, approximate **and exact** single-resonance centres and widths (`exact_resonances`), and overlap-chain boundaries | S15 Figs. 2, 3 and 5 to about 1% (findings 02) |
| `lkspin/secular.py` | full secular orbit (vector $\mathbf L$, $\mathbf e$), short-range forces, weak-friction tides, planet spin, stellar spins. Opt-in: eccentricity-extremum events, the winding diagnostic $W$, back-reaction, octupole, spin-down. | analytic orbit, A16's $e_{\max}$, conservation laws, REBOUND (findings 03; 06 when written) |
| `lkspin/regime.py` | step 5 features and the v2 outcome classifier | 50 hand-labelled cases, 17/18 on the held-out set (findings 04) |

Runs are cheap: one migrating system takes 0.5–10 s. The full step 5 sample of 3000 systems with 11 spins each took
170 s on 10 cores.

## Findings that matter for steps 7–9

- **SQ4: a third parameter controls chaos.** Among systems S16's map calls advective, the chaotic fraction rises from 7% to 100% as $\hat\Omega_L=N_{\max,0}/(2\bar{\mathcal A}_0)$ goes from 0.005 to 0.32. At fixed $N_{\max,0}$ it is 2–24% below the median and 81–87% above it. The split is near $\hat\Omega_L\approx0.13$–0.15, and $\hat\Omega_L$ tracks $\theta_{lb}^0$ with correlation −0.98 (findings 04). This is a concrete target: derive why, for example from resonance widths growing with $\epsilon\lvert C_N\rvert$ against spacing $1/N_{\max}$, which together set an overlap condition in $\hat\Omega_L$.
- **SQ4: the overlap theory doesn't reproduce its own published curve.** S15's analytic chaos boundary (Fig. 11) doesn't reproduce under any reading of its overlap rule (contiguous chain or anywhere; approximate or exact widths; edge at centre or separatrix). The best reading fails on the lower edge for ε ≈ 0.03–0.1 (findings 02). Overlaid on the numerical bifurcation diagram, though, it predicts the aligned trajectory's fate well. Any SQ4 claim built on overlap must say which reading it uses.
- **SQ2: get the phase sampling right.** A ±2% change in $P_\star$ is often not enough to spread class I ensembles over both branches; many land wholly on one (findings 04). Use wider $P_\star$ jitter or another phase variable, and show the branch probability has converged. Branch values to match: at 89° about 133–142° and 32–42°; at 87° about 119° and 52° (findings 03).
- **SQ3: N labels exist.** The winding $W$ per LK cycle locks to the resonance order during advection. `regime.best_lock_fraction` gives the $N$ and its locked fraction per system. Classes IV.1, IV.2, IV.3 and IV.4+ have 68, 51, 22 and 21 systems in step 5.
- **The analytic tools assume the orbit shape is frozen and has no SRFs.** `shape.py` builds shape functions from S15's analytic orbit, with no SRFs and fixed $a$. SQ2–SQ4 in the realistic setting need them for the SRF-modified orbit, and as they drift while $a$ decays. The brief's "usual failure" warns that freezing the shape functions while $a$ decays overstates advection. One route is to build $\beta,\gamma,\psi$ numerically from `secular.py`'s first-cycle (or per-epoch) orbit. Validate against `shape.py` in the no-SRF limit first.
- **Use the right paper's normalisations.** Our $N_{\max,0}$ and $\bar{\mathcal A}_0$ coefficients are about 5% below S16's quoted values, while their ratio agrees to 1% (findings 03). Use our own consistently, and don't mix in S16's numbers.
- **A16 has three equation-level errors.** Its tidal-bulge precession is inconsistent by a factor of 2 (findings 03). Its printed octupole $d\mathbf L/dt$ has $(\mathbf e\cdot\hat L_b)$ where $(\mathbf e\cdot\hat e_b)$ belongs, in two places (findings 05). Neither changes steps 7–9 directly. They're a reminder to check every transcribed equation against a conservation law.

## Suggested approach

1. **Step 9 first.** It is the best defined and has the strongest lead, $\hat\Omega_L$. In the non-dissipative setting (`spin.py` plus an orbit), build sections and FTLE maps across the $P_\star$ (or ε) range where S16 sees unexplained chaos with $N_{\max}\lesssim1$. Then test the candidate mechanisms in turn.
2. **Step 7 next.** The adiabatic invariant $J=\oint p\,d\phi$ of the $N=0$ separatrix at $t=0$, tracked as $a$ decays, against phase ensembles of about 200. This needs the drifting shape functions from the last bullet above.
3. **Step 8 last.** Neishtadt capture probability from the rates of change of the separatrix areas, against ensemble capture fractions (the $W$ lock), with striation labels on a $\theta_{sl,f}$-against-$P_\star$ plot.

## House rules

- **Write up each step.** Each one gets a findings file (`docs/findings/NN-*.md`) that opens with the brief's hash and covers method, validation numbers, comparison with the papers, corrections, assumptions and what was not done.
- **Sign before committing.** Run `git add`, then `scripts/sign.sh`, then `git add docs/SHA256SUMS`, and commit with the brief hash in the message.
- **Settings and outputs.** Settings go in `config.toml` with a source or ASSUMED for each value; outputs go in `results/stepN/`; never commit `sources/`.
- **Report honestly.** Never tune physics to match a figure. Hold out validation cases before looking at them. When a reproduction fails, say so and keep both readings.

## State at handoff (2026-11-17)

- **Steps 1–6, 10 and 11 are complete and committed** (findings 01–06). Step 11 shows the secular equations hold in this parameter range, matching direct N-body to under 1% in $1-e_{\max}$. Long N-body runs are checkpointed and resumable (`results/step11/state/`).
- **SQ1 is a negative result under fuller physics.** With back-reaction on, S16's map is no better than $\theta_{lb}^0$ alone in all three cases (findings 05). Advection (IV) all but disappears, from 162 systems to 0–2, and a third of the systems stop migrating in every ensemble member. This changes the stakes for SQ3: capture into resonance matters in S16's idealised model and at slow spins, but with back-reaction on, the advected population is tiny. State that up front in SQ3.
- **SQ5:** with back-reaction on, the outcome class changes for more than 10% of systems above $S/L\approx0.11$ and more than 50% above 0.59, with $S/L$ at the first $e_{\max}$.
- **A16's printed octupole $d\mathbf L/dt$ (Eq. Lvec, l.1763) has $(\mathbf e\cdot\hat L_b)$ in place of $(\mathbf e\cdot\hat e_b)$ twice.** The corrected form is implemented and checked against conservation laws and REBOUND (findings 05).
- **The classifier was never re-validated with back-reaction.** v2 was hand-validated only without it. Check it before relying on IV versus II there.
- **Large intermediates are not in git, because they regenerate.** These are `results/**/series*.npz` and the `chunks/` directories. `scripts/step5_run.py` takes about 3 minutes on an idle machine; `scripts/step6_run.py --case {b,c,d}` takes about 45 minutes under load. Classification (`step5_classify.py`, `step6_analyze.py`) needs the series.

## Prepared data (2026-11-18, findings 07)

Compute-only groundwork is ready:
- **SQ2.** Converged branch probabilities and peaks for the four fiducials and 16 class I systems, and the misaligned-spin peak shifts (`results/prep_sq2/`). Use $P_\star$ jitter, not $\omega_0$, as the phase variable. Low-$\bar{\mathcal A}_0$ systems need wider sampling.
- **SQ4.** FTLE scans, sections and an $\hat\Omega_L$ test in the frozen-shape system (`results/prep_sq4/`). S16's 87° chaos is already present at $t=0$ up to 9.7 d. At fixed $N_{\max,0}$, $\hat\Omega_L\ge0.21$ is a sharp chaos threshold. Use FTLE, not $\theta_{sl}$ spread, to identify chaos.
- **Sensitivity.** Our 5% coefficient offset from S16 is $R_\star$ ($\propto R_\star^3$). The 30 d branch is generic phase sensitivity.
