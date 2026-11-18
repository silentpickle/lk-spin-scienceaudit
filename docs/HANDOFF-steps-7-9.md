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
- **Large intermediates are not in git.** These are `results/**/series*.npz` and the `chunks/` directories. Since 2026-11-18 they are on the `v0.1-first-results` release: fetch them with `scripts/fetch_data.sh`, as described in `docs/DATA.md`, rather than regenerating them. Classification (`step5_classify.py`, `step6_analyze.py`) needs the series.
- **Kickoff prompt for the Fable session:** [`HANDOFF-fable-prompt.md`](HANDOFF-fable-prompt.md).
- **2026-11-18, Fable session, paused.** Done: Stage 0 (the "Fable plan" below) and unit T0, `lkspin/numshape.py`, which builds the shape functions from any LK cycle of a `secular.py` orbit and is validated (`results/step9/t0_validation.json`). Learned: first-cycle $\hat\Omega_L$ and $\hat\psi_0/\hat\Omega_L$ are 0.190 and 0.189 at 87°, 0.067 and 0.066 at 89° (SRFs on). Not done: 9a onward. Two responses were stopped by an automated safety classifier; what that cost and the open questions are in [`PROCESS-LOG-steps-7-9.md`](PROCESS-LOG-steps-7-9.md). Next: 9a.

## Prepared data (2026-11-18, findings 07)

Compute-only groundwork is ready:
- **SQ2.** Converged branch probabilities and peaks for the four fiducials and 16 class I systems, and the misaligned-spin peak shifts (`results/prep_sq2/`). Use $P_\star$ jitter, not $\omega_0$, as the phase variable. Low-$\bar{\mathcal A}_0$ systems need wider sampling.
- **SQ4.** FTLE scans, sections and an $\hat\Omega_L$ test in the frozen-shape system (`results/prep_sq4/`). S16's 87° chaos is already present at $t=0$ up to 9.7 d. At fixed $N_{\max,0}$, $\hat\Omega_L\ge0.21$ is a sharp chaos threshold. Use FTLE, not $\theta_{sl}$ spread, to identify chaos.
- **Sensitivity.** Our 5% coefficient offset from S16 is $R_\star$ ($\propto R_\star^3$). The 30 d branch is generic phase sensitivity.

- **Feedback channels (findings 08).** Apsidal precession from the stellar quadrupole is what stops migration. The rotation of the orbit about $\hat S$ is what changes outcome classes; on its own it reproduces the full SQ5 curve. When splitting feedback terms, keep $\mathbf L$ and $\mathbf e$ rotating together (cases b_P / b_F), or $\mathbf e$ leaves the orbit plane.

## Fable plan (2026-11-18, Stage 0)

Written by the Claude Fable 5.2 session before any new computation. What was run, why, and what each result decided
is logged unit by unit in [`PROCESS-LOG-steps-7-9.md`](PROCESS-LOG-steps-7-9.md); the reasoning for each result goes in
that step's findings file. Order: 9, 7, 8. Every run is minutes on 14 cores unless marked **long**; long runs are
chunked and resumable and will be handed over as a terminal command.

**Two ideas the plan rests on (both are leads until the tests below pass):**
- S16's $N=0$ Hamiltonian is exactly Colombo's top, with $\eta\cos I=\epsilon\psi_0=1/\bar{\mathcal B}$ and $\eta\sin I=\epsilon\beta_0=1/(2\bar{\mathcal A})$. Its separatrix splits the sphere into three zones with areas $A_I+A_{II}+A_{III}=4\pi$. The aligned spin sits exactly on the separatrix at $\bar{\mathcal A}=1$ when $\psi_0=0$.
- Averaging over an LK cycle is valid only while $\hat L$ moves little per cycle, and $\hat L$ moves $2\pi\hat\Omega_L$ rad per cycle. So $\hat\Omega_L$ measures how badly the $N=0$ picture fails.

**Deviation from the kickoff prompt.** 9a asks for `shape.exact_resonances` with no SRFs, compared against `results/prep_sq4/`. That data has SRFs on (at 87° the SRFs change $1-e_{\max}$ sevenfold), so the comparison would be between two different orbits. The numerical shape functions of 9d are therefore built first, as unit T0, and 9a uses them.

| Unit | Decisive computation | Prediction fixed in advance | Kill criterion |
|---|---|---|---|
| **T0** numerical shapes [W] | `lkspin/numshape.py`: $\beta,\gamma,\psi(\tau)$, $N_{\max}$, $\bar{\mathcal A}$, $\bar{\mathcal B}$, $C_N$ from one LK cycle of a `secular.py` orbit | — (validation gate) | Gate: SRFs off at 85° must match `shape.py` ($\beta_0,\psi_0,\mathcal N$, $C_N$ for $\lvert N\rvert\le10$) to $10^{-4}$; integrating S16 Eq. kozaiham with the numerical shapes must reproduce `secular.py` $(p,\phi)$ at minima over 50 cycles. Nothing downstream runs until it passes. |
| **9a** overlap test [W] | Exact first-order resonances from T0 shapes against existing FTLE: 87° and 89° scans (600 points) and the $\hat\Omega_L$ test (120 points). No new integrations. | Aligned spin is chaotic iff the resonance attached to $p=1$ overlaps a neighbour, under each of the four readings of findings 02 | Balanced accuracy below 0.75 on the $N_{\max,0}<1$ points for every reading: first-order overlap is not the mechanism; go to 9b |
| **9b.0** scoreboard map [W] **long?** | Frozen FTLE of the aligned spin over $\hat\Omega_L$ (40 cached orbits, 85–89.9°) × $N_{\max,0}\in[0.15,3]$ (about 120 $P_\star$). 3-orbit pilot first. Odd-numbered orbits are held out before anything is plotted. | — (this is what every candidate is scored on) | — |
| **9b.1** $N=0$ separatrix layer [W] | Where is the pole relative to the Colombo separatrix: $\bar{\mathcal A}_c$ (separatrix born) and $\bar{\mathcal A}_{\rm det}$ (pole on it), from T0's $\beta_0,\psi_0$ | Chaos at $N_{\max,0}<1$ only in a window around $\bar{\mathcal A}_{\rm det}$, widening steeply with $\hat\Omega_L$; none where no separatrix exists | Chaos where $\bar{\mathcal A}_0<\bar{\mathcal A}_c$, or under 80% of chaotic points inside the window. The existing $N_{\max,0}=0.3$ column (chaos at $\bar{\mathcal A}_0$ = 0.53–0.72) may already kill it there. |
| **9b.2** Mathieu / half-integer [W] | Newton-solve the fixed point of the one-LK-cycle map that continues from the pole; Greene residue $R$ and rotation number | Chaos tongues at rotation number ½, i.e. $R>1$ (period doubling) near $N_{\max}\approx\tfrac12,\tfrac32$ | $0<R<1$ throughout the chaotic $N_{\max,0}<1$ points |
| **9b.3** secondary resonances [W] | Winding number of the aligned orbit against $P_\star$ on the 9b.0 map | Each chaotic band is bracketed by locked plateaus $W=m/k$ whose measured widths overlap | No plateaus next to the bands |
| 9b fallback | Empirical boundary $\hat\Omega_{L,\rm crit}(N_{\max,0})$ with bootstrap intervals | — | — |
| **9c** threshold [W] | From the surviving mechanism, or the map: a migrating system is chaotic if its track ($N_{\max}$ rising as $a$ decays) crosses any chaotic band | $\hat\Omega_{L,\rm crit}^{\rm tides}=\min_{N\ge N_{\max,0}}\hat\Omega_{L,\rm crit}^{\rm frozen}(N)$, computed before comparing with step 5's 0.13–0.15 on a held-out half of the predicted-IV systems | Off by more than 0.03: report the gap and test the drift of $\hat\Omega_L$ itself in 9d |
| **9d** realistic [W] | T0 shapes at about 6 epochs along the 87° and 89° migration tracks; frozen FTLE at each | Union of chaotic $P_\star$ bands over epochs matches the scatter bands in `results/step4/scan_{87,89}.npz` | Bands disagree by more than 20% in $P_\star$ |
| **7a** drifting areas [W] | $(N_{\max},\bar{\mathcal A},\bar{\mathcal B})(t)$ per LK cycle along each of the 20 tracks (orbit only), and $A_I,A_{II},A_{III}(t)$; $J_0$ = area of the orbit through the pole at $t=0$ | — (gate: areas sum to $4\pi$; closed form at $\psi_0=0$; measured area of the true spin's stroboscopic orbit is constant to a few % until crossing) | Gate fails: the invariant is not conserved, SQ2's premise is wrong, report as negative |
| **7b** peak formula [W] | $A_{II}(t_c)=J_0$ gives $t_c$; $\cos\theta^{\rm low}=1-A_I(t_c)/2\pi$, $\cos\theta^{\rm high}=-1+A_{III}(t_c)/2\pi$; reduces to $\pm J_0/4\pi$ when $\psi_0=0$ | Hand check at Stage 0: 89°/30 d gives ±0.72 against a measured half-difference of 0.715; 27 d gives 0.756 against 0.748. Development set: the 4 fiducials and even-numbered `v2_*`. Held out: odd-numbered `v2_*`. | Median error above 5° on the held-out set |
| **7c** branch weights [W] | Henrard: $p_{\rm high}=\dot A_{III}/(\dot A_I+\dot A_{III})$ at $t_c$ | Inside ensemble A's 95% interval | Outside it for more than 1 of the 4 fiducials |
| **7d** misaligned start [W] | $J_0(\theta_{sl,0},\text{azimuth})$ distribution pushed through 7b | Shift of the high peak at 30°: large at 89°, small at 87° | Wrong inclination trend, or shift off by more than 5° |
| **8a** capture fractions [W] | Locked fraction per IV.$N$ class from `results/step5/series.npz`; new ensembles only where under 20 members | — | — |
| **8b** capture probability [W] | Henrard/Neishtadt ratio of separatrix-area rates for the single-resonance $H_N$ with drifting $C_N(t),N_{\max}(t)$; exit when $A_{\rm res}(t)=J$ | Formula first checked on a pendulum with slowly varying length (textbook case), then compared with 8a on a held-out half | Toy check fails: implementation wrong. Held-out fractions outside binomial 95% intervals for more than a third of bins: report as negative. |
| **8c** striations [W] | $\theta_{sl,f}$ against $P_\star$ with predicted $N$ per striation, against measured $W$ locks | Labels agree | More than 10% mislabelled |

Pauses as in the kickoff prompt: after 9e, 7e, 8d, and before the release.
