# Preparation data for SQ2 and SQ4, and sensitivity of step 4 to its assumed parameters

| | |
|---|---|
| Answers brief | `docs/brief/planet-formation-and-spin-orbit-dynamics--how-does-stellar-spin-evolve-during-lidov-kozai-migration.md` |
| Brief SHA-256 | `11966760d9b6fc5e008207e46cd21e6b5dfb3c2160a9453c8ac546143dd1c3c3` |
| What this is | Compute-only groundwork for steps 7 (SQ2) and 9 (SQ4), plus a sensitivity study that closes open items in findings 03. This file provides data; it contains no theory. |
| Date | 2026-11-18 |
| Done by | Claude: a Sonnet 5.5 subagent (code and runs) supervised by Opus 5.5 (spec, review, this write-up). Not checked by a physicist. |
| Reproduce | `uv run python scripts/prep_sq2.py`, `prep_sq4.py`, `prep_sens4.py` (settings in `config.toml` `[prep_sq2]`, `[prep_sq4]`, `[sens4]`). Each takes minutes on 14 cores. Tests in `tests/test_prep.py`. |

## 1. SQ2: bimodal branch probabilities (`results/prep_sq2/`)

The ensembles cover 20 class I systems: the four S16 fiducials, plus 16 step 5 class I systems stratified over $\bar{\mathcal A}_0$ = 0.05–0.96 (seed 12). Without feedback, every ensemble except B shares one orbit integration.

| Ensemble | Phase variable | Verdict |
|---|---|---|
| A | $P_\star$ jitter ±5%, 200 members | Converged for the fiducials: 50, 100 and 200 members agree to about 0.03 in probability and about 1° in the peaks. At 87°/20 d, 50 members are not enough. |
| B | initial $\omega_0$ uniform, 200 members | **Poor sampler.** At 89° all members land on one branch. Don't use it for branch probabilities. |
| C | spin tilted 10°, 20° or 30° from $\hat L$, random azimuth (fiducials) | Gives the initially misaligned peak shift. |
| W | $P_\star$ jitter ±30% (extra) | Finds both branches in 12 of the 13 one-branch step 5 systems. It also averages over $N_{\max,0}$, so it is not a fixed-$P_\star$ probability. |

**Fiducials (ensemble A):**

| System | $p_{\rm high}$ [95% CI] | Low peak | High peak |
|---|---|---|---|
| 89°, 30 d | 0.47 [0.40, 0.54] | 41.6° | 133.1° |
| 89°, 27 d | 0.57 [0.51, 0.64] | 38.9° | 135.8° |
| 87°, 30 d | 0.345 [0.28, 0.41] | 53.0° | 119.2° |
| 87°, 20 d | 0.43 [0.365, 0.50] | 51.5° | 120.6° |

There is no third mode. The branch probabilities are not 0.5, so a theory of the peaks also has to predict the branch weights.

**Initially misaligned spins (ensemble C).**
- At $\theta_{sl,0}=10°$ the peaks don't move beyond the confidence intervals.
- At 30° the **high peak moves 13–14° toward 90° at 89°** (133° → 120° at 30 d; 136° → 122° at 27 d). At 87° it moves only 0–6°.
- The low peak moves about ±5°, and $p_{\rm high}$ stays within about 0.2.

S16 said only that the peaks "move closer to 90°". The size of that shift depends on the inclination.

**Low-$\bar{\mathcal A}_0$ systems.** At ±5% jitter, 13 of the 16 step 5 class I systems sit on one branch. Five of them (v2_00, 01, 02, 05, 10) stay on one branch even at ±30%. Their exit phase barely varies, so a different phase variable is needed. For each system, `summary.json` records $\bar{\mathcal A}_0$, $N_{\max,0}$, $\hat\Omega_L$, first-cycle $e_{\min}$ and $e_{\max}$, the LK period and $t_{hj}$.

## 2. SQ4: chaos maps in the frozen-shape system (`results/prep_sq4/`)

**Setting.** Short-range forces on, tides off, $a$ fixed at 1.5 AU, no feedback. This is S16's Hamiltonian at $t=0$.

**Chaos measure.** The finite-time Lyapunov exponent (FTLE) is the slope of ln(separation) from a twin offset by $10^{-9}$ rad, measured until the separation saturates. The data split cleanly: regular cases fall below $1.1\times10^{-3}$ per LK cycle and chaotic ones above 0.03, with the threshold set at 0.02.

**Validation.** With short-range forces off, the scan reproduces the chaotic ε ranges of the step 3 integrator exactly, and matches S15 Fig. 12's spread in 93% of ε values.

**FTLE scan, 300 $P_\star$ values in 1–40 d:**

| | Chaotic $P_\star$ (FTLE > 0.02) |
|---|---|
| 87° | 71 of 300: bands at 1.28–1.38, 2.07–3.04, 4.6–4.9 and 7.4–9.7 d, plus isolated points at 1.83, 6.2 and 7.0 d. Regular strands at 3.07–4.56 d, 4.97–6.1 d, and everything above 9.8 d. |
| 89° | 3 of 300 (2.78–2.85 d) |

**This bears directly on S16's "unexplained chaos".** In the frozen-shape system the chaos at 87° extends to 9.7 d, where $\bar{\mathcal A}_0\lesssim1$ and $N_{\max,0}\lesssim1$. It is therefore already present at $t=0$, not created by tidal drift. Above 9.8 d the frozen system is regular. S16's scatter there (in its Fig. 8 class I/III region) must then come from separatrix crossing during migration.

**Spread is not chaos.** Regular orbits can swing $\theta_{sl}$ by about 50° (87°, $P_\star>10$ d). Use FTLE, not spread, to identify chaos.

**Sections.** There are sections at 8 periods for 87° (`prep_sq4_sections.png`):
- 3.0–4.7 d: a chaotic sea around a central island;
- 5.5 and 7.5 d: a shrinking sea;
- 9 d: island chains with a thin separatrix layer;
- 12 and 30 d: regular.

**$\hat\Omega_L$ test, holding $N_{\max,0}$ fixed while varying $\theta_{lb}^0$ (`prep_sq4_omega_test.png`):**

| $N_{\max,0}$ | Result |
|---|---|
| 1.5 | sharp threshold: chaotic exactly for $\hat\Omega_L\ge0.21$ |
| 0.7 | the same threshold, with two exceptions at the largest $\hat\Omega_L$ |
| 0.3 | chaos intermittent and weak |

So $\hat\Omega_L$ controls chaos at fixed $N_{\max,0}$ in the frozen system too. The threshold here (about 0.21) differs from findings 04's 0.13–0.15. That earlier value came from ensemble scatter, with tides on and a drifting orbit. $\hat\Omega_L$ depends only on the orbit shape, set by $\theta_{lb}^0$: it is 0.28 at 85° and 0.007 at 89.9°.

## 3. Sensitivity of step 4 to the assumed parameters (`results/sens4/`)

Each assumed parameter in `config.toml` `[step4]` was varied one at a time. The table shows 89°, with coefficients relative to S16's (3.34 and 25.5):

| Variant | $N_{\max,0}$ coeff. | $\bar{\mathcal A}_0$ coeff. | $\theta_{sl,f}$ at 30 d | $t_{hj}$ (Myr) |
|---|---|---|---|---|
| baseline | 0.943 | 0.952 | 41.6° | 1442 |
| χ = 5 / 20 | 0.946 / 0.941 | 0.955 / 0.947 | 41.9° / 41.5° | 2882 / 723 |
| $R_p$ = 0.9 / 1.2 $R_J$ | 1.077 / 0.714 | 0.999 / 0.863 | 134.6° / 131.6° | 1119 / 2959 |
| $e_0$ = 0.001 / 0.05 | 0.977 / 0.920 | 0.949 / 0.958 | 41.6° / 42.5° | 1446 / 1423 |
| bulge factor 15 | 0.900 | 0.944 | 42.2° | 2226 |
| planet spin 20 h / 5 h | 1.132 / 0.619 | 1.036 / 0.819 | 135.1° / 96° | 416 / no hot Jupiter in 15 Gyr |
| $R_\star$ = 0.9 / 1.1 $R_\odot$ | 0.687 / 1.254 | 0.694 / 1.268 | 128.9° / 27.9° | 1442 / 1442 |

1. **The coefficient offset is explained by $R_\star$.** Both coefficients scale as $R_\star^3$, so $R_\star\approx1.02R_\odot$ matches S16 for both while keeping their ratio. $R_p$ moves the two coefficients unequally and would break the ratio. χ and $e_0$ have no effect, and a bulge factor of 15 makes it worse. This also supports the $\tfrac{15}{2}$ choice in findings 03.
2. **The 30 d branch.** $R_p$, $R_\star$ and the initial planet spin all flip 30 d onto the high branch (129–135°). But 29.87 d flips with it, whereas in S16 30 d is high (134°) and 29.87 d low (42°). The branch at a given $P_\star$ is generic phase sensitivity, so no single assumed value reproduces S16's pair. This matches ensemble A's $p_{\rm high}\approx0.5$.
3. **Migration time.** χ, $R_p$ and the planet spin change the migration time several-fold. A 5 h initial planet spin stops migration within 15 Gyr.

## Caveats

- **Grouped runs.** At 87°, grouped integrations (many spins on one orbit) differ from step 4's one-spin runs by up to 116° in the bimodal and striated regions. Outcomes there depend on the step sequence, which is another face of the phase sensitivity. Compare distributions there, not individual values.
- **Ensemble W** mixes phase with $N_{\max,0}$ variation.
- **Settings.** The FTLE threshold (0.02 per cycle) and the jitter widths are ASSUMED, as noted in `config.toml`.
