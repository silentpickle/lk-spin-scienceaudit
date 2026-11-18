# Which part of spin-orbit feedback does what (extends SQ1 case b and SQ5)

| | |
|---|---|
| Answers brief | `docs/brief/planet-formation-and-spin-orbit-dynamics--how-does-stellar-spin-evolve-during-lidov-kozai-migration.md` |
| Brief SHA-256 | `11966760d9b6fc5e008207e46cd21e6b5dfb3c2160a9453c8ac546143dd1c3c3` |
| Extends | findings 05 (SQ1 case b, SQ5); closes its "which channel" item |
| Date | 2026-11-18 |
| Done by | Claude: a Sonnet 5.5 subagent (code, runs, analysis, and the physically consistent split) supervised by Opus 5.5 (spec, review, this write-up). Not checked by a physicist. |
| Reproduce | `uv run python scripts/step6_run.py --case {b_L,b_E,b_P,b_F}`, then `scripts/job4_validate.py` and `scripts/job4_analyze.py`. Each case takes 6–10 minutes on 14 cores. |

## Result

Feedback from the oblate star does two different jobs, and different terms do them:

1. **Lost migration comes from apsidal precession.** The precession of the planet's pericentre caused by the star's quadrupole lowers
   $e_{\max}$ and stops migration. 347 systems migrate without feedback (case a) but not with it (case b). With only
   this term on, 315 of them (91%) already fail to migrate.
2. **Changes of outcome class come from the orbit rotating about the stellar spin.** This is where the orbital plane, carrying both $\mathbf L$ and $\mathbf e$,
   precesses about $\hat S$. On its own it reproduces case b's SQ5 curve almost exactly: the class-change fraction crosses 10%
   at $S/L=0.104$ and 50% at 0.59, against case b's 0.106 and 0.59. Apsidal precession alone changes only about 11% of classes, and never reaches
   50%.
3. **Advection disappears through both together.** Apsidal precession alone leaves no class IV systems, mostly by stopping their
   migration. Orbit rotation alone leaves 27 of the original 162, most of the rest turning chaotic. Full feedback leaves 1.

## Splitting the terms consistently

The A16 feedback terms (findings 05) are $d\mathbf L/dt=\Omega_{ps}\mathbf S\times\hat L$ and
$d\mathbf e/dt=-\dot\omega_\star[\cos\theta_{sl}\hat S\times\mathbf e+\tfrac12(1-5\cos^2\theta_{sl})\hat L\times\mathbf e]$.

**The first split was wrong.** It was specified (by Opus) as "$d\mathbf L/dt$ only" (b_L) against "$d\mathbf e/dt$ only" (b_E). That is inconsistent. The
$\cos\theta_{sl}\hat S\times\mathbf e$ part of $d\mathbf e/dt$ is the partner of $d\mathbf L/dt$, since together they rotate the
whole orbit about $\hat S$. Separating them lets $\mathbf e$ leave the orbital plane: the final $\lvert\cos(\mathbf l,\mathbf e)\rvert$ exceeds 0.1
in 55–58% of runs, against $10^{-10}$ in case b. The subagent caught this and added a consistent split, and the conclusions above rest on it:

| Case | Contains | $\mathbf l\cdot\mathbf e$ | Conserved |
|---|---|---|---|
| b_P | apsidal precession only: $-\dot\omega_\star\tfrac12(1-5\cos^2\theta_{sl})\hat L\times\mathbf e$ | 0 (below $6\times10^{-8}$) | $\mathbf L$ is untouched by the star |
| b_F | orbit rotation about $\hat S$: $d\mathbf L/dt$ plus $-\dot\omega_\star\cos\theta_{sl}\hat S\times\mathbf e$ | 0 | $\mathbf S+\mathbf L$ (to $4\times10^{-15}$ with LK off) |

The b_L and b_E results are kept for completeness but flagged as unphysical.

**Checks.**
- New multipliers default to 1.
- Three stored case b pilot systems are bit-identical.
- The step 4 pilot is unchanged, and all 57 tests pass.
- Conservation for each variant: `results/job4/validation.json`.

## Numbers (`results/job4/summary.json`, `class_change_vs_sl.png`)

All 1028 systems form hot Jupiters in case (a). Classes are the v2 classes of findings 04, applied unchanged. Accuracy is for S16's
two-parameter map on the systems that still classify. "θ only" is the $\theta_{lb}^0$-only baseline on the held-out half.

| Case | Hot Jupiter / timeout / mixed | I / II / III / IV / V | Map accuracy (95% CI) | θ only | 10% / 50% class-change $S/L$ |
|---|---|---|---|---|---|
| b (both) | 681 / 204 / 143 | 445 / 45 / 190 / 1 / 0 | 57.4% (53.7–60.9) | 64.0% | 0.106 / 0.59 |
| **b_P** (apsidal) | 588 / 397 / 43 | 416 / 64 / 108 / 0 / 0 | 76.9% (73.3–80.3) | 69.7% | 0.20 / never (max 0.48) |
| **b_F** (rotation) | 807 / 79 / 142 | 418 / 99 / 263 / 27 / 0 | 55.1% (51.7–58.5) | 49.2% | 0.104 / 0.59 |
| b_L (unphysical) | 908 / 7 / 113 | 323 / 96 / 472 / 14 / 3 | 40.9% | 49.1% | 0.05 / 0.79 |
| b_E (unphysical) | 741 / 121 / 166 | 432 / 26 / 283 / 0 / 0 | 47.2% | 56.4% | 0.064 / 0.52 |

**Attribution, consistent split.**
- Of the 116 systems that change class from a to b (among 551 classified in a, b, b_P and b_F), orbit rotation alone reproduces 70 (60%) and apsidal precession alone reproduces 25 (22%).
- Of the 347 lost migrators, 315 are lost under apsidal precession alone and 187 under rotation alone.

**Not a fair comparison.** Under apsidal precession alone the map scores 77%. That is higher than in case (a), but it is measured on the 588 systems that still migrate, a heavily selected population. It does not mean the map works better with that physics.

## Caveats

- **Classifier.** The v2 classifier was hand-validated without feedback only (findings 05). Classes with feedback carry that caveat.
- **Additivity.** The two terms are not additive. Neither alone reproduces case b's near-total loss of class IV.
- **Time horizon.** Systems that never migrate are capped at 5 Gyr, as in steps 4–6.
