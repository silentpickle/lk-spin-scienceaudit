# Step 11: N-body spot check of the secular equations (REBOUND)

| | |
|---|---|
| Answers brief | `docs/brief/planet-formation-and-spin-orbit-dynamics--how-does-stellar-spin-evolve-during-lidov-kozai-migration.md` |
| Brief SHA-256 | `11966760d9b6fc5e008207e46cd21e6b5dfb3c2160a9453c8ac546143dd1c3c3` |
| Work plan step | 11 (check with N-body that the secular approximation holds at the $e_{\max}$ reached; flag systems where the orbital period approaches the precession timescales) |
| Date | 2026-11-17 |
| Done by | Claude: a Sonnet 5.5 subagent (setup, runs, comparison) supervised by Opus 5.5 (spec, checkpointing, completion check, this write-up). Not checked by a physicist. |
| Verify | `shasum -a 256 -c docs/SHA256SUMS` |

## Result

In the parameter range of steps 4–6, the double-averaged secular equations reproduce the direct three-body evolution. The range is $a_b$ = 200–500 AU, $e_b=0$, $M_p$ = 1–5 $M_J$, $\theta_{lb}^0$ = 85–89.9°, and $a_0$ = 1.5 AU. Across 20 systems, with $e_{\max}$ up to 0.9992:

| Quantity | N-body against secular (GR only, no tides) |
|---|---|
| first $e_{\max}$ | $\lvert\Delta e\rvert\le4.6\times10^{-5}$; $1-e$ differs by a median 0.09% and at most 0.77% (system X2) |
| time of first $e_{\max}$ | within 0.075% |
| LK period (11 systems with two peaks) | within 0.05% |
| second $e_{\max}$ | within 0.9% in $1-e$ |
| $\theta_{lb}$ at $e_{\max}$ | within 0.33°, median about 0.02° |

The N-body $e_{\max}$ comes out slightly lower in every system, by under 1% in $1-e$. The gap widens as the time spent near $e_{\max}$
shrinks relative to the outer period, which fits the block averaging smoothing a sharp peak. The osculating envelope
agrees with the secular value to about $10^{-5}$ in $e$. None of the analytic breakdown indicators trip:

| Indicator | Range | Flag |
|---|---|---|
| time near $e_{\max}$ / outer period, $t_k\sqrt{1-e_{\max}^2}/P_{\rm out}$ | 40–642 (S16 fiducials 108 and 116) | none; double averaging fails as this approaches 1 |
| $P_{\rm orb}\,\alpha$ at $e_{\max}$, $P_\star$ = 0.3, 3, 30 d | at most $8.7\times10^{-4}$ | none; orbit-averaging the stellar torque holds |
| pericentre / disruption radius $2.7R_p(M_\star/M_p)^{1/3}$, full SRFs | at least 1.21 | none |

The planet's tidal and rotational bulge terms matter much more than any secular error. They raise $1-e_{\max}$ by
factors of 1.09 to 8.3 over GR alone. With GR only, five systems (R02, R10 and X0–X2) would pass inside the disruption radius.
The bulge terms are tested here only against their formula (step 4), not in the N-body.

## Method

- **Settings.** Stored in `config.toml` `[step11]`. Scripts: `scripts/step11_common.py`, `step11_nbody.py`, `step11_analyze.py`. Tests: `tests/test_step11.py` (5 pass).
- **Dependencies.** REBOUND 5.2.2 and REBOUNDx 5.1.0. Both were checked on PyPI before install: authors Rein and Tamayo, links to github.com/hannorein/rebound and github.com/dtamayo/reboundx.
- **Systems.** The two S16 fiducials (F89, F87); 15 draws (R00–R14, seed 11) from the step 5 ranges; and the 3 systems with the highest secular $e_{\max}$ among 2000 further draws (X0–X2). The planet's and companion's mean anomalies are random and recorded.
- **N-body.** Star, planet and companion integrated with IAS15, plus REBOUNDx `gr_full`. `gr_full` reproduces the 1PN apsidal precession $3G(M_\star+M_p)n/(c^2aj^2)$ to $10^{-5}$; `gr` comes out 0.3% low and `gr_potential` 0.5% low. The planet's orbit is taken in Jacobi coordinates. $e$ is averaged over blocks of one outer-binary period, and $e_{\max}$ is the peak of that block mean.
- **Secular side.** `lkspin.secular` with tides off. "GR only" sets $k_{2p}=k_{qp}=0$; "full SRFs" uses the step 4 values.
- **Run length.** Up to two $e_{\max}$. Systems whose second peak falls after $2.2\times10^7$ yr (R01, R06, R08, R09, R14) were run through the first only.
- **Checkpointing.** Full REBOUND and REBOUNDx state (about 16 KB) is saved with every output chunk, written atomically. A run resumed from a checkpoint gives output bit-identical to an uninterrupted run (test). REBOUND's warning about resetting function pointers when it loads a binary is expected: `reboundx.Extras(sim, file)` re-attaches the force on the next line.
- **Cost.** From 0.6 to 4 h of wall time per system on a shared, heavily loaded 16-core machine. Each system is one serial integration.

## Caveats

- **Range tested.** Only $e_b=0$ and $a_b\ge200$ AU. Eccentric or closer companions (octupole terms, smaller time near $e_{\max}$ relative to $P_{\rm out}$) are untested here. The octupole was checked separately for one system in findings 05.
- **Physics left out.** No dissipative tides and no planetary bulges in the N-body. This tests the LK-plus-GR secular dynamics only.
- **Phases.** One planet and companion phase per system.

## Process note

The first launch had no state checkpoints and a hard 4-hour cap. Five runs (R01, R06, R08, R09, R14) hit the cap and could
not be resumed. Their partial outputs are kept in `results/step11/legacy_capped/`. All five were rerun from scratch with
checkpointing and completed. Four of the reruns were themselves interrupted (by the session's 2-hour limit on background commands) and finished from their
checkpoints in the owner's terminal, losing no work.
