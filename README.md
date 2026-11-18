# Lidov-Kozai stellar spin: working one ScienceAudit open question

Casey Handmer's [ScienceAudit physics map](https://science.caseyhandmer.com/physics/) lists open questions in
theoretical astrophysics that "anyone with spare Claude Max tokens could answer today", each with a `.md` research
brief. This repo works one of them, from the subfield *Planet formation & spin-orbit dynamics*:

> **How does the host star's spin evolve during Lidov-Kozai migration of a hot Jupiter, and what sets the final
> spin-orbit misalignment?**

The question as asked is the brief in [`docs/brief/`](docs/brief/). Everything here was done by Claude (Opus 5.5,
with Sonnet 5.5 subagents) working with the repo owner. **No physicist has checked it.** Read each result as a
lead to verify.

## Status

The brief's work plan, step by step:

| Step | What | Status | Write-up |
|---|---|---|---|
| 1 | Recover the definitions marked garbled | done; also found 8 errors in the brief | [01](docs/findings/01-definitions-and-calibration.md) |
| 2 | Recompute the published prefactors | done; all 5 asked for reproduce | [01](docs/findings/01-definitions-and-calibration.md) |
| 3 | Reproduce Storch & Lai 2015's bifurcation diagram and sections | done; analytic chaos-boundary curve does not reproduce | [02](docs/findings/02-dynamics-calibration.md) |
| 4 | Reproduce Storch et al. 2016's $\theta_{sl,f}$ vs $P_\star$ with short-range forces and tides | done; all five outcome classes at the published periods; minutes, not CPU-hours | [03](docs/findings/03-tidal-migration-reproduction.md) |
| 5 | Regime map from 1028 systems; outcome classifier checked by hand | done; S16's map 64% accurate vs 47% for $\theta_{lb}^0$ alone; chaos set by a third parameter $\hat\Omega_L$ | [04](docs/findings/04-regime-map-and-classifier.md) |
| 6 (SQ1) | The map under back-reaction, octupole and spin-down | done; **negative result**: with back-reaction the map is no better than $\theta_{lb}^0$ alone, and advection all but vanishes | [05](docs/findings/05-sq1-full-physics-and-sq5-feedback.md) |
| 10 (SQ5) | Back-reaction threshold | done; class changes exceed 10% at $S/L\approx0.11$ and 50% at 0.59 | [05](docs/findings/05-sq1-full-physics-and-sq5-feedback.md) |
| 11 | N-body spot check (REBOUND), 20 systems | done; secular equations match direct N-body to under 1% in $1-e_{\max}$ and 0.1% in timing; no breakdown flags | [06](docs/findings/06-nbody-spot-check.md) |
| prep | Compute groundwork for SQ2 and SQ4; sensitivity to assumed parameters | done | [07](docs/findings/07-prep-data-and-sensitivity.md) |
| job 4 | Which part of feedback does what | done; apsidal precession stops migration, rotation of the orbit about the spin changes classes | [08](docs/findings/08-feedback-channels.md) |
| 7–9 (SQ2–SQ4) | Theory: bimodal peaks, capture probabilities, chaos mechanism | **paused** after planning and one tool (numerical shape functions, validated); no physics result yet; see [handoff](docs/HANDOFF-steps-7-9.md) | [process log](docs/PROCESS-LOG-steps-7-9.md) |

## Quick start

    uv sync
    scripts/fetch_sources.sh        # pinned arXiv LaTeX sources, hash-checked (kept out of git)
    uv run pytest                   # fast checks
    shasum -a 256 -c docs/SHA256SUMS
    scripts/fetch_data.sh           # optional: the raw series, chunks and N-body checkpoints (about 4.8 GB), checksum-verified

**Data:** every result in the findings is backed by data in git or in the
[`v0.1-first-results` release](https://github.com/silentpickle/lk-spin-scienceaudit/releases/tag/v0.1-first-results).
[`docs/DATA.md`](docs/DATA.md) lists each archive: what it holds, which script wrote it, which scripts read it, its
checksum, and how to take the work in a different direction without rerunning anything.

Each findings file lists the scripts that produce its results; every script writes to `results/`.
Tools and versions: [`docs/TOOLING.md`](docs/TOOLING.md).

## Layout

| Path | What |
|---|---|
| `docs/brief/` | the research brief, unchanged |
| `docs/findings/NN-*.md` | one write-up per step: method, validation numbers, comparison with the papers, corrections, what was not done |
| `docs/SHA256SUMS` | hashes of every tracked file plus the pinned sources |
| `lkspin/` | the physics: `lk_orbit.py` (analytic LK orbit), `spin.py` (spin integrator), `shape.py` (Fourier and resonance analysis), `secular.py` (full secular orbit + tides + stellar spins, step 4 onward) |
| `scripts/` | one script per result; `calibrate.py`, `step3_*.py`, `step4_*.py`, `fetch_sources.sh`, `sign.sh` |
| `results/` | data (`.npz`, `.json`) and plots |
| `tests/` | `uv run pytest` |
| `config.toml` | settings (no secrets are used) |

## For agents picking this up

- **Read first:** the brief, then the findings files in order. Findings 01 lists where the brief is wrong. Don't
  repeat those errors, for example the regime-map thresholds or the $N_{\max}$ formula.
- **Name the papers** S15 (1411.0739), S16 (1607.03937), A16 (1510.08918), A17 (1610.02626), SAL14 (1409.3247).
  Cite equations by the source's `\label` name and `.tex` line number.
- **Four adiabaticity parameters** exist and run in different directions: $\epsilon$ (S15), $\bar{\mathcal A}$
  (S16), $\mathcal A$ (A16), $\mathscr A$ (A17). Keep them as separate symbols in code.
- **Stellar constants differ by paper:** S15, S16 and A16 use $k_q=0.05$, $k_\star=0.1$; A17 uses 0.01 and 0.06.
- **Units:** `lk_orbit`, `spin` and `shape` use time in $t_k$, and the non-dissipative spin problem depends only on
  $\epsilon$. `secular` uses AU, $M_\odot$ and years, and reads its physical parameters from `config.toml`
  `[step4]`, where every assumed value is marked.
- **Tidal bulge:** A16's tidal-bulge precession is inconsistent by a factor of 2; we use $\tfrac{15}{2}k_{2p}$
  (findings 03).
- **Sign each result:**
  1. Open every new findings file with the brief's SHA-256
     (`11966760d9b6fc5e008207e46cd21e6b5dfb3c2160a9453c8ac546143dd1c3c3`).
  2. Run `git add`, then `scripts/sign.sh`, then `git add docs/SHA256SUMS` before committing.
  3. Put the brief hash in the commit message.
- **Don't commit** `sources/`, which holds the arXiv tarballs and renders of the publishers' figures. They are
  copyrighted, and `fetch_sources.sh` recreates them.
- **Report honestly:**
  - Never tune physics to match a figure.
  - State what didn't reproduce, and every assumption you made.
  - Validate each integrator against an independent method before running it in production.
- **Code style:** small functions, named constants, no try/except. Lint strictness is not a goal here.

## Licence

Everything written for this repo is dedicated to the public domain under [CC0 1.0](LICENSE). The brief in
`docs/brief/` belongs to its author and is included unchanged. The papers are cited, not redistributed.
