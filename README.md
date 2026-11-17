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
| 4 | Reproduce Storch et al. 2016's $\theta_{sl,f}$ vs $P_\star$ with short-range forces and tides | in progress | — |
| 5–11 | Regime map, sub-questions SQ1–SQ5, N-body spot check | not started | — |

## Quick start

    uv sync
    scripts/fetch_sources.sh        # pinned arXiv LaTeX sources, hash-checked (kept out of git)
    uv run pytest                   # fast checks
    shasum -a 256 -c docs/SHA256SUMS

Each findings file lists the scripts that produce its results; every script writes to `results/`.
Tools and versions: [`docs/TOOLING.md`](docs/TOOLING.md).

## Layout

| Path | What |
|---|---|
| `docs/brief/` | the research brief, unchanged |
| `docs/findings/NN-*.md` | one write-up per step: method, validation numbers, comparison with the papers, corrections, what was not done |
| `docs/SHA256SUMS` | hashes of every tracked file plus the pinned sources |
| `lkspin/` | the physics: `lk_orbit.py` (analytic LK orbit), `spin.py` (spin integrator), `shape.py` (Fourier and resonance analysis) |
| `scripts/` | one script per result; `calibrate.py`, `step3_*.py`, `fetch_sources.sh`, `sign.sh` |
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
- **Units:** `lkspin` uses time in $t_k$, and the non-dissipative spin problem depends only on $\epsilon$.
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
