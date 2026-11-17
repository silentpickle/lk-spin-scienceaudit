# Lidov-Kozai stellar spin: working the ScienceAudit brief

Work on one open question from the ScienceAudit physics map (science.caseyhandmer.com/physics), subfield
"Planet formation & spin-orbit dynamics". The brief is in `docs/brief/`. Each findings file in `docs/findings/`
names the brief's SHA-256, so results stay tied to the exact version of the question they answer.

## Findings

- `docs/findings/01-definitions-and-calibration.md`: work plan steps 1–2. Definitions recovered from the papers'
  LaTeX; published prefactors recomputed; eight corrections to the brief.
- `docs/findings/02-dynamics-calibration.md`: step 3. Storch & Lai 2015's bifurcation diagram, surfaces of section,
  Fourier coefficients and resonance widths reproduced; their analytic chaos-boundary curve does not reproduce.

## Layout

- `lkspin/`: the physics. `lk_orbit.py` (analytic Lidov-Kozai orbit), `spin.py` (spin integrator), `shape.py`
  (shape functions, Fourier coefficients, resonance widths, chaos boundary), `spin_check.py` (scipy reference).
- `scripts/`: one script per result. `calibrate.py` (step 2), `step3_*.py` (step 3), `fetch_sources.sh`
  (re-downloads the pinned arXiv sources, kept out of git), `sign.sh` (rebuilds the hash manifest).
- `results/`: data (`.npz`, `.json`) and plots written by the scripts.
- `tests/`: fast checks (`uv run pytest`).

## Reproduce

    uv sync
    scripts/fetch_sources.sh                     # papers, checked against docs/SHA256SUMS
    uv run pytest
    uv run python scripts/calibrate.py           # step 2
    uv run python scripts/step3_bifurcation.py   # about 6 min on 16 cores
    uv run python scripts/step3_sections.py
    uv run python scripts/step3_fourier.py
    uv run python scripts/step3_exact.py
    uv run python scripts/step3_validate.py

## Verify

`shasum -a 256 -c docs/SHA256SUMS` checks the brief, findings, code, results and the pinned sources together. Each
findings file and commit message names the brief's SHA-256. The published figures used for comparison are rendered
locally from the arXiv EPS files (needs Ghostscript) and are not committed, because of publisher copyright.
