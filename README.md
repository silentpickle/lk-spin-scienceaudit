# Lidov-Kozai stellar spin: working the ScienceAudit brief

Work on one open question from the ScienceAudit physics map (science.caseyhandmer.com/physics), subfield
"Planet formation & spin-orbit dynamics". The brief is in `docs/brief/`. Each findings file in `docs/findings/`
names the brief's SHA-256, so results stay tied to the exact version of the question they answer.

- `docs/findings/01-definitions-and-calibration.md`: work plan steps 1–2
- `scripts/calibrate.py`: step 2 checks (`uv run python scripts/calibrate.py`)
- `scripts/fetch_sources.sh`: re-downloads the pinned arXiv sources (kept out of git)

Verify everything: `shasum -a 256 -c docs/SHA256SUMS` (the source lines need `scripts/fetch_sources.sh` first on a fresh clone).
