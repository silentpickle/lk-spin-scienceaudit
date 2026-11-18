# Data: everything the results rest on

The git repo holds the code, the findings, every summary result (`results/**/*.json`, `features.npz`, plots,
`comparison.csv`) and the step 11 N-body output. The large raw outputs are attached to the GitHub release
[`v0.1-first-results`](https://github.com/silentpickle/lk-spin-scienceaudit/releases/tag/v0.1-first-results)
as compressed archives: the per-cycle series, the per-system chunks and the N-body checkpoints. Nothing needs to be regenerated.

**Get everything:**

    scripts/fetch_data.sh          # downloads every archive, checks docs/SHA256SUMS-data.txt, unpacks in place
    scripts/fetch_sources.sh       # the papers' arXiv LaTeX (not redistributed here; see below)
    shasum -a 256 -c docs/SHA256SUMS

The four largest archives (`step6-b_E`, `step6-b_L`, `step6-b_P`, `step6-d`) are stored on the release as 200 MB parts (`<name>.part00`, `.part01`, ...), because GitHub uploads of about 700 MB kept failing. `fetch_data.sh` reassembles them, and the checksum covers the whole archive. To do it by hand, run `cat <name>.part* > <name>`.

Each archive unpacks at the repo root and puts its files back exactly where the scripts expect them. Checksums are
in [`SHA256SUMS-data.txt`](SHA256SUMS-data.txt), which is also attached to the release. The format is
`.tar.zst` (`zstd -d -c X.tar.zst | tar -xf -`).

## Archives

The base URL is `https://github.com/silentpickle/lk-spin-scienceaudit/releases/download/v0.1-first-results/`.

| Archive | Size | Files | Contents | Written by | Read by | Findings |
|---|---|---|---|---|---|---|
| `lkspin-data-step4.tar.zst` | 0.3 MB | 2012 | per-chunk scan results (θ_sl,f, outcome, t_end, steps) for the 2 × 1000 $P_\star$ runs; the merged `scan_*.npz` is in git | `scripts/step4_scan.py` | merged into `results/step4/scan_{89,87}.npz` | 03 |
| `lkspin-data-step5.tar.zst` | 134 MB | 3001 | `results/step5/series.npz` (per-LK-cycle $t,a,e_{\min}$, amplitude, $(p,\Psi)$ for member 0, the middle member and the twin, for every hot-Jupiter system) and `chunks/` (per-system raw results for all 3000 draws) | `scripts/step5_run.py` | `step5_classify.py` (v2 labels), `step5_plots.py` (montages), `prep_sq2.py` | 04 |
| `lkspin-data-step6-b.tar.zst` | 648 MB | 1030 | case (b), back-reaction: `series.npz` (member 0's own orbit), `series_members.npz` (every member's series), `chunks/` | `scripts/step6_run.py --case b` | `step6_analyze.py`, `step10_sq5.py`, `job4_analyze.py` | 05 |
| `lkspin-data-step6-c.tar.zst` | 577 MB | 1030 | case (c), back-reaction + octupole; same layout | `step6_run.py --case c` | `step6_analyze.py` | 05 |
| `lkspin-data-step6-d.tar.zst` | 728 MB | 1030 | case (d), back-reaction + octupole + spin-down; same layout | `step6_run.py --case d` | `step6_analyze.py` | 05 |
| `lkspin-data-step6-b_P.tar.zst` | 671 MB | 1030 | job 4, apsidal precession only (physically consistent) | `step6_run.py --case b_P` | `job4_analyze.py` | 08 |
| `lkspin-data-step6-b_F.tar.zst` | 612 MB | 1030 | job 4, orbit rotation about $\hat S$ only (physically consistent) | `step6_run.py --case b_F` | `job4_analyze.py` | 08 |
| `lkspin-data-step6-b_L.tar.zst` | 755 MB | 1030 | job 4, $d\mathbf L/dt$ only. **Unphysical split**: $\mathbf e$ leaves the orbital plane | `step6_run.py --case b_L` | `job4_analyze.py` | 08 |
| `lkspin-data-step6-b_E.tar.zst` | 582 MB | 1030 | job 4, $d\mathbf e/dt$ only. **Unphysical split** | `step6_run.py --case b_E` | `job4_analyze.py` | 08 |
| `lkspin-data-step6-pilots.tar.zst` | 29 MB | 66 | 20-system pilots of cases b, c and d | `step6_run.py --case X --pilot 20` | pilot analysis | 05 |
| `lkspin-data-step11.tar.zst` | 0.06 MB | 32 | full REBOUND/REBOUNDx checkpoints (`state/*.bin`, `*.rebx`) for every N-body run, so any run can be extended; the N-body output itself (`chunks/`) is in git | `scripts/step11_nbody.py` | `step11_nbody.py` (resume) | 06 |
| `lkspin-data-prep_sq2.tar.zst` | 0.3 MB | 216 | per-system ensemble chunks (A, B, C, W) | `scripts/prep_sq2.py` | `prep_sq2.py` summary | 07 |
| `lkspin-data-prep_sq4.tar.zst` | 10 MB | 130 | FTLE separation series (`series_ftle_*.npz`) and chunks for the scans, sections and $\hat\Omega_L$ test | `scripts/prep_sq4.py` | `prep_sq4.py` plots | 07 |

## What lives where

- **In git.** All code (`lkspin/`, `scripts/`, `tests/`), all findings, `config.toml` (every setting, each with its source or
  marked ASSUMED), and every summary product. The summaries are `results/step5/features.npz`, `labels_v2.npz`,
  `handlabels.csv` (our hand labels, with reasons), `results/step6/*/features.npz` and `summary.json`, the `results/step11/` chunks and
  `comparison.csv`, `results/prep_*/summary.json`, `results/sens4/`, and the plots. `docs/SHA256SUMS` hashes all of these.
- **In the release.** The raw series, chunks and checkpoints listed above, hashed in `docs/SHA256SUMS-data.txt`.
- **Not redistributed.** `sources/`, the arXiv LaTeX of the five papers and the renders of their figures, which remain
  the publishers' copyright. `scripts/fetch_sources.sh` downloads the exact pinned versions and checks them against
  `docs/SHA256SUMS`. Rendering the comparison figures needs Ghostscript (`docs/TOOLING.md`).

## Going a different direction

- **Reclassify** without rerunning anything. `lkspin.regime.v2_features` and `classify_v2` work on the step 5 and 6
  series; change the thresholds in `config.toml` `[step5]` (`v2_*`) and rerun `step5_classify.py` and `step6_analyze.py`.
- **Extend an N-body run** from its checkpoint. Raise its target in `[step11]` and rerun `step11_nbody.py --only NAME`.
  Resumes are bit-identical (`tests/test_step11.py`).
- **Change physics.** The integrators take every physical parameter from `config.toml`. `scripts/prep_sens4.py` shows how to vary one
  at a time. Back-reaction channels can be switched separately through the opt-in multipliers in `lkspin/secular.py`
  (findings 08).
- **Choose differently.** Every assumption is marked ASSUMED in `config.toml`. Every judgement call (hand labels, thresholds, the inconsistent
  feedback split) is recorded in the findings, with the reasoning, so it can be disputed.
