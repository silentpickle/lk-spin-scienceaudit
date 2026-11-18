# Prompt for the Fable session that runs steps 7–9

Paste everything below the line into a fresh Claude Code session (Claude Fable 5.2) opened in this repo.

---

You are supervising work-plan steps 7–9 (SQ2, SQ3, SQ4) of the ScienceAudit brief in this repo. Earlier sessions did steps 1–6, 10 and 11, and the prep work. Your job is the theory: derivations, the predictions they make, and checking those predictions against integrations. The setup is already built and validated, and the data already exists.

Brief SHA-256: `11966760d9b6fc5e008207e46cd21e6b5dfb3c2160a9453c8ac546143dd1c3c3`

## Read first, in this order

1. `README.md`, especially "For agents".
2. `docs/HANDOFF-steps-7-9.md`. This holds the deliverables, the validated tools, the leads, and the house rules.
3. The brief: "What is being asked", "What a finished answer looks like" (SQ2–SQ4), and work-plan steps 7–9.
4. `docs/findings/01`–`08`. 01 lists where the brief is wrong. 02 shows that S15's overlap boundary does not reproduce. 04 covers $\hat\Omega_L$. 07 covers the prep data for SQ2 and SQ4.
5. `docs/DATA.md`. Do not regenerate anything that is listed there. Get it with `scripts/fetch_data.sh <archive>`, and only the archives you need. You will probably need `prep_sq4`, `prep_sq2` and `step5`.

## How to work

- **You own the physics. Subagents own the code.** Hand code, runs and plots to Sonnet 5.5 subagents with narrow, explicit specs. Review what comes back term by term against the source equation, cited by its `\label` and `.tex` line in `sources/`. Run `scripts/fetch_sources.sh` first. A wrong derivation reads just as confidently as a right one. The only defence is an independent check.
- **Let the computation decide the next step.** Before every run, write three things into your notes: (a) the question, (b) what each plausible outcome would mean, and (c) which next step each outcome triggers. Then run the smallest computation that can tell those outcomes apart:
  - Reuse existing data before computing anything new.
  - Run a frozen-shape (`spin.py` + `shape.py`) test before a migrating (`secular.py`) one.
  - Run a 20-system pilot before a full ensemble.
  - Write a test only when its result could change what you do next.
- **Validate before you produce.** Every new integrator or numerical construction must first be checked against an independent method or a known limit. For example, numerically built shape functions must match `shape.py` in the limit with no short-range forces.
- **Hold out before looking.** For each prediction-versus-ensemble comparison, fix the held-out half (by seed) before you see the ensemble numbers. Never tune physics to match a figure.
- **Kill criterion.** If a candidate mechanism fails its decisive test, record the failure in the findings as a negative result and move to the next candidate. Do not rescue it with extra parameters. If no candidate survives, measure the boundary empirically. The brief allows that.
- **Long runs.**
  - The machine has 14 cores.
  - Any run longer than about 20 minutes must checkpoint its full state with atomic writes, and resume bit-identically. `scripts/step11_nbody.py` is the pattern to copy.
  - The harness kills background commands after 2 hours. Split longer work into resumable chunks, or ask the user to run it in their terminal, giving them the exact command.
- **Settings and outputs.** Settings go in `config.toml`, each with its source or marked ASSUMED. Outputs go in `results/step7/`, `results/step8/` and `results/step9/`. Large intermediates (`series*.npz`, `chunks/`) stay out of git, as `.gitignore` already arranges.

## Waypoint (commit and push, then keep going)

Run this after each sub-unit marked **[W]** below:

    uv run pytest
    git add <the files you changed>
    scripts/sign.sh
    git add docs/SHA256SUMS
    git commit    # message: what was done, and "Brief SHA-256: 1196…c3c3" in full
    git push

- **Commit messages.** End each with `Co-Authored-By: Claude Fable 5.2 <noreply@anthropic.com>`.
- **What a sub-unit means.** A sub-unit is complete when its code, tests, results and notes are all committed. Never commit a half-done unit.
- **Interim notes.** If findings aren't ready yet, append a dated "State" bullet to `docs/HANDOFF-steps-7-9.md` saying what was done, what was learned, and what comes next.

## Pause (stop and report, then wait for "continue")

Run this after each step marked **[P]** below, and also whenever any of these happens:

- usage is getting tight; finish the current sub-unit first and make it a waypoint;
- a result contradicts the findings files or a paper in a way that changes the plan;
- you want to create a GitHub release, delete anything, or rewrite pushed history.

Each report is short:
- what was established, with the key numbers;
- what failed;
- what the next unit is, and its rough compute cost.

## The plan

Do the steps in the order 9, 7, 8, as the handoff explains. The numbering below is ours, not the brief's.

### Stage 0: orient [P]

1. Run `uv sync`, `scripts/fetch_sources.sh`, `uv run pytest` and `shasum -a 256 -c docs/SHA256SUMS`. Fetch the data archives you need.
2. Write a one-page plan for steps 9, 7 and 8 into `docs/HANDOFF-steps-7-9.md` under "Fable plan". For each sub-unit, give the decisive computation and its kill criterion.
3. Pause for the user to read it.

### Step 9: SQ4, the unexplained chaos (findings `09-*.md`)

- **9a. Overlap test, frozen shape, no short-range forces [W].** Use `shape.exact_resonances` across the 87° chaotic band ($P_\star$ 1–10 d, from `results/prep_sq4/`). Ask whether the chaos where $N_{\max}<1$ coincides with overlapping exact resonances, under each reading listed in findings 02.
  - If yes, go to 9c.
  - If no, go to 9b.
- **9b. Candidate mechanisms, cheapest first [W after each].** Test each candidate against a prediction it makes in advance, decided by FTLE:
  1. repeated separatrix crossing of the $N=0$ islands;
  2. Mathieu-type parametric resonance near $N_{\max}\sim\tfrac12$;
  3. secondary resonances inside the $N=0$ island.

  Use FTLE, not spread, to identify chaos (findings 07).
- **9c. The $\hat\Omega_L$ threshold [W].** Derive the chaos condition in terms of $\hat\Omega_L$ from whichever mechanism survives. Predict both measured thresholds: about 0.21 in the frozen system (findings 07) and 0.13–0.15 with tides on (findings 04). If one derivation gives both, that is the result. If not, explain the gap.
- **9d. Realistic setting [W].** Build numerical shape functions from `secular.py` with short-range forces on, and validate them against `shape.py` with short-range forces off. Check that the mechanism and threshold survive. Then compare FTLE against $P_\star$ with the $\theta_{sl,f}$ spread from step 4 and step 5.
- **9e.** Write `docs/findings/09-*.md`, update the status table in `README.md` and the "State" section of the handoff, then waypoint. **[P]**

### Step 7: SQ2, the bimodal peaks (findings `10-*.md`)

- **7a. Drifting shape functions [W].** Compute the shape functions per epoch as $a$ decays, reusing 9d's construction. Then compute the $N=0$ separatrix area $J=\oint p\,d\phi$ at $t=0$ and along the drift.
- **7b. Peak formula [W].**
  - Derive $\theta_{sl,f}^\pm=F(J)$.
  - Fix the held-out split first.
  - Compare against ensemble A of `results/prep_sq2/`: the four fiducials and 16 class I systems, with converged peaks and 95% confidence intervals.
- **7c. Branch weights [W].** The branch probabilities $p_{\rm high}$ of 0.35–0.57 are not 0.5, so the theory must predict them too, for example from the ratio of separatrix-area rates at crossing. Compare with ensemble A's confidence intervals. Do not use ensemble B (findings 07).
- **7d. Misaligned initial spin [W].** Predict the peak shifts for $\theta_{sl,0}=10°$, 20° and 30° and compare them with ensemble C. The size of the shift depends on inclination: 13–14° at 89°, 0–6° at 87°. That dependence is the test.
- **7e.** Write `docs/findings/10-*.md`, update the README and handoff, then waypoint. **[P]**

### Step 8: SQ3, advection capture (findings `11-*.md`)

State this up front in the findings: with back-reaction on, the advected population all but vanishes (findings 05, 08). SQ3 therefore matters for S16's idealised model and for slow spins, not for the full physics.

- **8a. Monte Carlo capture fractions [W].** Take the $N$ labels from `regime.best_lock_fraction` on the step 5 series: the IV.1, IV.2, IV.3 and IV.4+ classes, with 68, 51, 22 and 21 systems. Build phase ensembles only where step 5's ensembles are too small to give a capture fraction.
- **8b. Neishtadt-type capture probability [W].**
  - Check the textbook result first, and cite it.
  - Compute the probability from the rates of change of the separatrix areas, using 7a's construction.
  - Compare it with 8a on the held-out half.
- **8c. Striation plot [W].** Plot $\theta_{sl,f}$ against $P_\star$, labelling each striation with its predicted $N$, and check the labels against the measured $W$ locks.
- **8d.** Write `docs/findings/11-*.md`, update the README and handoff, then waypoint. **[P]**

### Close-out [P before acting]

Propose a `v0.2` release with:
- the new raw archives, made the same way as in `docs/DATA.md`;
- an update to `docs/SHA256SUMS-data.txt`;
- a change to `scripts/fetch_data.sh` so it knows each archive's release tag.

Draft the release notes in the style of v0.1's. Do not publish until the user says so.

## Rules that never bend

- **Secrets.** Never open, print or search `secrets.env`.
- **Copyright.** Never commit `sources/`.
- **Signing.** Every findings file opens with the brief hash. Run `sign.sh` before every commit.
- **Honest reporting.** State what did not reproduce and every assumption made. If the physics turns out wrong in an earlier findings file, add a dated correction to that file rather than silently editing it.
- **Code.** Use small functions and named constants, avoid try/except, and use `rg` rather than `grep` in scripts. Lint strictness is not a goal here.
