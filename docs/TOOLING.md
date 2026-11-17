# Tooling index

What this project uses, what later work-plan steps will need, and what is worth evaluating. Versions are the ones
installed when this was last updated (2026-11-17, macOS arm64). Anything marked *to evaluate* has not been installed
or verified here yet; check its current release before relying on it.

## In use

| Tool | Version | Used for | Install |
|---|---|---|---|
| Python | 3.12.13 | everything | managed by uv |
| uv | 0.11.29 | environments, dependencies, `uv run` | `brew install uv` |
| numpy | 2.5.5 | arrays | `uv sync` |
| scipy | 1.18.1 | elliptic functions, quadrature, reference ODE solver (`solve_ivp` DOP853) | `uv sync` |
| numba | 0.68.1 | compiled integrators (`lkspin/spin.py`, `lkspin/lk_orbit.py`) | `uv sync` |
| matplotlib | 3.11.2 | all plots | `uv sync` |
| pytest | dev group | `tests/` | `uv sync` |
| Ghostscript | 10.08.0 | rendering the papers' EPS figures for comparison | `brew install ghostscript` |
| ImageMagick | 7.1.2 | `magick -density 110 fig.eps fig.png`; cropping figures to read values off by pixel | `brew install imagemagick` |
| ripgrep | 15.2.0 | searching the LaTeX sources | `brew install ripgrep` |
| shasum | macOS built-in | `docs/SHA256SUMS` signing | — |
| gh | 2.96.0 | GitHub; also `gh api /licenses/...` | `brew install gh` |
| arXiv e-print endpoint | — | `https://arxiv.org/e-print/<id>v<n>` returns the LaTeX source tarball, which is far better than PDF text for recovering equations | `scripts/fetch_sources.sh` |
| arXiv API | — | `https://export.arxiv.org/api/query?id_list=...` for versions and authors. It rate-limits and sometimes returns empty bodies; wait and retry, or read `\author` from the source instead. | — |

## Needed by later steps

| Tool | Step | Why | Status |
|---|---|---|---|
| REBOUND + REBOUNDx | 11 (N-body spot check) | direct N-body check that the secular (double-averaged) equations hold at the $e_{\max}$ reached. REBOUNDx has spin and tide effects; check which ones the current release implements. | to evaluate (`uv add rebound reboundx`) |
| MESA, or the Dartmouth Stellar Evolution Database (arXiv:0804.4473) | 6 (SQ1 part d, spin-down) | $R_\star(t)$. The Dartmouth tables lack $k_q$ and $k_\star$. | to evaluate |
| A TeX engine (tectonic is light) | any | compiling the papers to get **printed** equation numbers; findings currently cite `\label` names and `.tex` line numbers | to evaluate (`brew install tectonic`) |
| multiprocessing (stdlib) | 4–10 | parallel runs over parameter grids; 16 cores here | in use |

## Worth evaluating

| Tool | Why |
|---|---|
| JAX (`jax`, `diffrax`) | the brief suggests Numba **or** JAX. JAX would give batched integrators across a grid on one device, plus gradients (for example $d\theta_{sl,f}/dP_\star$ as a chaos indicator). Not needed while numba keeps up. |
| NASA Exoplanet Archive (TAP API), Southworth's TEPCat | optional comparison with observed obliquities. These are projected angles $\lambda$, not $\psi$, so simulated $\theta_{sl}$ has to be sky-projected first. |
| WebPlotDigitizer | reading values off published plots more systematically than by hand with ImageMagick crops (used for S15 Fig. 5) |
| An image-difference metric (`scikit-image`) | quantitative overlay of reproduced vs published figures instead of comparison by eye |
