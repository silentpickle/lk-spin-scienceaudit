# Steps 1–2: recovered definitions and prefactor calibration

| | |
|---|---|
| Answers brief | `docs/brief/planet-formation-and-spin-orbit-dynamics--how-does-stellar-spin-evolve-during-lidov-kozai-migration.md` |
| Brief SHA-256 | `11966760d9b6fc5e008207e46cd21e6b5dfb3c2160a9453c8ac546143dd1c3c3` |
| Work plan steps | 1 (recover definitions) and 2 (calibrate on numbers) |
| Date | 2026-11-17 |
| Done by | Claude (Opus 5.5), working with the repo owner. Not checked by a physicist. |
| Verify | `shasum -a 256 -c docs/SHA256SUMS` checks this file, the brief, the sources and the results together |

## Sources

LaTeX sources downloaded from arXiv, pinned by version (`scripts/fetch_sources.sh` re-fetches and checks them).

| Short name | arXiv | Paper | e-print SHA-256 |
|---|---|---|---|
| SAL14 | 1409.3247v1 | Storch, Anderson & Lai 2014 (Science) | `f352bb86…92de5` |
| S15 | 1411.0739v1 | Storch & Lai 2015 | `9f2fa443…a16097` |
| A16 | 1510.08918v2 | Anderson, Storch & Lai 2016 | `61c2e4df…7f096f` |
| S16 | 1607.03937v1 | Storch, Lai & Anderson 2016 | `b1733e7f…21aa9` |
| A17 | 1610.02626v2 | Anderson, Lai & Storch 2017 | `3b91b5e7…e10c` |

Equation references below are the LaTeX `\label` names, plus line numbers in the `.tex` source. I did not compile
the PDFs, so printed equation numbers are not given. The SAL14 source is a 7-line stub with no equations.

## Definitions (work plan step 1)

| Item | Recovered form | Source |
|---|---|---|
| $t_k$ | $t_k^{-1}=n\,(1-e_b^2)^{-3/2}(M_b/M_\star)(a/a_b)^3$, $n=\sqrt{GM_\star/a^3}$. S16 drops the $e_b$ factor (it sets $e_b=0$). A16 and A17 write $a_{b,\rm eff}=a_b\sqrt{1-e_b^2}$, which is the same thing. The brief's $-3/2$ is right. | S15 l.288; S16 l.213; A16 Eq. tk; A17 Eq. tk |
| $n_e$ | $n_e=t_k^{-1}\,\frac{6\pi\sqrt6}{8K}\sqrt{x_2-x_1}$, with $x=x_0+(x_1-x_0)\,\mathrm{cn}^2(\theta,k^2)$ and $\theta=\frac{K}{\pi}(n_et+\pi)$ (no SRFs, $\omega=0$ at $t=0$). | S15 Eq. x, l.322–341 |
| $\Omega_{pl,0}$ | $\frac{3}{4t_k}\sqrt h\simeq\frac{3}{4t_k}\lvert\cos\theta_{lb}^0\rvert$, $h=x_0\cos^2\theta_{lb}^0$ | S15 l.412 |
| $\alpha_0$ (the brief's $\Omega_{ps,0}$) | $\frac{3GM_p(I_3-I_1)}{2a^3I_3\Omega_\star}=1.19\times10^{-8}\frac{2\pi}{\rm yr}\left(\frac{2k_q}{k_\star}\right)\left(\frac{10^3M_p}{M_\star}\right)\left(\frac{\hat\Omega_\star}{0.05}\right)\left(\frac{a}{\rm AU}\right)^{-3}\left(\frac{M_\star}{M_\odot}\right)^{1/2}\left(\frac{R_\star}{R_\odot}\right)^{3/2}$ | S15 Eq. alpha0 l.379; S16 Eq. alpha0 l.268 |
| $k_q,\ k_\star$ | S15, S16 and A16: $k_q=0.05$, $k_\star=0.1$ (Claret 1992; Claret & Gimenez 1992). **A17 uses $k_{q,0}=0.01$, $k_\star=0.06$.** | S15 l.390; S16 l.280; A16 l.416; A17 l.177, l.424 |
| $\epsilon$ | $\Omega_{pl,0}/\alpha_0=1.17\left(\frac{k_\star}{2k_q}\right)\left(\frac{R_\star}{R_\odot}\right)^{-3/2}\left(\frac{\hat\Omega_\star}{0.1}\right)^{-1}\left(\frac{M_b}{10^3M_p}\right)\left(\frac{a}{\rm AU}\right)^{9/2}\left(\frac{a_b}{300\,\rm AU}\right)^{-3}\lvert\cos\theta_{lb}^0\rvert$ | S15 Eq. epsilon l.165 |
| Shape functions | $\epsilon\beta=-\frac{\Omega_{pl}}{\alpha}\sin\theta_{lb}$, $\epsilon\gamma=\frac{\dot\theta_{lb}}{\alpha}$, $\epsilon\psi=-\frac{\Omega_{pl}}{\alpha}\cos\theta_{lb}$; $H=\frac{\bar\alpha}{n_e}\{-\tfrac12p^2+\epsilon\psi p-\epsilon\sqrt{1-p^2}[\beta\cos\phi+\gamma\sin\phi]\}$. $\beta,\psi$ are cosine series, $\gamma$ a sine series; S16 sets $\gamma_0\equiv\beta_0$ "for convenience". | S15 Eqs. beta–psi l.549; S16 Eq. kozaiham l.319, l.369 |
| $\bar\alpha/n_e$ | $\frac{1}{2\pi}\int_0^{2\pi/n_e}\alpha\,dt$ | S16 Eq. alphabar l.334 |
| $N_{\max}$ | S15: $\lfloor\bar\alpha/n_e\rfloor=\lfloor\mathcal N(\cos\theta_{lb}^0;e_0)/\epsilon\rfloor$. S16: $\bar\alpha/n_e$, not rounded down. | S15 Eq. nmax l.621; S16 Eq. nmax l.413 |
| $\mathcal N$ | $\mathcal N=0.98$ at $\theta_{lb}^0=85^\circ$, $e_0=0.01$. No closed form is printed. From the definitions, $\mathcal N=\epsilon\bar\alpha/n_e=\frac{3\sqrt h}{8\pi}\int_{\rm 1\,LK\,period}x^{-3/2}\,d(t/t_k)$ (my derivation; check C8 reproduces 0.98). | S15 l.631–634 |
| $w_N$ | $w_N\simeq2\left[2\epsilon\lvert\beta_N+\gamma_N\rvert\sqrt{1-p_N^2}\right]^{1/2}$, $p_N\simeq-N/N_{\max}$ | S15 Eq. appwidth l.773 |
| $\bar{\mathcal A}$ (S16) | $\bar{\mathcal A}=(2\epsilon\beta_0)^{-1}=\frac{\bar\alpha}{2\langle\Omega_{pl}\sin\theta_{lb}\rangle}=\frac{N_{\max}}{2\hat\Omega_L}$, $\hat\Omega_L=\langle\Omega_{pl}\sin\theta_{lb}\rangle/n_e$ | S16 Eqs. abar, abar2 l.442–454 |
| $\mathcal A$ (A16) | $\lvert\Omega_{ps}/\Omega_L\rvert$, $\Omega_L=[(\Omega_{pl}\sin\theta_{lb})^2+\dot\theta_{lb}^2]^{1/2}$; $\mathcal A_0=0.07\,\frac{\bar k_q\bar M_p\bar M_{\rm tot}^{1/2}\bar R_\star^3\bar a_b^3}{\bar k_\star\bar M_\star\bar M_b\bar a^{9/2}\bar P_\star}\left\lvert\frac{\cos\theta_{sl,0}}{\sin2\theta_{lb,0}}\right\rvert$ | A16 Eqs. OmegaPL, Aparameter, adiabatic l.398–450 |
| $\mathscr A$ (A17) | $\left\lvert\frac{\Omega_{ps}t_kj^3}{\cos\theta_{sl}}\right\rvert\simeq58\frac{\bar k_{q,0}\bar m_1\bar m_{01}^{1/2}\bar R_0^3}{\bar k_\star\bar m_0\bar m_2}\left(\frac{P_\star}{5\,\rm d}\right)^{-1}\left(\frac{a_{\rm in}}{0.1\,\rm AU}\right)^{-9/2}\left(\frac{a_{\rm out,eff}}{10\,\rm AU}\right)^3$ | A17 Eq. Aparameter l.420 |
| $\mathscr A_{\rm crit}$ | $\simeq3$ | A17 l.428 |
| $S/L$ | $0.12\left(\frac{k_\star}{0.1}\right)\left(\frac{M_\star}{M_\odot}\right)^{1/2}\left(\frac{R_\star}{R_\odot}\right)^2\left(\frac{M_p}{M_J}\right)^{-1}\left[\frac{a(1-e^2)}{0.05\,\rm AU}\right]^{-1/2}\left(\frac{P_\star}{30\,\rm d}\right)^{-1}$. **The missing factor is $[a(1-e^2)/0.05\,\mathrm{AU}]^{-1/2}$.** | S16 Eq. SoverL l.601 |
| LK window | $(\cos I_0)_\pm=\frac{1}{10}\left(-\eta\pm\sqrt{\eta^2+60-20\varepsilon_{\rm SRF}}\right)$, real only if $\eta^2+60-20\varepsilon_{\rm SRF}\ge0$; $\varepsilon_{\rm SRF}=\tfrac43\varepsilon_{\rm GR}+2\varepsilon_{\rm Rot}+\tfrac43\varepsilon_{\rm Tide}$. The brief's reconstruction is right. | A17 Eqs. Icrit, kozai_condition l.229–244 |
| $I_{0,\lim}$ | $\cos I_{0,\lim}=\frac\eta2\left(\tfrac45j_{\lim}^2-1\right)$ (brief right) | A17 Eq. Ilim l.282 |
| Max perturber distance | $a_{\rm out,eff}<19.6\,{\rm AU}\,(\bar m_2/\bar m_{01}^2)^{1/3}(a_{\rm in}/0.1\,{\rm AU})^{4/3}(1-\tfrac53\cos^2I_0)^{1/3}$, from $\varepsilon_{\rm GR}>\tfrac94(1-\tfrac53\cos^2I_0)$ | A17 Eq. maxaout_zero_ecc l.326–332 |
| SRF strengths | $\varepsilon_{\rm GR}\simeq3\times10^{-2}$, $\varepsilon_{\rm Tide}\simeq9.1\times10^{-7}$, $\varepsilon_{\rm Rot}\simeq2.9\times10^{-5}(P_\star/10\,\rm d)^{-2}$, with $\bar a_{\rm out}=a_{\rm out,eff}/100$ AU and $\bar a_{\rm in}=a_{\rm in}/1$ AU | A17 Eq. epsilon_SRFs l.171 |
| Spin-down | Skumanich: $\dot{\bf S}_\star\rvert_{\rm MB}=-\alpha_{\rm MB}I_\star\Omega_\star^3\hat S$, $\alpha_{\rm MB}=1.5\times10^{-14}$ yr (G stars), $1.5\times10^{-15}$ yr (F stars, $1.4M_\odot$) | A16 l.420–423 |
| Chaos onset (S15 example) | widespread chaos for $\epsilon\gtrsim0.14$ | S15 l.184 |

### Fiducial parameters

- **S15** (non-dissipative, no SRFs): $e_0=0.01$, $\theta_{lb}^0=85^\circ$, $\omega$ circulating ($\omega=0$ at $t=0$). Fig. 1 varies $\epsilon$ through $M_p$. The $M_p$ range is in the figure, not the text. (S15 l.348, l.180)
- **S16** (SRFs and weak-friction tides; no feedback; $e_b=0$; $\theta_{sl,0}=0$): $M_p=5M_J$, $a_0=1.5$ AU,
  $a_b=300$ AU, $\theta_{lb}^0=89^\circ$ or $87^\circ$, and $P_\star$ from 50 d down to about 0.3 d (l.631–636, l.961–971).
  For $89^\circ$ with SRFs: $\bar{\mathcal A}\simeq25.5/P_\star$ and $N_{\max}\simeq3.34/P_\star$, with $P_\star$ in days (l.673).
  Outcome bands from Fig. variableP (l.976–994):

  | Class | $\theta_{lb}^0=89^\circ$ | $\theta_{lb}^0=87^\circ$ |
  |---|---|---|
  | I bimodal | $P_\star\gtrsim23$ d | $P_\star\gtrsim9$ d |
  | II stationary | $7\lesssim P_\star\lesssim20$ d | — |
  | III chaos | — | $0.5\lesssim P_\star\lesssim9$ d |
  | IV advection | $0.5\lesssim P_\star\lesssim7$ d | — |
  | V fully adiabatic | $P_\star\lesssim0.4$ d | $P_\star\lesssim0.3$ d |

  Regime map (l.1373–1377): $\bar{\mathcal A}_0\lesssim1$ gives I; $\bar{\mathcal A}_0\gtrsim1$ with $N_{\max,0}\lesssim0.5$ gives II;
  $\bar{\mathcal A}_0\gtrsim1$ with $N_{\max,0}\gtrsim0.5$ gives IV, except for chaos (III) at the lower values of $\bar{\mathcal A}_0$.
  **The parameter ranges of the Fig. 14 suite are not in the text.**
- **A17**: $m_0=m_1=1M_\odot$, $m_2=1$ or $0.1M_\odot$, $a_{\rm in}=0.17$ AU in the figures, $k_{2,0}=0.03$.

## Calibration (work plan step 2)

`uv run python scripts/calibrate.py` recomputes each prefactor from the definitions above with IAU 2015 nominal
constants and the Julian year. The pass test is the brief's: both numbers agree when rounded to 2 significant figures.
C8–C10 are extra checks I added that the brief did not ask for. Full output: `results/calibration.json`.

| Check | Published | Computed | Rel. diff | 2 s.f. |
|---|---|---|---|---|
| C1_spin_precession ($\alpha_0$, 2π/yr) | 1.19e-8 | 1.18924e-8 | −0.06% | ok |
| C2_epsilon_prefactor | 1.17 | 1.16786 | −0.18% | ok |
| C3_triple_adiabaticity ($\mathscr A$) | 58 | 58.0842 | +0.15% | ok |
| C4_gr_strength ($\varepsilon_{\rm GR}$) | 3e-2 | 2.96119e-2 | −1.29% | ok |
| C5_perturber_distance (AU, from constants) | 19.6 | 19.6595 | +0.30% | ok |
| C6_perturber_rounded (AU, using 3e-2) | 19.6 | 19.5743 | −0.13% | ok |
| C7_kozai_window | $\sqrt{3/5}$ | 0.774597 | 0.00% | ok |
| C8_shape_factor ($\mathcal N$, 85°, $e_0=0.01$) | 0.98 | 0.97922 | −0.08% | ok |
| C9_anderson_precession (A16 $\Omega_{ps}$, 1/yr) | 1.64e-7 | 1.65249e-7 | +0.76% | **differs** |
| C10_anderson_adiabaticity (A16 $\mathcal A_0$) | 0.07 | 0.0701354 | +0.19% | ok |

All five prefactors the brief asked for reproduce (C1–C3, C5–C7). C8 also checks the analytic LK solution that step 3 will
use as forcing. C9 is 0.76% off. That is smaller than C4's gap, but it falls across the 1.65 rounding boundary. My guess is
that A16 used slightly different constants; I haven't checked. It does not affect anything downstream, since C10 (built on
C9) passes.

## Where the brief is wrong

Each item is checked against the LaTeX source listed.

1. **Chaos threshold.** The brief says widespread chaos in S15's example sets in at $\epsilon$ "of order unity". The source says $\epsilon\gtrsim0.14$ (S15 l.184).
2. **$N_{\max}$ formula.** The brief writes $N_{\max}=1/[\epsilon\,\mathcal N]$. The source has $\lfloor\mathcal N/\epsilon\rfloor$ (S15 l.621). S16 drops the rounding.
3. **Shape-function names.** The brief calls them $\mathcal A,\mathcal B,\mathcal C$. The papers use $\psi,\beta,\gamma$, with $\epsilon$ factored out of the perturbation terms. The placeholder the brief mentions is $\gamma_0\equiv\beta_0$ (S16 l.369), not "$\mathcal A_0$".
4. **There are four adiabaticity parameters, not three.** The definition the brief gives for A17, $\lvert\Omega_{ps}j^3/\Omega_L\rvert$ at $e=0$, is really A16's ($\lvert\Omega_{ps}/\Omega_L\rvert$). A17's is $\lvert\Omega_{ps}t_kj^3/\cos\theta_{sl}\rvert$. Code needs four symbols: $\epsilon$, $\bar{\mathcal A}_{\rm S16}$, $\mathcal A_{\rm A16}$ and $\mathscr A_{\rm A17}$.
5. **Regime-map thresholds.** The brief puts the stationary/advection boundary at $N_{\max,0}\approx1$ and chaos at lower $N_{\max,0}$. The paper's text gives $N_{\max,0}\approx0.5$, with chaos at lower $\bar{\mathcal A}_0$ (S16 l.1373–1377).
6. **Unexplained chaos.** S16 says its theory fails there because $\bar{\mathcal A}_0\lesssim1$ for $P_\star>5$ d (l.1194–1200). The brief frames the region as $N_{\max}\lesssim1$. Both may hold in that region, but the paper's stated reason is $\bar{\mathcal A}_0$.
7. **Stellar constants differ between papers.** A17 uses $k_q=0.01$ and $k_\star=0.06$, not 0.05 and 0.1. Any cross-paper comparison has to carry this.
8. **"Garbled" items were extraction artefacts.** The $(1-e_b^2)^{-3/2}$ factor is clean in the S15 source. The garbling was only in the text the brief's author was given.

## Assumptions

- $\varepsilon_{\rm GR}=3Gm_{01}^2a_{\rm out,eff}^3/(c^2a_{\rm in}^4m_2)$ is taken from Liu, Muñoz & Lai 2015 from memory, not read. A17 prints only the prefactor. C4 reproduces that prefactor to 1.3%, which supports it.
- For prefactors, every normalised quantity is set to 1, with all masses $1M_\odot$ where the paper's barred units allow it.
- $\mathcal N$ uses $\Omega_{pl,0}=3\sqrt h/(4t_k)$ exactly. The $\lvert\cos\theta_{lb}^0\rvert$ approximation differs by $5\times10^{-5}$ here.

## Not done

- A16's population-synthesis parameters (Section 6, l.1071 onward) and its appendix equations (l.1758 onward). SQ1 needs both.
- The parameter ranges of the S16 Fig. 14 suite. They are not in the text, so they have to be read off the figure or obtained from the authors.
- The SAL14 weak-friction tide equations. The SAL14 source is a stub, so use the A16 appendix (l.1861) instead.
- The values attributed to Claret 1992 are confirmed as quoted in the papers, not checked against Claret 1992 itself.
- Work plan steps 3 onward.
