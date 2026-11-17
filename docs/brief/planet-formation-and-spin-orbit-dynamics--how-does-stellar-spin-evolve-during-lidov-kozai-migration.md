<!--
Research brief from the ScienceAudit physics map. Subfield: Planet formation & spin-orbit dynamics.
Written by a language model (Claude Opus) on 2026-11-11 from the papers listed inside. No person has checked it.
The question was judged "answerable now" by a model; nobody has attempted the work. Treat every equation, number and
resource named here as a lead to verify, starting with the calibration step in the work plan.
-->

# Lidov-Kozai migration: how does the host star's spin evolve, and what sets the final spin-orbit misalignment?

## What is being asked

A giant planet starts at a few AU around a star of mass $M_\star$. A distant stellar companion of mass $M_b$ drives Lidov-Kozai (LK) oscillations in the planet's eccentricity $e$ and in its inclination $\theta_{lb}$ to the binary axis. Tidal dissipation at periastron then shrinks the orbit until the planet becomes a hot Jupiter (HJ). The star is rotating and therefore oblate, so the planet torques it and the spin axis $\hat S$ precesses about the planet's orbital axis $\hat L$. The quantity of interest is the final spin-orbit angle $\theta_{sl,f}$, with $\cos\theta_{sl}=\hat S\cdot\hat L$, and its dependence on planet mass, stellar spin, orbital geometry and initial conditions.

Storch et al. 2016 (arXiv:1607.03937) claim that two dimensionless numbers evaluated at $t=0$, the adiabaticity $\mathcal{A}_0$ and the maximum resonance order $N_{\max,0}$, fully determine which of five outcome classes a system falls into. That claim rests on an analysis that leaves out the star's back-reaction torque on the orbit, octupole terms, stellar spin-down and initial misalignment. The question is therefore made precise as the following sub-questions.

1. **Classification under full physics.** Integrate systems with spin-orbit feedback, octupole terms, short-range forces, tides and magnetic spin-down switched on. What fraction land in the outcome class (I to V) that $(\mathcal{A}_0, N_{\max,0})$ predicts, and where in that plane does the prediction fail?
2. **Bimodal regime (I).** Can the two peak values of $\theta_{sl,f}$ be predicted quantitatively from the adiabatic invariant $J=\oint p\,d\phi$ of the $N=0$ separatrix at $t=0$, and does the predicted peak location match integrations to within a stated tolerance (for example $5^\circ$)?
3. **Adiabatic advection (IV).** When the evolving shape functions are included, what is the probability of capture into resonance $N$, and what final angle results along each "diagonal striation"? Is advection a significant contributor of misalignment, or is it confined to unphysically short spin periods?
4. **Unexplained chaos.** Storch et al. 2016 report widespread chaos in a spin-period range where $N_{\max}\lesssim 1$, which their perturbative overlap theory cannot explain. Which resonances, or which non-perturbative mechanism, produce it, and is there a criterion that predicts its extent?
5. **Feedback threshold.** At what value of the spin-to-orbit angular momentum ratio $S/L$ does the star's back-reaction on $\hat L$ change the outcome class or shift $\theta_{sl,f}$ by more than the class-internal spread?

The following are out of scope: the HJ formation fraction, which belongs to a separate open question; coagulation; and disc processes.

## What a finished answer looks like

- **SQ1:** a confusion matrix of predicted class against observed class over at least several thousand integrated systems that sample $M_p$, $P_\star$, $a_0$, $a_b$ and $\theta_{lb}^0$. Report it separately for (a) the idealized Storch et al. 2016 model, (b) the same model plus feedback, (c) plus octupole with $e_b\neq 0$, and (d) plus spin-down. Plot the misclassified systems in the $(\mathcal{A}_0, N_{\max,0})$ plane. The classification counts as settled if (a) reproduces the published regime map and the accuracy loss in (b) to (d) is quantified with bootstrap errors.
- **SQ2:** an analytic or semi-analytic formula, $\theta_{sl,f}^{\pm}=F(\text{separatrix area at }t=0)$, overlaid on histograms of $\theta_{sl,f}$ from integrations with randomized phase. It must also predict the peak shift for initially misaligned systems; Storch et al. 2016 state qualitatively that the peaks move closer to $90^\circ$.
- **SQ3:** capture probability as a function of $N$ and of the drift rate of the Hamiltonian coefficients, from a separatrix-crossing (Neishtadt-type) calculation (standard result, check against a dynamics textbook), compared with Monte Carlo fractions. Include a plot of $\theta_{sl,f}$ against $P_\star$ that shows the striations with the predicted $N$ labels.
- **SQ4:** surfaces of section at fixed $\mathcal{A}$ and $N_{\max}<1$ in the period range where chaos appears, together with an identified set of overlapping resonances (or a demonstration that no overlap exists and some other mechanism, such as separatrix crossing of the $N=0$ islands, is responsible). Finite-time Lyapunov exponents as a function of $P_\star$ should match the region of spread in $\theta_{sl,f}$.
- **SQ5:** curves of the class-change fraction and of the median $|\Delta\theta_{sl,f}|$ against $S/L$ evaluated at periastron.
- **Negative result.** The question counts as "not answerable this way" if the idealized model (a) cannot reproduce the published Storch et al. 2016 regime map with the published parameters, if the outcome under full physics is dominated by phase sensitivity so strongly that the class is undefined for most of parameter space, or if the classification accuracy under (b) to (d) is no better than a classifier that uses $\theta_{lb}^0$ alone.

## Background the assistant needs

**LK oscillations (quadrupole, test-particle planet).** The planet has semimajor axis $a$, eccentricity $e$, argument of pericentre $\omega$, node $\Omega$, and inclination $\theta_{lb}$ between $\hat L$ and the binary axis $\hat L_b$. The binary has semimajor axis $a_b$ and eccentricity $e_b$. The characteristic rate is
$$t_k^{-1}=n\,\frac{M_b}{M_\star}\left(\frac{a}{a_b}\right)^3(1-e_b^2)^{-3/2},$$
where $n$ is the planet's mean motion. Storch et al. 2016 set $e_b=0$. The sign of the $(1-e_b^2)$ exponent is garbled in the supplied text of Storch & Lai 2015, and $-3/2$ is the standard result (check against a textbook).

With $e_0\approx 0$ and no short-range forces, the maximum eccentricity is
$$e_{\max}=\left(1-\tfrac53\cos^2\theta_{lb}^0\right)^{1/2}.$$
The conserved quantity is $\sqrt{1-e^2}\cos\theta_{lb}$ (standard result). Storch & Lai 2015 write $x\equiv 1-e^2$ and give an explicit Jacobi-elliptic solution $x=x_0+(x_1-x_0)\,\mathrm{cn}^2(\psi,k^2)$ with $\psi=K(n_e t+\dots)$. Here $n_e=2\pi/P_e$ is the eccentricity-oscillation frequency, $K$ is the complete elliptic integral of the first kind, and $k^2=(x_0-x_1)/(x_2-x_1)$, where $x_{1,2}$ are the roots of $x^2-\tfrac13(5+5h-2x_0)x+\tfrac53h=0$. Storch et al. 2016 write $n_e=\mathcal{K}\,t_k^{-1}$, where $\mathcal{K}$ is of order unity and depends on $e_0$ and on the short-range forces. Two derived rates matter: the nodal precession rate $\Omega_{pl}=d\Omega/dt$, whose value at $e=0$ is $\Omega_{pl,0}\approx \tfrac{3}{4t_k}|\cos\theta_{lb}^0|$, and the nutation rate $d\theta_{lb}/dt$. Storch & Lai 2015 consider only circulating $\omega$.

**Stellar spin precession.**
$$\frac{d\hat S}{dt}=\Omega_{ps}\,\hat L\times\hat S,\qquad \Omega_{ps}=-\frac{3GM_p(I_3-I_1)}{2a^3(1-e^2)^{3/2}S}\cos\theta_{sl}\equiv-\alpha\cos\theta_{sl}.$$
Here $I_3,I_1$ are the star's principal moments, $S$ is the spin angular momentum, and $\alpha=\Omega_{ps,0}/(1-e^2)^{3/2}$, with
$$\Omega_{ps,0}=1.19\times10^{-8}\,\frac{2\pi}{\rm yr}\left(\frac{2k_q}{k_\star}\right)\left(\frac{10^3M_p}{M_\star}\right)\left(\frac{\hat\Omega_\star}{0.05}\right)\left(\frac{a}{1\,\rm AU}\right)^{-3}\left(\frac{M_\star}{M_\odot}\right)^{1/2}\left(\frac{R_\star}{R_\odot}\right)^{3/2}.$$
Both papers use $I_3-I_1=k_qM_\star R_\star^2\hat\Omega_\star^2$ and $S=k_\star M_\star R_\star^2\Omega_\star$, with $\hat\Omega_\star=\Omega_\star/(GM_\star/R_\star^3)^{1/2}$ (the last relation is a standard result, check it). The normalization implies solar values $k_q\approx0.05$ and $k_\star\approx0.1$. The supplied text drops the numbers it attributes to Claret, so confirm them from the papers.

Storch et al. 2016 give the spin-to-orbit angular momentum ratio as $S/L\approx0.12\,(k_\star/0.1)(M_\star/M_\odot)^{1/2}(R_\star/R_\odot)^2(M_p/M_J)^{-1}(P_\star/30\,{\rm d})^{-1}$ times a factor in $a$ (and possibly $e$) that is garbled in the supplied text. Recover that factor.

**Adiabaticity parameters.** The conventions differ between papers.
- Storch et al. 2014 and Storch & Lai 2015 define $\epsilon=|\Omega_{pl}/\Omega_{ps}|$ evaluated at $e=0$ and $\theta_{sl}=0$, so $\epsilon=\Omega_{pl,0}/\Omega_{ps,0}$, which is constant through an LK cycle. Numerically, $\epsilon=1.17\,(\dots)(R_\star/R_\odot)^{-3/2}(\hat\Omega_\star/0.1)^{-1}(M_b/10^3M_p)(a/1\,{\rm AU})^{9/2}(a_b/300\,{\rm AU})^{-3}|\cos\theta_{lb}^0|$. Large $\epsilon$ is non-adiabatic.
- Storch et al. 2016 define $\bar\alpha=(n_e/2\pi)\int_0^{2\pi/n_e}\alpha\,dt$ and $N_{\max}=\bar\alpha/n_e$, which is the maximum number of spin precessions per LK cycle. They also define a "more precise" adiabaticity $\mathcal{A}$, which measures the LK-averaged maximum spin precession rate against the average rate of change of $\hat L$. Large $\mathcal{A}$ is adiabatic. The exact expression is garbled in the supplied text (it involves $\epsilon$, $N_{\max}$ and $\langle\Omega_L\rangle$ with $\Omega_L=|d\hat L/dt|$). Recover it from the paper before doing anything else.
- Anderson et al. 2017 define a third parameter, also written $\mathcal{A}$, as $|\Omega_{ps}j^3/\Omega_L|$ at $e=0$, with $j=\sqrt{1-e^2}$. They state explicitly that it differs from the earlier definitions.

Keep three distinct symbols in code ($\epsilon$, $\mathcal{A}_{\rm S16}$, $\mathcal{A}_{\rm A17}$).

**Hamiltonian formulation (Storch & Lai 2015).** In the frame rotating with $\hat L$, take $p=\cos\theta_{sl}$ conjugate to the spin precession phase $\phi$, and rescaled time $\tau=(n_e/\bar\alpha)\int_0^t\alpha\,dt'$. Then
$$H'(p,\phi,\tau)=\frac{\bar\alpha}{n_e}\left\{-\tfrac12p^2+\mathcal{A}(\tau)p-\sqrt{1-p^2}\,[\mathcal{B}(\tau)\cos\phi+\mathcal{C}(\tau)\sin\phi]\right\}.$$
The shape functions $\mathcal{A},\mathcal{B},\mathcal{C}$ are built from $\Omega_{pl}/\alpha$, $\theta_{lb}$ and $\dot\theta_{lb}/\alpha$ (the exact forms are in Storch & Lai 2015, Eqs. near the "Renormalized Hamiltonian" section). They depend only on the orbit shape, set by $\theta_{lb}^0$ and $e_0$. The function $\mathcal{A}(\tau)$ here is a shape function and must not be confused with the adiabaticity $\mathcal{A}$. The Fourier expansions are $\mathcal{A}=\sum\mathcal{A}_M\cos M\tau$, $\mathcal{B}=\sum\mathcal{B}_M\cos M\tau$ and $\mathcal{C}=\sum\mathcal{C}_M\sin M\tau$.

Resonance occurs when the time-averaged precession frequency equals an integer multiple of the eccentricity frequency, $\bar\Omega_{ps}=-\bar\alpha\,p=Nn_e$. The resonance centres are therefore at $p_N=-N/N_{\max}$, neighbouring resonances are spaced by $1/N_{\max}$, and only $|N|\le N_{\max}$ exist. Storch & Lai 2015 also write $N_{\max}=1/[\epsilon\,\mathcal{N}(\theta_{lb}^0;e_0)]$. A single resonance reduces to a pendulum with width $w_N\approx2[\dots]^{1/2}$, which depends on a sum of Fourier coefficients $\mathcal{B}_{|N|}$ and $\mathcal{C}_{|N|}$. The width is non-monotonic in $N$ (it rings down) and differs between $\pm N$.

Chirikov overlap: chaos sets in when $\Delta p\lesssim\tfrac12(w_N+w_M)$ (Chirikov criterion, standard result). With slow tidal evolution, $J=\oint p\,d\phi$ is an adiabatic invariant.

**Short-range forces (SRFs).** Pericentre precession from GR and from the planet's tidal and rotational bulges lowers $e_{\max}$ and smooths the shape functions. It lowers both $N_{\max}$ and $\mathcal{A}$ and reduces the number of resonances (Storch et al. 2016). In Anderson et al. 2017 notation, $\varepsilon_{\rm GR}$, $\varepsilon_{\rm Tide}$ and $\varepsilon_{\rm Rot}$ have prefactors $3\times10^{-2}$, $9.1\times10^{-7}$ and $2.9\times10^{-5}(P_\star/10\,{\rm d})^{-2}$ (for the inner-binary scalings given there).

**Typical magnitudes.** Stellar spin periods are a few to 30 days. LK periods are of order $10^6$ yr scale times the geometry factors in the $\Omega_{pl,0}$ formula. In the Storch & Lai 2015 example, widespread chaos sets in once $\epsilon$ exceeds a value of order unity (the exact value is garbled in the supplied text).

## What the literature has established

**Storch et al. 2014 (arXiv:1409.3247).** Only the abstract and the summary in later papers are available. The paper integrated LK cycles together with stellar spin precession, short-range forces and tides. It found that the stellar spin can evolve chaotically, with "periodic islands in an ocean of chaos" as $M_p$ and $P_\star$ vary, and that the complex evolution imprints on $\theta_{sl,f}$. It introduced the non-adiabatic, adiabatic and trans-adiabatic regimes and the parameter $\epsilon$. It left open the origin of the chaos and any predictive theory.

**Storch & Lai 2015 (arXiv:1411.0739).** The paper took a non-dissipative system with strictly periodic quadrupole LK cycles and no SRFs, octupole or back-reaction. Using Hamiltonian perturbation theory, it derived the resonance condition $\bar\Omega_{ps}=Nn_e$, the widths $w_N$ and the overlap criterion. It checked these against surfaces of section taken once per cycle and against a bifurcation diagram of $\theta_{sl}$ at eccentricity minima as a function of $\epsilon$, starting from $\theta_{sl}=0$. Truncating the forcing to harmonics $|M|\le N_{\max}$ reproduces the full dynamics, and the agreement improves at small $\epsilon$. The theory predicts the extent of the chaotic zone against $\epsilon$. Narrow resonances produce quasi-periodic bands that split phase space. The resonance nearest $p=1$ controls initially aligned systems. The paper also showed restricted chaos, caused by trajectories riding a separatrix even without strict overlap, and proposed "adiabatic resonance advection".

Its own stated limits: the theory is valid only near the adiabatic regime (small $\epsilon$), and the advection demonstration keeps the orbit shape fixed while the real $e_{\min}$ drifts. The specific $\theta_{lb}^0$, the values of $\epsilon$ at the regime transitions and the numerical value of $\mathcal{N}$ are lost from the supplied text.

**Anderson et al. 2016 (arXiv:1510.08918).** This paper is known here only through Storch et al. 2016. It is a population synthesis of HJ formation by LK migration in stellar binaries, including the octupole term, SRFs, tides, magnetic-braking spin-down and spin-orbit feedback. It found that octupole terms mainly raise the tidal-disruption fraction and that most HJs form through the quadrupole channel. Bimodality of $\theta_{sl,f}$ remains "nearly ubiquitous" for Jupiter-mass planets despite the feedback. One of its figures shows the bimodal peaks moving toward $90^\circ$ for initially misaligned systems. Feedback grows in importance with $M_p$ and is discussed in its Section 4.3.

**Storch et al. 2016 (arXiv:1607.03937).** The paper extended the theory to the non-adiabatic regime with SRFs and weak-friction tides, but no feedback, no octupole ($e_b=0$, quadrupole only), no spin-down, and $\theta_{sl,0}=0$. It identified five paths:
- I, bimodal, when the initial governing separatrix is the $N=0$ central island and $\mathcal{A}_0\lesssim1$. The trajectory is forced across the shrinking separatrix and exits through the top or bottom depending on phase.
- II, stationary adiabatic, governed by a top side island. The final angle decreases smoothly as $P_\star$ decreases, and only modest prograde misalignment results.
- III, widespread chaos from resonance overlap.
- IV, adiabatic advection by resonance $N$, which produces striations.
- V, fully adiabatic, with essentially zero misalignment.

Only I and III can give retrograde orbits. The regime map in the $(\mathcal{A}_0,N_{\max,0})$ plane, drawn from a suite varying $M_p$, $P_\star$, $a_b$ and $\theta_{lb}^0$, is: bimodal at $\mathcal{A}_0\lesssim1$; stationary adiabatic at $\mathcal{A}_0\gtrsim1$ with $N_{\max,0}\lesssim1$; advection or chaos at both $\gtrsim1$, with chaos at the lower $N_{\max,0}$. Check these thresholds against the figure, because the inequalities are garbled in the supplied text.

The paper itself flags three open problems: widespread chaos in one $P_\star$ range where $N_{\max}\lesssim1$ is unexplained; the onset of regime V cannot be pinpointed; and neglecting feedback is "the single biggest simplifying assumption". It covers a deliberately unphysical range of $P_\star$. It asserts, without computing it, that spin-down leaves the classification unchanged.

**Anderson et al. 2017 (arXiv:1610.02626).** The same physics applied to stellar triples. The paper ran secular octupole integrations with SRFs on both stars, both spins and back-reaction, but no tides. It sampled solar-mass and brown-dwarf tertiaries, giving 1779 and 1742 retained equal-mass systems. Rejection cuts were the stability criterion, $e_{\lim}$ below a cut, and $a(1-e_{\max})<0.08$ AU.

Analytic results: an LK window with arbitrary angular momentum ratio $\eta=(L_{\rm in}/L_{\rm out})_{e=0}$ and SRFs, in which the window closes when $\eta^2+60-20\varepsilon_{\rm SRF}<0$ (reconstructed from garbled text; it reduces to $|\cos I_0|<\sqrt{3/5}$ as $\eta,\varepsilon_{\rm SRF}\to0$); a limiting eccentricity at $\cos I_{0,\lim}=\tfrac{\eta}{2}(\tfrac45j_{\lim}^2-1)$; $e_{\max}\le e_{\lim}$ even with the octupole term; and a maximum perturber distance $a_{\rm out,eff}<19.6\,{\rm AU}(\dots)(a_{\rm in}/0.1\,{\rm AU})^{4/3}(1-\tfrac53\cos^2I_0)^{1/3}$.

Spin result: appreciable obliquity requires $\mathcal{A}_{\rm A17}<\mathcal{A}_{\rm crit}$, where $\mathcal{A}_{\rm A17}\approx58(\dots)(P_\star/5\,{\rm d})^{-1}(a_{\rm in}/0.1\,{\rm AU})^{-9/2}(a_{\rm out,eff}/10\,{\rm AU})^3$. The transition is abrupt. The numerical $\mathcal{A}_{\rm crit}$ is lost from the supplied text. The paper applied this to DI Herculis and found that a tertiary must be very close, near the stability limit.

**Disagreements.** The papers use different adiabaticity definitions, and Anderson et al. 2017 notes the "fuzziness". Storch et al. 2016 calls its theory "somewhat pedagogical" for low $M_p$ and fast spin, while Anderson et al. 2016 finds that bimodality survives in that regime. No paper reconciles the two statements quantitatively.

## The gap

No one has tested the $(\mathcal{A}_0,N_{\max,0})$ classification against integrations that include spin-orbit feedback, the octupole term and spin-down, and reported a misclassification rate. The bimodal peak locations have not been derived quantitatively from the adiabatic invariant. Capture probabilities for advection with evolving shape functions have not been computed. The chaos at $N_{\max}\lesssim1$ is explicitly unexplained, and the boundary of regime V is unknown.

The pieces now exist. The Hamiltonian, its Fourier decomposition and the overlap criterion are written down. The secular equations with SRFs and tides are standard and cheap: one system takes seconds to minutes in compiled code. The separatrix-crossing probability theory is classical. Storch et al. 2016 supplies a reproducible baseline with a stated regime map. The whole programme needs only computation.

## Data, catalogues and codes

- **arXiv full texts** of the four core papers (https://arxiv.org/abs/1411.0739, https://arxiv.org/abs/1607.03937, https://arxiv.org/abs/1610.02626, https://arxiv.org/abs/1409.3247). Take the exact fiducial parameters, the garbled definitions (of $\mathcal{A}$, $\mathcal{N}$, $w_N$, $\mathcal{A}_{\rm crit}$ and $S/L$) and the figure axes from them. Their appendices and the references there (Liu, Muñoz & Lai 2015 is cited for the SRF equations) give the full vector equations of motion.
- **Anderson et al. 2016 (arXiv:1510.08918)** for the complete equation set with feedback, tides and spin-down, and for its population parameters.
- **REBOUND and REBOUNDx** (https://rebound.readthedocs.io). Use them as an independent direct N-body check of secular breakdown at high $e$. REBOUNDx has spin and tide effects; locate the current release and verify which effects it implements before relying on them.
- **SciPy** (`solve_ivp` with DOP853 or Radau) plus Numba or JAX for the secular vector equations. Write the integrator directly; a secular code is about 200 lines.
- **MESA** (locate the current release) or the **Dartmouth Stellar Evolution Database** (arXiv:0804.4473) for $R_\star(t)$. Apsidal and gyration constants $k_q$ and $k_\star$ are not in the Dartmouth tables, so watch that point.
- **Optional observational comparison:** NASA Exoplanet Archive (https://exoplanetarchive.ipac.caltech.edu) and Southworth's TEPCat obliquity compilation (locate the current release). These give projected obliquities $\lambda$, not $\psi$, so project the simulated $\theta_{sl}$ before comparing.

## Work plan

1. **Recover the definitions.** Download the full texts. Extract the exact forms of $\mathcal{A}_{\rm S16}$, $\mathcal{N}$, the shape functions, $w_N$, the fiducial parameters ($M_p$, $a_0$, $a_b$, $\theta_{lb}^0$, $P_\star$ ranges), $\mathcal{A}_{\rm crit}$ and the $S/L$ factor. Size: about one hour. Check: write each formula into a table with its source equation number. If anything is still ambiguous, derive it from the stated Hamiltonian and record the derivation.
2. **Calibrate on numbers.** From the definitions, reproduce the prefactor 1.19e-8 in $\Omega_{ps,0}$, the prefactor 1.17 in $\epsilon$, the prefactor 58 in $\mathcal{A}_{\rm A17}$, and $19.6$ AU. Reproduce the LK-window limit $\sqrt{3/5}$. Check: agreement to 2 significant figures. If a prefactor fails, look for a convention difference (for example $\hat\Omega_\star$ normalization or the $2k_q/k_\star$ factor) before going on.
3. **Calibrate on dynamics: reproduce the Storch & Lai 2015 bifurcation diagram** ($\theta_{sl}$ at eccentricity minima against $\epsilon$, initial $\theta_{sl}=0$) and one of its surfaces of section. Method: integrate the non-dissipative Hamiltonian equations with the analytic LK solution as forcing; use 500 values of $\epsilon$ and 2000 LK cycles each. Size: minutes. Check: the chaos onset and the periodic islands visually match the published figure. Also compute $w_N$ from the Fourier coefficients and overlay the analytic chaotic-zone boundary. If it fails, check the sign conventions of $N$ and of $\mathcal{C}$ (antisymmetric) and the definition of $\tau$.
4. **Reproduce the Storch et al. 2016 $\theta_{sl,f}$ against $P_\star$ figure** (both orbit shapes) with SRFs and weak-friction tides, no feedback, $e_b=0$. Use about 300 periods times 10 random initial phases, integrated to circularization. Size: CPU-hours. Check: regimes I to V appear at the published $P_\star$ ranges. If they shift, compare $e_{\max}$ with SRFs against the published shape-function figure.
5. **Reproduce the regime map** in $(\mathcal{A}_0,N_{\max,0})$ from about 2000 systems. This is the baseline classifier for SQ1. Write an automatic class labeller (bimodality test across phase ensembles; finite-time Lyapunov exponent for chaos; libration of $\bar\Omega_{ps}/n_e$ about $N$ for advection) and validate it by hand on 50 cases.
6. **SQ1.** Re-run the step 5 sample with feedback on ($d\hat L/dt$ term with $S/L$), then with octupole terms ($e_b>0$), then with spin-down. Produce the confusion matrices. Size: tens of CPU-hours. Check: with feedback artificially scaled to zero, step 5 is recovered exactly.
7. **SQ2.** Compute the area of the $N=0$ separatrix at $t=0$ from $H_0$. Track the separatrix area along the tidal track and find the crossing time when it equals $J$. Predict the exit $p$ for top and bottom crossings, and compare with phase ensembles of 200 per system. Also test initial $\theta_{sl,0}\neq0$.
8. **SQ3.** For regime IV systems, compute the analytic capture probability from the rates of change of the separatrix areas (standard separatrix-crossing theory, check against a textbook). Compare with phase-ensemble capture fractions and $\theta_{sl,f}$ striation locations.
9. **SQ4.** At the unexplained-chaos periods, build surfaces of section at frozen $\mathcal{A}$ and $N_{\max}$, compute Lyapunov exponents, and test candidate mechanisms: overlap of the $N=0$ side islands with $N=\pm1$; Mathieu-type parametric resonance at $N_{\max}\sim1/2$; repeated separatrix crossing. If no mechanism fits, report the measured boundary empirically.
10. **SQ5.** Scan $M_p$ and $P_\star$ to vary $S/L$ at fixed $(\mathcal{A}_0,N_{\max,0})$ where possible, then plot the class-change fraction.
11. **Spot-check with N-body.** Run 20 systems with REBOUND to confirm that the secular approximation holds at the $e_{\max}$ reached. Flag systems where $P_{\rm orb}$ at periastron is comparable to the precession timescales.

## Checks and pitfalls

- **Limiting cases.** As $\epsilon\to\infty$, $\theta_{sb}$ must be conserved. As $\mathcal{A}\to\infty$, $\theta_{sl}$ must be conserved. With feedback, total angular momentum $\vec S+\vec L+\vec L_b$ must be conserved to integrator tolerance, and energy must be conserved in the non-dissipative runs.
- **Convergence.** Halve the tolerance and confirm that $\theta_{sl,f}$ is unchanged in the regular regimes. In chaotic regimes only the distribution should converge, not individual outcomes. Report the number of LK cycles used for surfaces of section.
- **Phase sensitivity.** Classes I and III must be judged over phase ensembles, never from single runs.
- **Convention traps.** Three different "adiabaticity" parameters exist, and they run in opposite directions. Two angles share notation: $\theta_{lb}$ in the planet papers corresponds to $I$ in the triple paper. The planet papers use $x=1-e^2$, the triple paper $j=\sqrt{1-e^2}$. Storch & Lai 2015 sets $\mathcal{A}_0$ to a placeholder (the true mean of $\mathcal{C}$ is zero). The papers use $a_{\rm out,eff}=a_{\rm out}\sqrt{1-e_{\rm out}^2}$ or $e_b=0$. Storch et al. 2016 covers unphysical spin periods, so restrict any population claim to 1 to 30 days.
- **Selection effects.** The triple-paper rejection cuts on $e_{\lim}$ and pericentre shape the reported distributions. Any comparison with observed obliquities needs a sky-projection step and a model of the HJ selection.
- **Secular breakdown.** At $e\to1$ the double-averaging can fail (step 11). Tidal disruption within $\sim2.5R_p(M_\star/M_p)^{1/3}$ removes systems (standard result, check), so do not count them.
- **Usual failure.** Leaving the shape functions frozen while $a$ decays overstates advection. Storch & Lai 2015 itself warns that $e_{\min}$ drifts.

## Known problems in the nearby literature

The audit found no verified errors in this subfield.

## What to report

1. The calibration table (step 2) and the overlays of the reproduced figures (steps 3 and 4).
2. The confusion matrices and regime map for SQ1, the bimodal-peak formula and comparison for SQ2, the capture probabilities for SQ3, the chaos mechanism or empirical boundary for SQ4, and the $S/L$ threshold for SQ5.
3. All code and the random seeds.

The report must list explicitly what was not done, every assumption made (including every formula reconstructed because the supplied text was garbled), and every place where results disagree with Storch & Lai 2015, Storch et al. 2016, Anderson et al. 2016 or Anderson et al. 2017, with the size of each disagreement.

## Reading list

1. Storch & Lai 2015 (arXiv:1411.0739): the Hamiltonian, the resonances and the overlap theory; everything else builds on it.
2. Storch et al. 2016 (arXiv:1607.03937): the five paths and the two-parameter classification to be tested.
3. Anderson et al. 2016 (arXiv:1510.08918): the full equations with feedback and spin-down, and the population context.
4. Anderson et al. 2017 (arXiv:1610.02626): the third adiabaticity convention, the LK window with SRFs, and the stellar-triple analogue.
5. Storch et al. 2014 (arXiv:1409.3247): the original numerical discovery of the chaos.
6. Kozai 1962 (doi:10.1086/108790): the origin of the LK quadrupole secular equations.
7. Storch et al. 2015 (arXiv:1502.06550): an alternative planetary tidal model, if the weak-friction assumption is questioned.

## First message to your assistant

```
Attached is a research brief on stellar spin evolution during Lidov-Kozai
migration of hot Jupiters. Read all of it first. Then do Work plan steps 1
and 2: download the full texts of arXiv:1411.0739, 1607.03937, 1610.02626
and 1510.08918; recover every definition the brief marks as garbled
(the adiabaticity A of Storch et al. 2016, the function N, the resonance
width w_N, A_crit, the S/L factor, and the fiducial parameters); and
reproduce the calibration prefactors (1.19e-8, 1.17, 58, 19.6 AU, and
sqrt(3/5)). Report a table of the formulas with their source equation
numbers and the result of each check before running any integrations.
Flag every assumption you make and every point where you disagree with
the papers.
```
