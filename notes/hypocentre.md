# Hypocentre inversion for H5 (Our Position – Where Is Nemo): direction and distance of the seed (2026-09-05)

Script `hypocentre.py`, output `hypocentre_results.txt`. Companion to the row "Hypocentre inversion with distance" in `notes/epicentre_techniques.md`. Nothing in the paper was changed.

## Question

Geiger's method locates an earthquake by guessing a hypocentre, predicting what each station should see, and iterating. Here the "stations" are the observed large-scale axes and amplitudes, the source is the rebound of Kernel 0,0,0 treated as a spherically symmetric explosion, and the hypocentre is the pair (direction n̂, distance D) of the seed from us. `fit_slope.py` fitted n̂ only. This adds D, using the fact that a super-horizon gradient's observable amplitudes fall off with D in known ways, and asks whether the data fix D or only the Johnson–Mehl prior does.

## Physics of the distance dependence

Profile A(r) = A_D (r/D)^s around the seed; observer at r = D; horizon radius d_LS = 14 Gpc (comoving), Hubble radius d_H = c/H0 = 4.4 Gpc, so (d_LS/d_H)² ≈ 10. The linear gradient across the horizon is

    eps = (d_LS/D) · s,   s = d ln A / d ln r at r = D.

Two modes, because "the profile" can be a gravitating or a non-gravitating quantity, and the CMB treats them differently.

**Adiabatic (potential) slope, A = Φ.** By Poisson, the density contrast is δ = (2/3) s(s+1) Φ_D (d_H/D)² (r/D)^(s−2): super-horizon potentials carry a density contrast suppressed by (d_H/D)². Velocities follow v = −(2/3H)∇Φ. Expanding the velocity field to second order about the observer (the observer's own velocity is removed by the CMB frame):

- bulk flow of a sphere of radius R: v_bulk = (R²H/10)|∇δ|, i.e. v/(HR) = (R/10)|dδ/dr| ∝ R/D³;
- dipole of the apparent Hubble rate at depth R: ΔH/H = (3/10) R |dδ/dr|, three times the bulk-flow ratio, also ∝ R;
- Grishchuk–Zel'dovich quadrupole (Sachs–Wolfe 1/3, quadratic term 1/2, Y₂₀ rms 0.3): Q ≈ 0.05 |s(s−2)| Φ_D (d_LS/D)²;
- no intrinsic dipole (Turner 1991; Erickcek, Kamionkowski & Carroll 2008).

Both flow-type signals grow linearly with depth R and are proportional to Φ_D/D³, the quadrupole to Φ_D/D². Their ratio is ∝ 1/D: in principle this ratio fixes D, exactly as the S–P lag does. That is the Geiger content of the test.

**Isocurvature (composition) slope, A = S** (relic fraction or primordial temperature: P12). Then δ_m = S directly, no Poisson suppression: v/(HR) = (R/10) s S_D/D, ΔH/H three times that, intrinsic CMB dipole (1/3) s S_D (d_LS/D) bounded near 10⁻⁴ (Sec. 5.4), quadrupole 0.05|s(s−2)| S_D (d_LS/D)².

**What the LTB toys give for s.** `ltb_inversion.py` uses shell-speed profiles v ∝ (1−x)^(−α), α = 0.05–0.3, which vary by factors 1.4–8 across the ejecta mass coordinate; that is a slope of order unity in the rings, outside the kernel. Inside the kernel the paper assumes homogeneity (Sec. 5, P12), i.e. s ≈ 0 to the flatness limits of Sec. 5.4. So the physically motivated scan s = 0.1–3 covers "kernel interior slightly tilted" to "we sit on the rings' slope"; the amplitude A_D is left free (log-flat over 10⁻⁹–1) and also fixed at 1 (the literal eps = s d_LS/D of the task).

## Data used (amplitudes; direction data as in fit_slope.py)

| Datum | Depth R/d_LS | Value | Error | ΛCDM local floor |
|---|---|---|---|---|
| CF4 bulk flow, 150 h⁻¹ Mpc (Watkins et al. 2023) | 0.0153 | v/(HR) = 0.0258 | 0.0019 | 0.0113 (170 km/s) |
| CF4 bulk flow, 200 h⁻¹ Mpc (paper Table 3) | 0.0204 | 0.0214 | 0.0018 | 0.0075 (150 km/s) |
| Cluster H₀ dipole (Migkas et al. 2021; 9% peak-to-trough) | 0.030 | ΔH/H = 0.045 | 0.015 | 0.02 |
| SN H₀ dipole, Pantheon+/DES-SN | 0.036 | 0.01 | 0.01 | 0.01 |
| Planck kSZ bulk flow < 254 km/s at 1 Gpc | 0.071 | one-sided | | |
| GZ quadrupole ΔT/T < 10⁻⁵ | 1 | one-sided | | |
| Intrinsic dipole < 10⁻⁴ (isocurvature only) | 1 | one-sided | | |

Two treatments: "signal" (the observed flows are entirely the slope; measurement errors only) and "floor" (ΛCDM local structure added in quadrature, the honest one). "Merged" drops the cluster dipole as a separate datum (it traces the same structure as the CF4 flow). One-sided bounds enter as a Gaussian CDF with the 95% limit at 1.645σ.

Direction likelihood: the p-values of fit_slope.py are scores, not likelihoods (small p = good fit), so each observed axis gets a proper density given n̂: von Mises–Fisher exp(κ cos θ) for the along-n̂ tests (CF4 flow, hemispherical power asymmetry, cluster H₀ dipole), a Watson girdle exp(−κ (â·n̂)²) for the planarity ring, κ = 1/σ² with σ = 30° (sensitivity 20°–45° reported). Fisher's combined p is still computed and reproduces fit_slope.py: best (262°, −30°) at nside 16, p = 5 × 10⁻⁵.

Prior on D: an offset-only Johnson–Mehl Monte Carlo (same process as `jm_cells.py`, 4000 realizations) gives offset/L 10/50/90% = 0.34/0.66/1.05, matching the paper's 0.35/0.67/1.05. D = offset × L with L = 8 d_LS (the paper's lower limit), L log-flat over 8–32 d_LS, and a flat-in-log-D prior for comparison. Grid: nside 16 (3072 directions), D ∈ [0.5, 50] d_LS (80 log steps), s ∈ {0.1, 0.3, 1, 3}, A_D 91 log steps. Runtime 2 minutes, memory well under 1 GB.

## Results

**Direction.** MAP n̂ = (270°, −22°) for σ = 20°, 30° and 45° alike; 68% sky region 1900 deg² (l 242°–301°, b −42°…−2°), 95% region 5100 deg². The void centroid (305°, −30°) is outside the 68% region and inside the 95%. Merging the cluster dipole with the flow moves the MAP to (264°, −27°) and widens the 68% region to 2800 deg². This is fit_slope.py's answer with a proper likelihood, and it is independent of D: the four axis tests carry no distance information, so the joint posterior factorises.

**Distance: the amplitudes do not constrain D.** With A_D marginalised, the amplitude likelihood varies by less than 0.7 in ln L over the whole grid 0.5–50 d_LS in every mode and treatment; the Kullback–Leibler information of the amplitude-only posterior against a flat-log prior is 0.01 nats. The posterior on D is the prior:

| Prior | MAP D/d_LS | median | 68% | 95% | P(2.8 < D/d_LS < 8.4) |
|---|---|---|---|---|---|
| JM, L = 8 d_LS | 6.1 | 5.2 | 3.2–7.3 | 1.6–9.4 | 0.81 |
| JM, L log-flat 8–32 d_LS | 9.8 | 10.1 | 5.3–18.7 | 2.6–29 | 0.36 |
| flat in log D | 50 (edge) | 6.8 | 1.3–28 | 0.6–46 | 0.24 |

Identical to two decimals for adiabatic vs isocurvature, signal vs floor, merged vs not. The inferred D is consistent with the tessellation's 0.7 L (5.6 d_LS at L = 8 d_LS) because that is what the prior says, not because the data say it.

**What the amplitudes do constrain: the flatness of the profile at our position.** The horizon-scale gradient eps = s A_D d_LS/D is bounded at 95%:

| Mode | s = 0.3 | s = 1 | set by |
|---|---|---|---|
| adiabatic | 7.6 × 10⁻³ | 1.6 × 10⁻² | GZ quadrupole |
| isocurvature | 4.9 × 10⁻⁴ | 4.9 × 10⁻⁴ | intrinsic dipole |

At D = 5.6 d_LS this means s·A_D < 0.04–0.09 (adiabatic) or < 3 × 10⁻³ (isocurvature): the rebound profile must be flat to a few per cent (potential) or a few per mille (composition) at our position, which is P12's flatness statement made quantitative for an observer 0.7 L from the seed. With an order-unity profile (A_D = 1), the quadrupole alone requires D > 71 √|s(s−2)| d_LS (25 d_LS for s = 0.1, 40 d_LS for s = 0.3, beyond the grid for s ≥ 1) and the intrinsic dipole, for isocurvature, D > 3300 s d_LS. An order-unity slope at 0.7 L is excluded; the seed at 0.7 L is allowed only if the profile is nearly flat there.

**The flows are not the slope.** If the observed flow-type amplitudes were the slope with A_D = 1, they would need D/(s S_D) = 0.06–0.2 d_LS (CF4 at two depths, cluster dipole) and 1.1 d_LS (SN dipole) in the isocurvature mode, or D³ ∝ 0.4–7 d_LS³ in the adiabatic mode: a source inside the horizon, contradicted by the quadrupole and dipole bounds by 3–5 orders of magnitude. Treating the flows as pure slope signal costs Δ ln L = −163 relative to the floor treatment at every D and A_D, and the maximum-likelihood A_D is zero; the amplitude machinery simply says "these flows come from something else", which is the mapped attractors of Sec. 7.

**Wadati variant.** Converting the four flow-type data to a common Hubble-dipole equivalent (bulk flow × 3) and fitting amplitude versus depth: a slope predicts a line through the origin (signal ∝ R). The weighted linear fit gives intercept 0.127 ± 0.012 and slope −3.1 ± 0.5 per d_LS (χ² = 1.0 for 2 dof), while the through-origin fit has χ² = 121 for 3 dof (Δχ² = 120). The signal is the intercept, and it falls with depth, reaching zero at R ≈ 0.04 d_LS ≈ 570 Mpc. The amplitudes extrapolate to a local origin at a few hundred Mpc, not to a super-horizon one. This is the depth-persistence criterion of P4 (Off-Centre Flow) read the other way: the direction persists (Sec. 7) but the amplitude decays.

## Verdict

- Direction: (270°, −22°), 68% region ≈ 1900 deg² toward the southern void excess, same as fit_slope.py and subject to its 2–3σ global significance and its dependence caveats.
- Distance: unconstrained by the data; posterior = JM prior, D = 5.2 (3.2–7.3) d_LS for L = 8 d_LS, or 10 (5–19) d_LS if L is let float to 32 d_LS. Consistent with 0.7 L by construction.
- The amplitude data constrain the combination s·A_D/D (flatness at our position), not D. They exclude the observed flows as a slope signal at high odds and require the profile to be flat to 10⁻² (adiabatic) or 10⁻³ (isocurvature) across our horizon.

## Caveats

- The direction tests carry no D information at all; the "joint" posterior is a product of a sky map and a one-dimensional prior. A true Geiger inversion needs the ratio of a quadrupole-type to a flow-type signal from the same source, and at present the quadrupole is a bound and the flows are local.
- Dependence between tests: the CF4 flow, the cluster dipole and the 2M++ void centroid trace one structure; the merged variant (cluster dipole dropped) is the more honest direction posterior. The two CF4 depths share a catalogue and are correlated; their errors were taken as independent, which overstates the Wadati Δχ² but not its sign.
- Look-elsewhere: the direction MAP inherits fit_slope.py's global p = 0.003 (2.8σ) from 400 random-axis nulls; nothing here adds significance.
- The amplitude model is leading-order in R/D and (d_LS/D) with order-unity coefficients (1/10, 3/10, 0.05) derived for a matter-dominated growing mode; a factor of two either way changes the flatness limits, not the conclusions. The ΛCDM floors (150–170 km/s bulk flow, 2% cluster dipole, 1% SN dipole) are round numbers.
- Watkins et al. (2023) values at 150 h⁻¹ Mpc (387 ± 28 km/s) are from memory of that paper and should be checked; the paper's own Table 3 value (428 ± 36 at 200 h⁻¹ Mpc) is used as stated.
- The L ≥ 8 d_LS bound sets the D scale entirely; the prior's spread in D is dominated by the unknown L, not by the offset distribution.
