# Second ringing in the galaxy power spectrum: pre-registered search

Script `ring_echo.py`; results `ring_echo_results.txt`; figure `figures/ring_echo.png`. Run 2026-09-05.

## Question

Sec. 3.1 (acoustic structure of the kernel, "no-echo" falsifier of H8) allows recompression waves crossing the rings to imprint a ringing in the matter spectrum at a ring-crossing scale, not at the drag-epoch sound horizon. The BAO is a damped sinusoid in k with scale r_d ~ 147 Mpc = 100 Mpc/h. A second system would be a sinusoid sin(k L2 + phi) with a different scale L2 and its own damping. This note records a search for it.

## Data

| Sample | Bins | k range [h/Mpc] | Covariance | Window |
|---|---|---|---|---|
| eBOSS DR16 LRG NGC, SGC, 0.6 < z < 1.0, pre-recon P0 (Gil-Marín et al. 2020) | 32 each | 0.0075–0.315 | 1000 EZmocks (Hartlap 0.967) | W_l²(s) table |
| BOSS DR12 z1, z2, z3 × NGC, SGC, pre-recon P0 (Beutler et al. 2017) | 14 each | 0.0156–0.145 | 2045/2048 Patchy (Hartlap 0.993); published only to k = 0.15 | RR pair counts (smoothed, normalised) |

148 points total. Fiducial cosmology of both releases: flat ΛCDM, Ω_m = 0.31, h = 0.676. URLs in the README data table. DESI DR1 full-shape P(k) with covariance was not used (not found as a plain public file at the time of the run).

## Model and procedure (fixed before the scan)

P0_model(k) = W[ (1 + A2 S(k)) B² P_nw(k) {1 + (O_lin(k/α_z) − 1) e^{−k²Σ_nl²/2}} ] + (1 + A2 S(k)) Σ_i a_i k^i

* S(k) = sin(k L2 + φ) e^{−k²Σ2²/2}; L2 scanned on 61 log-spaced values 20–600 Mpc/h, Σ2 ∈ {0, 5, 10, 20} Mpc/h, φ free.
* BAO template: CAMB linear P(k) for the fiducial cosmology (r_drag = 147.8 Mpc); smooth part = Eisenstein–Hu no-wiggle × degree-7 polynomial in ln k; one α per redshift bin (four), one Σ_nl.
* Broadband: 5 terms (k⁻² … k²) for the 32-bin eBOSS samples, 3 terms (k⁻¹, 1, k) for the 14-bin BOSS samples; B² per sample. 36 linear nuisance parameters, solved exactly at every step (variable projection).
* Window: monopole convolution with the survey window plus the Kaiser quadrupole-leakage term (β = 0.35); bin averaging over ±dk/2; operator verified against a direct k-space convolution with an analytic window to 0.2%.
* Scan: α_z, Σ_nl fixed at the null best fit; with A2 sin(x+φ) = a_c sin x + a_s cos x the model is linear in (a_c, a_s) to first order about the null fit, so each grid point is one weighted linear solve; the best points are re-fitted exactly with everything free.
* Thresholds: detection requires Δχ² > 25 globally; Δχ² > 9 is a hint, not a detection. Otherwise the 95% upper limit on A2 is the worst-phase one-sided limit max_φ [â·u(φ) + 1.645 σ(φ)].
* Look-elsewhere: analytic (χ² with 2 dof, N_eff = 40) and 500 Gaussian mocks drawn from the covariances around the null best fit, each put through the same null refit and scan.

## Results

**BAO check.** α(z1) = 1.013 ± 0.020, α(z2) = 1.019 ± 0.018, α(z3) = 0.987 ± 0.018, α(eBOSS) = 0.997 ± 0.012; all within 1.1σ of unity and consistent with the published values (which are within 1–2% of 1 for these fiducials). Removing the wiggles raises χ² by 94 (9.7σ BAO). χ²_null = 110.0 for 107 dof. Σ_nl = 3.1 ± 2.5 Mpc/h is lower than the expected 6–9 and poorly constrained: a monopole-only pre-recon fit with a free broadband per sample and BOSS truncated at k = 0.15 does not pin the wiggle damping (the mocks return 3.0 ± 1.8, so the low value is a property of this setup, not a feature of the data).

**Scan.** No detection. Maximum Δχ² = 13.1 at L2 = 204 Mpc/h (Σ2 = 0, A2 = 1.2 ± 0.3%; exact refit with all nuisance free gives 13.4 and the same amplitude). Grid points above the hint threshold (Δχ² > 9):

| L2 [Mpc/h] | Δχ² | A2 |
|---|---|---|
| 30–37 (Σ2 = 0) | 9.2–9.4 | 2–4% |
| 204 (Σ2 = 0, 5, 10) | 10–13 | 1.2–2% |
| 271–287 (Σ2 = 5–20) | 9–12 | 1.5–6% |
| 322–430 (various Σ2) | 9–12.5 | 1.4–45% (broadband-degenerate) |
| 600 (Σ2 = 0) | 10 | 20% |

**Global significance.** Analytic: local p = 0.0014 (2 dof), p_global = 0.055 with N_eff = 40. Mocks (500): the maximum Δχ² over the same grid has median 6.7, 95th percentile 12.6, 99th 15.9, largest 18.5; 4.0% of noise-only mocks exceed the data's 13.1, 22% exceed 9, none exceeds 25. The hints are what noise does in a scan of this size. The mock noise peaks cluster at L2 > 150 Mpc/h (see below).

**Upper limits (95%, worst phase, envelope over Σ2).**

| L2 [Mpc/h] | A2 < (median over the band) | range |
|---|---|---|
| 20–30 | 20% | 16–25% |
| 30–60 | 10% | 7–15% |
| 60–100 | 5% | 3.8–6.7% |
| 100–200 | 4.5% | 3.8–5.4% |
| 200–300 | 6.4% | 6–10% |
| 300–600 | 14% | 10–74% |

The envelope is set by the most damped case (Σ2 = 20 Mpc/h); for an undamped or lightly damped ringing (Σ2 ≤ 5) the limit is A2 < 1–2% over 60–250 Mpc/h and < 4% over 30–60 Mpc/h (per-Σ2 columns in `ring_echo_results.txt`, lower-right panel of the figure). For comparison the BAO amplitude in the same units is O_lin − 1 ≈ 5–7% at k ~ 0.05–0.1 h/Mpc before damping. A second ringing at ring-crossing scales of 60–300 Mpc/h with more than about 4–6% amplitude (roughly the BAO's own amplitude) is excluded; for 30–60 Mpc/h the bound is ~10%, and at 20–30 Mpc/h the data have no power (fewer than two bins per oscillation are needed only at the top of the k range, but the eBOSS errors there are 2% and BOSS does not reach it).

Per-dataset scans (eBOSS alone, BOSS alone) give maximum Δχ² of 9.8 each at different scales (600 and 287 Mpc/h), i.e. the joint hint at 204 Mpc/h is not shared by the two surveys individually.

## Caveats

* **Broadband degeneracy at long scales.** For L2 > 300 Mpc/h the k-period 2π/L2 is below twice the bin width (0.01 h/Mpc), so the sinusoid is undersampled and partly aliased into the polynomial broadband; the amplitudes and limits there are unreliable (the 45% at 322 Mpc/h is a degenerate solution). 71% of the mock noise peaks land at L2 > 150 Mpc/h for this reason. The reliable range of the scan is 30–300 Mpc/h.
* **Window function.** Only the monopole convolution with a fixed Kaiser leakage term was modelled; the integral constraint and the exact multipole mixing were not. Residual window error is a smooth function of k and is absorbed by the broadband, but a 1% wobble at the survey scale (L ~ 1/Δk_window ~ 300–1000 Mpc/h) would land in the aliased region anyway.
* **Fibre collisions and non-linear shape at high k.** The eBOSS points at k > 0.2 carry fibre-collision and non-linear-bias shape corrections that the 5-term polynomial only approximates; a second ringing at L2 < 40 Mpc/h lives there and its limit is correspondingly weak and template-dependent.
* **Σ_nl.** Poorly constrained here (see above); the BAO-scale part of the scan (L2 ~ 80–120 Mpc/h) is entangled with the BAO template, and a second system with L2 within ~10% of r_d would be absorbed by α and Σ_nl.
* **Redshift dependence.** L2, φ, Σ2 are shared across the four redshift bins (comoving ring-crossing scale). A ringing that evolves in phase or scale with redshift would be diluted.
* **Gaussian mocks.** The look-elsewhere calibration assumes the published covariances are right and the residuals Gaussian; it does not include template or window systematics.

## Verdict

No second ringing. The pre-registered detection threshold (Δχ² > 25) is not met; the maximum Δχ² of 13 is reproduced by 4% of noise-only mocks and every hint is at the level 22% of mocks produce somewhere in the scan. Any acoustic-structure echo of the kernel in the matter spectrum at scales 60–300 Mpc/h has amplitude below about 4–6% (95%), below or comparable to the BAO's own 5–7%.
