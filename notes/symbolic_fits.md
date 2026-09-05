# Symbolic regression on the two toy grids (2026-09-05)

Driver: `symbolic_fits.py` (rerun with `--pysr` for the PySR front); output `symbolic_fits_results.txt`.
Samples: `jm_cells_dump.py` (2 seeds x 300 realizations) and `jm_sizes_dump.py` (3 seeds), copies of the
originals that also save per-realization / per-cell arrays into `symbolic_runs/*.npz`.

**Engine.** Two engines. (1) A hand-written enumeration over a library of simple forms (logistic, erfc, power
law with threshold, saturating exponential, rational, quadratic, the tautology 1 - (f_c/f_eff)^3, linear in
KE/|U|), each fitted by least squares from several starts, ranked by BIC, with leave-one-f-out (LOO) cross
validation. (2) PySR 2.2.1 (Julia backend, installed under `~/.julia`, 60 iterations, 4 procs) as an
independent check on Grid 1. PySR's Pareto front never beat the library forms and its ballistic front is
junk, so the library results are the ones reported. Distributions in Grid 2 were fitted by maximum likelihood
and scored by the Kolmogorov-Smirnov (KS) statistic.

## Grid 1: energy-conserving rebound toy (`contact_grid_results.txt`)

28 grid points (f = 0.3..0.9, sigma = 0, 0.3, 0.5, 1.0). The three contact thresholds give identical
numbers and were pooled. Seed scatter (3 seeds) is sigma_seed ~ 0.03 (contact) and ~ 0.04 (ballistic).

### (a) Contact runs

| form | k | RMS | RMS/sigma_seed | LOO RMS | BIC |
|---|---|---|---|---|---|
| quadratic polynomial (6 coeff) | 6 | 0.019 | 0.60 | 0.025 | -203 |
| **1 - exp(-1.23 (f sqrt(1+sigma^2) - 0.51)_+^1.09)** | 3 | 0.024 | 0.78 | 0.031 | -199 |
| 0.98 (f - 0.58 + 0.23 sigma)_+^0.86 | 4 | 0.023 | 0.74 | 0.028 | -198 |
| **x^2/(x^2 + 0.70^2), x = (f sqrt(1+sigma^2) - 0.38)_+** | 2 | 0.026 | 0.84 | 0.027 | -198 |
| 0.77 (f sqrt(1+sigma^2) - 0.545)_+^0.755 | 3 | 0.026 | 0.85 | 0.035 | -194 |
| 0.5 erfc((0.94 - f)/sqrt(2(0.20^2 + sigma^2 f^2))) (single-shell escape) | 2 | 0.035 | 1.1 | 0.037 | -182 |
| logistic in f, sigma-dependent width (baseline) | 3 | 0.050 | 1.6 | 0.076 | -158 |
| clip(0.24 KE/|U| - 0.07) (baseline) | 2 | 0.045 | 1.5 | 0.051 | -168 |
| max(0, 1 - (0.75/f_eff)^3) (the tautology) | 1 | 0.109 | 3.5 | 0.113 | -121 |

Reading: the unbound fraction of the contact runs is a function of the single combination
f_eff = f sqrt(1 + sigma^2) = sqrt(KE/2|U|) to within the seed scatter, with a threshold at f_eff ~ 0.5-0.55
and an almost linear rise above it (exponent 0.76-1.09 depending on the form). The empirical U = 0.1 crossing
sits at f_eff = 0.60-0.66 for every sigma (KE/|U| = 0.73-0.88), which is the compact statement. LOO-RMS is
close to the in-sample RMS for every form, so these are not grid interpolations. The two bolded forms are
the ones worth quoting (2-3 parameters, residual below the seed scatter). The quadratic polynomial wins BIC
only by using six coefficients and is not worth quoting. The tautology 1 - f_eff^-3 is 3.5 sigma_seed off;
the logistic baseline is 1.6 sigma_seed off. PySR's best compact expression, (f - 0.308^(1.865^sigma))^2 with
loss 0.0013 (RMS 0.036), is worse than the library forms at the same complexity.

Caveat: the grid stops at f = 0.9 and U ~ 0.5, so the saturation to U = 1 is unconstrained; the forms differ
there (rational -> 1 slowly, exp-sat -> 1 at f_eff ~ 2.5, power law -> clipped). Quote the threshold and the
slope, not an extrapolation.

### (b) Ballistic null

No compact form. Best library RMS is 0.13-0.15 (3-4 sigma_seed) even for the six-coefficient polynomial,
LOO 0.18; PySR's front reaches only loss 0.016 (RMS 0.13). Cause: U_ballistic is non-monotone in f
(0.36, 0.35, 0.36, 0.30, 0.14, 1.00, 1.00 at sigma = 0) because shell crossing through the caustic decides
which shells end up unbound; the sigma = 0 column at f <= 0.6 is a fixed 0.30-0.36 that no smooth function
of the kick energy reproduces. Do not quote a formula; quote the table row or the two facts "about 1/3
unbound for f <= 0.6, all unbound at f >= 0.8".

### (c) Contact minus ballistic

Difference is -0.36 to 0 for f <= 0.7 and -0.71 to -0.10 at f = 0.8-0.9, roughly -A exp(-b sigma) with
A = 0.17-0.4, b = 1.25, but RMS 0.12-0.15 (2.3-2.9 sigma_seed). The ratio is undefined where contact U = 0.
Not worth quoting beyond "contact suppresses the unbound fraction by 0.1-0.7, most strongly at small sigma".

## Grid 2: Johnson-Mehl tessellation (I = c = 1, L = 1)

### (a) Nearest-wall distance w/L (600 realizations)

Weibull CDF **1 - exp(-(w/0.172)^1.26)**, KS 0.047, or Gamma(shape 1.42, scale 0.112), KS 0.048. The pure
exponential 1 - exp(-w/0.160) has KS 0.079 and over-predicts the small-w tail (P(w < 0.05): 0.20 measured,
0.27 exponential); the shape parameter above 1 is real, but the ray grid (200 directions, step 0.0075 L)
also suppresses w below ~0.02 L, so quote the Weibull as descriptive only. Mean 0.160 L, median 0.125 L.

### (b) Observer offset r/L (600 realizations) -- exact result, zero parameters

The arrival time A of the first front at a random point has P(A > a) = exp(-pi I c^3 a^4/3); given A the
winning seed is uniform in the past light-cone, so

    p(r) = 4 pi r^2 int_r^inf exp(-pi a^4/3) da = pi (3/pi)^(1/4) r^2 Gamma(1/4, pi r^4/3)      (L = 1).

KS 0.026 against the samples with nothing fitted; mean 0.672 (measured 0.671), median 0.664 (measured 0.655).
Fitted 1-parameter competitors: r^2 exp(-(r/0.747)^3) KS 0.027, Maxwell KS 0.044; a 3-parameter generalized
gamma reaches KS 0.019 only by fitting. **This is the one to quote in the paper**: it replaces the
10/50/90 percentile line for "where is Nemo" with a formula that follows from the model.

Side numbers: neighbours touched 14.9 +/- 4.6; age differences 10/50/90 = -0.22/0.11/0.56 L/c; fraction of
neighbours older 0.33.

### (c) Cell volume V/L^3 (552 interior cells, 3 seeds)

Gamma(shape 0.64, scale 1.50), KS 0.063; Weibull(k 0.76, scale 0.84), KS 0.070; exponential KS 0.11;
lognormal KS 0.13. The Kiang gamma (shape 5.6) for Poisson-Voronoi is rejected outright, KS 0.35: Johnson-Mehl
cells are far broader (sd/mean = 1.04 versus 0.42 for Poisson-Voronoi) because late seeds fill small gaps.
Measured mean 0.964 L^3 against the analytic Johnson-Mehl 1/N = 4 (pi/3)^(1/4) / Gamma(1/4) = 1.116 L^3;
the 14 % deficit is the interior-box truncation in `jm_sizes.py` (cells of interior seeds are only counted
inside the sampling box), not a discrepancy. Quote "Gamma with shape ~0.64, i.e. sd ~ mean, 20 % of cells
below 0.1 L^3, 7 % below 0.01 L^3" and the analytic mean 1.12 L^3.

### (d) Birth time versus final volume

Binned medians (t = 0.09 .. 1.1): 2.12, 1.08, 0.50, 0.14, 0.05, 0.009 L^3. Best fit
**median V(t) = 2.38 exp(-(t/0.361)^1.53)** L^3, RMS 0.026 dex over the bins (BIC -28.5); the alternative
3.13 (1 - t/1.76)^5.9 is nearly as good (0.032 dex). Simple exp(-t/tau) (0.12 dex) and the naive
exp(-pi t^4/3) untransformed-fraction scaling (0.59 dex) fail. Per-cell scatter about the relation is 0.47 dex,
so this is a median relation, not a one-to-one map. Worth quoting as the stretched exponential with the
0.47 dex scatter. Birth times of real seeds: 10/50/90 = 0.09/0.40/0.88 L/c; the density exp(-pi t^4/3) predicts median 0.45.

## Artefact warnings

- Grid 1 has 7 f values; LOO over f was run for every form and reported. The ballistic and difference fits
  are interpolations of noise and are not offered as formulae.
- Grid 1 contact: seeds only 3, min-max spread used as the scatter estimate (range/1.69).
- Grid 2a: ray resolution floor at w ~ 0.02 L; 2c: volume truncation lowers the mean by ~14 %.
- Runtime: JM cells 2 x 300 realizations ~ 45 min wall each (parallel); sizes 3 seeds ~ 3 min; PySR ~ 10 min.
