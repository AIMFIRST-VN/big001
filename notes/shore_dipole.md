# Beach-shore test: cluster-count dipole along the slope axis (2026-09-05)

Script `shore_dipole.py`, output `shore_dipole_results.txt`. Nothing in the paper was changed.

## Idea

A super-horizon slope of amplitude eps along n boosts the abundance of massive clusters on the dense side by
R eps with R = delta_c/sigma^2(M) (peak-background split), i.e. a count dipole of a few percent for eps ~ 1e-2
with R = 3-7. The slope predicts the same sign and amplitude at every depth; local structure (2M++) decays with
depth; the kinematic (Ellis-Baldwin) term (2+x) beta = 3.7e-3 along (264,48) projects to +1.3e-3 on n.
Axis n = (l,b) = (270,-22) from `hypocentre.py` (68% region l 242-301, b -42..-2). Sign: d_n > 0 means more
clusters toward (270,-22); the slope model is d_n = R eps + local + kinematic.

## Method (fixed before looking)

* Fit pixels HEALPix nside 16 (13 deg^2). Footprint: |b| > 20 (area fraction of each pixel from nside-64
  subsampling); PSZ2 analytic only (all sky); eRASS1 additionally l > 180 and an nside-8 occupancy map
  (~30 objects per nside-8 pixel, so occupancy is safe); ACT DR5 and SPT-SZ nside-16 occupancy (footprint-limited,
  checks only).
* Model N_i = n0 A_i (1 + d . r_i); Poisson-weighted linear least squares (weights 1/(nbar A_i)).
  Fixed-axis statistic: two-parameter fit (n0, d_n). Free axis: four-parameter fit.
* Null: 2000 multinomial shuffles of the N objects over footprint pixels with probability proportional to A_i;
  p_fixed = P(|d_n,null| >= |d_n|), p_free = P(|d_null| >= |d|).
* Redshift split z < 0.05 (D < ~215 Mpc), 0.05-0.2, > 0.2 (eRASS1 also 0.2-0.5, > 0.5); the |b| > 30 variant
  as a Galactic-latitude check. R from an EH98 no-wiggle sigma(M) with sigma_8 = 0.81 at the catalogue's typical
  M500 (stated per catalogue); the task's R = 3-7 range is also quoted.
* LCDM local expectation: 2M++ delta_g (Carrick et al. 2015) integrated 10-135 Mpc/h with r^2 dr weighting
  (volume-limited local sample), cluster/galaxy bias ratio 2.5, dipole fitted over the same footprint.

Revised after looking (stated honestly): the eRASS1 flux column is in 1e-14 cgs, so the pre-registered
F500 > 3e-13 cut selected everything; after seeing the density-vs-ecliptic-latitude table (10x more clusters per
deg^2 at the south ecliptic pole than at the ecliptic equator) the cuts were reset to 4e-13, 8e-13, 1.5e-12
before their dipoles were looked at, and an exposure model N_i ~ A_i E_i^gamma (E_i = pixel geometric-mean
catalogue EXP) was added. Subsamples with N < 50 are skipped.

## Results

| Catalogue | N in footprint | fsky | M_typ, R | fixed-axis d_n | null sigma, p_fixed | free dipole | p_free |
|---|---|---|---|---|---|---|---|
| Planck PSZ2 (all) | 1392 | 0.661 | 4.5e14, 3.2 | +0.026 +- 0.050 | 0.050, 0.60 | 0.092 toward (156,-52) | 0.29 |
| PSZ2 SNR > 6 | 590 | 0.661 | | +0.041 +- 0.077 | 0.077, 0.61 | 0.070 toward (22,-82) | 0.83 |
| PSZ2 COSMO | 485 | 0.661 | | +0.002 +- 0.084 | 0.084, 0.98 | 0.050 | 0.94 |
| PSZ2 |b| > 30 | 1077 | 0.500 | | -0.004 +- 0.060 | 0.060, 0.95 | 0.129 toward (117,-37) | 0.20 |
| eRASS1 raw (all) | 12107 | 0.330 | 2.0e14, 2.2 | +2.14 +- 0.05 | 0.026, 0 | 2.29 toward (270,-20) | 0 |
| eRASS1 exposure-corrected | 12107 | 0.330 | | +0.168 +- 0.026 | 0.024, 0 | 0.155 toward (199,-39) | 0 |
| eRASS1 F500 > 4e-13 (raw A) | 2631 | | | +0.231 +- 0.060 | 0.055, 0 | 0.281 toward (252,-15) | 0.001 |
| eRASS1 F500 > 8e-13 | 1238 | | | +0.194 +- 0.087 | 0.082, 0.023 | 0.195 toward (274,-22) | 0.11 |
| eRASS1 F500 > 1.5e-12 | 576 | | | +0.167 +- 0.126 | 0.117, 0.16 | 0.341 toward (301,-4) | 0.04 |
| ACT DR5 (all) | 4195 | 0.301 | 3.0e14, 2.7 | +0.016 +- 0.035 | 0.035, 0.66 | 0.250 toward (289,40) | 0 |
| SPT-SZ (all) | 674 | 0.064 | 3.5e14, 2.9 | +0.61 +- 0.25 | 0.19, 0.009 | 0.69 toward (264,-34) | 0.22 |

Redshift splits (fixed-axis d_n):

| Catalogue | z < 0.05 | 0.05-0.2 | 0.2-0.5 | > 0.5 | z >= 0.05 (far) | no z |
|---|---|---|---|---|---|---|
| PSZ2 | +0.02 +- 0.26 (N=50) | -0.02 +- 0.10 (368) | -0.13 +- 0.08 (572, z>0.2) | | -0.085 +- 0.061, p=0.16 | +0.29 +- 0.09, p=0.002 (402) |
| eRASS1 exposure-corrected | -0.34 +- 0.13 (298) | +0.07 +- 0.05 (3403) | +0.12 +- 0.04 (5304) | +0.46 +- 0.06 (3102) | +0.185 +- 0.027 | |
| ACT DR5 | (4) | -0.20 +- 0.14 (265) | +0.03 +- 0.04 (3926, z>0.2) | | +0.016 +- 0.036 | |
| SPT-SZ | (28) | (33) | +0.39 +- 0.23 (613) | | +0.35 +- 0.22 | |

2M++ local prediction (volume-limited 10-135 Mpc/h, bias ratio 2.5): d_n = -0.09 (PSZ2 footprint),
-0.18 (eRASS1), -0.21 (ACT), +0.15 (SPT). The local field is slightly underdense toward (270,-22) on the
full-sky footprint. Measured z < 0.05 samples are too small to test this (PSZ2 N=50: +0.02 +- 0.26).

Implied eps (d_n minus kinematic, divided by R; error = fit and null in quadrature):

| Sample | eps | 95% |eps| bound (R as stated) | R = 3-7 range |
|---|---|---|---|---|
| PSZ2 all | +0.008 +- 0.022 | < 0.044 | < 0.047 - 0.020 |
| PSZ2 z >= 0.05 | -0.027 +- 0.027 | < 0.071 | < 0.076 - 0.032 |
| ACT DR5 all | +0.006 +- 0.019 | < 0.036 | < 0.032 - 0.014 |
| eRASS1 exposure-corrected | +0.075 +- 0.016 | < 0.10 | < 0.075 - 0.032 |
| SPT-SZ all | +0.21 +- 0.11 | < 0.39 | (check only) |

## Reading

* **PSZ2, the only genuinely full-sky sample, is null on the fixed axis**: d_n = +0.026 +- 0.050 (p = 0.60),
  stable under SNR > 6, the COSMO subsample, |b| > 30 and excluding z < 0.05. The free-axis dipole
  (0.09 toward (156,-52)) is consistent with the shuffle null (p = 0.29). The "no z" subsample carries a
  strong dipole toward (283,-62): redshift follow-up is anisotropic (SDSS in the north), not a signal, and
  it makes the "with z" split (-0.08 +- 0.06) a selection effect rather than a depth test.
* **ACT DR5 fixed-axis is null** (+0.016 +- 0.035, p = 0.66) but its free dipole (0.25 toward (289,40),
  p = 0) is the known depth variation between the deep (D56, BOSS-N) and wide fields; it happens to be nearly
  orthogonal to n. Check only.
* **eRASS1 is dominated by the exposure pattern.** The south ecliptic pole (276,-30) is 8 deg from the axis,
  so the raw dipole (2.1, an artefact: density 3x higher at the hemisphere centre than at its edge) is pure
  depth. After the exposure model (gamma = 0.82, expected ~0.5-1 for a counts-limited survey) the residual
  is +0.17 +- 0.03; the flux-limited samples still have a residual exposure index gamma = 0.09-0.15 and give
  +0.17 to +0.23 +- 0.06-0.13 with the axis fit uncorrected. These residuals GROW with redshift (0.07,
  0.12, 0.46 for 0.05-0.2, 0.2-0.5, > 0.5): the faintest, most distant objects are the most
  exposure-dependent. That is the signature of residual incompleteness, not of a slope (which predicts a
  depth-independent d_n) nor of local structure (which decays with depth). An eRASS1 bound needs the real
  sensitivity map, not a catalogue-derived exposure proxy; the number above should be read as a systematics
  floor of ~0.2 in d_n, i.e. no bound better than |eps| ~ 0.1 from eRASS1 here.
* **SPT-SZ** (2500 deg^2, 19 fields of different depth, patch centred near the axis) gives +0.61 +- 0.25
  (p_fixed = 0.009, p_free = 0.22): a gradient across a small patch of unequal field depths. Check only.
* **Bound.** From PSZ2 (and independently ACT DR5) |eps| < 0.04-0.05 at 95% for R = 3, about 0.02 for
  R = 7. A slope large enough to produce the Sec. 7 amplitudes (eps ~ 1e-2 in the paper's own accounting) is
  allowed; the test is not yet at the interesting level. To reach eps ~ 0.01 one needs a count dipole error
  ~0.03 R^-1 ~ 0.01 on a full-sky mass-limited sample (~10^4 clusters with a modelled selection: eRASS1 with
  its sensitivity map, or PSZ2 + ACT + SPT joined with a common mass limit).

## Caveats

* Footprint systematics dominate low-resolution dipoles: any survey depth pattern that is not in A_i lands in
  d. For eRASS1 the depth pattern is aligned with the test axis, the worst case.
* The |b| > 20 cut removes the plane symmetrically, but Galactic absorption / dust residuals and the CMB noise
  map (SZ selection) are not symmetric; the |b| > 30 variant moves PSZ2 by 0.03, within the null sigma.
* PSZ2 completeness is not uniform (the union catalogue's completeness proxy is the SNR threshold and the
  local noise); the SNR > 6 and COSMO subsamples agree with the full sample, which is the best available check
  short of using the PSZ2 completeness masks.
* The 2M++ prediction assumes a volume-limited local sample and a bias ratio 2.5; it is a scale, not a
  measurement, and the local subsamples are too small to test it.
* R = delta_c/sigma^2 is the high-peak limit of the Lagrangian bias; at M500 ~ 2-4.5e14 (M_vir ~ 1.4x) it is
  2.2-3.2, below the task's 3-7 range which applies to higher masses; both are quoted.
* The pre-registered eRASS1 flux threshold was mis-scaled and reset after a look at the depth pattern (not at
  the dipoles); the exposure model was added after seeing the raw dipole. Those eRASS1 numbers are post hoc.
