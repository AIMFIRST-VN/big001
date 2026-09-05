# Archival dark-lens search for the Kernel 0,0,0 remnant (Sec. 3.1, "A remnant?")

Script `core_search.py`; full output `core_search_results.txt`; data URLs in the README data table.

## What was searched
DES Y3 Kaiser-Squires E-mode convergence map (Jeffrey et al. 2021; full source sample, HEALPix nside 1024,
4742 deg^2 inside the released mask). KiDS-1000 and HSC Y3 were not used: neither collaboration serves a
ready-made convergence map at a public URL (only shear catalogues, far above the download budget), so this
is a DES-footprint-only search (about 11% of the sky).

Counterpart catalogues (entries inside the DES mask in brackets): DES Y3 redMaPPer lambda>=20 (20840),
DES Y1 redMaPPer (6476), Planck PSZ2 union (248), ACT DR5 (1775), SPT-SZ 2500d Bocquet+19 (550),
eROSITA eRASS1 primary clusters (4604; covers the western Galactic hemisphere, i.e. only part of DES),
MCXC (277), 2MRS groups Tully 2015 (2870, z<0.03).

## Thresholds fixed before looking at the map
- Gaussian smoothing FWHM 5, 10, 20 arcmin (normalised convolution with the smoothed mask).
- Noise: 1.4826 x MAD of the smoothed map over the survey interior (smoothed-mask fraction > 0.95).
  This includes real structure, so S/N is conservative. No noise-rotation maps are in the release.
  Note: sigma barely changes with smoothing (0.00248, 0.00245, 0.00234), so the released map is already
  band-limited/smoothed at roughly the 10-20 arcmin scale; the 5 arcmin pass adds little.
- Peaks: local maxima (8 neighbours) with S/N >= 4 listed, S/N >= 5 is the headline threshold.
- Match radius 5 arcmin + FWHM to any catalogue entry.
- Edge flag: smoothed mask fraction < 0.95 at the peak or peak within one FWHM of the mask boundary.
- Enclosed mass: M_ap = Sigma_crit(z_l, z_s = 0.63) x sum over the 2 FWHM disc of (kappa - annulus mean at
  2-3 FWHM) x pixel area at z_l, tabulated for z_l = 0.02, 0.05, 0.1, 0.2, 0.4. Order-of-magnitude only:
  the KS map is noise-dominated, the smoothing spreads a compact lens, and the annulus subtraction
  removes part of an extended profile (so the 5 arcmin masses are low by design).

## Result
| FWHM | peaks S/N>=4 | S/N>=5 | matched | unmatched |
|---|---|---|---|---|
| 5'  | 10 | 2 | 7 | 3 |
| 10' | 11 | 3 | 11 | 0 |
| 20' | 10 | 2 | 10 | 0 |

Every S/N >= 5 peak is a known cluster (Abell 133 region at RA 14, Dec -1.4; the A3128/A3125 complex at
RA 49.6, Dec -44.4; A3827/A3822 field at RA 326.6, Dec -57.1; etc.), each matched by three or more
catalogues.

Three peaks are unmatched, all at the 5 arcmin scale only, all S/N 4.1-4.4, none edge-flagged:

| RA | Dec | (l,b) | S/N | M_ap (z_l=0.1) | nearest counterpart just outside the 10' match radius | 2M++ LOS 20-200 Mpc/h |
|---|---|---|---|---|---|---|
| 49.552 | -45.537 | (255.1,-55.6) | 4.42 | 1.3e13 | redMaPPer z=0.100 lambda 63 at 12.1'; eRASS1 z=0.077 at 13.1'; 2MRS group 4.9e15 Msun at 12.5' | underdense, mean delta -0.28 |
| 28.477 | -6.843 | (161.9,-64.9) | 4.37 | 1.1e13 | redMaPPer z=0.632 lambda 20 at 14.2' | overdense, +0.24 |
| 71.297 | -44.994 | (250.3,-40.7) | 4.09 | 0.8e13 | eRASS1 z=0.279 at 11.3'; redMaPPer z=0.237/0.671 at 17-19' | overdense, +0.13 |

At the 10 and 20 arcmin scales the same three positions reappear as peaks and are matched. The first is
an outskirt of the A3128 complex (68' from the S/N 5.4 peak) and sits in a 2M++ underdense line of sight
only because the 2M++ box ends at z~0.07 and the structure is behind/at its edge. The expected number of
pure-noise local maxima above 4 sigma in ~4400 deg^2 at 5 arcmin resolution is of order 10-20, so three
unmatched 4.1-4.4 sigma peaks are what noise alone gives. 2M++ predicted LOS velocities toward them are
+36, -130, +34 km/s (nothing coherent); angles to the CF4 bulk-flow axis are 58, 100, 53 degrees.
CF4 per-group velocities were not used: `table3.dat` has no header in this checkout and a trial parse
gave a direction-independent zero-point offset.

## Verdict
Null. No convergence peak with S/N >= 5 lacks a cluster counterpart in the DES Y3 footprint; the three
S/N ~4 unmatched peaks are noise-level and lie within ~15' of catalogued clusters. Mass sensitivity:
a compact lens with M_ap(2 FWHM) >~ 1e14 Msun at z_l ~ 0.1 (>~ 1e13 at z ~ 0.02) would have reached
S/N 5 in the 10' map; the bound-table remnant mass at 100 Mpc (2e46 kg = 1e16 Msun) is excluded in
this 11% of the sky, and a remnant at 10 Mpc (2e43 kg = 1e13 Msun) would be below the detection threshold
at its expected 10-arcmin-scale convergence. Nothing here bears on the other 89% of the sky.
