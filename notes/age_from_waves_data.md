# Verified inputs for age_from_waves.py (checked 2026-09-05 against the arXiv full texts)

Format: memory value -> published value; note.

## f sigma_8(z) compilation

| z | memory | published | source (arXiv) | note |
|---|---|---|---|---|
| 0.067 | 0.423 +- 0.055 | 0.423 +- 0.055 | 6dFGS, Beutler+ 2012, 1204.4725 | unchanged |
| 0.15 | 0.49 +- 0.145 | 0.49 +0.15/-0.14 | SDSS MGS, Howlett+ 2015, 1409.3238 | unchanged (symmetrised); eBOSS final table (2007.08991, Tab. 3) quotes MGS 0.53 +- 0.16 from a later reanalysis |
| 0.38 | 0.497 +- 0.045 | 0.497 +- 0.039 (stat) +- 0.024 (sys) | BOSS DR12 consensus BAO+FS, Alam+ 2017, 1607.03155, Tab. 7 | quadrature total 0.046; eBOSS final table re-quotes 0.500 +- 0.047 |
| 0.51 | 0.459 +- 0.038 | 0.458 +- 0.035 +- 0.015 | same | total 0.038; eBOSS table: 0.455 +- 0.039 |
| 0.61 | 0.436 +- 0.034 | 0.436 +- 0.034 +- 0.009 | same | total 0.035; eBOSS table: 0.448 +- 0.043 |
| 0.70 | 0.473 +- 0.041 | 0.473 +- 0.041 | eBOSS LRG consensus (Bautista+ 2021 2007.08993; Gil-Marin+ 2020 2007.08994), Alam+ 2021 2007.08991 Tab. 3 | unchanged |
| 0.85 | 0.315 +- 0.095 | 0.315 +- 0.095 | eBOSS ELG (de Mattia+ 2021 2007.09008; Tamone+ 2020 2007.09009) | unchanged |
| 1.48 | 0.462 +- 0.045 | 0.462 +- 0.045 | eBOSS QSO (Hou+ 2021 2007.08998; Neveux+ 2020 2007.08999) | unchanged |
| 0.44 | 0.413 +- 0.080 | 0.413 +- 0.080 | WiggleZ, Blake+ 2012, 1204.3674, Tab. 1 | unchanged; three slices are correlated (paper gives covariance) and overlap BOSS volume |
| 0.60 | 0.390 +- 0.063 | 0.390 +- 0.063 | same | unchanged |
| 0.73 | 0.437 +- 0.072 | 0.437 +- 0.072 | same | unchanged. (Blake+ 2011, 1104.2948, four-slice version: 0.42+-0.07, 0.45+-0.04, 0.43+-0.04, 0.38+-0.04 at z=0.22,0.41,0.60,0.78 -- not used, superseded by the AP-marginalised 2012 values) |
| 0.60 | 0.55 +- 0.12 | 0.55 +- 0.12 | VIPERS, Pezzotta+ 2017, 1612.05645 | unchanged |
| 0.86 | 0.40 +- 0.11 | 0.40 +- 0.11 | same | unchanged |
| 1.40 | 0.482 +- 0.116 | 0.482 +- 0.116 (fixed sigma_v); 0.494 +0.126/-0.120 (sigma_v marginalised) | FastSound, Okumura+ 2016, 1511.08083 | unchanged |

### DESI DR1 full-shape (added)
DESI 2024 V (2411.12021), Tab. 9 (ShapeFit+BAO, MAP) reports f sigma_s8 / (f sigma_s8)_fid per bin; fiducial f sigma_s8 from Tab. 11.
f sigma_8 below = ratio x fiducial (symmetrised errors). This ignores the small sigma_s8 -> sigma_8 rescaling
(alpha_iso and m-correction, Sec. 4.2.3), so treat as approximate (few-per-cent level); DESI does not publish a plain f sigma_8 table.

| z_eff | ratio (Tab. 9) | fid f sigma_s8 (Tab. 11) | f sigma_8 used |
|---|---|---|---|
| 0.295 BGS | 0.84 +- 0.19 | 0.4723 | 0.397 +- 0.090 |
| 0.510 LRG1 | 1.16 +- 0.13 | 0.4733 | 0.549 +- 0.062 |
| 0.706 LRG2 | 1.04 +0.11/-0.092 | 0.4608 | 0.479 +- 0.047 |
| 0.919 LRG3 | 0.997 +0.10/-0.084 | 0.4398 | 0.438 +- 0.040 |
| 1.317 ELG2 | 0.945 +0.097/-0.077 | 0.3944 | 0.373 +- 0.034 |
| 1.491 QSO | 1.16 +- 0.12 | 0.3750 | 0.435 +- 0.045 |

Caveat: DESI DR1 LRGs overlap the BOSS footprint and redshift range, so BOSS+DESI are not independent; the script
has a switch (USE_DESI) and both results are reported.

## S8
memory: 0.776 +- 0.017, attributed to DES Y3 + KiDS-1000 joint cosmic shear.
published: DES Y3 + KiDS-1000 joint cosmic shear (DES & KiDS Collaborations 2023, 2305.17173): S8 = 0.790 +0.018/-0.014 (MAP 0.801).
The memory number 0.776 +- 0.017 is the DES Y3 3x2pt value (Abbott+ 2022, 2105.13549), misattributed. Corrected to 0.790 +- 0.016 (symmetrised).

## Omega_m h^2, CMB-independent (galaxy full-shape)
memory: 0.142 +- 0.006 "BOSS/eBOSS full-shape".
- Ivanov, Simonovic & Zaldarriaga 2020 (1909.05277), BOSS P(k), BBN prior on omega_b: Omega_m = 0.295 +- 0.010, H0 = 67.9 +- 1.1 -> Omega_m h^2 ~ 0.136.
- d'Amico+ 2020 (1909.05271), BOSS P(k), Omega_b/Omega_c held fixed (no free omega_b): Omega_m = 0.309 +- 0.010, H0 = 68.5 +- 2.2 -> ~0.145.
- Philcox & Ivanov 2022 (2112.04515), BOSS P+Q0+BAO+B0, BBN prior, Tab. III: omega_cdm = 0.141 +0.011/-0.013, h = 0.696 +- 0.011, Omega_m = 0.338 +0.016/-0.017 -> Omega_m h^2 = 0.164 +- 0.012 (n_s free; with Planck n_s, H0 = 68.31).
- DESI 2024 VII (2411.12022), DESI DR1 FS+BAO + BBN prior + loose n_s prior, Eq. 3.1: Omega_m = 0.2962 +- 0.0095, H0 = 68.56 +- 0.75, sigma_8 = 0.842 +- 0.034 -> Omega_m h^2 = 0.139 +- ~0.005 (error propagated ignoring the Omega_m-H0 correlation).
All of these use a BBN (or fixed) baryon density; none is a pure k_eq measurement. The spread 0.136-0.164 exceeds the memory error.
Adopted: 0.139 +- 0.006 (DESI DR1, the tightest and newest; error rounded up).
