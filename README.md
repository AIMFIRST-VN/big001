# The Big Pool Bangs Theory

Speculative cosmology paper: a jammed hard-sphere Planck-relic core, a loop-quantum-cosmology-type bounce, chaotic ejection as the origin of dark matter, and a set of hypotheses about where we sit inside our own Big Pool Bang and where a second one might be (the Cold Spot direction). Every hypothesis carries explicit falsifiers and a pre-registered test list.

The framework's decidable questions, with status, are collected in [CELESTIAL_QUESTIONS.md](CELESTIAL_QUESTIONS.md) (twenty-three, four groups).

## Paper

Single current version at the head of the repo: `big-pool-bangs.md` (source), `big-pool-bangs.tex`, `big-pool-bangs.pdf`.

Versioning uses git tags and GitHub Releases: every released version is a tag (`v5-2`, ...) with the PDF attached under https://github.com/AIMFIRST-VN/big001/releases. Older drafts (v1–v5-1) remain in the git history.

Build: `pandoc big-pool-bangs.md -s -o big-pool-bangs.pdf --pdf-engine=pdflatex -V geometry:margin=2.5cm -V fontsize=11pt  (title and subtitle come from the YAML header)`

## Analysis scripts

| Script | Purpose |
|---|---|
| `bounce_model.py` | LQC-form modified Friedmann bounce |
| `nbody_ejection.py`, `nbody_turbulent.py`, `ejecta_analysis.py` | Toy N-body ejection, ejected fraction, velocity sorting |
| `nucleation_kicks.py` | Poisson / cascade nucleation kick spectra |
| `sn_hubble.py`, `sn_des.py` | Pantheon+ and DES-SN5YR distance–redshift comparisons |
| `wz_effective.py`, `ltb_relativistic.py`, `ltb_inversion.py` | Toy w(z), LTB lightcone, profile inversion |
| `coldspot_profile.py`, `spot_locus.py` | Cold Spot ring profile; extreme-spot locus test |
| `curl_qe.py` | Temperature-only quadratic-estimator curl search at the spots |
| `neighbour_tail.py` | Point-mass Sachs–Wolfe tail (low-ℓ template) and tidal bulk-flow direction |
| `rings_hydro.py` | First 1-D hydro toy (RETRACTED numbers: kick profile and band-finder artefacts, see paper Sec. 3.2); kept for the record with `rings_hydro_scan.txt` |
| `rings_hydro2.py` | Controlled 1-D hydro toy: kick at a fixed fraction of the local escape speed, pressure-to-gravity ratio beta, fixed-contrast edge finder, pressureless run as ballistic null |
| `jm_cells.py` | Johnson–Mehl nucleation-and-growth Monte Carlo for the pool of bangs: nearest-wall distance, neighbour count, age differences, observer offset; results in `jm_cells_results.txt` |
| `rings_hybrid.py` | Hybrid 1-D toy: shells are fluid where the local packing fraction exceeds a threshold and billiard balls elsewhere |
| `jm_spots.py` | Same tessellation: can eight compact extremes be eight neighbouring walls? (no; `jm_spots_results.txt`) |
| `jm_sizes.py` | Cell-size distribution of the pool (small bangs, observer weighting); ten seeds pooled in `jm_sizes_results.txt` |
| `rings_contact.py` | Energy-conserving rebound toy (Hénon shells, elastic contact); grid results in `contact_grid_results.txt` |
| `spot_population.py`, `run_spots.sh` | Spot population as a clock: width–amplitude track vs Gaussian mocks (`spot_population_results.txt`) |
| `slope_axes.py`, `fit_slope.py` | Slope / Axis-of-Evil axis comparison and joint direction fit (Sec. 7) |
| `contact_annihilation.py`, `grav_baryo_sign.py` | Baryogenesis checks: wall-annihilation percolation; sign of dR/dt after the bounce (Sec. 13, 4b) |
| `frag_ps.py`, `frag_ps2.py`, `frag_blue.py` | Fragment-and-detonate toy (road (iii), Sec. 13): horizon-crossing Press–Schechter of the kernel, Planck-star clock, survivor/tail/injection bookkeeping; results in `frag_*_results.txt`; referee reports in `notes/` |
| `contact_annihilation_long.py` | Relic/anti-relic annihilation with many reshuffles: survivor law ~1/rounds; `contact_annihilation_long_results.txt`, `notes/annihilation_radiation.md` |
| `pool_spectrum.py` | Energy–frequency spectrum of the universe's mass (3-D, projected on a plane, one line cut; Planck lensing convergence as measured projection) against ocean spectra; `figures/pool_spectrum.png`, `pool_spectrum_results.txt` |
| `hypocentre.py` | Hypocentre inversion for H5: joint posterior over seed direction and distance from the four axis tests plus flow/dipole/quadrupole amplitudes; `hypocentre_results.txt`, `notes/hypocentre.md` |
| `shore_dipole.py` | Shore test: cluster-count dipole along the fitted slope axis in PSZ2, ACT DR5, eRASS1, SPT-SZ with footprint correction, shuffles and kinematic subtraction; `shore_dipole_results.txt`, `notes/shore_dipole.md` |
| `pool_spectrum_3d.py` | 3-D spectrum: P(k,z) growing with time vs JONSWAP developing with fetch; `figures/pool_spectrum_3d.png` |
| `extreme_events.py` | Extreme-value test: most extreme CMB peak and peak counts vs Gaussian skies with the map's own spectrum; nucleation action as a sigma excursion; `extreme_events_results.txt` |
| `pool_spectrum_shells.py` | Projected-mass spectrum of thin shells at increasing radius (Limber): the peak walks from l~2 nearby to l~160 at last scattering, amplitude falls as D(z)^2; `figures/pool_spectrum_shells.png` |
| `cascade_decay.py` | Decay law of the kernel's turbulence: turnovers to recombination, residual velocity under t^-10/7 and t^-6/5, what free decay would need; `cascade_decay_results.txt` |
| `wave_components.py` | Wave systems in the matter spectrum: smooth growing mode, BAO acoustic component, isocurvature and tensor bounds, with energy share and epoch; `figures/wave_components.png` |
| `age_from_waves.py` | Age of the universe from the waves alone: spectrum peak (Omega_m h^2) + growth rates f sigma_8(z) + lensing amplitude S8, no CMB peaks, no distance ladder; `age_from_waves_results.txt` |
| `symbolic_fits.py`, `jm_cells_dump.py`, `jm_sizes_dump.py` | Symbolic regression on the contact-toy grid and the Johnson–Mehl samples (closed forms for the unbound fraction, the exact observer-offset density, cell-volume gamma law); `symbolic_fits_results.txt`, `notes/symbolic_fits.md` |
| `void_dipole.py` | Void-fraction dipole compass for H5 in WISE×SuperCOSMOS (2-D) and 2M++ (3-D) |
| `core_search.py` | Pre-registered dark-lens search for the Kernel 0,0,0 remnant in the DES Y3 convergence map vs. cluster/group/SZ/X-ray catalogues; `core_search_results.txt`, `notes/core_search.md` |
| `ring_echo.py` | Pre-registered search for a second ringing in the galaxy P0(k) (eBOSS DR16 LRG + BOSS DR12), distinct from the BAO: scan of scales 20–600 Mpc/h with covariance, window, look-elsewhere mocks; `ring_echo_results.txt`, `notes/ring_echo.md`, `figures/ring_echo.png` |
| `energy_budget.py` | Bounce-to-equality energy budget: relic kinetic energy, Kernel 0,0,0 conversion radiation, latent heat; ejecta momentum at equality |

Figures: `bounce.png`, `ejecta_kinematics.png`, `sn_hubble.png`, `wz_effective.png`, `ltb_wz.png`, `coldspot_profile.png`.

## Data (not in the repo, download links)

Large datasets are excluded via `.gitignore`. Place them in the repo root with the file names below to rerun the scripts.

| File(s) used | Source |
|---|---|
| `smica.fits` (Planck 2018 SMICA I/Q/U + masks, 2 GB) | Planck Legacy Archive: https://pla.esac.esa.int/ (COM_CMB_IQU-smica_2048_R3.00_full.fits) |
| `PR4_lensing_maps.tar`, `PR4_variations/` (PR4 lensing κ alms, mean field, mask) | https://github.com/carronj/planck_PR4_lensing |
| `pantheon.dat`, `STAT+SYS.npz` (Pantheon+) | https://github.com/PantheonPlusSH0ES/DataRelease |
| `DES-Dovekie_HD.csv` (DES-SN5YR Hubble diagram) | https://github.com/des-science/DES-SN5YR |
| `twompp_density.npy`, `twompp_velocity.npy` (2M++ density and velocity fields, Carrick et al. 2015) | https://cosmicflows.iap.fr/ (2M++ reconstruction) |
| `quaia_G20.0.fits`, `quaia_selfunc_ns64.fits` (Quaia quasar catalogue + selection function) | https://zenodo.org/records/8060755 |
| `wisescos_counts.csv` (WISE×SuperCOSMOS photo-z counts, z 0.1–0.4) | SSA server-side SQL: http://ssa.roe.ac.uk/WISExSCOS.html |
| `spincat/` (Shamir SDSS galaxy spin catalogues) | https://people.cs.ksu.edu/~lshamir/data/ |
| `table2.dat`, `table3.dat` (CosmicFlows-4 groups/velocities) | https://edd.ifa.hawaii.edu/ (CF4) |
| `data/KS_full.fits`, `data/wiener_full.fits`, `data/glimpse_mask.fits` (DES Y3 convergence maps + mask, Jeffrey et al. 2021, nside 1024) | https://desdr-server.ncsa.illinois.edu/despublic/y3a2_files/massmaps/ |
| `data/y3_redmapper_v6.4.22+2_release.h5` (DES Y3 redMaPPer, 344 MB) | https://desdr-server.ncsa.illinois.edu/despublic/y3a2_files/y3kp_clusters/data/y3_redmapper_v6.4.22+2_release.h5 |
| `data/redmapper_y1a1_public_v6.4_catalog.fits.gz` (DES Y1 redMaPPer) | https://desdr-server.ncsa.illinois.edu/despublic/y1a1_files/redmapper/ |
| `data/HFI_PCCS_SZ-union_R2.08.fits` (Planck PSZ2 union) | https://cdsarc.cds.unistra.fr/ftp/J/A+A/594/A27/fits/ |
| `data/DR5_cluster-catalog_v1.1.fits` (ACT DR5 SZ clusters, Hilton et al. 2021) | https://lambda.gsfc.nasa.gov/data/suborbital/ACT/ACT_dr5/DR5_cluster-catalog_v1.1.fits |
| `data/spt_sz_bocquet2019_table5.dat` (+ReadMe) (SPT-SZ 2500d, Bocquet et al. 2019) | https://cdsarc.cds.unistra.fr/ftp/J/ApJ/878/55/ |
| `data/erass1cl_primary_v3.2.fits` (eROSITA eRASS1 clusters, Bulbul et al. 2024) | https://erosita.mpe.mpg.de/dr1/AllSkySurveyData_dr1/Catalogues_dr1/BulbulE_DR1/erass1cl_primary_v3.2.fits.tgz |
| `data/mcxc_piffaretti2011.dat` (+ReadMe) (MCXC X-ray clusters) | https://cdsarc.cds.unistra.fr/ftp/J/A+A/534/A109/ |
| `data/tully2015_table3.dat` (+ReadMe) (2MRS galaxy groups, Tully 2015) | https://cdsarc.cds.unistra.fr/ftp/J/AJ/149/171/ |
| `data/ring_echo/Data_LRGPk_{NGC,SGC}_0.6z1.0_prerecon.txt`, `Covariance_LRGPk_*`, `Window_LRGPk_*`, `README_LRG-QSOPk.txt` (eBOSS DR16 LRG P(k) multipoles, Gil-Marín et al. 2020) | https://svn.sdss.org/public/data/eboss/DR16cosmo/tags/v1_0_1/dataveccov/lrg_elg_qso/LRG_Pk/ |
| `data/ring_echo/Beutler_fs.tar.gz` → `beutler_fs/public_material_RSD/` (BOSS DR12 P(k) multipoles, Patchy covariances, RR windows, Beutler et al. 2017) | https://data.sdss.org/sas/dr12/boss/papers/clustering/Beutler_etal_DR12COMBINED_fullshape_powspec.tar.gz |

Derived intermediate maps (`*_ns64.npy`, `*_ns256.npy`, `curl_omega_map.npy`) are regenerated by the scripts.

## Licence

Text and figures CC BY 4.0; code MIT.
