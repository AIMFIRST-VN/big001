# The Big Pool Bangs Theory

Speculative cosmology paper: a jammed hard-sphere Planck-relic core, a loop-quantum-cosmology-type bounce, chaotic ejection as the origin of dark matter, and a set of hypotheses about where we sit inside our own Big Pool Bang and where a second one might be (the Cold Spot direction). Every hypothesis carries explicit falsifiers and a pre-registered test list.

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
| `void_dipole.py` | Void-fraction dipole compass for H5 in WISE×SuperCOSMOS (2-D) and 2M++ (3-D) |
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

Derived intermediate maps (`*_ns64.npy`, `*_ns256.npy`, `curl_omega_map.npy`) are regenerated by the scripts.

## Licence

Text and figures CC BY 4.0; code MIT.
