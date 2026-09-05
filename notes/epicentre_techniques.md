# Earthquake-location techniques and their analogues for H5 (Our Position – Where Is Nemo) (saved 2026-09-05)

## Classical
| Seismic tool | What it does | Our analogue | Status |
|---|---|---|---|
| S–P triangulation | S trails P by a distance-proportional lag; circles from >=3 stations intersect at the epicentre | two signals with different speeds: sound (c/sqrt3) in the plasma vs gravity/light (c); no transient in our case, but lookback depth supplies time | not used |
| Wadati diagram | plot S–P lag vs P arrival; the intercept gives the ORIGIN TIME without knowing the location | extrapolate depth-dependent signals (flow, cluster dipole, CMB) back to the rebound's origin time | not done |
| Single-station back-azimuth | P-wave particle motion points along the ray | P1 (Polarization Ring): radial polarization around a source gives its direction | proposed test |
| First-motion polarity (beach ball) | up/down first motions map compression/dilatation quadrants; explosion vs double-couple | isotropic rebound = explosion source (all first motions outward); displaced observer sees a dipole in the flow; hot/cold extremes as quadrants | P4 done |
| Geiger's method / grid search | guess hypocentre, predict arrivals, least-squares iterate | fit_slope.py joint fit of four axes (direction only) | done |
| Hypocentre inversion with distance | amplitude/attenuation and lag fix distance | void-dipole amplitude ~ d_LS/D, flow-vs-depth ~ d_LS/D, GZ quadrupole ~ (d_LS/D)^2: fit direction AND distance; compare with 0.7 L of Sec 11f | not done (proposed) |
| Double-difference relocation | relative arrival times between event pairs cancel path structure | relative geometry of the eight extreme spots instead of absolute positions | not done |
| Array back-projection / beamforming | delay-and-stack toward trial source points | stack CMB / mass maps along trial directions with the expected profile | not done |
| Macroseismic intensity map | contour damage; maximum = epicentre | void compass (void abundance as intensity) | done, hint |
| Template matching / matched filter | cross-correlate with known waveforms to detect weak repeats | P6 (Spot Population) sub-threshold matched filter for Cold-Spot-class spots | proposed |
| Coda / scattering | late scattered energy constrains medium | P9 (Low-ℓ Tail) | pre-registered |
| Moment-tensor inversion | full source mechanism from waveforms | monopole/dipole/quadrupole decomposition of the flow and count fields (explosion + displacement + tide) | partial |
| Ambient-noise interferometry | cross-correlate noise between stations to get Green's functions without a source | cross-correlate CMB/lensing/galaxy fields to extract the transfer function of the rebound without assuming a source | not done |
| Probabilistic location (NonLinLoc, Bayesian) | full posterior over hypocentre, not a point | posterior over (n, D) with all four tests and the JM prior (offset 0.67 L, 0.35–1.05) | not done |

## Machine-learning / AI in seismology (and what they suggest)
- Deep-learning phase pickers (PhaseNet, Zhu & Beroza 2019; EQTransformer, Mousavi et al. 2020; GPD): pick arrivals from raw waveforms -> analogue: learned spot/feature detectors on CMB maps trained on Gaussian mocks (only as detectors, with pre-registered thresholds).
- Direct location networks (ConvNetQuake, Perol et al. 2018; graph neural nets on arrays, e.g. GENIE, McBrearty & Beroza 2023): map waveforms straight to location -> analogue: a network trained on simulated displaced-observer skies to regress (n, D); risky (trained on our own model), useful only as a fast likelihood emulator.
- Aftershock-pattern learning (DeVries et al. 2018; and its critique, Mignan & Broccardo 2019: a one-neuron model did as well): a caution about ML-found "relations".
- Symbolic regression / AI-discovered formulae: AI Feynman (Udrescu & Tegmark 2020), PySR (Cranmer 2023), Eureqa: fit compact closed-form relations to data. Uses here: (i) rings_contact.py grid -> unbound fraction as a closed form in (f, sigma); (ii) jm_cells.py / jm_sizes.py outputs -> closed forms for nearest-wall, offset and size distributions vs (L, S); (iii) the energy-budget relation f = T_eq/T_conv checked against the depletion integrals; (iv) NOT for the observed axes (four numbers, any formula fits).
- Physics-informed neural nets (Raissi et al. 2019) for wave propagation in the rings: an emulator for the acoustic structure of Sec 3.1.
- Empirical "pleasant" laws seismology found by data before theory: Gutenberg–Richter (b ~ 1), Omori (aftershock decay 1/t), Båth (largest aftershock 1.2 magnitudes below), Wadati. Our nearest: survivors ~ 1/rounds (annihilation), unbound = 1 - fbar^-3 (the retracted kick tautology: a warning), cell size L = e^{S/4} l_P.

## More techniques (added on request)
- Prompt elasto-gravity signals (PEGS; Vallée et al. 2017; Montagner et al. 2016): the gravity change of a rupture travels at c and arrives before the P wave; used for fast magnitude/location -> the sharpest analogue we have: gravity (potential, GZ quadrupole, ISW) versus sound (BAO, c/sqrt3) as two 'phases' of one source; a source-linked offset between the potential feature and the acoustic feature at last scattering is an S–P lag.
- ETAS / epidemic aftershock models (Ogata 1988): cascades of triggered events with Omori decay and Båth's law -> our nucleation cascade and Johnson–Mehl population: seeds triggering neighbours (the packet-triggered chain tested in Sec 11f) is an ETAS process; its branching ratio is the quantity to bound.
- Coulomb stress transfer (King et al. 1994): where one rupture loads the next -> which neighbouring cell is triggered by ours (the wall stress of a growing bubble on the false vacuum around it).
- Receiver functions and tomography (travel-time inversion for Earth structure; adjoint/full-waveform inversion): use many sources to image the medium -> use many spots/neighbours to image the tessellation (wall positions) rather than one source at a time.
- Rayleigh/Love surface-wave polarization and dispersion: the direction and the depth structure from a single station -> E/B decomposition (P1, P2) is our polarization; dispersion (scale-dependent speed) is the acoustic-peak structure.
- Hydroacoustic T-phase and infrasound arrays (CTBTO): locate sources in a medium with a different wave speed and a waveguide (SOFAR channel) -> propagation in the mantle ring with its own sound speed (Sec 3.1 acoustic structure) would carry a source's signal differently from the kernel; a waveguide is a ring.
- InSAR and GNSS geodesy (static displacement field of the rupture): the permanent deformation, not the waves -> the static displacement field of the rebound is the large-scale velocity/void field itself; its divergence and curl maps are the geodetic inversion (P4, P8).
- Tsunami inversion (arrival times at tide gauges give the source area): a slow, long-range messenger -> the kSZ / bulk-flow field as the slow messenger versus the CMB as the fast one; consistency of the two locations is a test.
- Distributed acoustic sensing (DAS on fibres): dense, cheap, low-quality sensors beat few good ones for location -> many low-significance probes (voids, dipoles, flows) combined beat one clean test; that is Fisher's method in Sec 7.
- Early-warning triage (ShakeAlert, JMA EEW): decide in seconds from the first arrivals -> nothing here; a reminder that the first-arriving signal (the largest-scale one) is the least precise.
