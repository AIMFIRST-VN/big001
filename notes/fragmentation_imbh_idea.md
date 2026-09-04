# Candidate mechanism for open problem 4(a): bounce, fall back, fragment, detonate — and its black-hole tail

Context: The Big Pool Bangs Theory (big-pool-bangs.md). The universe is a rebound of a jammed ball of Planck-mass relics
(hard spheres, Compton = Schwarzschild at m_P). Outer ~88% free-streams away (collisionless; it is the dark matter as
individual relics). The bound interior (Kernel 0,0,0) is held at jamming density by contact pressure. The hot Big Bang
is the conversion of relic mass to radiation inside the kernel. Open problem 4(a): two-body contact mergers need a
cross-section 10^28–10^44 l_P^2 c to convert enough mass before equality (geometric l_P^2 freezes out at the bounce).

The idea (from the author, 2026-09-04):
1. The first bounce throws relics apart; the bound part falls back (recompression; the energy-conserving toy shows
   a jammed interior throughout).
2. Where fall-back re-jams at Planck density it fragments gravitationally into clumps. Any jammed clump of
   macroscopic mass is inside its Schwarzschild radius (Planck density: 1e17 kg has R = 2e-27 m vs r_s = 1.5e-10 m).
3. Under the loop bounce each clump is a Planck star (Rovelli & Vidotto 2014) that explodes after t ~ (M/m_P)^2 t_P,
   returning its whole mass as radiation into all species. Assembly is gravitational so the relic cross-section is
   irrelevant; conversion is complete; survivor fraction f = mass that did not fragment.
4. Clock: explode before BBN (1 s) needs M < 1e14 kg; before recombination M < 3e20 kg; before today M < 6e22 kg.
   Heavier fragments are permanent primordial black holes (PBHs) inside our universe.
5. Fragment mass function is NOT derived. At Planck density with c_s ~ c the Jeans mass is ~m_P, so gravity is
   unstable at all scales and Jeans does not set the size; dynamics of the fall-back must.
6. Heavy tail: if the mass function reaches 10–100 M_sun at f_PBH ~ 1e-3 of DM it could supply part of the LIGO rate;
   if it reaches 1e2–1e5 M_sun (intermediate-mass, IMBH) at the allowed f_PBH ~ 1e-4–1e-3 it puts ~1e4–1e5 holes of
   1e4 M_sun per Milky-Way halo, i.e. one per dwarf/globular, gives heavy seeds for z>7 quasars and JWST overmassive
   holes, and predicts dwarf-nucleus occupation near unity, LISA IMBH merger rate rising past z=10.
7. Later clumping: hierarchical mergers of stellar-mass PBHs stall at a few hundred M_sun (GW recoil ejects products
   from small clusters); Eddington accretion from 1e2–1e4 M_sun seeds reaches 1e9 M_sun by z~7.
8. Claimed status: LIGO and IMBH are constraints on the (unknown) mass function, not successes; sub-solar mergers or a
   merger rate rising at high z would favour it.

## Toy results (2026-09-04; frag_ps.py, frag_ps2.py, frag_blue.py; results *.txt alongside)
Horizon-crossing Press-Schechter of the re-jammed kernel: M_H = m_P t/t_P, collapse when smoothed contrast > delta_c,
sequential depletion dU/dlnM = -U beta, beta = erfc(dc/(sqrt2 sigma)). Planck-star clock sorts fragments by mass.
Bookkeeping: survivors U(M_L) = f are the dark matter (Sec 5), so f_PBH = tail(>6e22 kg)/f; radiation at time t is what
exploded before t, so injection in each era is (mass exploding in era)/(mass exploded before it); limits: BBN hadronic
1e-6, late BBN 1e-4, mu 1.2e-4, y 6e-5, post-recombination 1e-7 (rough).
- Kolmogorov (sigma rising with M, the framework's own cascade): NO fragmentation below the kernel scale for any A,
  delta_c; collapse only at the outer scale. The mechanism does not run on the paper's spectrum.
- Flat (n~1) and shot-noise+flat: every configuration killed. Structural reason: tail/f = U(6e22)/f - 1 = exp(int beta
  dlnM above 6e22 kg) - 1, so the permanent tail exceeds the survivors unless collapse stops above 6e22 kg; and
  f in 1e-28..1e-13 needs beta ~ 0.3 per e-fold below 1e14 kg. Pure shot noise (A=0) gives f = 0.2-0.8: no radiation era.
- Blue power law sigma = (dc/2)(M/M_piv)^(-q), i.e. n = 1+6q: ALIVE for q = 0.1-0.5 with M_piv = 1e12-1e13 kg
  (q=0.1 needs M_piv >= 1e12; q >= 0.3 also at 1e8-1e10): f = 1e-14..1e-20, tail/f < 1e-80, all injection limits met
  by many orders. The surviving version predicts ZERO primordial black holes above 6e22 kg from the kernel; the LIGO
  and IMBH tails are gone. Extrapolating q=0.1, pivot 1e12 kg to the kernel scale gives sigma(M_L) ~ 2e-5 (CMB-like
  amplitude) but with n = 1.6, so the spectrum must break to the observed n_s = 0.965 at large scales.
- Consequence: 'bounce, fall back, fragment, detonate' is a live candidate for 4(a) only with a blue small-scale contrast
  spectrum pivoting near 1e12 kg, and then it forbids kernel-made IMBHs; a confirmed primordial IMBH or sub-solar merger
  would falsify the channel, not support it.
