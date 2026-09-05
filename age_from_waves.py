"""An age of the universe from the wave distribution alone: no CMB peaks, no distance ladder.
Ingredient 1 (spectrum shape): the peak of P(k) is the equality scale, k_eq ~ 0.073 Omega_m h^2 /Mpc; galaxy-survey
  full-shape fits give Omega_m h^2 = 0.139 +- 0.006 (DESI DR1 full-shape + BBN baryon prior, ladder- and CMB-independent; verified,
  see notes/age_from_waves_data.md; BOSS-only full-shape fits span 0.136-0.164).
Ingredient 2 (growth with radius): f sigma_8(z) from redshift-space distortions; the fall of amplitude with distance
  fixes Omega_m through the growth factor D(z) (GR, LCDM expansion form). Compilation verified against the papers (6dFGS, SDSS MGS,
  BOSS DR12, eBOSS, WiggleZ, VIPERS, FastSound, DESI DR1) -- see notes/age_from_waves_data.md.
Then h = sqrt(Omega_m h^2 / Omega_m), and t0 = int da/(a H) with H = 100 h sqrt(Omega_m a^-3 + 1 - Omega_m).
Grid posterior over (Omega_m, sigma_8); age posterior by propagating Omega_m h^2 error. Usage: python3 age_from_waves.py"""
import numpy as np
from scipy.integrate import quad
fs8 = np.array([  # z, f sigma8, err  -- verified 2026-09-05, see notes/age_from_waves_data.md
 [0.067,0.423,0.055],   # 6dFGS, Beutler+12 (1204.4725)
 [0.15,0.49,0.145],     # SDSS MGS, Howlett+15 (1409.3238), +0.15/-0.14 symmetrised
 [0.38,0.497,0.046],[0.51,0.458,0.038],[0.61,0.436,0.035],   # BOSS DR12 consensus, Alam+17 (1607.03155) Tab.7, stat+sys in quadrature
 [0.70,0.473,0.041],[0.85,0.315,0.095],[1.48,0.462,0.045],   # eBOSS LRG/ELG/QSO, Alam+21 (2007.08991) Tab.3
 [0.44,0.413,0.080],[0.60,0.390,0.063],[0.73,0.437,0.072],   # WiggleZ, Blake+12 (1204.3674); correlated slices, overlap BOSS
 [0.60,0.55,0.12],[0.86,0.40,0.11],                          # VIPERS, Pezzotta+17 (1612.05645)
 [1.40,0.482,0.116]])                                        # FastSound, Okumura+16 (1511.08083)
USE_DESI = True   # DESI DR1 full-shape (2411.12021 Tab.9 ratio x Tab.11 fiducial f sigma_s8; approximate, overlaps BOSS)
desi = np.array([[0.295,0.397,0.090],[0.510,0.549,0.062],[0.706,0.479,0.047],[0.919,0.438,0.040],[1.317,0.373,0.034],[1.491,0.435,0.045]])
if USE_DESI: fs8 = np.vstack([fs8, desi])
omh2, omh2_err = 0.139, 0.006   # DESI DR1 FS+BAO+BBN (2411.12022 eq.3.1: Om=0.2962+-0.0095, H0=68.56+-0.75); BOSS full-shape spans 0.136-0.164 (1909.05277, 1909.05271, 2112.04515); all use a BBN omega_b prior
S8, S8_err = 0.790, 0.016   # weak-lensing amplitude of the mass map, DES Y3 + KiDS-1000 joint cosmic shear (2305.17173: 0.790 +0.018/-0.014); breaks the f sigma8 degeneracy
def growth(z, Om):
    a = 1/(1+z); E = lambda x: np.sqrt(Om/x**3+1-Om); return 2.5*Om*E(a)*quad(lambda x: 1/(x*E(x))**3, 1e-6, a)[0]
def f_of(z, Om):
    a = 1/(1+z); return (Om/(a**3*(Om/a**3+1-Om)))**0.55   # gamma = 0.55 (GR)
def age(Om, h):
    return quad(lambda a: 1/(a*100*h*np.sqrt(Om/a**3+1-Om)), 1e-8, 1)[0]*977.8   # Gyr (1/(km/s/Mpc) = 977.8 Gyr)
Oms = np.linspace(0.15, 0.5, 71); s8s = np.linspace(0.6, 1.0, 81)
chi2 = np.zeros((len(Oms), len(s8s)))
for i, Om in enumerate(Oms):
    D0 = growth(0, Om)
    pred = np.array([f_of(z, Om)*growth(z, Om)/D0 for z in fs8[:,0]])
    for j, s8 in enumerate(s8s):
        chi2[i, j] = np.sum(((fs8[:,1]-s8*pred)/fs8[:,2])**2) + ((s8*np.sqrt(Om/0.3)-S8)/S8_err)**2
L = np.exp(-0.5*(chi2-chi2.min())); pOm = L.sum(axis=1); pOm /= pOm.sum()
Om_mean = (Oms*pOm).sum(); Om_sd = np.sqrt(((Oms-Om_mean)**2*pOm).sum())
print(f"growth + lensing-amplitude fit (gamma=0.55): Omega_m = {Om_mean:.3f} +- {Om_sd:.3f}; sigma_8 = {s8s[np.unravel_index(chi2.argmin(), chi2.shape)[1]]:.2f}; chi2_min = {chi2.min():.1f} for {len(fs8)-1} dof")
# propagate: sample Om from pOm, omh2 from gaussian
rng = np.random.default_rng(0); n = 20000
Om_s = rng.choice(Oms, size=n, p=pOm) + rng.uniform(-0.0025, 0.0025, n); oh = rng.normal(omh2, omh2_err, n)
h_s = np.sqrt(oh/Om_s); t_s = np.array([age(o, hh) for o, hh in zip(Om_s[:4000], h_s[:4000])])
q = lambda a: np.percentile(a, [16, 50, 84])
print(f"h = {np.median(h_s):.3f} (16-84%: {q(h_s)[0]:.3f}-{q(h_s)[2]:.3f})")
print(f"AGE from the wave distribution: t0 = {q(t_s)[1]:.2f} Gyr (16-84%: {q(t_s)[0]:.2f}-{q(t_s)[2]:.2f})")
print(f"for reference: Omega_m=0.315, h=0.6736 gives {age(0.315,0.6736):.2f} Gyr; h=0.73 gives {age(0.315*(0.6736/0.73)**2,0.73):.2f} Gyr at fixed Omega_m h^2")
print("assumptions: GR growth (gamma=0.55), LCDM expansion form (w=-1), flat; data values verified 2026-09-05 (notes/age_from_waves_data.md); overlapping surveys treated as independent")
