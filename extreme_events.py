"""Frequencies and extreme events: for a Gaussian field the spectrum fixes the statistics of its extremes
(Rice 1944; Longuet-Higgins 1957; BBKS 1986). Test on the CMB: take the SMICA map's own pseudo-spectrum, make
Gaussian skies with the same spectrum, mask and smoothing, and compare the distribution of the most extreme peak
(in sigma units) and the number of |nu| > 3.5 peaks with the data. Also: the ocean's Rayleigh law for wave heights
(rogue wave H > 2 H_s: p = exp(-8) = 3.4e-4 per wave) and the nucleation action as a sigma-excursion: S = nu^2/2.
Usage: python3 extreme_events.py [nmock]"""
import numpy as np, healpy as hp, os, sys
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(1)
nmock = int(sys.argv[1]) if len(sys.argv) > 1 else 300
m = np.load(os.path.join(HERE, 'smica_nside128.npy')); nside = 128; npix = hp.nside2npix(nside)
theta, phi = hp.pix2ang(nside, np.arange(npix)); b = 90 - np.degrees(theta); mask = np.abs(b) > 20
fwhm = np.radians(1.0)
def peaks(x):
    """local maxima of |x| among the 8 neighbours, masked; return sigma-units of |x| at peaks"""
    xm = np.where(mask, x, 0.0); sig = xm[mask].std()
    nb = hp.get_all_neighbours(nside, np.arange(npix)); nbv = np.where(nb >= 0, np.abs(xm)[np.clip(nb, 0, npix-1)], -np.inf)
    ismax = (np.abs(xm) > nbv.max(axis=0)) & mask
    return np.abs(xm)[ismax] / sig
ms = hp.smoothing(m, fwhm=fwhm, verbose=False); cl = hp.anafast(np.where(mask, m - m[mask].mean(), 0), lmax=3*nside-1) / mask.mean()
pd = peaks(ms); dmax = pd.max(); dn35 = (pd > 3.5).sum(); dn4 = (pd > 4).sum()
mx, n35, n4 = [], [], []
for i in range(nmock):
    g = hp.smoothing(hp.synfast(cl, nside, verbose=False), fwhm=fwhm, verbose=False); pg = peaks(g)
    mx.append(pg.max()); n35.append((pg > 3.5).sum()); n4.append((pg > 4).sum())
mx, n35, n4 = map(np.array, (mx, n35, n4))
print(f"SMICA, 1 deg smoothing, |b|>20, nside {nside}: most extreme peak = {dmax:.2f} sigma; peaks >3.5 sigma: {dn35}; >4 sigma: {dn4}")
print(f"Gaussian mocks ({nmock}, same spectrum/mask/smoothing): max peak median {np.median(mx):.2f} (5-95%: {np.percentile(mx,5):.2f}-{np.percentile(mx,95):.2f}); P(max >= data) = {np.mean(mx >= dmax):.3f}")
print(f"  peaks >3.5 sigma: median {np.median(n35):.0f} (5-95%: {np.percentile(n35,5):.0f}-{np.percentile(n35,95):.0f}), P(>= data) = {np.mean(n35 >= dn35):.3f};  >4 sigma: median {np.median(n4):.0f}, P(>= data) = {np.mean(n4 >= dn4):.3f}")
# Gumbel fit to the mock maxima: location/scale
from scipy.stats import gumbel_r
loc, sc = gumbel_r.fit(mx); print(f"  Gumbel fit to mock maxima: mode {loc:.2f} sigma, scale {sc:.2f}; the spectrum (through its moments) fixes these two numbers")
print("Ocean (Rayleigh, narrow band): P(H > 2 H_s) = exp(-8) = 3.4e-4 per wave, one in 3000; observed freak-wave rates exceed this by 2-10x in steep seas (Benjamin-Feir instability; Dysthe et al. 2008)")
for S in (184, 280, 560): print(f"Nucleation action S = {S}: as a Gaussian excursion nu = sqrt(2S) = {np.sqrt(2*S):.1f} sigma")
