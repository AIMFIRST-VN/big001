"""P6 (Spot Population) sub-test: population as a clock. Take the K most extreme peaks of the smoothed SMICA map
(|b| > BCUT), measure each peak's amplitude and half-max radius, and test the size-amplitude relation
(Spearman rho of half-width vs |dT|, and the slope of log width on |dT|) against Gaussian mocks with the map's own
pseudo-spectrum and the same selection. A decaying-vortex population would put the real spots on a track that the
Gaussian population does not follow. Usage: python3 spot_population.py <nmocks> <seed> [K=50] [fwhm_deg=1.0]
Prints one line per mock (rho, slope) and, for seed 0 only, the data values."""
import numpy as np, healpy as hp, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
nm = int(sys.argv[1]); seed = int(sys.argv[2]); K = int(sys.argv[3]) if len(sys.argv) > 3 else 50
fwhm = float(sys.argv[4]) if len(sys.argv) > 4 else 1.0
BCUT = 20.0; NS = 128; npix = hp.nside2npix(NS); LMAX = 3 * NS - 1
th, ph = hp.pix2ang(NS, np.arange(npix)); b = 90 - np.degrees(th)
good = np.abs(b) > BCUT
neigh = hp.get_all_neighbours(NS, np.arange(npix))          # (8, npix)
def analyse(m):
    ms = hp.smoothing(m, fwhm=np.radians(fwhm), verbose=False) if fwhm > 0 else m
    a = np.abs(ms); nb = np.where(neigh >= 0, a[np.maximum(neigh, 0)], -np.inf)
    peak = good & (a > nb.max(0))
    idx = np.where(peak)[0]; idx = idx[np.argsort(a[idx])[::-1][:K]]
    amps, widths = [], []
    for i in idx:
        v = hp.pix2vec(NS, i); half = a[i] / 2
        for r in np.arange(0.5, 15.0, 0.25):
            ring = hp.query_disc(NS, v, np.radians(r)); ring = ring[np.linalg.norm(np.array(hp.pix2vec(NS, ring)).T - v, axis=1) > np.sin(np.radians(r - 0.25))]
            prof = np.mean(ms[ring] * np.sign(ms[i]))
            if prof < half: break
        amps.append(a[i]); widths.append(r)
    amps, widths = np.array(amps), np.array(widths)
    from scipy.stats import spearmanr
    rho = spearmanr(widths, amps)[0]; slope = np.polyfit(amps, np.log(widths), 1)[0]
    return rho, slope, amps, widths
m = np.load(os.path.join(HERE, 'smica_nside128.npy'))
if np.std(m) < 1e-2: m *= 1e6
if seed == 0:
    rho, slope, A, W = analyse(m)
    print(f"DATA K={K} fwhm={fwhm}: rho={rho:+.3f} slope={slope:+.5f}  |dT| range {A.min():.0f}-{A.max():.0f} uK, width range {W.min():.2f}-{W.max():.2f} deg, median width {np.median(W):.2f}")
cl = hp.anafast(np.where(good, m - m[good].mean(), 0.0), lmax=LMAX) / good.mean()
rng = np.random.default_rng(seed)
for k in range(nm):
    np.random.seed(rng.integers(1 << 31))
    mm = hp.synfast(cl, NS, lmax=LMAX, verbose=False)
    rho, slope, A, W = analyse(mm)
    print(f"MOCK seed={seed} k={k}: rho={rho:+.3f} slope={slope:+.5f} medwidth={np.median(W):.2f}")
