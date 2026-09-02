"""Extreme large-scale CMB spots from SMICA (Nside 128): per-spot (DT, theta)
with theta = half-max radius of the radial profile around each peak.
Locus test: common-distance-shell bangs predict DT ~ theta^3 (slope +3);
BBKS Gaussian peaks predict weak amplitude-size correlation.
Null: 60 Gaussian realizations from the map's own pseudo-Cl, same pipeline."""
import numpy as np, healpy as hp

NSIDE = 128
m = np.load("/home/ubuntu/big-pool-bang/smica_nside128.npy").astype(float)
m -= m.mean()
if np.std(m) < 1e-2:      # map in K -> uK
    m *= 1e6
npix = hp.nside2npix(NSIDE)
l, b = hp.pix2ang(NSIDE, np.arange(npix), lonlat=True)
mask = np.abs(b) > 25
vecs = np.array(hp.pix2vec(NSIDE, np.arange(npix))).T

def catalog(mp, nspots=8):
    sm = hp.smoothing(mp, fwhm=np.radians(2.0))
    order = np.argsort(-np.abs(sm))
    used = np.zeros(npix, bool)
    cat = []
    for p in order:
        if not mask[p] or used[p]:
            continue
        d = np.degrees(np.arccos(np.clip(vecs @ vecs[p], -1, 1)))
        # radial profile in 1-deg rings out to 25 deg
        prof = np.array([sm[(d >= i) & (d < i + 1)].mean() for i in range(25)])
        dt = prof[0]
        half = np.where(np.sign(dt) * prof < 0.5 * abs(dt))[0]
        th = float(half[0]) + 0.5 if len(half) else 25.0
        cat.append((l[p], b[p], dt, th))
        used |= d < max(2.5 * th, 5.0)
        if len(cat) >= nspots:
            break
    return cat

def slope_r(cat):
    x = np.log([c[3] for c in cat])
    y = np.log([abs(c[2]) for c in cat])
    if np.std(x) < 1e-6:
        return np.nan, np.nan
    s = np.polyfit(x, y, 1)[0]
    return s, np.corrcoef(x, y)[0, 1]

cat = catalog(m)
print("DATA spots (l, b, DT_uK, theta_halfmax_deg):")
for c in cat:
    M = (abs(c[2]) / c[3]) ** 1.5
    D = abs(c[2]) ** 0.5 * c[3] ** -1.5
    print("  (%6.1f,%+6.1f)  DT=%+7.1f uK  th=%4.1f deg   M~%7.0f  D~%6.2f" %
          (c[0], c[1], c[2], c[3], M, D))
s, r = slope_r(cat)
print("DATA: slope dlnDT/dlnTheta = %+.2f   corr r = %+.2f" % (s, r))
print("  (locus on a common distance shell: +3; Gaussian peaks: ~0)")

cl = hp.anafast(m, lmax=3 * NSIDE - 1)
slopes, rs = [], []
for i in range(60):
    np.random.seed(i)
    g = hp.synfast(cl, NSIDE)
    sg, rg = slope_r(catalog(g))
    if np.isfinite(sg):
        slopes.append(sg); rs.append(rg)
slopes, rs = np.array(slopes), np.array(rs)
print("GAUSSIAN NULL (%d sims): slope = %+.2f +- %.2f   r = %+.2f +- %.2f" %
      (len(slopes), slopes.mean(), slopes.std(), rs.mean(), rs.std()))
print("percentile of data slope in null: %.1f%%" % (100 * np.mean(slopes < s)))
