"""Approximate TT quadratic-estimator curl (omega) reconstruction on SMICA.
Relative test only: curl-map variance stacked at the 8 extreme spots vs random
positions, mean-field subtracted via Gaussian sims through the same pipeline.
Normalization cancels in the spot-vs-null comparison."""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np, healpy as hp

NS, LMAX = 256, 512
T = hp.read_map(os.path.join(HERE, "smica.fits"), field=0)
TM = hp.read_map(os.path.join(HERE, "smica.fits"), field=3)
T = hp.ud_grade(T, NS); TM = hp.ud_grade(TM, NS) > 0.9
if np.std(T[TM]) < 1e-2: T *= 1e6
T = T - np.mean(T[TM])
Tm = T * TM

ctot = hp.anafast(Tm, lmax=LMAX) / np.mean(TM)   # pseudo-Cl as total power filter
ctot[ctot <= 0] = np.min(ctot[ctot > 0])
cl_th = ctot.copy()                              # use same for signal filter (Wiener-ish)
finv = np.zeros(LMAX + 1); finv[2:] = 1.0 / ctot[2:]
fwei = np.zeros(LMAX + 1); fwei[2:] = cl_th[2:] / ctot[2:]

def curl_map(tmap):
    """TT QE: product of inverse-variance-filtered T and gradient of Wiener T;
    curl part of the resulting vector field."""
    alm = hp.map2alm(tmap * TM, lmax=LMAX)
    t1 = hp.alm2map(hp.almxfl(alm, finv), NS)
    _, dth, dph = hp.alm2map_der1(hp.almxfl(alm, fwei), NS)
    # vector field V = t1 * grad(t2); decompose spin-1: E=gradient(phi), B=curl(omega)
    almE, almB = hp.map2alm_spin([t1 * dth, t1 * dph], 1, lmax=LMAX)
    w = hp.alm2map(almB, NS)                     # curl-type potential map (unnormalized)
    return hp.smoothing(w, fwhm=np.radians(2.0))

print("building data curl map...", flush=True)
wd = curl_map(T)

print("mean field from 40 sims...", flush=True)
rng_seeds = range(100, 140)
mf = np.zeros(hp.nside2npix(NS))
for i in rng_seeds:
    np.random.seed(i)
    mf += curl_map(hp.synfast(ctot, NS))
mf /= len(list(rng_seeds))
wd = wd - mf
np.save(os.path.join(HERE, "curl_omega_map.npy"), wd)

npix = hp.nside2npix(NS)
lp, bp = hp.pix2ang(NS, np.arange(npix), lonlat=True)
vecs = np.array(hp.pix2vec(NS, np.arange(npix))).T
spots = [(157.1,-70.5),(79.5,-33.2),(203.2,-56.3),(155.0,-29.3),
         (304.5,-29.0),(210.6,-35.0),(170.3,-46.6),(184.2,-54.3)]
def discvar(mp, l, b, r=6.0):
    lr, br = np.radians(l), np.radians(b)
    n = np.array([np.cos(br)*np.cos(lr), np.cos(br)*np.sin(lr), np.sin(br)])
    d = np.degrees(np.arccos(np.clip(vecs @ n, -1, 1)))
    s = (d < r) & TM
    return np.mean(mp[s]**2) if s.sum() > 100 else np.nan

print("per-spot curl variance (arb units):")
vals = [discvar(wd, l, b) for l, b in spots]
for (l, b), v in zip(spots, vals):
    print("  (%6.1f,%+6.1f)  varOmega=%10.4g%s" % (l, b, v, "   <-- ColdSpot" if abs(l-203.2)<1 else ""))
stat = np.nanmean(vals)
cs = vals[2]

print("null: 60 sims through identical pipeline (mean-field subtracted)...", flush=True)
stats_n, cs_n = [], []
for i in range(200, 260):
    np.random.seed(i)
    wn = curl_map(hp.synfast(ctot, NS)) - mf
    vv = [discvar(wn, l, b) for l, b in spots]
    stats_n.append(np.nanmean(vv)); cs_n.append(vv[2])
stats_n, cs_n = np.array(stats_n), np.array(cs_n)
print("8-spot stack: data %.4g  null %.4g +- %.4g  percentile %.1f%%" %
      (stat, stats_n.mean(), stats_n.std(), 100*np.mean(stats_n < stat)))
print("ColdSpot:     data %.4g  null %.4g +- %.4g  percentile %.1f%%" %
      (cs, cs_n.mean(), cs_n.std(), 100*np.mean(cs_n < cs)))
