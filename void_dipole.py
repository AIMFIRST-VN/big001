"""Void-metric compass for H5 (Our Position): dipole of void fraction in (a) WISE x SuperCOSMOS
2-D density (z 0.1-0.4, 10-degree high-passed, Nside 64) with a rotation null at three latitude cuts,
and (b) the 2M++ 3-D density field (Carrick et al. 2015) within 100-200 Mpc/h: centroid of void cells
(delta < -0.5, -0.7) and hemispheric void-fraction asymmetry along the centroid axis vs random axes.
Results (2026-09-02): (a) significant only at |b|>30 (p<0.01), gone at |b|>45 (p=0.28) and |b|>60 (p=0.43):
a low-latitude systematic, null. (b) centroid offsets 10-20 Mpc/h toward (l,b) ~ (295-322, -15 to -44),
asymmetry 1.5-2 sigma against random axes before any look-elsewhere correction: not a detection."""
import numpy as np, healpy as hp, os
HERE = os.path.dirname(os.path.abspath(__file__))
d = np.load(os.path.join(HERE, 'wisescos_delta_ns64.npy')); ns = 64
lp, bp = hp.pix2ang(ns, np.arange(hp.nside2npix(ns)), lonlat=True)
vec = np.array(hp.pix2vec(ns, np.arange(hp.nside2npix(ns)))).T
def dipole(mapv, mask):
    A = np.column_stack([np.ones(mask.sum()), vec[mask]]); x, *_ = np.linalg.lstsq(A, mapv[mask], rcond=None)
    dv = x[1:]; amp = np.linalg.norm(dv); return amp, np.degrees(np.arctan2(dv[1], dv[0])) % 360, np.degrees(np.arcsin(dv[2] / amp))
rng = np.random.default_rng(5)
for bcut in (30, 45, 60):
    good = np.isfinite(d) & (d != 0) & (np.abs(bp) > bcut); mp = (d < -0.3).astype(float)
    amp, l, b = dipole(mp, good); amps = []
    for _ in range(100):
        rot = hp.Rotator(rot=list(rng.uniform(0, 360, 3)), deg=True)
        mr = (rot.rotate_map_pixel(np.where(np.isfinite(d), d, 0.0)) < -0.3).astype(float); amps.append(dipole(mr, good)[0])
    amps = np.array(amps)
    print(f"WISExSCOS |b|>{bcut}: void-fraction dipole {amp:.4f} toward ({l:.0f},{b:.0f}); rotation null {amps.mean():.4f}+-{amps.std():.4f} p={np.mean(amps >= amp):.2f}")
den = np.load(os.path.join(HERE, 'twompp_density.npy')); ax = np.linspace(-200, 200, 257)
X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij'); R = np.sqrt(X**2 + Y**2 + Z**2)
for rmax in (100, 150, 200):
    sel = (R < rmax) & (R > 10); pos = np.stack([X[sel], Y[sel], Z[sel]], 1)
    for thr in (-0.5, -0.7):
        isv = den[sel] < thr; vx = pos[isv].mean(0); amp = np.linalg.norm(vx); u = vx / amp
        h = (pos @ u) > 0; a_c = isv[h].mean() - isv[~h].mean()
        dirs = rng.normal(size=(200, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
        asym = np.array([isv[(pos @ w) > 0].mean() - isv[(pos @ w) <= 0].mean() for w in dirs])
        print(f"2M++ R<{rmax} delta<{thr}: void centroid {amp:.1f} Mpc/h toward ({np.degrees(np.arctan2(u[1],u[0]))%360:.0f},{np.degrees(np.arcsin(u[2])):.0f}); asymmetry {a_c:+.3f} vs random-axis rms {asym.std():.3f}")
