"""Pre-registered archival 'dark lens' search (paper Sec. 3.1, 'A remnant?').

Question: is there a single compact mass concentration (>= 1e13 Msun) in a public
weak-lensing convergence map that has NO galaxy cluster / group / SZ / X-ray
counterpart?  The paper predicts at most one such object (the cold remnant of
Kernel 0,0,0) and a priori puts it within reach of a lensing survey with ~1% chance,
so a null is expected.

Thresholds fixed before looking at the map (see notes/core_search.md):
  smoothing FWHM      : 5, 10, 20 arcmin
  peak thresholds     : S/N >= 4 (list), S/N >= 5 (headline)
  noise               : robust (MAD) sigma of the smoothed kappa map in the survey
                        interior (mask fraction > 0.95); this INCLUDES signal, so the
                        S/N is conservative
  matching radius     : 5 arcmin + FWHM
  counterpart catalogs: DES Y3 redMaPPer (lambda>=5, full sample), DES Y1 redMaPPer,
                        Planck PSZ2 union, ACT DR5, SPT-SZ 2500d (Bocquet+19),
                        eROSITA eRASS1 primary clusters, MCXC, 2MRS groups (Tully 2015)
  edge flag           : smoothed mask fraction < 0.95 at the peak, or peak within
                        1 FWHM of the mask boundary

Map: DES Y3 Kaiser-Squires E-mode convergence (Jeffrey et al. 2021), HEALPix
nside 1024 RING, full-sample (4 tomographic bins combined), with the glimpse mask.
Data files live in data/ (see README data table for URLs).
"""
import os, sys, numpy as np, healpy as hp
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from astropy.cosmology import Planck18 as cosmo

HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, 'data')
FWHMS = [5.0, 10.0, 20.0]          # arcmin, fixed before looking
SN_LIST, SN_HEAD = 4.0, 5.0
MATCH_BASE = 5.0                   # arcmin, plus FWHM
NSIDE = 1024
out = open(os.path.join(HERE, 'core_search_results.txt'), 'w')
def log(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()

# ---------------------------------------------------------------- map
kap = hp.read_map(os.path.join(D, 'KS_full.fits'))
mask = hp.read_map(os.path.join(D, 'glimpse_mask.fits')).astype(float)
log(f"DES Y3 KS map: nside {NSIDE}, mask area {mask.sum()*hp.nside2pixarea(NSIDE, degrees=True):.0f} deg^2, "
    f"raw pixel std inside mask {np.std(kap[mask > 0]):.5f}")
kap = np.where(mask > 0, kap, 0.0)

# ---------------------------------------------------------------- catalogues
cats = {}
def add(name, ra, dec, z=None, extra=None):
    ra = np.asarray(ra, float); dec = np.asarray(dec, float); ok = np.isfinite(ra) & np.isfinite(dec)
    cats[name] = dict(ra=ra[ok], dec=dec[ok], z=(None if z is None else np.asarray(z, float)[ok]),
                      extra=(None if extra is None else np.asarray(extra)[ok]))
    log(f"  catalogue {name:22s} {ok.sum():7d} entries")
log("Counterpart catalogues:")
try:
    import h5py
    with h5py.File(os.path.join(D, 'y3_redmapper_v6.4.22+2_release.h5'), 'r') as f:
        g = None
        def find(n, o):
            global g
            if isinstance(o, h5py.Dataset) and n.lower().endswith('/ra') and g is None: g = n.rsplit('/', 1)[0]
        f.visititems(find)
        grp = f[g]; keys = list(grp.keys())
        lam = grp['lambda_chisq'][:] if 'lambda_chisq' in keys else grp['lambda'][:]
        zc = grp['z_lambda'][:] if 'z_lambda' in keys else np.full(len(lam), np.nan)
        add('DESY3_redMaPPer', grp['ra'][:], grp['dec'][:], zc, lam)
        log(f"    (group {g}, lambda range {lam.min():.1f}-{lam.max():.1f})")
except Exception as e:
    log(f"  DES Y3 redMaPPer NOT loaded: {e}")
t = fits.open(os.path.join(D, 'redmapper_y1a1_public_v6.4_catalog.fits.gz'))[1].data
add('DESY1_redMaPPer', t['RA'], t['DEC'], t['Z_LAMBDA'], t['LAMBDA'])
t = fits.open(os.path.join(D, 'HFI_PCCS_SZ-union_R2.08.fits'))[1].data
add('Planck_PSZ2', t['RA'], t['DEC'], t['REDSHIFT'], t['SNR'])
t = fits.open(os.path.join(D, 'DR5_cluster-catalog_v1.1.fits'))[1].data
add('ACT_DR5', t['RADeg'], t['decDeg'], t['redshift'], t['SNR'])
rows = [l for l in open(os.path.join(D, 'spt_sz_bocquet2019_table5.dat'))]
add('SPT-SZ', [float(l[22:30]) for l in rows], [float(l[31:39]) for l in rows],
    [float(l[69:74] or 'nan') for l in rows], [float(l[40:45]) for l in rows])
t = fits.open(os.path.join(D, 'erass1cl_primary_v3.2.fits'))[1].data
add('eRASS1', t['RA'], t['DEC'], t['BEST_Z'], t['EXT_LIKE'])
rows = [l for l in open(os.path.join(D, 'mcxc_piffaretti2011.dat'), encoding='latin-1')]
add('MCXC', [float(l[108:115]) for l in rows], [float(l[116:123]) for l in rows], [float(l[140:146]) for l in rows])
rows = [l for l in open(os.path.join(D, 'tully2015_table3.dat'))]
sgl = np.array([float(l[22:30]) for l in rows]); sgb = np.array([float(l[31:39]) for l in rows])
vcmb = np.array([float(l[46:51]) for l in rows]); mlum = np.array([float(l[74:83] or 'nan') for l in rows])
c = SkyCoord(sgl=sgl * u.deg, sgb=sgb * u.deg, frame='supergalactic').icrs
add('2MRS_groups', c.ra.deg, c.dec.deg, vcmb / 299792.458, mlum)

def vec(ra, dec): return hp.ang2vec(ra, dec, lonlat=True)
cat_vec = {k: vec(v['ra'], v['dec']) for k, v in cats.items()}
# coverage check: is each catalogue's footprint inside the DES mask? (informational)
for k, v in cats.items():
    pix = hp.ang2pix(NSIDE, v['ra'], v['dec'], lonlat=True); log(f"  {k:22s} entries inside DES mask: {int(mask[pix].sum())}")

# ---------------------------------------------------------------- 2M++ and CF4 (step 3)
den = np.load(os.path.join(HERE, 'twompp_density.npy')); vel = np.load(os.path.join(HERE, 'twompp_velocity.npy'))
def los_density(l, b, rmin=20, rmax=200, n=37):
    """Mean and min 2M++ galaxy density contrast along the line of sight (l,b) at 20-200 Mpc/h,
    and mean line-of-sight peculiar velocity (km/s, +ve = away from us, i.e. flowing toward the direction)."""
    r = np.linspace(rmin, rmax, n); u_ = np.array(hp.ang2vec(l, b, lonlat=True))
    idx = np.rint(np.outer(r, u_) / (400. / 256.) + 128).astype(int)
    ok = np.all((idx >= 0) & (idx <= 256), 1); idx = idx[ok]
    dd = den[idx[:, 0], idx[:, 1], idx[:, 2]]; vv = vel[:, idx[:, 0], idx[:, 1], idx[:, 2]]
    return dd.mean(), dd.min(), (vv.T @ u_).mean()
# CosmicFlows-4 group tables (table2/3.dat) have no column header in this checkout and a trial parse gave a
# direction-independent ~+800 km/s 'outflow' in every cone, i.e. a zero-point artefact; the per-group flow is
# therefore NOT used.  We report the 2M++ predicted LOS velocity and the angle to the published CF4 bulk-flow axis.
cf4 = None
log("CosmicFlows-4 per-group velocities: skipped (table3.dat column layout unverified); using 2M++ velocity field and the published CF4 bulk-flow axis (298,-8) instead.")

# ---------------------------------------------------------------- lensing mass scale
# Sigma_crit for the DES Y3 full source sample. Jeffrey+21 use the four MagLim/metacal bins;
# effective mean source redshift ~0.63.  We tabulate M(<theta_ap) for several lens redshifts.
ZS = 0.63
def sigma_crit(zl):  # Msun / Mpc^2 (physical)
    from astropy.constants import c as cc, G
    Ds = cosmo.angular_diameter_distance(ZS); Dl = cosmo.angular_diameter_distance(zl)
    Dls = cosmo.angular_diameter_distance_z1z2(zl, ZS)
    return (cc ** 2 / (4 * np.pi * G) * Ds / (Dl * Dls)).to(u.Msun / u.Mpc ** 2).value
ZLS = [0.02, 0.05, 0.1, 0.2, 0.4]
log("Mass convention: M_ap = Sigma_crit(z_l, z_s=0.63) * sum_pix (kappa - kappa_annulus) * A_pix(z_l) inside theta_ap = 2 FWHM,"
    " background from the 2-3 FWHM annulus; KS maps are noise-dominated and biased low by the smoothing, so treat as order of magnitude.")

# ---------------------------------------------------------------- peak finding
lmax = 3 * NSIDE - 1
alm_k = hp.map2alm(kap, lmax=lmax, iter=0); alm_m = hp.map2alm(mask, lmax=lmax, iter=0)
neigh = None
results = []
for fwhm in FWHMS:
    bl = hp.gauss_beam(np.radians(fwhm / 60), lmax=lmax)
    ks = hp.alm2map(hp.almxfl(alm_k, bl), NSIDE)
    ms = hp.alm2map(hp.almxfl(alm_m, bl), NSIDE)
    ks = np.where(ms > 0.5, ks / np.clip(ms, 0.5, None), 0.0)     # normalised convolution
    interior = ms > 0.95
    sig = 1.4826 * np.median(np.abs(ks[interior] - np.median(ks[interior])))
    sn = np.where(interior, ks / sig, 0.0)
    # local maxima
    if neigh is None: neigh = hp.get_all_neighbours(NSIDE, np.arange(hp.nside2npix(NSIDE)))
    cand = np.where(sn >= SN_LIST)[0]
    nb = neigh[:, cand]; nbv = np.where(nb >= 0, sn[np.clip(nb, 0, None)], -np.inf)
    peaks = cand[np.all(sn[cand][None, :] > nbv, 0)]
    # distance to mask edge
    edge = np.where((ms > 0.05) & (ms < 0.95))[0]; edge_vec = np.array(hp.pix2vec(NSIDE, edge)).T
    log(f"\n=== FWHM {fwhm:.0f} arcmin: sigma_MAD = {sig:.5f}, interior area {interior.sum()*hp.nside2pixarea(NSIDE,degrees=True):.0f} deg^2, "
        f"peaks S/N>={SN_LIST}: {len(peaks)}, S/N>={SN_HEAD}: {(sn[peaks] >= SN_HEAD).sum()}")
    # Gaussian expectation for comparison (Bond & Efstathiou peak counts are not needed; report simple pixel-count rate)
    log(f"    interior pixels with S/N>=4: {(sn >= 4).sum()}, >=5: {(sn >= 5).sum()} (of {interior.sum()})")
    rmatch = np.radians((MATCH_BASE + fwhm) / 60)
    for p in peaks[np.argsort(-sn[peaks])]:
        ra, dec = hp.pix2ang(NSIDE, p, lonlat=True); l, b = hp.pix2ang(NSIDE, p, lonlat=True)
        g_ = SkyCoord(ra=ra * u.deg, dec=dec * u.deg).galactic; gl, gb = g_.l.deg, g_.b.deg
        v = hp.pix2vec(NSIDE, p)
        matches = []
        for k, cv in cat_vec.items():
            dist = np.arccos(np.clip(cv @ v, -1, 1)); j = np.argmin(dist)
            if dist[j] < rmatch:
                z = cats[k]['z']; ex = cats[k]['extra']
                matches.append(f"{k}(sep {np.degrees(dist[j])*60:.1f}', z={'' if z is None else f'{z[j]:.3f}'}"
                               f"{'' if ex is None else f', {ex[j]:.3g}'})")
        d_edge = np.degrees(np.arccos(np.clip(edge_vec @ v, -1, 1)).min()) * 60 if len(edge) else 999
        flags = []
        if ms[p] < 0.95 or d_edge < fwhm: flags.append('EDGE')
        if abs(gb) < 20: flags.append('LOWLAT')
        # aperture mass
        disc = hp.query_disc(NSIDE, v, np.radians(2 * fwhm / 60)); ann = hp.query_disc(NSIDE, v, np.radians(3 * fwhm / 60))
        ann = np.setdiff1d(ann, disc); disc = disc[ms[disc] > 0.5]; ann = ann[ms[ann] > 0.5]
        ksum = (ks[disc] - ks[ann].mean()).sum() if len(ann) else ks[disc].sum()
        masses = {}
        for zl in ZLS:
            apix = (cosmo.angular_diameter_distance(zl).value * np.sqrt(hp.nside2pixarea(NSIDE))) ** 2  # Mpc^2 physical
            masses[zl] = sigma_crit(zl) * ksum * apix
        results.append(dict(fwhm=fwhm, pix=p, sn=sn[p], kappa=ks[p], ra=ra, dec=dec, l=gl, b=gb, matches=matches,
                            flags=flags, d_edge=d_edge, masses=masses))

# ---------------------------------------------------------------- report
log("\n\n================ ALL PEAKS (S/N >= 4) ================")
log("FWHM  S/N   kappa    RA       Dec      l       b     edge'  flags   M_ap(z_l=0.05,0.1,0.2) [Msun]   counterparts")
unmatched = []
for r in sorted(results, key=lambda r: (r['fwhm'], -r['sn'])):
    m = r['masses']
    log(f"{r['fwhm']:4.0f} {r['sn']:5.2f} {r['kappa']:+.4f} {r['ra']:8.3f} {r['dec']:+8.3f} {r['l']:7.2f} {r['b']:+7.2f} {r['d_edge']:5.0f}  "
        f"{','.join(r['flags']) or '-':7s} {m[0.05]:.1e} {m[0.1]:.1e} {m[0.2]:.1e}   {'; '.join(r['matches']) or 'NONE'}")
    if not r['matches']: unmatched.append(r)

log("\n\n================ UNMATCHED PEAKS ================")
log(f"{len(unmatched)} unmatched of {len(results)} peaks (all scales pooled); "
    f"unmatched with S/N>=5: {sum(r['sn'] >= SN_HEAD for r in unmatched)}; "
    f"unmatched, S/N>=5, no EDGE flag: {sum(r['sn'] >= SN_HEAD and 'EDGE' not in r['flags'] for r in unmatched)}")
for r in sorted(unmatched, key=lambda r: -r['sn']):
    dmean, dmin, vlos = los_density(r['l'], r['b'])
    s = (f"FWHM {r['fwhm']:.0f}' S/N {r['sn']:.2f} kappa {r['kappa']:+.4f} RA {r['ra']:.3f} Dec {r['dec']:+.3f} (l,b)=({r['l']:.1f},{r['b']:+.1f}) "
         f"edge {r['d_edge']:.0f}' flags [{','.join(r['flags']) or '-'}]\n"
         f"     M_ap: " + ' '.join(f"z_l={z}: {r['masses'][z]:.1e}" for z in ZLS) + "\n"
         f"     2M++ LOS 20-200 Mpc/h: mean delta {dmean:+.2f}, min delta {dmin:+.2f} ({'UNDERDENSE' if dmean < 0 else 'overdense'}); "
         f"2M++ mean LOS peculiar velocity {vlos:+.0f} km/s (+ = flowing toward this direction)")
    v = hp.ang2vec(r['ra'], r['dec'], lonlat=True)
    bf_ = SkyCoord(l=298 * u.deg, b=-8 * u.deg, frame='galactic').icrs; bfc = hp.ang2vec(bf_.ra.deg, bf_.dec.deg, lonlat=True)
    s += f"\n     angle between peak and CF4 bulk-flow axis (l,b)=(298,-8) (Watkins+23): {np.degrees(np.arccos(np.clip(np.dot(v, bfc), -1, 1))):.0f} deg"
    log(s)
out.close()
