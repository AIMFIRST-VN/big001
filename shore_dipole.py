"""Beach-shore test: cluster-count dipole along the fitted slope axis (l,b) = (270,-22).

A super-horizon slope eps boosts massive-cluster abundance on the dense side by R*eps with
R = delta_c/sigma^2(M) (peak-background split), i.e. a count dipole of a few percent for eps ~ 1e-2.
We measure the count dipole in PSZ2, eRASS1, ACT DR5 and SPT-SZ over modelled footprints.

Pre-registered choices (fixed before looking at any dipole):
  * fit pixels nside 16; footprint from analytic cuts (|b| > 20; eRASS1: l > 180) times an
    nside-8 occupancy map (eRASS1 only; PSZ2 analytic; ACT/SPT nside-16 occupancy); pixel area fractions from nside-64 subsampling of the analytic cut
  * model N_i = n0 A_i (1 + d . r_i), Poisson-weighted linear least squares (weights 1/(nbar A_i))
  * fixed-axis statistic: 2-parameter fit (n0, d_n) with n = (270,-22); free-axis: 4-parameter fit
  * null: 2000 multinomial shuffles of the N objects over footprint pixels with probability prop. to A_i
  * kinematic term (2 + x) beta, x = 1, beta = 370/299792, along (264,48), subtracted from the free vector
  * response R = delta_c / sigma^2(M_typ), EH98 no-wiggle, sigma_8 = 0.81, M_typ per catalogue
  * redshift split: z < 0.05 ("local", D < ~215 Mpc), 0.05-0.2, > 0.2 (eRASS1 also 0.2-0.5, > 0.5)
Runtime ~1 minute, memory < 1 GB.
"""
import os, numpy as np, healpy as hp
from astropy.io import fits
from astropy.coordinates import SkyCoord
import astropy.units as u
from scipy import integrate

HERE = os.path.dirname(os.path.abspath(__file__)); D = os.path.join(HERE, 'data')
NS, NS_OCC, NS_SUB = 16, 8, 64
BCUT = 20.0
AXIS = (270.0, -22.0)
APEX = (264.0, 48.0)
NSHUF = 2000
BETA = 370.0 / 299792.458
XKIN = 1.0
D_KIN = (2 + XKIN) * BETA
rng = np.random.default_rng(20260905)
out = open(os.path.join(HERE, 'shore_dipole_results.txt'), 'w')
def log(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n'); out.flush()

def gal(ra, dec):
    c = SkyCoord(ra=np.asarray(ra) * u.deg, dec=np.asarray(dec) * u.deg).galactic
    return c.l.deg, c.b.deg
def unit(l, b): return np.array(hp.ang2vec(l, b, lonlat=True))
NHAT = unit(*AXIS); AHAT = unit(*APEX)
def lb_of(v): return np.degrees(np.arctan2(v[1], v[0])) % 360, np.degrees(np.arcsin(v[2] / np.linalg.norm(v)))

# ---------------------------------------------------------------- catalogues (Galactic l, b, z, extra)
cats = {}
t = fits.open(os.path.join(D, 'HFI_PCCS_SZ-union_R2.08.fits'))[1].data
z = np.where(t['REDSHIFT'] > 0, t['REDSHIFT'], np.nan)
cats['PSZ2'] = dict(l=t['GLON'], b=t['GLAT'], z=z, snr=t['SNR'], cosmo=t['COSMO'].astype(bool), region='allsky', M=4.5e14)
t = fits.open(os.path.join(D, 'erass1cl_primary_v3.2.fits'))[1].data
l, b = gal(t['RA'], t['DEC'])
cats['eRASS1'] = dict(l=l, b=b, z=np.asarray(t['BEST_Z'], float), f500=np.asarray(t['F500'], float) * 1e-14,  # catalogue unit 1e-14 cgs
                      exp=np.asarray(t['EXP'], float),
                      detlike=np.asarray(t['DET_LIKE_0'], float), region='west', M=2.0e14)
t = fits.open(os.path.join(D, 'DR5_cluster-catalog_v1.1.fits'))[1].data
l, b = gal(t['RADeg'], t['decDeg'])
cats['ACT_DR5'] = dict(l=l, b=b, z=np.asarray(t['redshift'], float), region='occ', M=3.0e14)
rows = [r for r in open(os.path.join(D, 'spt_sz_bocquet2019_table5.dat'))]
l, b = gal([float(r[22:30]) for r in rows], [float(r[31:39]) for r in rows])
cats['SPT-SZ'] = dict(l=l, b=b, z=np.array([float(r[69:74] or 'nan') for r in rows]), region='occ', M=3.5e14)
for k, v in cats.items():
    v['l'] = np.asarray(v['l'], float); v['b'] = np.asarray(v['b'], float)

# ---------------------------------------------------------------- footprints
NPIX = hp.nside2npix(NS)
pix_vec = np.array(hp.pix2vec(NS, np.arange(NPIX))).T           # (NPIX,3)
sub = np.arange(hp.nside2npix(NS_SUB)); sub_l, sub_b = hp.pix2ang(NS_SUB, sub, lonlat=True)
sub_parent = hp.ring2nest(NS_SUB, sub) // (NS_SUB // NS) ** 2
sub_parent = hp.nest2ring(NS, sub_parent)

def footprint(cat, bcut=BCUT):
    """Return area fraction A_i (nside NS, ring) of each fit pixel inside the modelled footprint."""
    l, b = cat['l'], cat['b']
    good_sub = np.abs(sub_b) > bcut
    if cat['region'] == 'west': good_sub &= (sub_l > 180.0)
    A = np.bincount(sub_parent, weights=good_sub.astype(float), minlength=NPIX) / (NS_SUB // NS) ** 2
    # occupancy at nside 8 (drops unobserved regions, e.g. eRASS1 holes, and defines ACT/SPT footprints)
    occ8 = np.zeros(hp.nside2npix(NS_OCC), bool); occ8[hp.ang2pix(NS_OCC, l, b, lonlat=True)] = True
    p8 = hp.ang2pix(NS_OCC, *hp.pix2ang(NS, np.arange(NPIX), lonlat=True), lonlat=True)
    if cat['region'] == 'occ':
        # footprint-limited surveys: occupancy at nside 16 itself (parent-8 occupancy overestimates edges)
        occ16 = np.zeros(NPIX, bool); occ16[hp.ang2pix(NS, l, b, lonlat=True)] = True
        A = A * occ16
    elif cat['region'] == 'west':
        A = A * occ8[p8]          # eRASS1: ~30 objects per nside-8 pixel, occupancy safe; PSZ2 (allsky) uses the analytic cut only
    return A

# ---------------------------------------------------------------- dipole fits
def design(A, fixed):
    if fixed: X = np.column_stack([A, A * (pix_vec @ NHAT)])
    else:     X = np.column_stack([A, A[:, None] * pix_vec])
    return X
def fit(N, A, fixed, nbar=None):
    """Poisson-weighted linear LS. Returns (params, cov) with params = (n0, n0*d...)."""
    sel = A > 0; X = design(A, fixed)[sel]; y = N[sel]
    if nbar is None: nbar = y.sum() / A[sel].sum()
    w = 1.0 / np.maximum(nbar * A[sel], 1e-3)
    XtW = X.T * w; C = np.linalg.inv(XtW @ X); p = C @ (XtW @ y)
    return p, C
def dipole_from(p, C):
    n0 = p[0]; d = p[1:] / n0
    err = np.sqrt(np.diag(C)[1:]) / n0
    return d, err
def shuffles(A, Ntot, nshuf=NSHUF):
    sel = A > 0; prob = A[sel] / A[sel].sum()
    Xf = design(A, True)[sel]; Xv = design(A, False)[sel]
    nbar = Ntot / A[sel].sum(); w = 1.0 / np.maximum(nbar * A[sel], 1e-3)
    Pf = np.linalg.solve((Xf.T * w) @ Xf, Xf.T * w); Pv = np.linalg.solve((Xv.T * w) @ Xv, Xv.T * w)
    counts = rng.multinomial(Ntot, prob, size=nshuf).astype(float)   # (nshuf, npixsel)
    pf = counts @ Pf.T; pv = counts @ Pv.T
    dn = pf[:, 1] / pf[:, 0]; dv = pv[:, 1:] / pv[:, [0]]
    return dn, np.linalg.norm(dv, axis=1)

def analyse(name, A, l, b, label, do_null=True):
    N = np.bincount(hp.ang2pix(NS, l, b, lonlat=True), minlength=NPIX).astype(float)
    inside = A[hp.ang2pix(NS, l, b, lonlat=True)] > 0
    N = np.bincount(hp.ang2pix(NS, l[inside], b[inside], lonlat=True), minlength=NPIX).astype(float)
    Ntot = int(N.sum()); fsky = A.sum() / NPIX
    if Ntot < 50:
        log(f"  [{label}] N={Ntot} too few, skipped"); return None
    pf, Cf = fit(N, A, True); dn, dn_err = dipole_from(pf, Cf); dn, dn_err = dn[0], dn_err[0]
    pv, Cv = fit(N, A, False); dv, dv_err = dipole_from(pv, Cv)
    amp = np.linalg.norm(dv); lv, bv = lb_of(dv)
    dv_corr = dv - D_KIN * AHAT; dn_corr = dv_corr @ NHAT
    res = dict(N=Ntot, fsky=fsky, dn=dn, dn_err=dn_err, amp=amp, dir=(lv, bv), dn_free=dv @ NHAT, dn_kincorr=dn_corr,
               amp_err=np.sqrt(np.mean(dv_err ** 2)))
    s = (f"  [{label:28s}] N={Ntot:5d} fsky={fsky:.3f} | fixed-axis d_n = {dn:+.4f} +- {dn_err:.4f}"
         f" | free |d| = {amp:.4f} toward ({lv:.0f},{bv:.0f}), proj on n {dv @ NHAT:+.4f}, kin-corr {dn_corr:+.4f}")
    if do_null:
        sn, sa = shuffles(A, Ntot)
        p_fixed = np.mean(np.abs(sn) >= abs(dn)); p_free = np.mean(sa >= amp)
        res.update(p_fixed=p_fixed, p_free=p_free, null_dn_std=sn.std(), null_amp_med=np.median(sa))
        s += f" | null: sigma(d_n)={sn.std():.4f} p_fixed={p_fixed:.3f}, median|d|={np.median(sa):.4f} p_free={p_free:.3f}"
    log(s); return res

# ---------------------------------------------------------------- response factor R = delta_c / sigma^2(M)
def sigma_M(M, h=0.674, Om=0.315, Ob=0.049, ns=0.965, s8=0.81):
    """EH98 no-wiggle sigma(M) (M in Msun, top-hat), normalised to sigma_8."""
    rho_m = 2.775e11 * Om * h ** 2           # Msun/Mpc^3 (comoving, h-free units: uses h^2 explicitly)
    th = 2.725 / 2.7; wm = Om * h ** 2; wb = Ob * h ** 2
    s = 44.5 * np.log(9.83 / wm) / np.sqrt(1 + 10 * wb ** 0.75); ag = 1 - 0.328 * np.log(431 * wm) * wb / wm + 0.38 * np.log(22.3 * wm) * (wb / wm) ** 2
    def T(k):   # k in 1/Mpc
        G = ag + (1 - ag) / (1 + (0.43 * k * s) ** 4); q = k * th ** 2 / (wm * G)
        L = np.log(2 * np.e + 1.8 * q); Cc = 14.2 + 731 / (1 + 62.5 * q); return L / (L + Cc * q * q)
    def W(x): return 3 * (np.sin(x) - x * np.cos(x)) / x ** 3
    def s2(R):
        f = lambda lk: (lambda k: k ** (3 + ns) * T(k) ** 2 * W(k * R) ** 2)(np.exp(lk))
        return integrate.quad(f, np.log(1e-4), np.log(50), limit=200)[0]
    R8 = 8 / h; norm = s8 ** 2 / s2(R8)
    R = (3 * M / (4 * np.pi * rho_m)) ** (1 / 3)
    return np.sqrt(norm * s2(R))

# ---------------------------------------------------------------- LCDM local expectation from 2M++
den = np.load(os.path.join(HERE, 'twompp_density.npy')); ax = np.linspace(-200, 200, 257)
def local_dipole(A, rmin=10.0, rmax=135.0, bias_ratio=2.5, nr=50):
    """Predicted count map for a volume-limited local sample (r^2 dr weighting) from 2M++ delta_g,
    cluster/galaxy bias ratio bias_ratio; dipole fitted over the same footprint."""
    r = np.linspace(rmin, rmax, nr); wr = r ** 2; wr /= wr.sum()
    from scipy.ndimage import map_coordinates
    pts = pix_vec[:, None, :] * r[None, :, None]                 # (NPIX, nr, 3) Mpc/h
    idx = (pts + 200) / (400 / 256)
    dg = map_coordinates(den, [idx[..., 0].ravel(), idx[..., 1].ravel(), idx[..., 2].ravel()], order=1, mode='nearest').reshape(NPIX, nr)
    dens = 1 + bias_ratio * dg
    Nmap = (dens * wr).sum(1) * A * 1e4
    pf, Cf = fit(Nmap, A, True); pv, Cv = fit(Nmap, A, False)
    return dipole_from(pf, Cf)[0][0], dipole_from(pv, Cv)[0]

# ================================================================ run
log("Beach-shore test: cluster-count dipole along the slope axis (l,b)=(270,-22)")
log(f"nside {NS}, |b|>{BCUT}, {NSHUF} shuffles, kinematic (2+x)beta = {D_KIN:.4f} along (264,48); "
    f"cos(apex,n) = {AHAT @ NHAT:+.3f} -> kinematic term along n = {D_KIN * (AHAT @ NHAT):+.4f}")
log("Sign convention: d_n > 0 means more clusters toward (270,-22). Slope model: d_n = R*eps (+ local + kinematic).")
results = {}
for name, cat in cats.items():
    A = footprint(cat); l, b, z = cat['l'], cat['b'], cat['z']
    sig = sigma_M(cat['M']); R = 1.686 / sig ** 2
    log(f"\n=== {name}: N_total {len(l)}, footprint {A.sum() * hp.nside2pixarea(NS, degrees=True):.0f} deg^2 "
        f"(fsky {A.sum() / NPIX:.3f}); M_typ = {cat['M']:.1e} Msun, sigma(M) = {sig:.3f}, R = delta_c/sigma^2 = {R:.2f}")
    r_all = analyse(name, A, l, b, 'all')
    results[name] = dict(all=r_all, R=R)
    hz = np.isfinite(z)
    if name == 'PSZ2':
        analyse(name, A, l[cat['snr'] > 6], b[cat['snr'] > 6], 'SNR>6')
        analyse(name, A, l[cat['cosmo']], b[cat['cosmo']], 'COSMO sample')
        analyse(name, A, l[hz], b[hz], 'with z')
        analyse(name, A, l[~hz], b[~hz], 'no z')
    if name == 'eRASS1':
        # depth variation: eRASS1 exposure rises ~10x toward the south ecliptic pole (276,-30), 8 deg from the axis.
        # (a) flux-limited subsamples (the pre-registered 3e-13 cut was mis-scaled: F500 is in 1e-14 cgs; cuts revised
        #     after seeing the density-vs-ecliptic-latitude table, before seeing their dipoles)
        for fcut in (4e-13, 8e-13, 1.5e-12):
            fl = cat['f500'] > fcut
            analyse(name, A, l[fl], b[fl], f'F500>{fcut:.0e}')
        analyse(name, A, l[cat['detlike'] > 40], b[cat['detlike'] > 40], 'DET_LIKE>40')
        # (b) exposure-modelled effective area A'_i = A_i (E_i/<E>)^gamma, E_i = pixel geometric-mean catalogue EXP
        p16 = hp.ang2pix(NS, l, b, lonlat=True); p8 = hp.ang2pix(NS_OCC, l, b, lonlat=True)
        lnE16 = np.bincount(p16, weights=np.log(cat['exp']), minlength=NPIX) / np.maximum(np.bincount(p16, minlength=NPIX), 1)
        lnE8 = np.bincount(p8, weights=np.log(cat['exp']), minlength=hp.nside2npix(NS_OCC)) / np.maximum(np.bincount(p8, minlength=hp.nside2npix(NS_OCC)), 1)
        par8 = hp.ang2pix(NS_OCC, *hp.pix2ang(NS, np.arange(NPIX), lonlat=True), lonlat=True)
        lnE = np.where(np.bincount(p16, minlength=NPIX) > 0, lnE16, lnE8[par8]); lnE -= lnE[A > 0].mean()
        N = np.bincount(p16, minlength=NPIX).astype(float)
        from scipy.optimize import minimize
        sel = A > 0
        def nll(th, dip):
            lnn0, g = th[0], th[1]; mu = np.exp(lnn0 + g * lnE[sel]) * A[sel]
            if dip: mu = mu * (1 + pix_vec[sel] @ th[2:5])
            mu = np.maximum(mu, 1e-9); return (mu - N[sel] * np.log(mu)).sum()
        r0 = minimize(nll, [np.log(N[sel].sum() / A[sel].sum()), 1.0], args=(False,), method='Nelder-Mead')
        gam = r0.x[1]
        r1 = minimize(nll, np.r_[r0.x, 0, 0, 0], args=(True,), method='Nelder-Mead', options=dict(maxiter=4000, xatol=1e-6, fatol=1e-6))
        dj = r1.x[2:5]; lj, bj = lb_of(dj)
        log(f"  exposure model: N_i ~ A_i E_i^gamma; gamma (no dipole) = {gam:.3f}; joint fit gamma = {r1.x[1]:.3f}, "
            f"dipole |d| = {np.linalg.norm(dj):.4f} toward ({lj:.0f},{bj:.0f}), d_n = {dj @ NHAT:+.4f}")
        Aeff = A * np.exp(gam * lnE); Aeff *= A.sum() / Aeff.sum()
        r_exp = analyse(name, Aeff, l, b, 'exposure-corrected')
        # flux-limited subsamples: is the residual exposure dependence gone?  (gamma ~ 0 means yes)
        for fcut in (4e-13, 8e-13, 1.5e-12):
            fl = cat['f500'] > fcut; N = np.bincount(p16[fl], minlength=NPIX).astype(float)
            rf = minimize(nll, [np.log(N[sel].sum() / A[sel].sum()), 0.3], args=(False,), method='Nelder-Mead')
            log(f"  flux-limited F500>{fcut:.0e}: residual exposure index gamma = {rf.x[1]:+.3f} (0 = complete everywhere)")
        N = np.bincount(p16, minlength=NPIX).astype(float)
        for z0, z1 in [(0, 0.05), (0.05, 0.2), (0.2, 0.5), (0.5, 5)]:
            m = hz & (z >= z0) & (z < z1); analyse(name, Aeff, l[m], b[m], f'exp-corr {z0:.2f}<=z<{z1:.2f}')
        m = hz & (z >= 0.05); r_expfar = analyse(name, Aeff, l[m], b[m], 'exp-corr z>=0.05 (far)')
        results[name]['exp'] = r_exp; results[name]['expfar'] = r_expfar
    # redshift split
    bins = [(0, 0.05), (0.05, 0.2), (0.2, 0.5), (0.5, 5)] if name == 'eRASS1' else [(0, 0.05), (0.05, 0.2), (0.2, 5)]
    for z0, z1 in bins:
        m = hz & (z >= z0) & (z < z1); analyse(name, A, l[m], b[m], f'{z0:.2f}<=z<{z1:.2f}')
    m = hz & (z >= 0.05); r_far = analyse(name, A, l[m], b[m], 'z>=0.05 (far)')
    m = ~(hz & (z < 0.05)); analyse(name, A, l[m], b[m], 'excluding z<0.05')
    results[name]['far'] = r_far
    # Galactic-latitude systematic check
    if name in ('PSZ2', 'eRASS1'):
        A30 = footprint(cat, 30.0); analyse(name, A30, l, b, 'all, |b|>30')
    # LCDM local expectation
    dloc_n, dloc_v = local_dipole(A)
    lv, bv = lb_of(dloc_v)
    log(f"  2M++ local prediction (10-135 Mpc/h volume-limited, b_cl/b_g=2.5): d_n = {dloc_n:+.4f}, "
        f"free |d| = {np.linalg.norm(dloc_v):.4f} toward ({lv:.0f},{bv:.0f})")
    # eps bound
    pairs = [('all', r_all), ('far', r_far)]
    if name == 'eRASS1': pairs += [('exp', results[name]['exp']), ('expfar', results[name]['expfar'])]
    for tag, r in pairs:
        if r is None: continue
        dn_c = r['dn'] - D_KIN * (AHAT @ NHAT)
        eps = dn_c / R; eps_err = np.hypot(r['dn_err'], r.get('null_dn_std', r['dn_err'])) / R
        log(f"  eps from {tag:6s}: (d_n - kin)/R = {eps:+.4f} +- {eps_err:.4f}; 95% |eps| < {abs(eps) + 1.645 * eps_err:.4f} "
            f"(R = {R:.1f}; R = 3-7 range -> |eps| < {(abs(dn_c) + 1.645 * eps_err * R) / 3:.4f} - {(abs(dn_c) + 1.645 * eps_err * R) / 7:.4f})")
log("\nDone.")
out.close()
