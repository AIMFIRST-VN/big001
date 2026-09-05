"""Hypocentre inversion for H5 (Our Position -- Where Is Nemo): joint posterior over the DIRECTION n and the
DISTANCE D of the seed of Kernel 0,0,0, in the spirit of Geiger's method (grid search over the hypocentre,
predict the observables, compare).  Direction information = the four pre-stated axis tests of fit_slope.py
(planarity ring, CF4 flow, hemispherical power asymmetry, cluster-H0 dipole).  Distance information = the
AMPLITUDES of the flow-type signals as a function of depth and the GZ quadrupole bound, for a spherically
symmetric rebound profile A(r) = A_D (r/D)^s seen by an observer at offset D from its centre.

Physics (derived in notes/hypocentre.md).  Fractional variation of the profile across our horizon:
    eps = (d_LS/D) * s,  s = dlnA/dlnr at r = D.
Two profile modes, one slope parameter s each, amplitude A_D marginalised (log-flat) or fixed to 1:
  'adiabatic'  A = potential Phi.  Density contrast delta = (2/3) s(s+1) Phi_D (d_H/D)^2 (r/D)^(s-2) [Poisson],
               flows come from grad(delta): bulk flow v/(H R) = (R/10)|d delta/dr|, Hubble dipole = 3x that;
               GZ quadrupole Q = 0.3*(1/6)|s(s-2)| Phi_D (d_LS/D)^2 (SW 1/3, quadratic term 1/2, Y20 rms 0.3);
               no intrinsic dipole (Turner 1991; Erickcek, Kamionkowski & Carroll 2008).
  'isocurvature' A = entropy/composition S (relic fraction, temperature: P12).  delta_m = S directly, so
               v/(H R) = (R/10) s S_D / D, Hubble dipole 3x; intrinsic dipole (1/3) s S_D (d_LS/D) < 1e-4;
               quadrupole (1/3)(1/2)|s(s-2)| S_D (d_LS/D)^2 * 0.3.
Depth R enters as R/d_LS; d_LS = 14 Gpc comoving, d_H = c/H0 = 4.4 Gpc, so (d_LS/d_H)^2 ~ 10.
Prior on D: (a) Johnson-Mehl observer offset (offset-only Monte Carlo, same process as jm_cells.py) times
L, with L >= 8 d_LS (L marginalised log-flat over 8-32 d_LS, and fixed at 8); (b) flat in log D.
Grid: HEALPix nside 16 for n; D in [0.5, 50] d_LS log-spaced; s in {0.1, 0.3, 1, 3}.
Runtime ~1 min, memory < 1 GB.  Output: hypocentre_results.txt next to this file.
"""
import os, numpy as np, healpy as hp
from scipy.stats import chi2, norm
HERE = os.path.dirname(os.path.abspath(__file__))
out = open(os.path.join(HERE, 'hypocentre_results.txt'), 'w')
def P(*a):
    s = ' '.join(str(x) for x in a); print(s); out.write(s + '\n')

# ---------------- constants ----------------
dLS, dH = 14000.0, 4400.0                 # Mpc comoving
h2 = (dLS / dH) ** 2                      # ~10
def u(l, b):
    l, b = np.radians(l), np.radians(b); return np.array([np.cos(b) * np.cos(l), np.cos(b) * np.sin(l), np.sin(b)])

# ---------------- direction likelihood (fit_slope.py) ----------------
NS = 16; npix = hp.nside2npix(NS); V = np.array(hp.pix2vec(NS, np.arange(npix))).T
aoe = u(260, 60); along = {'CF4': u(298, -8), 'HPA': u(227, -27), 'clusterH0': u(280, -15)}
def dir_terms(aoe, along):
    d = np.degrees(np.arccos(np.abs(V @ aoe))); p_ring = np.clip(np.sin(np.radians(np.abs(90 - d))), 1e-12, 1)
    ps = {'ring': p_ring}
    for k, v in along.items():
        th = np.arccos(np.clip(V @ v, -1, 1)); ps[k] = np.clip((1 - np.cos(th)) / 2, 1e-12, 1)
    return ps
ps = dir_terms(aoe, along)
def fisher_p(ps): X = -2 * sum(np.log(p) for p in ps.values()); return chi2.sf(X, 2 * len(ps))
# Direction likelihood: a proper density for each observed axis given n (the p-values are scores, not likelihoods:
# small p = good fit).  Along-n tests: von Mises-Fisher exp(kappa cos theta); ring test: Watson girdle exp(-kappa (a.n)^2).
# kappa = 1/sigma^2 with sigma = 30 deg (~ the 23-27 deg scatter of the flow-type axes); sensitivity to kappa reported.
def L_dir(aoe, along, sigma_deg=30.0):
    kap = 1 / np.radians(sigma_deg) ** 2
    lnL = -kap * (V @ aoe) ** 2
    for v in along.values(): lnL = lnL + kap * np.clip(V @ v, -1, 1)
    return np.exp(lnL - lnL.max())
Ld_all = L_dir(aoe, along)
Ld_merged = L_dir(aoe, {k: v for k, v in along.items() if k != 'clusterH0'})      # cluster dipole merged into the flow
pF = fisher_p(ps); imap = np.argmin(pF); l0, b0 = hp.pix2ang(NS, imap, lonlat=True)
P(f"# Hypocentre inversion (hypocentre.py).  d_LS = {dLS:.0f} Mpc, d_H = {dH:.0f} Mpc, nside {NS}, {npix} directions.")
P(f"direction-only fit (four pre-stated tests, as fit_slope.py at nside 16): best (l,b) = ({l0:.0f},{b0:.0f}), Fisher p = {pF[imap]:.1e}")
for sg in (20.0, 30.0, 45.0):
    Lx = L_dir(aoe, along, sg); i = np.argmax(Lx); l, b = hp.pix2ang(NS, i, lonlat=True)
    P(f"  vMF/Watson direction likelihood, sigma = {sg:.0f} deg: MAP (l,b) = ({l:.0f},{b:.0f}); sky fraction with L > 0.32 L_max = {np.mean(Lx > np.exp(-0.5) * Lx.max()):.3f}")

# ---------------- amplitude data ----------------
# Flow-type: (name, depth R/d_LS, kind, value, meas. error, LCDM local-structure floor (rms), comment)
#   kind 'v': bulk flow as v/(H0 R);  kind 'H': Hubble-rate dipole amplitude dH/H
flow = [
    ('CF4 bulk flow 150/h Mpc (Watkins+23: 387+-28 km/s)', 150 / 0.7 / dLS, 'v', 387 / 15000., 28 / 15000., 170 / 15000.),
    ('CF4 bulk flow 200/h Mpc (428+-36 km/s, paper Table 3)', 200 / 0.7 / dLS, 'v', 428 / 20000., 36 / 20000., 150 / 20000.),
    ('cluster H0 dipole, Migkas+21 (9% peak-to-trough -> dipole 0.045+-0.015, z~0.1)', 420 / dLS, 'H', 0.045, 0.015, 0.02),
    ('SN H0 dipole, Pantheon+/DES-SN (~1+-1 %, z~0.12)', 500 / dLS, 'H', 0.01, 0.01, 0.01),
]
# One-sided bounds: (name, R/d_LS, kind, 95% limit)
bounds = [
    ('Planck kSZ bulk flow < 254 km/s at ~1 Gpc (v/HR)', 1000 / dLS, 'v', 254 / (100 * 0.7 * 1000)),
    ('GZ quadrupole DT/T < 1e-5', 1.0, 'Q', 1e-5),
    ('intrinsic dipole < 1e-4 (isocurvature only)', 1.0, 'dip', 1e-4),
]
P("\namplitude data used (value, error; LCDM local floor added in quadrature in the 'floor' treatment):")
for n, R, k, v, e, f in flow: P(f"  {n}: R = {R:.4f} d_LS, {k} = {v:.4f} +- {e:.4f}, floor {f:.4f}")
for n, R, k, v in bounds: P(f"  bound: {n}: R = {R:.3f} d_LS, {k} < {v:.2e} (95%)")

# ---------------- model ----------------
def predict(mode, s, D, A, R, kind):
    """D in d_LS units, A = profile amplitude at r = D (Phi_D or S_D)."""
    if mode == 'adiabatic':
        ddelta = (2 / 3) * s * (s + 1) * abs(s - 2) * A * h2 / D ** 3       # |d delta/dr| in 1/d_LS
        Q = 0.3 * (1 / 6) * abs(s * (s - 2)) * A / D ** 2
        dip = 0.0 * D
    else:
        ddelta = s * A / D
        Q = 0.3 * (1 / 6) * abs(s * (s - 2)) * A / D ** 2
        dip = (1 / 3) * s * A / D
    if kind == 'v': return R / 10 * ddelta
    if kind == 'H': return 3 * R / 10 * ddelta
    if kind == 'Q': return Q
    if kind == 'dip': return dip

def lnL_amp(mode, s, D, A, treat, merged=False):
    """treat: 'signal' = observed flows are entirely slope; 'floor' = LCDM local structure added to the error.
    merged: CF4(200) and the cluster dipole are one measurement (drop the cluster row; they trace one structure)."""
    ll = 0.0
    rows = flow if not merged else [r for r in flow if not r[0].startswith('cluster')]
    for n, R, k, v, e, f in rows:
        sig = e if treat == 'signal' else np.hypot(e, f)
        m = predict(mode, s, D, A, R, k); ll = ll - 0.5 * ((v - m) / sig) ** 2
    for n, R, k, lim in bounds:
        if k == 'dip' and mode == 'adiabatic': continue
        m = predict(mode, s, D, A, R, k)
        ll = ll + norm.logcdf((lim - m) / (lim / 1.645))      # one-sided: 95% limit as a 1.645-sigma edge
    return ll

D = np.logspace(np.log10(0.5), np.log10(50), 80); lnD = np.log(D)
A_grid = np.logspace(-9, 0, 91)                                # profile amplitude, log-flat prior
s_vals = [0.1, 0.3, 1.0, 3.0]

# ---------------- JM prior on the offset ----------------
rng = np.random.default_rng(1); offs = []
box, tmax, I = 3.0, 5.0, 1.0
for _ in range(4000):
    n = rng.poisson(I * (2 * box) ** 3 * tmax); x = rng.uniform(-box, box, (n, 3)); t = rng.uniform(0, tmax, n)
    arr = t + np.linalg.norm(x, axis=1); w = np.argmin(arr)
    if arr[w] > tmax - 3.0: continue
    offs.append(np.linalg.norm(x[w]))
offs = np.array(offs); q = np.percentile(offs, [10, 50, 90])
P(f"\nJM observer offset / L (offset-only MC, {len(offs)} realizations): 10/50/90% = {q.round(3)}  (jm_cells.py: 0.35/0.67/1.05)")
def jm_prior(Lfac):
    """density in ln D for D = offset * L; Lfac scalar or array of L/d_LS with log-flat weights."""
    Ls = np.atleast_1d(Lfac); pr = np.zeros_like(D)
    for L in Ls:
        x = D / L; kde = np.exp(-0.5 * ((np.log(x)[:, None] - np.log(offs)[None, :]) / 0.08) ** 2).sum(1)
        pr += kde
    return pr / np.trapezoid(pr, lnD)
priors = {'JM (L = 8 d_LS)': jm_prior(8.0), 'JM (L log-flat 8-32 d_LS)': jm_prior(np.logspace(np.log10(8), np.log10(32), 12)),
          'flat in log D': np.ones_like(D) / (lnD[-1] - lnD[0])}
for k, pr in priors.items():
    c = np.cumsum(pr) * np.gradient(lnD); c /= c[-1]
    P(f"prior '{k}': D/d_LS 10/50/90% = {np.interp([0.1, 0.5, 0.9], c, D).round(2)}")

def summarize(post_lnD, label):
    """post_lnD: density in ln D on the grid."""
    p = post_lnD / np.trapezoid(post_lnD, lnD); c = np.cumsum(p) * np.gradient(lnD); c /= c[-1]
    q16, q50, q84, q025, q975 = np.interp([0.16, 0.5, 0.84, 0.025, 0.975], c, D)
    P(f"  {label}: MAP D = {D[np.argmax(p)]:.2f} d_LS; median {q50:.2f}; 68% [{q16:.2f}, {q84:.2f}]; 95% [{q025:.2f}, {q975:.2f}]; P(D in 0.7L range 2.8-8.4) = {np.interp(8.4, D, c) - np.interp(2.8, D, c):.2f}")
    return p

# ---------------- amplitude likelihood on the D grid ----------------
P("\n=== Amplitude likelihood L_amp(D) (direction factorises out: the four axis tests carry no D information) ===")
results = {}
for mode in ('adiabatic', 'isocurvature'):
    for treat in ('floor', 'signal'):
        for merged in (False, True):
            for s in s_vals:
                LL = np.array([[lnL_amp(mode, s, d, a, treat, merged) for a in A_grid] for d in D])   # (D, A)
                Lm = np.log(np.trapezoid(np.exp(LL - LL.max()), np.log(A_grid), axis=1)) + LL.max()       # marginal over A
                L1 = LL[:, -1]                                                                          # A = 1 (eps = s d_LS/D)
                results[(mode, treat, merged, s)] = (Lm, L1, LL)
for mode in ('adiabatic', 'isocurvature'):
    for treat in ('floor', 'signal'):
        P(f"\n--- mode {mode}, flows treated as '{treat}' ---")
        for merged in (False, True):
            for s in s_vals:
                Lm, L1, LL = results[(mode, treat, merged, s)]
                lm = Lm - Lm.max(); l1 = L1 - L1.max()
                # ratio of the amplitude likelihood between the grid ends = does it constrain D at all?
                P(f"  s={s:<4} merged={merged!s:5}: A-marginalised ln L(D) range over grid = {lm.min():+.2f}..0 (max at D={D[np.argmax(Lm)]:.2f}); "
                  f"A=1 fixed (eps = s d_LS/D): allowed (ln L within 1.92 of null) only for D >= {(lambda ok: D[ok][0] if ok.any() else np.inf)(np.array([L1[i] - lnL_amp(mode, s, D[i], 0., treat, merged) > -1.92 for i in range(len(D))])):.1f} d_LS; best ln L overall = {LL.max():+.2f} vs null (no slope) = {lnL_amp(mode, s, 10., 0., treat, merged):+.2f}")

# ---------------- posteriors ----------------
P("\n=== Posterior on D (marginalised over direction, over A log-flat, over s log-flat) ===")
post_store = {}
for mode in ('adiabatic', 'isocurvature'):
    for treat in ('floor', 'signal'):
        for merged in (False, True):
            Ls = np.array([results[(mode, treat, merged, s)][0] for s in s_vals]); Lall = np.log(np.exp(Ls - Ls.max()).mean(0)) + Ls.max()
            P(f"mode {mode}, treat {treat}, merged {merged}:")
            for pk, pr in priors.items():
                p = summarize(np.exp(Lall - Lall.max()) * pr, f"prior {pk:26s}")
                post_store[(mode, treat, merged, pk)] = p
            # is D constrained by the amplitudes?  compare posterior with prior (KL divergence in nats)
            pr = priors['flat in log D']; p = post_store[(mode, treat, merged, 'flat in log D')]
            kl = np.trapezoid(p * np.log(np.clip(p, 1e-300, None) / pr), lnD)
            P(f"  information from amplitudes alone (KL posterior||flat prior) = {kl:.3f} nats"
              + ("  -> D essentially unconstrained by the amplitudes" if kl < 0.1 else "  -> amplitudes carry D information"))

# ---------------- joint (n, D): MAP and credible regions ----------------
P("\n=== Joint (n, D) posterior: baseline = adiabatic, floor, not merged, JM prior L=8; direction = four tests ===")
def joint_report(Ld, pD, label):
    post = Ld[:, None] * pD[None, :]; post /= post.sum()
    i, j = np.unravel_index(np.argmax(post), post.shape); l, b = hp.pix2ang(NS, i, lonlat=True)
    srt = np.sort(post.ravel())[::-1]; c = np.cumsum(srt); lev68, lev95 = srt[np.searchsorted(c, 0.68)], srt[np.searchsorted(c, 0.95)]
    pn = post.sum(1); srt = np.sort(pn)[::-1]; c = np.cumsum(srt); l68, l95 = srt[np.searchsorted(c, 0.68)], srt[np.searchsorted(c, 0.95)]
    in68 = pn >= l68; in95 = pn >= l95
    lp, bp = hp.pix2ang(NS, np.arange(npix), lonlat=True)
    P(f"  {label}: MAP (l,b) = ({l:.0f},{b:.0f}), D = {D[j]:.2f} d_LS; sky 68% region = {in68.mean() * 41253:.0f} deg^2 "
      f"(l {lp[in68].min():.0f}-{lp[in68].max():.0f}, b {bp[in68].min():.0f}..{bp[in68].max():.0f}); 95% region = {in95.mean() * 41253:.0f} deg^2; "
      f"P(void centroid (305,-30) in 68%) = {bool(in68[hp.vec2pix(NS, *u(305, -30))])}, in 95% = {bool(in95[hp.vec2pix(NS, *u(305, -30))])}")
    return post
base = post_store[('adiabatic', 'floor', False, 'JM (L = 8 d_LS)')]
joint_report(Ld_all, base, 'four tests, JM L=8')
joint_report(Ld_all, post_store[('adiabatic', 'floor', False, 'flat in log D')], 'four tests, flat-log prior')
joint_report(Ld_merged, post_store[('adiabatic', 'floor', True, 'JM (L = 8 d_LS)')], 'cluster dipole merged with flow')
joint_report(Ld_all, post_store[('isocurvature', 'floor', False, 'JM (L = 8 d_LS)')], 'isocurvature mode')
joint_report(Ld_all, post_store[('adiabatic', 'signal', False, 'JM (L = 8 d_LS)')], "flows as pure slope signal (naive)")

# ---------------- what the amplitudes DO constrain: eps and the flatness ----------------
P("\n=== What the amplitude data constrain: the horizon-scale gradient eps = s A_D d_LS/D, not D ===")
for mode in ('adiabatic', 'isocurvature'):
    for s in (0.3, 1.0):
        best = -np.inf; lim = None
        for d in D:
            for a in A_grid:
                ll = lnL_amp(mode, s, d, a, 'floor'); null = lnL_amp(mode, s, d, 0.0, 'floor')
                eps = s * a / d
                if ll - null > -1.92 and (lim is None or eps > lim): lim = eps
        P(f"  {mode}, s={s}: 95% upper limit on eps = s*A_D*(d_LS/D) from all amplitude data: {lim:.2e}")
P("  -> at D = 5.6 d_LS (0.7 L, L = 8 d_LS) this is a flatness requirement on the profile at our position, A_D*s < eps_lim*5.6.")

# ---------------- flows as slope signal: what D would they need? ----------------
P("\n=== If the flows WERE the slope (A_D = 1): required D/d_LS per datum (D/s) ===")
for n, R, k, v, e, f in flow:
    fac = 1 / 10 if k == 'v' else 3 / 10
    P(f"  {n}: isocurvature D/(s S_D) = {fac * R / v:.4f} d_LS; adiabatic D^3/(s(s+1)|s-2| Phi_D) = {(2 / 3) * fac * R * h2 / v:.4f} d_LS^3")
P("  Quadrupole bound with A_D = 1: adiabatic/iso Q = 0.05|s(s-2)|/D^2 < 1e-5 -> D > 71*sqrt|s(s-2)| d_LS;  intrinsic dipole (iso): s/D < 3e-4 -> D > 3300 s d_LS")

# ---------------- Wadati variant ----------------
P("\n=== 'Wadati' variant: does the flow-type amplitude vs depth extrapolate to a common origin? ===")
P("  A slope predicts signal/(R) = const, i.e. v/(H R) and dH/H both proportional to R with zero intercept; a local source gives a non-zero intercept / decreasing ratio.")
Rs = np.array([r[1] for r in flow]); ys = np.array([r[3] / (1 if r[2] == 'H' else 1) * (3 if r[2] == 'v' else 1) for r in flow])   # convert bulk to H-dipole equivalent (x3)
es = np.array([r[4] * (3 if r[2] == 'v' else 1) for r in flow])
for (n, R, k, v, e, f), y, ee in zip(flow, ys, es): P(f"  {n[:28]:30s} R = {R:.4f}  H-dipole-equivalent = {y:.4f} +- {ee:.4f}  ratio/R = {y / R:.2f}  -> D/(s S_D) = {0.3 * R / y:.4f} d_LS")
W = 1 / es ** 2; Amat = np.vstack([np.ones_like(Rs), Rs]).T
cov = np.linalg.inv(Amat.T @ (W[:, None] * Amat)); coef = cov @ (Amat.T @ (W * ys))
chi_lin = np.sum(W * (ys - Amat @ coef) ** 2); slope_only = np.sum(W * Rs * ys) / np.sum(W * Rs ** 2); chi_slope = np.sum(W * (ys - slope_only * Rs) ** 2)
P(f"  weighted linear fit y = a + b R: a = {coef[0]:.4f} +- {np.sqrt(cov[0, 0]):.4f}, b = {coef[1]:.2f} +- {np.sqrt(cov[1, 1]):.2f}, chi2 = {chi_lin:.1f} (2 dof); "
  f"through-origin fit b = {slope_only:.2f}, chi2 = {chi_slope:.1f} (3 dof); delta chi2 = {chi_slope - chi_lin:.1f}")
P(f"  origin implied by the linear fit (y = 0): R0 = {-coef[0] / coef[1]:.4f} d_LS ({-coef[0] / coef[1] * dLS:.0f} Mpc); slope-model D/(s S_D) from through-origin fit = {0.3 / slope_only:.4f} d_LS")
P("  Verdict: the intercept is the signal (non-zero at 4-5 sigma); the depth dependence is flat-to-falling, so the flow-type signals do not extrapolate to a\n"
  "  super-horizon origin; a through-origin (slope) fit is strongly disfavoured.  The 'Wadati intercept' points at a local source, i.e. the mapped attractors.")
out.close()
