#!/usr/bin/env python3
"""
ring_echo.py -- pre-registered search for a SECOND oscillatory component in the
galaxy power-spectrum monopole, distinct from the BAO (Sec. 3.1 of the paper:
recompression waves crossing the rings of H8 could imprint a ringing at a
ring-crossing scale, not at the drag-epoch sound horizon).

Data (all public, downloaded to data/ring_echo/, see README data table):
  * eBOSS DR16 LRG P0(k), 0.6<z<1.0, NGC and SGC, pre-recon, 32 bins
    0.0075 < k < 0.315 h/Mpc, EZmock covariance (1000 mocks), window W_l^2(s)
    (Gil-Marin et al. 2020, arXiv:2007.08994).
  * BOSS DR12 combined-sample P0(k), z1 (0.2-0.5), z2 (0.4-0.6), z3 (0.5-0.75),
    NGC and SGC, pre-recon, Patchy covariance (2045/2048 mocks) which is only
    published for 0.01 < k < 0.15 h/Mpc (14 bins), RR-pair-count window
    (Beutler et al. 2017, MNRAS 466, 2242).
  Fiducial cosmology of both: flat LCDM Omega_m = 0.31, h = 0.676.

Model, per sample s:
  P0_model(k) = W_s[ (1 + A2 S(k)) B_s^2 P_nw(k) {1 + (O_lin(k/alpha_z) - 1) exp(-k^2 Sig_nl^2/2)} ]
                + (1 + A2 S(k)) sum_i a_{s,i} k^i
  S(k) = sin(k L2 + phi) exp(-k^2 Sig2^2 / 2)     (the second ringing; L2 = 2pi/k2 is the
          real-space scale in Mpc/h, k2 the period in k-space; the BAO has L = r_d ~ 100 Mpc/h)
  P_nw, O_lin from the Eisenstein & Hu (1998) fitting formulae (no-wiggle and
  full transfer functions); W_s is the survey window convolution (monopole,
  with the Kaiser quadrupole leakage term); the additive polynomial has 5 terms
  (k^-2 .. k^2) for the 32-bin eBOSS samples and 3 terms (k^-1, 1, k) for the
  14-bin BOSS samples.  alpha per redshift bin (4), one Sig_nl, and the second
  component (A2, phi, Sig2, period 2pi/k2) shared by all eight samples.

Pre-registered procedure (fixed before looking at the scan):
  1. BAO check: null fit (A2 = 0) with alpha_z, Sig_nl free; expect
     alpha ~ 1 (within a few %) and Sig_nl ~ 6-9 Mpc/h.
  2. Scan: 61 log-spaced scales L2 = 2pi/k2 of 20-600 Mpc/h x damping Sig2 in {0,5,10,20}
     Mpc/h; alpha_z and Sig_nl fixed at their null best fit.  With
     A2 sin(x+phi) = a_c sin x + a_s cos x the model is linear in (a_c, a_s)
     to first order (linearised about the null best fit), so each grid point
     is one weighted linear solve; the best period is refined with the exact
     non-linear model, all nuisance parameters free.
  3. Detection threshold: delta chi2 (null - best) > 25 globally.  Otherwise
     the 95% one-sided upper limit on A2 at each period is the worst-phase
     limit max_phi [ a.u(phi) + 1.645 sqrt(u^T Cov_a u) ].
  4. Look-elsewhere: (a) analytic, p_global = 1-(1-p_local)^N_eff with
     p_local from chi2(2 dof) and N_eff = 40; (b) 500 Gaussian mock
     realisations drawn from the covariances around the null best fit, each
     put through the same null refit + scan, giving the distribution of the
     maximum delta chi2.
  A period with delta chi2 > 9 is flagged as a hint, not a detection.

  6. Crest ripples: harmonics locked to the BAO, at 2 k_BAO and 3 k_BAO (per-redshift
     k_BAO = alpha_z 2pi/r_d), free amplitude and phase, same damping grid; reported with
     errors, delta chi2 and 95% limits under the same thresholds.

Outputs: ring_echo_results.txt, figures/ring_echo.png.
Runtime ~ a few minutes on 16 cores; memory < 1.5 GB.
"""
import os, sys, time
import numpy as np
from scipy import optimize, special, stats
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
DD = os.path.join(HERE, "data", "ring_echo")
BEU = os.path.join(DD, "beutler_fs", "public_material_RSD")

NMOCK = int(sys.argv[1]) if len(sys.argv) > 1 else 500
SEED = 20260905

# ----------------------------------------------------------------- cosmology
OM, OB, H0, NS, TCMB = 0.31, 0.048, 0.676, 0.97, 2.7255


def eh98(k, wiggles=True):
    """Eisenstein & Hu 1998 transfer function, k in h/Mpc. Returns T(k)."""
    h = H0
    om, ob = OM * h * h, OB * h * h
    fb = OB / OM
    th = TCMB / 2.7
    zeq = 2.5e4 * om * th ** -4
    keq = 7.46e-2 * om * th ** -2  # 1/Mpc
    b1 = 0.313 * om ** -0.419 * (1 + 0.607 * om ** 0.674)
    b2 = 0.238 * om ** 0.223
    zd = 1291 * om ** 0.251 / (1 + 0.659 * om ** 0.828) * (1 + b1 * ob ** b2)
    Rd = 31.5 * ob * th ** -4 * (1e3 / zd)
    Req = 31.5 * ob * th ** -4 * (1e3 / zeq)
    s = 2 / (3 * keq) * np.sqrt(6 / Req) * np.log((np.sqrt(1 + Rd) + np.sqrt(Rd + Req)) / (1 + np.sqrt(Req)))
    kk = k * h  # 1/Mpc
    if not wiggles:
        a_g = 1 - 0.328 * np.log(431 * om) * fb + 0.38 * np.log(22.3 * om) * fb ** 2
        gam = OM * h * (a_g + (1 - a_g) / (1 + (0.43 * kk * s) ** 4))
        q = k * th ** 2 / gam
        L0 = np.log(2 * np.e + 1.8 * q)
        C0 = 14.2 + 731 / (1 + 62.5 * q)
        return L0 / (L0 + C0 * q * q)
    ksilk = 1.6 * ob ** 0.52 * om ** 0.73 * (1 + (10.4 * om) ** -0.95)
    a1 = (46.9 * om) ** 0.670 * (1 + (32.1 * om) ** -0.532)
    a2 = (12.0 * om) ** 0.424 * (1 + (45.0 * om) ** -0.582)
    ac = a1 ** (-fb) * a2 ** (-fb ** 3)
    bb1 = 0.944 / (1 + (458 * om) ** -0.708)
    bb2 = (0.395 * om) ** -0.0266
    bc = 1 / (1 + bb1 * ((1 - fb) ** bb2 - 1))
    q = kk / (13.41 * keq)

    def T0(q, a, b):
        C = 14.2 / a + 386 / (1 + 69.9 * q ** 1.08)
        L = np.log(np.e + 1.8 * b * q)
        return L / (L + C * q * q)
    f = 1 / (1 + (kk * s / 5.4) ** 4)
    Tc = f * T0(q, 1, bc) + (1 - f) * T0(q, ac, bc)
    y = (1 + zeq) / (1 + zd)
    G = y * (-6 * np.sqrt(1 + y) + (2 + 3 * y) * np.log((np.sqrt(1 + y) + 1) / (np.sqrt(1 + y) - 1)))
    ab = 2.07 * keq * s * (1 + Rd) ** -0.75 * G
    bnode = 8.41 * om ** 0.435
    st = s / (1 + (bnode / (kk * s)) ** 3) ** (1 / 3)
    bbb = 0.5 + fb + (3 - 2 * fb) * np.sqrt((17.2 * om) ** 2 + 1)
    x = kk * st
    j0 = np.sin(x) / x
    Tb = (T0(q, 1, 1) / (1 + (kk * s / 5.2) ** 2) + ab / (1 + (bbb / (kk * s)) ** 3) * np.exp(-(kk / ksilk) ** 1.4)) * j0
    return fb * Tb + (1 - fb) * Tc


def plin(k, wiggles=True):
    return k ** NS * eh98(k, wiggles) ** 2


# ------------------------------------------------------------------ k grids
KIN = np.logspace(-4, np.log10(3.0), 4000)       # model grid
DLNK = np.log(KIN[1] / KIN[0])


RD_H = 99.9   # sound horizon in Mpc/h (overwritten by CAMB)


def camb_plin(k):
    """Linear P(k) at z=0 from CAMB for the fiducial cosmology (falls back to EH98)."""
    try:
        import camb
    except ImportError:
        print("CAMB not available: using EH98 wiggles (sound horizon ~3% off)")
        return plin(k, True), "EH98"
    pars = camb.CAMBparams()
    pars.set_cosmology(H0=100 * H0, ombh2=OB * H0 ** 2, omch2=(OM - OB) * H0 ** 2 - 0.0006, mnu=0.06, tau=0.06)
    pars.InitPower.set_params(As=2.1e-9, ns=NS)
    pars.set_matter_power(redshifts=[0.0], kmax=5.0)
    res = camb.get_results(pars)
    kh, z, pk = res.get_matter_power_spectrum(minkh=1e-4, maxkh=3.0, npoints=1500)
    rd = res.get_derived_params()["rdrag"]
    global RD_H
    RD_H = rd * H0
    return np.exp(np.interp(np.log(k), np.log(kh), np.log(pk[0]))), f"CAMB (r_drag = {rd:.2f} Mpc = {rd * H0:.2f} Mpc/h)"


PNW_IN = plin(KIN, False)
_PL, TEMPLATE_SRC = camb_plin(KIN)
# smooth broadband: EH98 no-wiggle x low-order polynomial in ln k fitted to the ratio (removes residual shape)
_m = (KIN > 1e-3) & (KIN < 1.0)
_c = np.polyfit(np.log(KIN[_m]), np.log(_PL[_m] / PNW_IN[_m]), 7)
PNW_IN = PNW_IN * np.exp(np.polyval(_c, np.log(KIN)))
OLIN_IN = _PL / PNW_IN
OLIN_IN[KIN < 1e-3] = 1.0
BETA = 0.35                                       # Kaiser leakage P2/P0 ratio
R2 = (4 * BETA / 3 + 4 * BETA ** 2 / 7) / (1 + 2 * BETA / 3 + BETA ** 2 / 5)
SGRID = np.arange(0.5, 3200.0, 1.0)


def window_operator(kout, dk, s_w, W0, W2):
    """Linear operator M (nout x nin): P(KIN) -> bin-averaged, window-convolved P0(kout).
    P0_w(k) = P0(k) + 4pi int s^2 [ (W0-1) xi0 + W2 xi2 / 5 ] j0(ks) ds   (Wilson et al. 2017),
    xi2 from the Kaiser ratio R2 (fixed).  The correction term is a narrow convolution in k
    (width ~ 1/survey size), so it is evaluated on a fine linear k grid up to 0.6 h/Mpc with a
    smooth taper above 0.45 (modes above that cannot feed k_out < 0.32 through the window).
    Bin averaging: 5 sub-points over +-dk/2."""
    sub = np.linspace(-0.4, 0.4, 5) * dk
    ksub = (kout[:, None] + sub[None, :]).ravel()
    w0 = np.interp(SGRID, s_w, W0, left=W0[0], right=0.0) - 1.0
    w2 = np.interp(SGRID, s_w, W2, left=W2[0], right=0.0)
    klin = np.arange(2e-4, 0.6, 2e-4)
    taper = 0.5 * (1 - np.tanh((klin - 0.45) / 0.03))
    # xi_l(s) = int k^2 dk/(2pi^2) P(k) j_l(ks)
    x = np.outer(SGRID, klin)
    j0 = np.sin(x) / x
    j2 = (3 / x ** 2 - 1) * j0 - 3 * np.cos(x) / x ** 2
    wk = klin ** 2 * (klin[1] - klin[0]) / (2 * np.pi ** 2) * taper
    B = (w0[:, None] * j0 + (R2 / 5) * w2[:, None] * j2) * wk[None, :]      # ns x nlin
    del x, j0, j2
    xs = np.outer(ksub, SGRID)
    A = 4 * np.pi * SGRID[None, :] ** 2 * np.sin(xs) / xs                      # nsub x ns
    K = A @ B                                                                  # nsub x nlin
    del A, B, xs
    # interpolation matrices KIN -> klin (log-log linear) and KIN -> ksub
    def interp_matrix(kt):
        I = np.zeros((len(kt), len(KIN)))
        idx = np.clip(np.searchsorted(KIN, kt) - 1, 0, len(KIN) - 2)
        f = (np.log(kt) - np.log(KIN[idx])) / DLNK
        I[np.arange(len(kt)), idx] = 1 - f
        I[np.arange(len(kt)), idx + 1] = f
        return I
    M = interp_matrix(ksub) + K @ interp_matrix(klin)
    return M.reshape(len(kout), len(sub), len(KIN)).mean(axis=1)


# --------------------------------------------------------------- load data
class Sample:
    pass


def load_eboss(cap):
    S = Sample()
    S.name = f"eBOSS LRG {cap} z0.6-1.0"
    S.zbin = 3
    d = np.loadtxt(os.path.join(DD, f"Data_LRGPk_{cap}_0.6z1.0_prerecon.txt"))
    S.k, S.P = d[:, 0], d[:, 1]
    rows = []
    with open(os.path.join(DD, f"Covariance_LRGPk_{cap}_0.6z1.0_prerecon.txt")) as f:
        for line in f:
            p = line.split()
            if len(line) and not line.startswith("#") and p[2] == "0" and p[3] == "0":
                rows.append((float(p[4]), float(p[5]), float(p[6])))
    rows = np.array(rows)
    n = len(S.k)
    S.C = rows[:, 2].reshape(n, n)
    S.hartlap = (1000 - n - 2) / (1000 - 1)
    w = np.loadtxt(os.path.join(DD, f"Window_LRGPk_{cap}_0.6z1.0.txt"))
    S.s_w, S.W0, S.W2 = w[:, 0], w[:, 1] / w[0, 1], w[:, 2] / w[0, 1]   # normalise W0(s->0) = 1 (B^2 absorbs it)
    S.npoly = 5
    S.powers = np.array([-2, -1, 0, 1, 2])
    return S


def load_boss(z, cap):
    S = Sample()
    S.name = f"BOSS DR12 {cap} z{z}"
    S.zbin = z - 1
    nm = {("z1", "NGC"): "2045", ("z2", "NGC"): "2045", ("z3", "NGC"): "2045",
          ("z1", "SGC"): "2048", ("z2", "SGC"): "2048", ("z3", "SGC"): "2048"}
    zz = f"z{z}"
    cov = np.loadtxt(os.path.join(BEU, f"Beutleretal_cov_patchy_{zz}_{cap}_1_15_1_15_1_10_{nm[(zz, cap)]}_60.dat"), skiprows=4)
    i, j = cov[:, 0].astype(int), cov[:, 1].astype(int)
    m = (i <= 14) & (j <= 14)
    n = 14
    C = np.zeros((n, n))
    C[i[m] - 1, j[m] - 1] = cov[m, 4]
    S.C = C
    S.hartlap = (int(nm[(zz, cap)]) - n - 2) / (int(nm[(zz, cap)]) - 1)
    pk = np.array([[float(x) for x in l.split()] for l in open(os.path.join(BEU, f"Beutleretal_pk_monopole_DR12_{cap}_{zz}_prerecon_120.dat"))
                   if len(l.split()) == 4 and l.split()[0][0].isdigit()])
    kc = cov[m & (i == j), 2]
    sel = [np.argmin(np.abs(pk[:, 0] - kk)) for kk in kc]
    S.k = pk[sel, 0]          # mean k of the modes in the bin (bin centres in column 2, dk = 0.01)
    S.P = pk[sel, 2]
    assert np.allclose(pk[sel, 0], kc, rtol=1e-3)
    w = np.loadtxt(os.path.join(BEU, f"Beutleretal_window_{zz}_{cap}.dat"), skiprows=1)
    s = w[:, 0]
    ds = np.gradient(s)
    W = w[:, 2:5] / (s[:, None] ** 2 * ds[:, None])
    # pair-count shot noise dominates below s ~ 20 Mpc/h: smooth in ln s (boxcar, ~10% in s) and
    # normalise W0(s -> 0) = 1 on the 5-40 Mpc/h plateau
    from scipy.ndimage import uniform_filter1d
    W = uniform_filter1d(W, 81, axis=0, mode="nearest")
    norm = np.median(W[(s > 5) & (s < 40), 0])
    W /= norm
    W[s < 5, 0] = 1.0
    S.s_w, S.W0, S.W2 = s, W[:, 0], W[:, 1]
    S.npoly = 3
    S.powers = np.array([-1, 0, 1])
    return S


def prepare(S):
    dk = np.median(np.diff(S.k))
    S.M = window_operator(S.k, dk, S.s_w, S.W0, S.W2)
    Cinv = np.linalg.inv(S.C) * S.hartlap
    S.L = np.linalg.cholesky(Cinv)          # chi2 = |L^T r|^2
    S.Lc = np.linalg.cholesky(S.C)          # for mocks
    S.poly = np.stack([S.k ** p for p in S.powers], axis=1)
    return S


# --------------------------------------------------------------- the model
def template(alpha, signl):
    """B^2-free BAO template on KIN."""
    ol = np.interp(KIN / alpha, KIN, OLIN_IN)
    return PNW_IN * (1 + (ol - 1) * np.exp(-KIN ** 2 * signl ** 2 / 2))


def ring(k, period, sig2):
    """sin and cos basis of the second ringing (damped)."""
    x = k * period                       # real-space scale 'period' [Mpc/h] -> k-period 2pi/period
    d = np.exp(-k ** 2 * sig2 ** 2 / 2)
    return np.sin(x) * d, np.cos(x) * d


def columns(S, T, ac, as_, period, sig2):
    """Design columns for sample S given the template T on KIN and ring amplitudes."""
    if period is None:
        mod_in = 1.0
        mod_out = 1.0
    else:
        si, ci = ring(KIN, period, sig2)
        so, co = ring(S.k, period, sig2)
        mod_in = 1 + ac * si + as_ * ci
        mod_out = 1 + ac * so + as_ * co
    c0 = S.M @ (mod_in * T)
    return np.column_stack([c0, mod_out[:, None] * S.poly if np.ndim(mod_out) else S.poly])


def whitened_residual(samples, theta, ring_par):
    """Variable projection: solve linear params per sample, return whitened residuals."""
    alphas, signl = theta[:4], theta[4]
    if ring_par is None:
        ac = as_ = 0.0
        period = sig2 = None
    else:
        ac, as_, period, sig2 = ring_par
    res = []
    lin = []
    for S in samples:
        T = template(alphas[S.zbin], signl)
        X = columns(S, T, ac, as_, period, sig2)
        Xw = S.L.T @ X
        dw = S.L.T @ S.P
        b, *_ = np.linalg.lstsq(Xw, dw, rcond=None)
        res.append(dw - Xw @ b)
        lin.append(b)
    return np.concatenate(res), lin


def fit_null(samples, theta0=None):
    th0 = np.array([1, 1, 1, 1, 7.0]) if theta0 is None else theta0
    r = optimize.least_squares(lambda t: whitened_residual(samples, t, None)[0], th0,
                               bounds=([0.8] * 4 + [1.0], [1.2] * 4 + [20.0]), x_scale=[0.01] * 4 + [1.0])
    res, lin = whitened_residual(samples, r.x, None)
    J = r.jac
    try:
        cov = np.linalg.inv(J.T @ J)
        err = np.sqrt(np.abs(np.diag(cov)))
    except np.linalg.LinAlgError:
        err = np.full(len(r.x), np.nan)
    return r.x, err, float(res @ res), lin


def fit_ring_exact(samples, theta_null, lin_null, period, sig2_0, free_nuis=True):
    """Exact non-linear fit of (theta, ac, as, sig2) at fixed period."""
    best = None
    for phi in (0, np.pi / 2, np.pi, 3 * np.pi / 2):
        for a0 in (0.02,):
            p0 = np.concatenate([theta_null, [a0 * np.cos(phi), a0 * np.sin(phi), sig2_0]])
            lo = [0.8] * 4 + [1.0, -1, -1, 0.0]
            hi = [1.2] * 4 + [20.0, 1, 1, 40.0]
            if free_nuis:
                fun = lambda p: whitened_residual(samples, p[:5], (p[5], p[6], period, p[7]))[0]
                r = optimize.least_squares(fun, p0, bounds=(lo, hi), x_scale=[0.01] * 4 + [1, 0.01, 0.01, 2])
                x = r.x
            else:
                fun = lambda p: whitened_residual(samples, theta_null, (p[0], p[1], period, p[2]))[0]
                r = optimize.least_squares(fun, p0[5:], bounds=(lo[5:], hi[5:]), x_scale=[0.01, 0.01, 2])
                x = np.concatenate([theta_null, r.x])
            chi2 = float(r.fun @ r.fun)
            if best is None or chi2 < best[0]:
                best = (chi2, x, r.jac)
    chi2, x, J = best
    try:
        cov = np.linalg.inv(J.T @ J)
    except np.linalg.LinAlgError:
        cov = np.full((len(x), len(x)), np.nan)
    return chi2, x, cov


# ------------------------------------------------- linearised scan (fast)
def scan_linear(samples, theta, lin, periods, sig2s):
    """For each (period, sig2): joint linear solve for all per-sample nuisance
    plus shared (ac, as), linearised about the null best fit (lin = null b's).
    Returns chi2[np, ns], ac, as, cov(ac,as) [np, ns, 2, 2]."""
    alphas, signl = theta[:4], theta[4]
    # precompute per-sample pieces
    pre = []
    for S, b in zip(samples, lin):
        T = template(alphas[S.zbin], signl)
        base = np.column_stack([S.M @ T, S.poly])          # nk x (1+npoly)
        model0 = base @ b                                   # null model at kout
        pre.append((S, T, base, b, model0))
    nlin = sum(1 + S.npoly for S in samples)
    ntot = sum(len(S.k) for S in samples)
    chi2 = np.zeros((len(periods), len(sig2s)))
    AC = np.zeros_like(chi2)
    AS = np.zeros_like(chi2)
    COV = np.zeros(chi2.shape + (2, 2))
    Xw = np.zeros((ntot, nlin + 2))
    dw = np.zeros(ntot)
    r0 = 0
    c0 = 0
    for S, T, base, b, model0 in pre:
        n = len(S.k)
        Xw[r0:r0 + n, c0:c0 + base.shape[1]] = S.L.T @ base
        dw[r0:r0 + n] = S.L.T @ S.P
        r0 += n
        c0 += base.shape[1]
    for ip, per in enumerate(periods):
        for js, s2 in enumerate(sig2s):
            r0 = 0
            for S, T, base, b, model0 in pre:
                n = len(S.k)
                pS = per(S) if callable(per) else per
                si, ci = ring(KIN, pS, s2)
                so, co = ring(S.k, pS, s2)
                # d model / d ac = W[si T] b0 + so * poly b  (linearised)
                cs = b[0] * (S.M @ (si * T)) + so * (S.poly @ b[1:])
                cc = b[0] * (S.M @ (ci * T)) + co * (S.poly @ b[1:])
                Xw[r0:r0 + n, nlin] = S.L.T @ cs
                Xw[r0:r0 + n, nlin + 1] = S.L.T @ cc
                r0 += n
            beta, res, rk, sv = np.linalg.lstsq(Xw, dw, rcond=None)
            r = dw - Xw @ beta
            chi2[ip, js] = r @ r
            AC[ip, js], AS[ip, js] = beta[nlin], beta[nlin + 1]
            cov = np.linalg.inv(Xw.T @ Xw)
            COV[ip, js] = cov[nlin:, nlin:]
    return chi2, AC, AS, COV


def upper_limit(ac, as_, cov, z=1.645):
    """worst-phase one-sided 95% upper limit on A2 = sqrt(ac^2+as^2)."""
    phis = np.linspace(0, 2 * np.pi, 181)
    u = np.stack([np.cos(phis), np.sin(phis)], axis=1)
    proj = u @ np.array([ac, as_])
    sig = np.sqrt(np.einsum("ij,jk,ik->i", u, cov, u))
    return float(np.max(proj + z * sig))


# ----------------------------------------------------------------- mocks
_G = {}


def _mock_worker(seed):
    samples, theta_data, periods, sig2s = _G["samples"], _G["theta"], _G["periods"], _G["sig2s"]
    rng = np.random.default_rng(seed)
    # generate
    mocks = []
    for S, mod in zip(samples, _G["model0"]):
        Sm = Sample()
        Sm.__dict__.update(S.__dict__)
        Sm.P = mod + S.Lc @ rng.standard_normal(len(S.k))
        mocks.append(Sm)
    th, _, chi2n, lin = fit_null(mocks, theta_data)
    chi2, AC, AS, COV = scan_linear(mocks, th, lin, periods, sig2s)
    dchi = chi2n - chi2
    ip, js = np.unravel_index(np.argmax(dchi), dchi.shape)
    # exact refinement at the best grid point (nuisance fixed, as pre-registered for the scan)
    c2, x, cv = fit_ring_exact(mocks, th, lin, periods[ip], sig2s[js], free_nuis=False)
    return float(dchi.max()), float(chi2n - c2), float(periods[ip]), th[4], float(np.hypot(x[5], x[6]))


def main():
    t0 = time.time()
    out = []

    def P(*a):
        s = " ".join(str(x) for x in a)
        print(s, flush=True)
        out.append(s)

    P("# ring_echo.py -- second-ringing search in the galaxy P0(k)")
    P(f"# run {time.strftime('%Y-%m-%d %H:%M')}, NMOCK={NMOCK}, seed={SEED}")
    P(f"# BAO template: {TEMPLATE_SRC}; broadband EH98 no-wiggle x deg-7 polynomial in ln k")
    samples = [load_eboss("NGC"), load_eboss("SGC")] + [load_boss(z, c) for z in (1, 2, 3) for c in ("NGC", "SGC")]
    for S in samples:
        prepare(S)
        P(f"  loaded {S.name}: {len(S.k)} bins, k = {S.k[0]:.4f}..{S.k[-1]:.4f}, Hartlap {S.hartlap:.3f}")
    ntot = sum(len(S.k) for S in samples)
    nlin = sum(1 + S.npoly for S in samples)
    P(f"  total {ntot} data points, {nlin} linear nuisance + 5 non-linear (4 alpha, Sig_nl)")

    # ---- 1. BAO check
    theta, err, chi2_null, lin = fit_null(samples)
    P("\n## 1. Null fit (A2 = 0): BAO check")
    for i, lab in enumerate(["alpha(z1, zeff 0.38)", "alpha(z2, zeff 0.51)", "alpha(z3, zeff 0.61)", "alpha(eBOSS LRG, zeff 0.70)"]):
        P(f"  {lab:28s} = {theta[i]:.4f} +- {err[i]:.4f}")
    P(f"  Sigma_nl                     = {theta[4]:.2f} +- {err[4]:.2f} Mpc/h")
    P(f"  chi2_null = {chi2_null:.2f} for {ntot} points, {nlin + 5} parameters -> dof {ntot - nlin - 5}")
    # no-BAO comparison (wiggles fully damped) for the BAO significance
    thn = theta.copy()
    thn[4] = 20.0
    r_nobao = optimize.least_squares(lambda t: whitened_residual(samples, np.concatenate([t, [60.0]]), None)[0], theta[:4],
                                     bounds=([0.8] * 4, [1.2] * 4))
    chi2_nobao = float(r_nobao.fun @ r_nobao.fun)
    P(f"  chi2 with BAO wiggles removed (Sigma_nl -> 60) = {chi2_nobao:.2f}; delta chi2(BAO) = {chi2_nobao - chi2_null:.1f}")
    model0 = []
    for S, b in zip(samples, lin):
        T = template(theta[S.zbin], theta[4])
        model0.append(np.column_stack([S.M @ T, S.poly]) @ b)
        r = S.P - model0[-1]
        P(f"    {S.name:26s} chi2 = {float((S.L.T @ r) @ (S.L.T @ r)):.1f} / {len(S.k)} bins; B^2 = {b[0]:.3g}")

    # ---- 2. scan
    periods = np.logspace(np.log10(20), np.log10(600), 61)
    sig2s = np.array([0.0, 5.0, 10.0, 20.0])
    chi2, AC, AS, COV = scan_linear(samples, theta, lin, periods, sig2s)
    dchi = chi2_null - chi2
    A2 = np.hypot(AC, AS)
    UL = np.array([[upper_limit(AC[i, j], AS[i, j], COV[i, j]) for j in range(len(sig2s))] for i in range(len(periods))])
    sigA = np.array([[np.sqrt(np.einsum("i,ij,j", np.array([AC[i, j], AS[i, j]]) / max(A2[i, j], 1e-12), COV[i, j],
                                        np.array([AC[i, j], AS[i, j]]) / max(A2[i, j], 1e-12)))
                      for j in range(len(sig2s))] for i in range(len(periods))])
    P("\n## 2. Scan (linearised, alpha and Sigma_nl fixed at null best fit)")
    P("# scale L2[Mpc/h]  k-period 2pi/L2[h/Mpc] | Sig2=0: A2 +- err  dchi2  UL95 | Sig2=5: A2 +- err dchi2 UL95 | Sig2=10: ... | Sig2=20: ...")
    for i, per in enumerate(periods):
        row = f"{per:8.1f} {2 * np.pi / per:8.4f} |"
        for j in range(len(sig2s)):
            row += f" {A2[i, j]:.4f} +- {sigA[i, j]:.4f} {dchi[i, j]:5.1f} {UL[i, j]:.4f} |"
        P(row)
    ib, jb = np.unravel_index(np.argmax(dchi), dchi.shape)
    P(f"\n  best grid point: period {periods[ib]:.1f} Mpc/h, Sig2 = {sig2s[jb]:.0f}, A2 = {A2[ib, jb]:.4f} +- {sigA[ib, jb]:.4f}, delta chi2 = {dchi[ib, jb]:.2f}")
    hints = [(periods[i], sig2s[j], dchi[i, j], A2[i, j]) for i in range(len(periods)) for j in range(len(sig2s)) if dchi[i, j] > 9]
    if hints:
        P("  grid points with delta chi2 > 9 (hints, not detections):")
        for h in hints:
            P(f"    period {h[0]:.1f} Mpc/h, Sig2 {h[1]:.0f}: dchi2 {h[2]:.1f}, A2 {h[3]:.4f}")
    else:
        P("  no grid point with delta chi2 > 9")
    # exact refinements
    P("\n  exact non-linear refinement at the best grid point:")
    c2f, xf, cvf = fit_ring_exact(samples, theta, lin, periods[ib], sig2s[jb], free_nuis=False)
    P(f"    nuisance fixed : A2 = {np.hypot(xf[5], xf[6]):.4f}, Sig2 = {xf[7]:.1f}, delta chi2 = {chi2_null - c2f:.2f}")
    c2a, xa, cva = fit_ring_exact(samples, theta, lin, periods[ib], sig2s[jb], free_nuis=True)
    P(f"    all free       : A2 = {np.hypot(xa[5], xa[6]):.4f}, Sig2 = {xa[7]:.1f}, delta chi2 = {chi2_null - c2a:.2f}; "
      f"alphas {xa[0]:.3f} {xa[1]:.3f} {xa[2]:.3f} {xa[3]:.3f}, Sig_nl {xa[4]:.2f}")
    # also exact refinement for the top-5 distinct periods with all free
    order = np.argsort(dchi.max(axis=1))[::-1]
    P("  exact refinement (all nuisance free) at the five best periods:")
    top = []
    for i in order[:5]:
        j = np.argmax(dchi[i])
        c2, x, cv = fit_ring_exact(samples, theta, lin, periods[i], sig2s[j], free_nuis=True)
        top.append((periods[i], chi2_null - c2, np.hypot(x[5], x[6]), x[7]))
        P(f"    period {periods[i]:7.1f}: grid dchi2 {dchi[i, j]:5.1f} -> exact {chi2_null - c2:5.1f}, A2 {np.hypot(x[5], x[6]):.4f}, Sig2 {x[7]:.1f}")

    # ---- 3. per-dataset scans (Sig2 = 0 and 10 only, informational)
    P("\n## 3. Per-dataset scans (linearised; columns: period, A2, dchi2, UL95 at Sig2=0 | at Sig2=10)")
    for lab, sub in (("eBOSS LRG only", samples[:2]), ("BOSS DR12 only (k<0.15)", samples[2:])):
        th_s, err_s, c2n_s, lin_s = fit_null(sub, theta)
        c2s, ACs, ASs, COVs = scan_linear(sub, th_s, lin_s, periods, np.array([0.0, 10.0]))
        d_s = c2n_s - c2s
        P(f"  {lab}: alphas {np.round(th_s[:4], 4)}, Sig_nl {th_s[4]:.2f}, chi2_null {c2n_s:.1f}, max dchi2 {d_s.max():.1f} at period {periods[np.unravel_index(np.argmax(d_s), d_s.shape)[0]]:.1f}")
        for i in range(0, len(periods), 4):
            P(f"    {periods[i]:7.1f} | {np.hypot(ACs[i, 0], ASs[i, 0]):.4f} {d_s[i, 0]:5.1f} {upper_limit(ACs[i, 0], ASs[i, 0], COVs[i, 0]):.4f} | "
              f"{np.hypot(ACs[i, 1], ASs[i, 1]):.4f} {d_s[i, 1]:5.1f} {upper_limit(ACs[i, 1], ASs[i, 1], COVs[i, 1]):.4f}")

    # ---- 4. look-elsewhere
    P("\n## 4. Global significance")
    dmax = float(dchi.max())
    p_loc = stats.chi2.sf(dmax, 2)
    p_glob = 1 - (1 - p_loc) ** 40
    P(f"  data max delta chi2 = {dmax:.2f} (2 dof local p = {p_loc:.3g}; N_eff = 40 -> p_global = {p_glob:.3f})")
    P(f"  pre-registered detection threshold delta chi2 > 25: {'MET' if dmax > 25 else 'NOT met'}")
    if NMOCK > 0:
        _G.update(samples=samples, theta=theta, periods=periods, sig2s=sig2s, model0=model0)
        with Pool(min(16, os.cpu_count())) as pool:
            res = pool.map(_mock_worker, range(SEED, SEED + NMOCK))
        res = np.array(res)
        dm, dme, per_m, snl_m, a2m = res.T
        P(f"  {NMOCK} Gaussian mocks (covariance around the null best fit, null refit + full scan each):")
        P(f"    max delta chi2 (grid): median {np.median(dm):.1f}, 68% {np.percentile(dm, 16):.1f}-{np.percentile(dm, 84):.1f}, 95th pct {np.percentile(dm, 95):.1f}, 99th {np.percentile(dm, 99):.1f}, max {dm.max():.1f}")
        P(f"    max delta chi2 (exact refit at best grid point): 95th pct {np.percentile(dme, 95):.1f}, max {dme.max():.1f}")
        P(f"    fraction of mocks with max dchi2 >= data ({dmax:.1f}): {np.mean(dm >= dmax):.3f}  <-- global p-value")
        P(f"    fraction of mocks with max dchi2 > 9: {np.mean(dm > 9):.3f};  > 25: {np.mean(dm > 25):.3f}")
        P(f"    mock Sigma_nl refits: {np.mean(snl_m):.2f} +- {np.std(snl_m):.2f} Mpc/h; best-fit A2 of the noise peak: median {np.median(a2m):.4f}")
        h, e = np.histogram(per_m, bins=np.logspace(np.log10(20), np.log10(600), 11))
        P("    where the mock noise peaks land (period bins 20-600, log): " + " ".join(f"{int(x)}" for x in h))
        np.save(os.path.join(DD, "ring_echo_mocks.npy"), res)

    # ---- 6. crest ripples locked to the BAO frequency (harmonics 2 k_BAO, 3 k_BAO)
    P("\n## 6. BAO harmonics ('crest ripples'): sinusoids at n k_BAO, n = 2, 3, with k_BAO = alpha_z 2pi/r_d per redshift bin")
    P(f"  r_d (template) = {RD_H:.2f} Mpc/h; fitted BAO scales r_d/alpha_z: " + " ".join(f"{RD_H / a:.1f}" for a in theta[:4]) + " Mpc/h")
    P("  each harmonic: free amplitude and phase, damping exp(-k^2 Sig2^2/2) with Sig2 in {0,5,10,20} Mpc/h (and the fitted Sig_nl);")
    P("  linearised about the null fit (as in the scan), alpha_z and Sig_nl fixed.  Second-order (mode-coupling) expectation: well under 1%.")
    sig2h = np.array([0.0, 5.0, 10.0, 20.0, theta[4]])
    P("# n  Sig2 |   A_n +- err   dchi2   UL95 (worst phase)")
    harm = {}
    for nh in (2, 3):
        per_fn = (lambda nh: (lambda S: nh * RD_H / theta[S.zbin]))(nh)
        c2h, ACh, ASh, COVh = scan_linear(samples, theta, lin, [per_fn], sig2h)
        for j, s2 in enumerate(sig2h):
            a = np.hypot(ACh[0, j], ASh[0, j])
            u = np.array([ACh[0, j], ASh[0, j]]) / max(a, 1e-12)
            e = np.sqrt(u @ COVh[0, j] @ u)
            d = chi2_null - c2h[0, j]
            ul = upper_limit(ACh[0, j], ASh[0, j], COVh[0, j])
            harm[(nh, j)] = (a, e, d, ul)
            P(f"  {nh}  {s2:4.1f} | {a:.4f} +- {e:.4f}  {d:5.2f}  {ul:.4f}")
    # joint 2+3 fit (4 linear parameters) at each damping: exact chi2 via the linearised design
    P("  joint n=2 + n=3 fit (4 dof):")
    for j, s2 in enumerate(sig2h):
        # build design: reuse scan_linear machinery by stacking two ring columns -> do it directly
        alphas = theta[:4]
        rows, cols = [], []
        Xw_list, dw_list = [], []
        for S, b in zip(samples, lin):
            T = template(alphas[S.zbin], theta[4])
            base = np.column_stack([S.M @ T, S.poly])
            extra = []
            for nh in (2, 3):
                pS = nh * RD_H / alphas[S.zbin]
                si, ci = ring(KIN, pS, s2)
                so, co = ring(S.k, pS, s2)
                extra.append(b[0] * (S.M @ (si * T)) + so * (S.poly @ b[1:]))
                extra.append(b[0] * (S.M @ (ci * T)) + co * (S.poly @ b[1:]))
            Xw_list.append((S.L.T @ base, S.L.T @ np.column_stack(extra)))
            dw_list.append(S.L.T @ S.P)
        nl = sum(x[0].shape[1] for x in Xw_list)
        nt = sum(len(d) for d in dw_list)
        X = np.zeros((nt, nl + 4)); dvec = np.zeros(nt); r0 = c0 = 0
        for (xb, xe), d in zip(Xw_list, dw_list):
            X[r0:r0 + len(d), c0:c0 + xb.shape[1]] = xb
            X[r0:r0 + len(d), nl:] = xe
            dvec[r0:r0 + len(d)] = d
            r0 += len(d); c0 += xb.shape[1]
        beta = np.linalg.lstsq(X, dvec, rcond=None)[0]
        rr = dvec - X @ beta
        P(f"    Sig2 {s2:4.1f}: A2 = {np.hypot(beta[nl], beta[nl + 1]):.4f}, A3 = {np.hypot(beta[nl + 2], beta[nl + 3]):.4f}, delta chi2 = {chi2_null - rr @ rr:.2f} (4 dof, p = {stats.chi2.sf(chi2_null - rr @ rr, 4):.3f})")
    dh = max(v[2] for v in harm.values())
    P(f"  max delta chi2 over harmonics and dampings = {dh:.2f}: {'above' if dh > 25 else 'below'} the detection threshold (25); "
      f"2-dof local p = {stats.chi2.sf(dh, 2):.3f} (10 trials -> p_global ~ {1 - (1 - stats.chi2.sf(dh, 2)) ** 10:.3f})")

    # ---- 5. upper-limit summary
    P("\n## 5. 95% upper-limit summary (worst phase, envelope over Sig2 in {0,5,10,20})")
    ULenv = UL.max(axis=1)
    for lo, hi in ((20, 30), (30, 60), (60, 100), (100, 200), (200, 300), (300, 600)):
        m = (periods >= lo) & (periods < hi * 1.0001)
        P(f"  periods {lo:3d}-{hi:3d} Mpc/h: A2 < {ULenv[m].min():.4f} .. {ULenv[m].max():.4f} (median {np.median(ULenv[m]):.4f})")
    P(f"\n# runtime {time.time() - t0:.0f} s")
    with open(os.path.join(HERE, "ring_echo_results.txt"), "w") as f:
        f.write("\n".join(out) + "\n")

    # figure
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        fig, ax = plt.subplots(2, 2, figsize=(11, 8))
        for S, m0 in zip(samples, model0):
            e = np.sqrt(np.diag(S.C))
            ax[0, 0].errorbar(S.k, S.k * S.P, S.k * e, fmt=".", ms=3, label=S.name)
            ax[0, 0].plot(S.k, S.k * m0, "k-", lw=0.6)
        ax[0, 0].set_xlabel("k [h/Mpc]"); ax[0, 0].set_ylabel("k P0(k)"); ax[0, 0].legend(fontsize=6); ax[0, 0].set_title("data and null (BAO-only) fits")
        for S, m0 in zip(samples[:2], model0[:2]):
            e = np.sqrt(np.diag(S.C))
            ax[0, 1].errorbar(S.k, S.P / m0 - 1, e / m0, fmt=".", label=S.name)
        ax[0, 1].axhline(0, color="k", lw=0.5); ax[0, 1].set_ylim(-0.15, 0.15)
        ax[0, 1].set_xlabel("k [h/Mpc]"); ax[0, 1].set_ylabel("P/P_null - 1"); ax[0, 1].set_title("eBOSS residuals"); ax[0, 1].legend(fontsize=7)
        for j, s2 in enumerate(sig2s):
            ax[1, 0].plot(periods, dchi[:, j], label=f"Sig2 = {s2:.0f}")
        ax[1, 0].axhline(9, color="r", ls=":", lw=0.8); ax[1, 0].axhline(25, color="r", ls="--", lw=0.8)
        ax[1, 0].set_xscale("log"); ax[1, 0].set_xlabel("period 2pi/k2 [Mpc/h]"); ax[1, 0].set_ylabel("delta chi2 vs A2=0"); ax[1, 0].legend(fontsize=7)
        ax[1, 0].axvline(147 * 0.676, color="gray", lw=0.5)
        for j, s2 in enumerate(sig2s):
            ax[1, 1].plot(periods, UL[:, j], label=f"Sig2 = {s2:.0f}")
        ax[1, 1].set_xscale("log"); ax[1, 1].set_yscale("log"); ax[1, 1].set_xlabel("period [Mpc/h]"); ax[1, 1].set_ylabel("95% upper limit on A2"); ax[1, 1].legend(fontsize=7)
        ax[1, 1].axvline(147 * 0.676, color="gray", lw=0.5)
        fig.tight_layout()
        os.makedirs(os.path.join(HERE, "figures"), exist_ok=True)
        fig.savefig(os.path.join(HERE, "figures", "ring_echo.png"), dpi=130)
    except Exception as ex:
        print("figure failed:", ex)


if __name__ == "__main__":
    main()
