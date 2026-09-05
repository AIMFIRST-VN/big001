"""Symbolic regression / closed-form search on two toy-model outputs (see notes/symbolic_fits.md).
Grid 1: rings_contact.py grid (contact_grid_results.txt): unbound fraction U(f, sigma) for contact runs
        (eta_c = 0.03/0.1/0.3 pooled) and for the ballistic null (eta_c = 1e9), plus their difference.
Grid 2: Johnson-Mehl tessellation samples produced by jm_cells_dump.py / jm_sizes_dump.py
        (symbolic_runs/jm_cells_samples_*.npz, symbolic_runs/jm_sizes_samples_*.npz).
Engine: a hand-written enumeration over a library of simple forms fitted by least squares with a BIC penalty,
plus leave-one-f-out cross-validation; PySR is used in addition if it imports (--pysr), otherwise skipped.
Usage: python3 symbolic_fits.py [--pysr]   -> symbolic_fits_results.txt"""
import numpy as np, os, sys, re, itertools, warnings
from scipy.optimize import least_squares
from scipy.special import erfc, gamma as Gamma, gammainc
from scipy.stats import ks_1samp, gamma as gamma_dist, lognorm, weibull_min, kstest
from scipy import stats
HERE = os.path.dirname(os.path.abspath(__file__))
warnings.filterwarnings("ignore")
out_lines = []
def P(*a):
    s = " ".join(str(x) for x in a); print(s); out_lines.append(s)

# ---------------------------------------------------------------- Grid 1
rows = []
for line in open(os.path.join(HERE, "contact_grid_results.txt")):
    m = re.match(r"([\d.]+) ([\d.]+) ([\de+.]+) \| ([\d.]+) \(([\d.]+)-([\d.]+)\)", line)
    if m:
        f, s, e, u, lo, hi = map(float, m.groups()); rows.append((f, s, e, u, lo, hi))
rows = np.array(rows)
def pool(mask):
    keys = sorted(set(map(tuple, rows[mask][:, :2])))
    F, S, U, LO, HI = [], [], [], [], []
    for f, s in keys:
        r = rows[mask & (rows[:, 0] == f) & (rows[:, 1] == s)]
        F.append(f); S.append(s); U.append(r[:, 3].mean()); LO.append(r[:, 4].min()); HI.append(r[:, 5].max())
    return map(np.array, (F, S, U, LO, HI))
Fc, Sc, Uc, LOc, HIc = pool(rows[:, 2] < 1)
Fb, Sb, Ub, LOb, HIb = pool(rows[:, 2] > 1)
assert np.allclose(Fc, Fb) and np.allclose(Sc, Sb)
halfrange_c = 0.5 * (HIc - LOc); halfrange_b = 0.5 * (HIb - LOb)
# seed scatter: half of the min-max range over 3 seeds ~ 1.0 sigma_seed (E[range] for n=3 normals = 1.69 sigma)
scat_c = np.sqrt(np.mean(((HIc - LOc) / 1.69) ** 2)); scat_b = np.sqrt(np.mean(((HIb - LOb) / 1.69) ** 2))

clip = lambda u: np.clip(u, 0, 1)
feff = lambda f, s: f * np.sqrt(1 + s**2)               # sqrt(KE/|U| / 2)
# library: name -> (function(f, s, p), initial guesses list, complexity = #params + #nonlinear ops, formula string)
LIB = {
 "logistic_sigwidth":  (lambda f, s, p: 1 / (1 + np.exp(-(f - p[0]) / (np.abs(p[1]) + np.abs(p[2]) * s))),
                        [[0.7, 0.05, 0.1]], 5, "1/(1+exp(-(f-fc)/(w0+w1 sigma)))"),
 "logistic_feff":      (lambda f, s, p: 1 / (1 + np.exp(-(feff(f, s) - p[0]) / np.abs(p[1]))),
                        [[0.8, 0.1]], 4, "1/(1+exp(-(f sqrt(1+sigma^2)-fc)/w))"),
 "erfc_escape":        (lambda f, s, p: 0.5 * erfc((p[0] - f) / np.sqrt(2 * (p[1]**2 + (p[2] * s * f)**2))),
                        [[0.8, 0.05, 1.0], [1.0, 0.1, 0.5]], 6, "0.5 erfc((fc-f)/sqrt(2(w0^2+(k sigma f)^2)))"),
 "erfc_escape_k1":     (lambda f, s, p: 0.5 * erfc((p[0] - f) / np.sqrt(2 * (p[1]**2 + (s * f)**2))),
                        [[0.8, 0.05], [1.0, 0.1]], 5, "0.5 erfc((fc-f)/sqrt(2(w0^2+sigma^2 f^2)))   [single-shell escape, k=1]"),
 "tautology_1-fbar^-3": (lambda f, s, p: clip(1 - (p[0] / feff(f, s)) ** 3),
                        [[0.8]], 3, "max(0, 1-(fc/(f sqrt(1+sigma^2)))^3)"),
 "KE_over_U_linear":   (lambda f, s, p: clip(p[0] * 2 * f**2 * (1 + s**2) + p[1]),
                        [[0.5, -0.3]], 3, "clip(a*2f^2(1+sigma^2)+b)"),
 "powerlaw_feff":      (lambda f, s, p: clip(np.abs(p[0]) * np.maximum(feff(f, s) - p[1], 0) ** np.abs(p[2])),
                        [[1.0, 0.5, 1.5], [2.0, 0.6, 2.0]], 5, "a (f sqrt(1+sigma^2)-fc)_+^n"),
 "powerlaw_f_sigma":   (lambda f, s, p: clip(np.abs(p[0]) * np.maximum(f - p[1] + p[3] * s, 0) ** np.abs(p[2])),
                        [[1.0, 0.5, 1.5, 0.1]], 6, "a (f-fc+c sigma)_+^n"),
 "exp_sat_feff":       (lambda f, s, p: 1 - np.exp(-np.abs(p[0]) * np.maximum(feff(f, s) - p[1], 0) ** np.abs(p[2])),
                        [[1.0, 0.5, 1.5]], 6, "1-exp(-a (f sqrt(1+sigma^2)-fc)_+^n)"),
 "exp_sat_KE":         (lambda f, s, p: 1 - np.exp(-np.abs(p[0]) * np.maximum(2 * f**2 * (1 + s**2) - p[1], 0) ** np.abs(p[2])),
                        [[1.0, 0.5, 1.0]], 6, "1-exp(-a (KE/|U| - q0)_+^n)"),
 "rational_feff":      (lambda f, s, p: clip(np.maximum(feff(f, s) - p[0], 0) ** 2 / (np.maximum(feff(f, s) - p[0], 0) ** 2 + p[1]**2)),
                        [[0.5, 0.3]], 4, "x^2/(x^2+b^2), x=(f sqrt(1+sigma^2)-fc)_+"),
 "poly2_f_sigma":      (lambda f, s, p: clip(p[0] + p[1] * f + p[2] * s + p[3] * f * f + p[4] * f * s + p[5] * s * s),
                        [[0, 0, 0, 1, 0, 0]], 6, "quadratic polynomial in f, sigma (6 coeff)"),
}

def fit(func, guesses, f, s, u, w=None):
    best = None
    for g in guesses:
        for jit in range(6):
            g0 = np.array(g, float) * (1 + 0.3 * np.random.default_rng(jit).standard_normal(len(g)) * (jit > 0))
            try:
                r = least_squares(lambda p: (func(f, s, p) - u) * (1 if w is None else w), g0, max_nfev=4000)
            except Exception:
                continue
            if best is None or r.cost < best.cost: best = r
    return best

def bic(rss, n, k, sig):
    return n * np.log(max(rss / n, 1e-12)) + k * np.log(n)

def study(label, F, S, U, halfrange, scat, clipf=True):
    cl = clip if clipf else (lambda u: u)
    P(f"\n=== Grid 1: {label}  (n={len(U)} grid points; seed scatter sigma_seed ~ {scat:.3f}, mean half-range {halfrange.mean():.3f}) ===")
    P(f"{'form':24s} {'k':>2s} {'cx':>3s} {'RMS':>6s} {'RMS/sc':>6s} {'LOO-RMS':>7s} {'BIC':>8s}  parameters")
    res = {}
    fvals = np.unique(F)
    for name, (func, guesses, cx, formula) in LIB.items():
        r = fit(func, guesses, F, S, U)
        if r is None: continue
        pred = cl(func(F, S, r.x)); rms = np.sqrt(np.mean((pred - U) ** 2))
        # leave-one-f-out: refit without one f column, predict it
        loo = []
        for fv in fvals:
            m = F != fv
            rr = fit(func, [r.x] + guesses, F[m], S[m], U[m])
            loo.append(cl(func(F[~m], S[~m], rr.x)) - U[~m])
        loo = np.sqrt(np.mean(np.concatenate(loo) ** 2))
        k = len(r.x); b = bic(np.sum((pred - U) ** 2), len(U), k, scat)
        res[name] = dict(rms=rms, loo=loo, bic=b, p=r.x, cx=cx, formula=formula, pred=pred)
        P(f"{name:24s} {k:2d} {cx:3d} {rms:6.3f} {rms/scat:6.2f} {loo:7.3f} {b:8.1f}  {np.round(r.x, 3).tolist()}   {formula}")
    return res

resC = study("CONTACT (eta_c=0.03/0.1/0.3 pooled)", Fc, Sc, Uc, halfrange_c, scat_c)
resB = study("BALLISTIC null (eta_c=1e9)", Fb, Sb, Ub, halfrange_b, scat_b)
D = Uc - Ub
P(f"\n=== Grid 1: DIFFERENCE U_contact - U_ballistic ===")
P("grid (rows f, cols sigma=0,0.3,0.5,1.0):")
for fv in np.unique(Fc):
    m = Fc == fv; P(f"  f={fv}: contact {np.round(Uc[m],2).tolist()}  ballistic {np.round(Ub[m],2).tolist()}  diff {np.round(D[m],2).tolist()}")
P("Ratio contact/ballistic is undefined at f<=0.5 (contact U=0) and not smooth (ballistic non-monotone in f); the")
P("difference is fitted below with the same library (no clipping to [0,1] is appropriate for a difference, values in [-1,0]).")
LIBD = {
 "poly2_unclipped":   (lambda f, s, p: p[0] + p[1] * f + p[2] * s + p[3] * f * f + p[4] * f * s + p[5] * s * s,
                       [[0, 0, 0, 0, 0, 0]], 6, "quadratic polynomial in f, sigma (6 coeff), no clipping"),
 "sigma_damped_gap":  (lambda f, s, p: -np.abs(p[0]) * np.exp(-np.abs(p[1]) * s) * (1 + p[2] * f),
                       [[0.4, 1.0, 0.0]], 4, "-A exp(-b sigma) (1 + c f)"),
 "const":             (lambda f, s, p: p[0] + 0 * f, [[-0.3]], 1, "constant"),
}
LIB_saved = LIB; LIB = LIBD
resD = study("DIFFERENCE", Fc, Sc, D, np.sqrt(halfrange_c**2 + halfrange_b**2), np.hypot(scat_c, scat_b), clipf=False)
LIB = LIB_saved

# per-column check of the contact best form: threshold f_c(sigma) from linear interpolation of U=0.1 crossing
P("\nEmpirical crossing f where U_contact = 0.10 (linear interpolation per sigma):")
for sv in np.unique(Sc):
    m = Sc == sv; ff, uu = Fc[m], Uc[m]
    idx = np.where(uu >= 0.10)[0]
    if len(idx) and idx[0] > 0:
        i = idx[0]; fx = ff[i-1] + (0.10 - uu[i-1]) * (ff[i] - ff[i-1]) / (uu[i] - uu[i-1])
        P(f"  sigma={sv}: f_0.1 = {fx:.3f}   (f sqrt(1+sigma^2) = {fx*np.sqrt(1+sv**2):.3f}, KE/|U| = {2*fx**2*(1+sv**2):.3f})")

# optional PySR
if "--pysr" in sys.argv:
    try:
        from pysr import PySRRegressor
        for label, U in (("contact", Uc), ("ballistic", Ub)):
            X = np.column_stack([Fc, Sc])
            model = PySRRegressor(niterations=60, binary_operators=["+", "-", "*", "/", "^"],
                                  unary_operators=["exp", "square", "sqrt"], maxsize=14, populations=8, procs=4,
                                  progress=False, temp_equation_file=True, verbosity=0, random_state=0, deterministic=True, parallelism="serial")
            model.fit(X, U, variable_names=["f", "sigma"])
            P(f"\n--- PySR Pareto front, {label} ---"); P(str(model.equations_[["complexity", "loss", "equation"]]))
    except Exception as e:
        P(f"\nPySR not used: {e!r}")

# ---------------------------------------------------------------- Grid 2
def ks_report(name, x, cdf, formula, npar):
    ks = kstest(x, cdf).statistic
    P(f"  {name:34s} KS = {ks:.4f}  (n={len(x)}, {npar} fitted par)   {formula}")
    return ks
cells = sorted(f for f in os.listdir(os.path.join(HERE, "symbolic_runs")) if f.startswith("jm_cells_samples"))
if cells:
    d = [np.load(os.path.join(HERE, "symbolic_runs", f)) for f in cells]
    wall = np.concatenate([x["wall"] for x in d]); offs = np.concatenate([x["offs"] for x in d])
    nnb = np.concatenate([x["nnb"] for x in d]); dts = np.concatenate([x["dts"] for x in d])
    P(f"\n=== Grid 2a: nearest-wall distance w/L ({len(wall)} realizations from {cells}) ===")
    P(f"  percentiles 10/50/90: {np.percentile(wall,[10,50,90]).round(3).tolist()}; mean {wall.mean():.3f}; s_grid resolution 0.0075")
    # candidate CDFs for w: exponential 1-exp(-w/a); Weibull 1-exp(-(w/a)^k); half-normal; gamma; 'wall ~ min over directions'
    a = wall.mean(); ks_report("exponential (a=mean)", wall, lambda x: 1 - np.exp(-x / a), f"1-exp(-w/{a:.3f})", 1)
    k, loc, lam = weibull_min.fit(wall, floc=0); ks_report("Weibull", wall, lambda x: weibull_min.cdf(x, k, 0, lam), f"1-exp(-(w/{lam:.3f})^{k:.3f})", 2)
    ga, gl, gs = gamma_dist.fit(wall, floc=0); ks_report("gamma", wall, lambda x: gamma_dist.cdf(x, ga, 0, gs), f"Gamma(shape={ga:.3f}, scale={gs:.3f})", 2)
    # 1 - (1 - w/w0)^n form (Weibull-like linear near 0)  and  1-(1+w/b)^-n (Lomax)
    r = least_squares(lambda p: np.sort(wall) - stats.lomax.ppf((np.arange(len(wall)) + 0.5) / len(wall), np.abs(p[0]), 0, np.abs(p[1])), [3, 0.3])
    ks_report("Lomax 1-(1+w/b)^-n", wall, lambda x: stats.lomax.cdf(x, abs(r.x[0]), 0, abs(r.x[1])), f"1-(1+w/{abs(r.x[1]):.3f})^-{abs(r.x[0]):.3f}", 2)
    # linear-density near zero check: fraction below 0.05 vs exponential prediction
    P(f"  P(w<0.05)={np.mean(wall<0.05):.3f} (exponential predicts {1-np.exp(-0.05/a):.3f}); P(w<0.25)={np.mean(wall<0.25):.3f}; P(w<0.5)={np.mean(wall<0.5):.3f}")
    P(f"\n=== Grid 2b: observer offset r/L ({len(offs)} realizations) ===")
    P(f"  percentiles 10/50/90: {np.percentile(offs,[10,50,90]).round(3).tolist()}; mean {offs.mean():.3f}; rms {np.sqrt(np.mean(offs**2)):.3f}")
    # Maxwell (3-D Gaussian radius): pdf ~ r^2 exp(-r^2/2s^2); generalized gamma r^2 exp(-(r/b)^k) (k=4 is the JM space-time analogue:
    # P(no seed beats the winner within r) ~ exp(-const r^4) for a uniform rate in space-time)
    s = np.sqrt(np.mean(offs**2) / 3); ks_report("Maxwell r^2 exp(-r^2/2s^2)", offs, lambda x: stats.maxwell.cdf(x, 0, s), f"s={s:.3f} (from rms)", 1)
    gg = stats.gengamma.fit(offs, floc=0); ks_report("gen-gamma r^(d-1) exp(-(r/b)^k)", offs, lambda x: stats.gengamma.cdf(x, *gg), f"a={gg[0]:.3f} (d=a*c), c=k={gg[1]:.3f}, b={gg[3]:.3f}", 3)
    # fixed-shape r^2 exp(-(r/b)^4): gengamma with a=3/4, c=4
    r4 = least_squares(lambda p: stats.gengamma.cdf(np.sort(offs), 0.75, 4, 0, abs(p[0])) - (np.arange(len(offs)) + 0.5) / len(offs), [0.8])
    ks_report("r^2 exp(-(r/b)^4)  [JM ansatz]", offs, lambda x: stats.gengamma.cdf(x, 0.75, 4, 0, abs(r4.x[0])), f"b={abs(r4.x[0]):.3f}; mean = b Gamma(1)/Gamma(3/4) = {abs(r4.x[0])*Gamma(1)/Gamma(0.75):.3f}", 1)
    # exact JM prediction, zero parameters: arrival time A at the origin has P(A>a)=exp(-pi a^4/3) (I=c=1); given A the winner
    # sits uniformly in the past light-cone |x|<A, so p(r) = 4 pi r^2 int_r^inf exp(-pi a^4/3) da = pi r^2 (3/pi)^(1/4) Gamma(1/4, pi r^4/3).
    from scipy.special import gammaincc
    from scipy.integrate import cumulative_trapezoid
    rg = np.linspace(0, 4, 4001); pr = np.pi * rg**2 * (3 / np.pi) ** 0.25 * Gamma(0.25) * gammaincc(0.25, np.pi * rg**4 / 3)
    cdf_r = cumulative_trapezoid(pr, rg, initial=0); cdf_r /= cdf_r[-1]
    ks_report("JM exact p(r) = 4 pi r^2 int_r^inf exp(-pi a^4/3) da", offs, lambda x: np.interp(x, rg, cdf_r), f"0 par; mean = {np.trapezoid(rg*pr, rg):.3f}, median = {rg[np.searchsorted(cdf_r, 0.5)]:.3f}", 0)
    r3 = least_squares(lambda p: stats.gengamma.cdf(np.sort(offs), 1.0, 3, 0, abs(p[0])) - (np.arange(len(offs)) + 0.5) / len(offs), [0.8])
    ks_report("r^2 exp(-(r/b)^3)  [Poisson-Voronoi]", offs, lambda x: stats.gengamma.cdf(x, 1.0, 3, 0, abs(r3.x[0])), f"b={abs(r3.x[0]):.3f}", 1)
    P(f"  neighbours touched: mean {nnb.mean():.2f}, sd {nnb.std():.2f}, 10/50/90 {np.percentile(nnb,[10,50,90]).round(1).tolist()}")
    P(f"  age differences: 10/50/90 {np.percentile(dts,[10,50,90]).round(3).tolist()}, fraction older {np.mean(dts<0):.3f}")
sizes = sorted(f for f in os.listdir(os.path.join(HERE, "symbolic_runs")) if f.startswith("jm_sizes_samples"))
if sizes:
    d = [np.load(os.path.join(HERE, "symbolic_runs", f)) for f in sizes]
    vol = np.concatenate([x["vol"] for x in d]); age = np.concatenate([x["age"] for x in d])
    P(f"\n=== Grid 2c: cell volume V/L^3 ({len(vol)} interior cells from {sizes}) ===")
    P(f"  mean {vol.mean():.3f} L^3 (JM analytic mean volume per nucleus for I=c=1: 1/N, N = Gamma(1/4)/(4 (pi/3)^(1/4)) -> {1/(Gamma(0.25)/(4*(np.pi/3)**0.25)):.3f} L^3; interior-box truncation biases the MC low), median {np.median(vol):.3f}, sd {vol.std():.3f}, sd/mean {vol.std()/vol.mean():.3f}")
    P(f"  percentiles 5/25/50/75/95: {np.percentile(vol,[5,25,50,75,95]).round(3).tolist()}; P(V<0.1)={np.mean(vol<0.1):.3f}; P(V<0.01)={np.mean(vol<0.01):.3f}")
    ga, gl, gs = gamma_dist.fit(vol, floc=0); ks_report("gamma (2 par)", vol, lambda x: gamma_dist.cdf(x, ga, 0, gs), f"shape={ga:.3f}, scale={gs:.3f} (mean {ga*gs:.3f})", 2)
    ks_report("Kiang gamma shape 5.6 (Poisson-Voronoi)", vol, lambda x: gamma_dist.cdf(x, 5.6, 0, vol.mean() / 5.6), "shape 5.6, scale mean/5.6", 0)
    ks_report("exponential (gamma shape 1)", vol, lambda x: 1 - np.exp(-x / vol.mean()), f"1-exp(-V/{vol.mean():.3f})", 0)
    ls = lognorm.fit(vol, floc=0); ks_report("lognormal", vol, lambda x: lognorm.cdf(x, *ls), f"sigma_ln={ls[0]:.3f}, median={ls[2]:.3f}", 2)
    wb = weibull_min.fit(vol, floc=0); ks_report("Weibull", vol, lambda x: weibull_min.cdf(x, *wb), f"k={wb[0]:.3f}, scale={wb[2]:.3f}", 2)
    # mixture-free physical form: V distribution induced by uniform birth time t through V(t): if V ~ V0 exp(-b t)... test via 2d
    P(f"\n=== Grid 2d: birth time vs final volume ===")
    bins = np.array([0, 0.2, 0.4, 0.6, 0.8, 1.0, 1.3, 1.6, 2.0, 3.0])
    tm, med, lo, hi, cnt = [], [], [], [], []
    for a_, b_ in zip(bins[:-1], bins[1:]):
        m = (age >= a_) & (age < b_)
        if m.sum() >= 5:
            tm.append(age[m].mean()); med.append(np.median(vol[m])); lo.append(np.percentile(vol[m], 25)); hi.append(np.percentile(vol[m], 75)); cnt.append(m.sum())
    tm, med, lo, hi = map(np.array, (tm, med, lo, hi))
    for a, b, c, d_, n_ in zip(tm, med, lo, hi, cnt): P(f"  t_birth ~ {a:.2f}: n={n_:4d} median V = {b:.3f}  (IQR {c:.3f}-{d_:.3f})")
    # fits on log(median): (i) exp(-b t), (ii) exp(-b t^2), (iii) exp(-b t^4) (JM untransformed fraction), (iv) power (1 - t/t*)^n, (v) A exp(-(t/tau)^k)
    y = np.log(med)
    forms = {"V0 exp(-t/tau)": (lambda p, t: np.log(abs(p[0])) - t / abs(p[1]), [2, 0.5]),
             "V0 exp(-(t/tau)^2)": (lambda p, t: np.log(abs(p[0])) - (t / abs(p[1]))**2, [2, 0.7]),
             "V0 exp(-(t/tau)^4)": (lambda p, t: np.log(abs(p[0])) - (t / abs(p[1]))**4, [2, 1.0]),
             "V0 exp(-(t/tau)^k)": (lambda p, t: np.log(abs(p[0])) - (t / abs(p[1]))**abs(p[2]), [2, 0.7, 2]),
             "V0 (1-t/t*)^n": (lambda p, t: np.log(abs(p[0])) + abs(p[2]) * np.log(np.maximum(1 - t / abs(p[1]), 1e-6)), [2, 2.5, 3]),
             "V0 exp(-pi t^4/3) [JM untransformed]": (lambda p, t: np.log(abs(p[0])) - np.pi * t**4 / 3, [2])}
    P("  fits to log(median V) vs mean birth time in bin (RMS in dex, BIC):")
    for name, (fn, g) in forms.items():
        r = least_squares(lambda p: fn(p, tm) - y, g)
        rms = np.sqrt(np.mean(r.fun**2)) / np.log(10); b = bic(np.sum(r.fun**2), len(y), len(g), 1)
        P(f"    {name:38s} p={np.round(np.abs(r.x),3).tolist()}  RMS={rms:.3f} dex  BIC={b:.1f}")
    # per-cell scatter around the median relation (log space) for the best 2-par exp(-(t/tau)^k) form
    r = least_squares(lambda p: forms["V0 exp(-(t/tau)^k)"][0](p, tm) - y, [2, 0.7, 2])
    resid = np.log10(vol) - forms["V0 exp(-(t/tau)^k)"][0](r.x, age) / np.log(10)
    P(f"  per-cell scatter of log10 V about V0 exp(-(t/tau)^k): sd = {resid.std():.3f} dex (IQR {np.percentile(resid,75)-np.percentile(resid,25):.3f} dex)")
    # implied V-distribution from uniform births in transformed-fraction weight: check untransformed-fraction weighted birth-time density
    P(f"  birth-time distribution of real seeds: 10/50/90 = {np.percentile(age,[10,50,90]).round(3).tolist()}; JM predicts density ~ exp(-pi t^4/3), median t = 0.452 (numerical; MC lower because late seeds near the box edge are excluded)")
open(os.path.join(HERE, "symbolic_fits_results.txt"), "w").write("\n".join(out_lines) + "\n")
