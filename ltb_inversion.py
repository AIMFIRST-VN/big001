"""Inversion: which ejection profiles give an accelerating (q0<0) LTB lightcone?

Profile family for shell speeds v(x), x = mass fraction coordinate in (0,1):
  v(x) = v_min * [ (1-x)^(-alpha) ]            smooth power-law rise outward
         + tail: outermost f_tail of mass boosted by factor B
  lumpiness: multiply each shell's v by lognormal noise, sigma_lump
M(r): interior mass Mint + cumulative shell mass (shells re-sorted by v after
noise so no crossing at t=0).
LTB lightcone -> q0_eff. Grid over (alpha, f_tail, B, sigma_lump) x seeds.
"""
import numpy as np

NS, NT, TMAX, T0 = 300, 6000, 300.0, 280.0
Mint = 15/85 + 0.91
vmin = np.sqrt(2 * Mint) * 1.02

def q0_of_profile(vsh):
    vsh = np.sort(vsh)
    mk = 1.0 / NS * 0.09        # kept ejecta mass fraction (same units as before)
    M = Mint + mk * (np.arange(NS) + 1)
    E = 0.5 * vsh ** 2 - M
    if (E <= 0).any():
        sel = E > 0
        M, vsh, E = M[sel], vsh[sel], E[sel]
        if len(M) < 50: return np.nan
    n = len(M)
    tg = np.linspace(0, TMAX, NT); dt = tg[1] - tg[0]
    R = np.empty((NT, n)); Rd = np.empty((NT, n))
    Ri = np.full(n, 1.0) + 1e-9 * np.arange(n); Vi = vsh.copy()
    for i in range(NT):
        R[i] = Ri; Rd[i] = Vi
        Vi = Vi - M / Ri ** 2 * dt; Ri = Ri + Vi * dt
    sq = np.sqrt(1 + 2 * E)
    def it(A, t):
        i = np.clip(np.searchsorted(tg, t) - 1, 0, NT - 2); w = (t - tg[i]) / dt
        return A[i] * (1 - w) + A[i + 1] * w
    lnz, t = 0.0, T0
    zs, ds = [], []
    for j in range(1, n):
        Rr, Rdr = it(R, t), it(Rd, t)
        Rp = np.gradient(Rr); Rdp = np.gradient(Rdr)
        t -= Rp[j] / sq[j]; lnz += Rdp[j] / sq[j]
        if t <= 1: break
        z = np.exp(lnz) - 1
        zs.append(z); ds.append((1 + z) ** 2 * it(R, t)[j])
    zs, ds = np.array(zs), np.array(ds)
    lo = zs < 0.4
    if lo.sum() < 6: return np.nan
    A = np.vstack([zs[lo], zs[lo] ** 2]).T
    c1, c2 = np.linalg.lstsq(A, ds[lo], rcond=None)[0]
    return 1 - 2 * c2 / c1

rng = np.random.default_rng(7)
x = (np.arange(NS) + 0.5) / NS
rows = []
for alpha in (0.05, 0.15, 0.3):
    for f_tail, B in ((0.0, 1.0), (0.1, 1.5), (0.1, 2.5), (0.3, 2.0)):
        for slump in (0.0, 0.15, 0.4):
            q0s = []
            for seed in range(3 if slump > 0 else 1):
                r2 = np.random.default_rng(100 + seed)
                v = vmin * (1 - 0.999 * x) ** (-alpha)
                nt = int(f_tail * NS)
                if nt: v[-nt:] *= B
                if slump: v = v * r2.lognormal(0, slump, NS)
                q0s.append(q0_of_profile(v))
            q0s = [q for q in q0s if np.isfinite(q)]
            if q0s:
                rows.append((alpha, f_tail, B, slump, np.mean(q0s)))
                print(f"alpha={alpha:4.2f} tail={f_tail:3.1f}x{B:3.1f} lump={slump:4.2f}  q0={np.mean(q0s):+7.3f}")
best = min(rows, key=lambda r: r[4])
print("\nmost accelerating:", best)
