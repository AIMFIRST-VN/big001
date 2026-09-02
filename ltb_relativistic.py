"""Relativistic (LTB dust) recalculation of the observer's d_L(z) and w_eff.

LTB metric: ds^2 = -dt^2 + R'(r,t)^2/(1+2E(r)) dr^2 + R^2 dOmega^2, dust.
Shell eq: (dR/dt)^2 = 2M(r)/R + 2E(r)   (exact GR for dust; same form as Newton).
Seed: mass and energy profiles from the measured N-body ejecta speed
distribution (Sec 9.1), shells sorted by speed (no crossing), furnace core
mass at center. Units G=1, M_ejecta=1, initial radius R0=1.

Observables for the central observer at t0 (today, when mean R ~ 30 R0):
- inward radial null ray:      dt/dr = -R'(r,t)/sqrt(1+2E)
- redshift along the ray:      d ln(1+z)/dr = Rdot'(r,t)/sqrt(1+2E)
- luminosity distance:         d_L = (1+z)^2 R(r_emit, t_emit)
Then compare mu(z) shape against FLRW coasting and LCDM (offset-marginalized),
and extract effective q0 from the low-z expansion of d_L(z).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ---- seed profiles from measured ejecta ----
spd = []
for s in range(1, 16):
    d = np.load(f"ejecta_full_{s}.npz")
    u = d["unbound"]
    spd.append(np.linalg.norm(d["v"][u], axis=1))
v0 = np.sort(np.concatenate(spd))
NS = 400                                        # radial shells
q = np.linspace(0, 1, NS + 2)[1:-1]
vsh = np.quantile(v0, q)                        # shell speeds (sorted)
m = 1.0 / NS
Mcore = 15.0 / 85.0
M = Mcore + m * (np.arange(NS) + 1)             # M(r) enclosed
E = 0.5 * vsh ** 2 - M / 1.0                    # energy function at R0=1
# keep only unbound shells (same convention as the Newtonian toy: bound
# shells fall back and are treated as virialized interior, not added to M)
keep = E > 0
M, vsh, E = M[keep], vsh[keep], E[keep]
NS = len(M)
print(f"LTB shells kept: {NS} of 400 ({NS/4:.0f}% of ejecta mass)")

# ---- integrate shells R_i(t) on a time grid ----
NT = 4000
t_grid = np.linspace(0, 60.0, NT)
dt = t_grid[1] - t_grid[0]
R = np.empty((NT, NS)); Rd = np.empty((NT, NS))
Ri = np.full(NS, 1.0) + 1e-9 * np.arange(NS)
Vi = vsh.copy()
for i, t in enumerate(t_grid):
    R[i] = Ri; Rd[i] = Vi
    acc = -M / Ri ** 2
    Vi = Vi + acc * dt
    Ri = Ri + Vi * dt
# label shells by comoving coordinate r = index; derivatives w.r.t. r via gradient
def interp_t(A, t):  # A on t_grid -> row at time t
    i = np.clip(np.searchsorted(t_grid, t) - 1, 0, NT - 2)
    w = (t - t_grid[i]) / dt
    return A[i] * (1 - w) + A[i + 1] * w

# ---- null geodesic from center at t0 ----
t0 = 55.0
sq = np.sqrt(1 + 2 * E)
rs = np.arange(NS).astype(float)
zs, ts, ds = [0.0], [t0], [0.0]
lnz = 0.0; t = t0
for j in range(1, NS):
    Rrow0, Rdrow0 = interp_t(R, t), interp_t(Rd, t)
    Rp = np.gradient(Rrow0)       # dR/dr with dr=1
    Rdp = np.gradient(Rdrow0)
    dtdr = -Rp[j] / sq[j]
    t = t + dtdr                  # dr = 1 step outward into the past lightcone? inward ray: we go from r=0 outward in r, backwards in t
    lnz = lnz + Rdp[j] / sq[j]
    if t <= t_grid[0] + 1: break
    zs.append(np.exp(lnz) - 1.0); ts.append(t)
    ds.append((1 + zs[-1]) ** 2 * interp_t(R, t)[j])
zs, ds = np.array(zs[1:]), np.array(ds[1:])
sel = (zs > 1e-3) & (zs < 3.0) & (ds > 0)
zs, ds = zs[sel], ds[sel]
print(f"lightcone: {len(zs)} points, z up to {zs.max():.2f}")

# ---- effective q0 from low-z expansion: d_L = z/H0 (1 + (1-q0)/2 z + ...) ----
lo = zs < 0.35
A = np.vstack([zs[lo], zs[lo] ** 2]).T
c1, c2 = np.linalg.lstsq(A, ds[lo], rcond=None)[0]
q0_eff = 1 - 2 * c2 / c1
print(f"effective q0 (relativistic LTB lightcone) = {q0_eff:+.3f}")
print(f"  -> w_eff(z~0) = {(2*q0_eff-1)/3:+.3f}   [q0<0 = acceleration; LCDM ~ -0.55, coasting = 0]")

# ---- shape comparison vs FLRW models (offset-marginalized in mu) ----
mu_ltb = 5 * np.log10(ds)
def mu_flrw(model):
    out = []
    for zv in zs:
        zg = np.linspace(0, zv, 200)
        Ez = np.sqrt(0.334*(1+zg)**3 + 0.666) if model == "lcdm" else (1 + zg)
        out.append(5*np.log10((1+zv)*np.trapezoid(1/Ez, zg)))
    return np.array(out)
for model in ("lcdm", "coast"):
    mm = mu_flrw(model)
    off = np.mean(mu_ltb - mm)
    rms = np.sqrt(np.mean((mu_ltb - mm - off) ** 2))
    print(f"shape residual vs {model:6s}: rms = {rms:.4f} mag over z<{zs.max():.1f}")

plt.figure(figsize=(7, 4.5))
mmc = mu_flrw("coast"); offc = np.mean(mu_ltb - mmc)
mml = mu_flrw("lcdm"); offl = np.mean(mu_ltb - mml)
plt.plot(zs, mu_ltb - mmc - offc, "b.", ms=2, label="LTB lightcone (this work) − coasting")
plt.plot(zs, mml + offl - mmc - offc - np.mean(mml+offl-mmc-offc)*0, "r-", lw=1, label=r"$\Lambda$CDM − coasting")
plt.axhline(0, color="0.6", lw=0.5)
plt.xlabel("z"); plt.ylabel(r"$\Delta\mu$ [mag]")
plt.title("Relativistic LTB distance modulus vs FLRW references")
plt.legend(fontsize=8); plt.tight_layout(); plt.savefig("ltb_wz.png", dpi=150)
