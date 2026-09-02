"""Hybrid 1-D spherical toy for H8: the same shells obey different mechanics by local density.
 - every shell: gravity from the enclosed mass (radii sorted each step, so collisionless shells may cross)
 - shells whose local packing fraction eta > eta_c are FLUID: hard-sphere pressure (Carnahan-Starling with a
   jamming divergence) acts between them and an artificial viscosity captures shocks; they cannot cross
   (a jamming floor keeps eta < 0.64 and enforces ordering among fluid shells)
 - shells with eta < eta_c are BILLIARD BALLS: no pressure, free streaming through everything
Kick: v = f * v_esc,local(r) * (1 + sigma * gaussian per shell), so noise is allowed (collisionless shells cross).
Units G = M = R = 1.  Usage: python3 rings_hybrid.py f sigma eta_c N tend seed
Reports: unbound fraction (self-consistent final potential), fraction of shells ever fluid, final density
edges (fixed 0.5-dex contrast over 5 sorted shells), radial structure, and the pressureless (eta_c=1e9) null
must be run separately for comparison."""
import numpy as np, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
f = float(sys.argv[1]) if len(sys.argv) > 1 else 0.7
sig = float(sys.argv[2]) if len(sys.argv) > 2 else 0.5
eta_c = float(sys.argv[3]) if len(sys.argv) > 3 else 0.1
N = int(sys.argv[4]) if len(sys.argv) > 4 else 200
tend = float(sys.argv[5]) if len(sys.argv) > 5 else 20.0
seed = int(sys.argv[6]) if len(sys.argv) > 6 else 1
rng = np.random.default_rng(seed)
eta0 = 0.60
m = np.full(N, 1.0 / N)
r = (np.cumsum(m)) ** (1 / 3)
rho0 = 3 / (4 * np.pi); vol_per_mass = eta0 / rho0
e = np.full(N, 1e-3)                       # small thermal energy; pressure matters only via jamming
Phi0 = -(3 - r**2) / 2.0
v = f * np.sqrt(2 * np.abs(Phi0)) * (1 + sig * rng.standard_normal(N))
ever_fluid = np.zeros(N, bool)


def cs_pressure(rho, e):
    eta = rho * vol_per_mass; eta_j = 0.64
    eta_cl = np.clip(eta, 0, 0.55)
    Z = (1 + eta_cl + eta_cl**2 - eta_cl**3) / (1 - eta_cl) ** 3
    Z = np.where(eta > 0.55, Z * (eta_j - 0.55) / np.maximum(eta_j - eta, 1e-3), Z)
    return rho * (2.0 / 3.0) * e * Z


t, dt, step = 0.0, 1e-4, 0
while t < tend:
    order = np.argsort(r); rs = r[order]; vs = v[order]; ms = m[order]; es = e[order]
    Menc = np.cumsum(ms)
    r_in = np.concatenate(([0.0], rs[:-1]))
    vol = 4 / 3 * np.pi * np.maximum(rs**3 - r_in**3, 1e-12); rho = ms / vol
    eta = rho * vol_per_mass
    fluid = (eta > eta_c) & (vol > 1e-9)      # crossing shells (zero volume) are not 'fluid'
    ever_fluid[order[fluid]] = True
    P = np.where(fluid, cs_pressure(rho, es), 0.0)
    dv = vs - np.concatenate(([0.0], vs[:-1]))
    q = np.where(fluid & (dv < 0), 2.0 * rho * dv**2, 0.0)
    Pt = np.nan_to_num(P + q, nan=0.0, posinf=1e6)
    Pout = np.concatenate((Pt[1:], [0.0]))
    area = 4 * np.pi * rs**2
    mass_face = 0.5 * (ms + np.concatenate((ms[1:], [ms[-1]])))
    # pressure force only on fluid shells (and only from fluid neighbours, since P=0 elsewhere)
    acc = np.where(fluid, -(Pout - Pt) * area / mass_face, 0.0) - Menc / np.maximum(rs, 1e-4) ** 2
    vs = vs + acc * dt
    rs_new = np.maximum(rs + vs * dt, 1e-4)
    # jamming floor + ordering among fluid shells only (collisionless shells may cross freely)
    if fluid.any():
        vmin3 = ms * vol_per_mass / 0.639 * 3 / (4 * np.pi)
        idx = np.where(fluid)[0]
        C = np.cumsum(vmin3[idx]); r3 = C + np.maximum.accumulate(rs_new[idx]**3 - C)
        r_floor = r3 ** (1 / 3)
        stopped = r_floor > rs_new[idx] + 1e-12
        rs_new[idx] = np.maximum(rs_new[idx], r_floor)
        vs[idx] = np.where(stopped, np.maximum(vs[idx], 0.0), vs[idx])
    # thermal energy update for fluid shells (pdV)
    r_in_new = np.concatenate(([0.0], rs_new[:-1]))
    vol_new = 4 / 3 * np.pi * np.maximum(rs_new**3 - r_in_new**3, 1e-12)
    es = np.where(fluid, np.maximum(es - Pt * (vol_new - vol) / ms, 1e-8), es)
    # scatter back
    r[order] = rs_new; v[order] = vs; e[order] = es
    t += dt
    cs = np.sqrt(np.maximum(np.nan_to_num((5 / 3) * Pt / rho, nan=0.0, posinf=1e6), 1e-12))
    dr = np.maximum(rs_new - r_in_new, 1e-6)
    dt_c = 0.3 * np.nanmin(dr / (cs + np.abs(np.nan_to_num(vs)) + 1e-9))
    dt = min(dt_c if np.isfinite(dt_c) and dt_c > 0 else 1e-6, 5e-3)
    if not (np.isfinite(rs_new).all() and np.isfinite(vs).all()):
        print(f'NONFINITE at step {step} t {t:.4f}', file=sys.stderr); break
    step += 1

order = np.argsort(r); rs = r[order]; vs = v[order]; ms = m[order]
Menc = np.cumsum(ms)
Phi_f = -Menc / rs - (np.cumsum((ms / rs)[::-1])[::-1] - ms / rs)
Etot = 0.5 * vs**2 + e[order] + Phi_f
unbound = ms[Etot > 0].sum()
r_in = np.concatenate(([0.0], rs[:-1])); vol = 4 / 3 * np.pi * np.maximum(rs**3 - r_in**3, 1e-12)
eta = ms / vol * vol_per_mass
# binned density profile (20 equal-mass bins) so that shell crossing does not create spurious spikes
nb = 20; idx = (np.arange(N) * nb) // N
rb = np.array([rs[idx == k].max() for k in range(nb)]); rb_in = np.concatenate(([0.0], rb[:-1]))
eta_b = (np.array([ms[idx == k].sum() for k in range(nb)]) / (4 / 3 * np.pi * (rb**3 - rb_in**3))) * vol_per_mass
le = np.log10(np.maximum(eta_b, 1e-30)); grad = np.abs(np.diff(le))
edges = [int(i) for i in np.where(grad > 0.5)[0]]
Menc_b = np.cumsum([ms[idx == k].sum() for k in range(nb)])
print(f"f={f} sigma={sig} eta_c={eta_c} N={N} seed={seed} t_end={t:.1f} steps={step}  KE/|U| at start={2*f**2*(1+sig**2):.2f}")
print(f"unbound fraction = {unbound:.3f}   shells ever fluid = {ever_fluid.mean():.2f}")
print(f"edges (>0.5 dex between adjacent 5%-mass bins) = {len(edges)} at mass fractions {[round(float(Menc_b[i]),2) for i in edges]}")
for qf in (0.1, 0.3, 0.5, 0.7, 0.9, 1.0):
    i = min(int(qf * N) - 1, N - 1)
    print(f"  M<{qf:3.1f}: r={rs[i]:8.3f} eta={eta[i]:8.2e} E={Etot[i]:+7.3f} {'fluid-history' if ever_fluid[order[i]] else 'ballistic'}")
