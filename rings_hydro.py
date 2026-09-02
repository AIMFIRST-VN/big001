"""1-D spherical Lagrangian hydrodynamics toy for the rings (H8).
Medium: hard-sphere gas (Carnahan-Starling pressure) with self-gravity, adiabatic index 5/3,
von Neumann-Richtmyer artificial viscosity for shocks. Units G = M_total = R_initial = 1.
Initial state: jammed ball (packing eta0 = 0.60) with small random thermal energy and an outward
radial kick v = fbar * v_esc(r) * (r/R) (1 + sigma_rel * gaussian per shell).
Outputs: energy partition (bound core / unbound), core compression history, number of
shock/wave passages through the core, and the final radial structure (density, sound speed,
collisionality proxy) to see whether kernel / mantle / shell separate with sharp edges."""
import numpy as np, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
np.random.seed(int(sys.argv[3]) if len(sys.argv) > 3 else 1)
fbar = float(sys.argv[1]) if len(sys.argv) > 1 else 0.7
sig = float(sys.argv[2]) if len(sys.argv) > 2 else 1.0
N = int(os.environ.get('NSHELL', 300))
eta0 = 0.60
# shell masses equal; initial radii from uniform density
m = np.full(N, 1.0 / N)
Menc = np.cumsum(m)
r = (Menc) ** (1 / 3)              # outer radius of each shell, R=1
rho0 = 3 / (4 * np.pi)
d_sphere = (eta0 / rho0) ** (1 / 3)   # per-unit-mass sphere volume scale so that eta = rho * vol_per_mass
vol_per_mass = eta0 / rho0            # eta = rho * vol_per_mass


def cs_pressure(rho, e):
    """Carnahan-Starling: P = rho*T*(1+eta+eta^2-eta^3)/(1-eta)^3 with T = (gamma-1)e, gamma=5/3."""
    eta = rho * vol_per_mass
    eta_j = 0.64                      # random close packing: pressure diverges here (jamming)
    eta_c = np.clip(eta, 0, 0.55)
    Z = (1 + eta_c + eta_c**2 - eta_c**3) / (1 - eta_c) ** 3
    # free-volume divergence approaching jamming: Z -> Z(0.55) * (eta_j-0.55)/(eta_j-eta)
    near = eta > 0.55
    Z = np.where(near, Z * (eta_j - 0.55) / np.maximum(eta_j - eta, 1e-3), Z)
    return rho * (2.0 / 3.0) * e * Z, Z


# thermal energy: small
e = np.full(N, 1e-3)
# kick
r_mid = np.concatenate(([0.0], r[:-1])); r_mid = 0.5 * (r_mid + r)
vesc = np.sqrt(2 * Menc / r)
v = fbar * vesc * (r / 1.0) * (1 + sig * np.random.randn(N))
v = np.maximum(v, -0.5 * vesc)
E0 = 0.5 * np.sum(m * v**2) + np.sum(m * e) - np.sum(m * Menc / r) * 0.6  # rough
t, dt, tend = 0.0, 1e-4, float(sys.argv[4]) if len(sys.argv) > 4 else 60.0
core_hist, ring_hist, passages, last_sign = [], [], 0, 0
rho = m / (4 / 3 * np.pi * (r**3 - np.concatenate(([0.0], r[:-1]))**3))
step = 0
while t < tend:
    r_in = np.concatenate(([0.0], r[:-1]))
    vol = 4 / 3 * np.pi * (r**3 - r_in**3)
    rho = m / vol
    P, Z = cs_pressure(rho, e)
    # artificial viscosity on compression
    dv = v - np.concatenate(([0.0], v[:-1]))
    q = np.where(dv < 0, 2.0 * rho * dv**2, 0.0)
    Pt = P + q
    # pressure gradient across shell boundaries (outer boundary of shell i between shell i and i+1)
    Pout = np.concatenate((Pt[1:], [0.0]))
    area = 4 * np.pi * r**2
    mass_face = 0.5 * (m + np.concatenate((m[1:], [m[-1]])))
    acc = -(Pout - Pt) * area / mass_face - Menc / r**2
    v = v + acc * dt
    r_new = r + v * dt
    # enforce ordering (no shell crossing)
    r_new = np.maximum(r_new, 1e-4)
    vmin3 = m * vol_per_mass / 0.639 * 3 / (4 * np.pi)     # r^3 increment for a shell at the jamming floor
    r_prev = 0.0
    for k in range(N):                                    # sequential inside-out floor: ordering + minimum volume
        r_k = max(r_new[k], (r_prev**3 + vmin3[k]) ** (1 / 3))
        r_new[k] = r_k; r_prev = r_k
    v = np.where(r_new > r + v * dt + 1e-12, np.maximum(v, 0.0), v)   # a shell stopped by jamming loses inward velocity
    r = r_new
    r_in = np.concatenate(([0.0], r[:-1]))
    vol_new = 4 / 3 * np.pi * (r**3 - r_in**3)
    e = np.maximum(e - Pt * (vol_new - vol) / m, 1e-8)
    t += dt
    # adaptive dt: CFL
    cs = np.sqrt(np.maximum((5 / 3) * Pt / rho, 1e-12))
    dr = r - r_in
    dt = min(0.3 * np.min(dr / (cs + np.abs(v) + 1e-9)), 5e-3)
    step += 1
    if step % 5000 == 0:
        print(f'progress t={t:.3f} dt={dt:.2e} step={step}', file=sys.stderr, flush=True)
    if step % 200 == 0:
        # core = innermost 10% of mass
        core_r = r[int(0.1 * N)]
        core_hist.append((t, core_r, np.mean(rho[:int(0.1 * N)] * vol_per_mass), np.mean(e[:int(0.1 * N)])))
        jj = np.abs(np.diff(np.log10(np.maximum(rho * vol_per_mass, 1e-30))))
        ring_hist.append((t, int(np.sum(jj > 0.3)) + 1))
        # count reversals of the core boundary velocity (wave passages)
        s = np.sign(v[int(0.1 * N)])
        if s != 0 and s != last_sign and last_sign != 0:
            passages += 1
        last_sign = s

# final diagnostics
r_in = np.concatenate(([0.0], r[:-1]))
vol = 4 / 3 * np.pi * (r**3 - r_in**3)
rho = m / vol
eta = rho * vol_per_mass
P, Z = cs_pressure(rho, e)
cs = np.sqrt((5 / 3) * P / rho)
Etot = 0.5 * v**2 + e - Menc / r   # specific energy (potential approx by enclosed mass)
bound = Etot < 0
# collisionality proxy: mean free path / radius  (lambda ~ 1/(n sigma); with sigma ~ d^2, n ~ eta/d^3 -> lambda ~ d/eta)
lam_over_r = (d_sphere / np.maximum(eta, 1e-12)) / r
print(f"fbar={fbar} sigma_rel={sig}  t_end={t:.1f} dyn times  steps={step}")
print(f"bound mass fraction (kernel+mantle) = {m[bound].sum():.3f}   unbound (shell) = {m[~bound].sum():.3f}")
print(f"core-boundary velocity reversals (wave passages through core) = {passages}")
print(f"core packing history (t, r_core, eta_core, e_core): first/last: {core_hist[0]}  {core_hist[-1]}")
print("ring count over time (t, rings):", [ring_hist[k] for k in range(0, len(ring_hist), max(1, len(ring_hist)//8))])
print("radial structure (mass shells at 5%,10%,...,100%): eta, sound speed, mean-free-path/r, bound?")
for q in np.linspace(0.05, 1.0, 20):
    i = min(int(q * N) - 1, N - 1)
    print(f"  M<{q:4.2f}: r={r[i]:8.3f} eta={eta[i]:8.2e} cs={cs[i]:8.2e} lam/r={lam_over_r[i]:8.2e} {'BOUND' if bound[i] else 'free'}")
# sharpness of boundaries: largest log-density jumps between adjacent shells
jumps = np.abs(np.diff(np.log10(np.maximum(eta, 1e-30))))
sig_j = np.where(jumps > 0.3)[0]
print(f"number of sharp density steps (>0.3 dex between adjacent shells) = {len(sig_j)}  -> rings = {len(sig_j)+1}")
print("steps at mass fractions:", [round(float(Menc[j]), 3) for j in sig_j])
# coarse-grained ring count: smooth log-eta over 5 shells, count alternating dense/rarefied bands (contrast > 1 dex)
le = np.convolve(np.log10(np.maximum(eta, 1e-30)), np.ones(5)/5, mode='same')
med = np.median(le[:int(0.9*N)])
band = (le > med).astype(int)
changes = np.where(np.diff(band) != 0)[0]
bands = len(changes) + 1
print(f"coarse-grained bands (dense/rarefied alternation, 5-shell smoothing) = {bands}; edges at mass fractions {[round(float(Menc[c]),2) for c in changes]}")
top = np.argsort(jumps)[-5:][::-1]
print("largest five density jumps (mass fraction, dex):", [(round(float(Menc[j]), 3), round(float(jumps[j]), 2)) for j in top])
np.savez(os.path.join(HERE, f"rings_hydro_f{fbar}_s{sig}.npz"), r=r, v=v, e=e, eta=eta, cs=cs, bound=bound, core_hist=np.array(core_hist))
