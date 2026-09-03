"""Energy-conserving hybrid rebound toy (replaces rings_hybrid.py, withdrawn in Sec. 3.3 after an energy blow-up).
Spherical shells of equal mass with radial velocity v and tangential specific angular momentum j (Henon shells):
    r'' = -G M_enc(r)/r^2 + j^2/r^3            (softened at small r)
integrated with fixed-step leapfrog (symplectic; energy drift is reported and must stay < 1%).
Contact: a shell has a radial thickness d_c(r) = m*vol_per_mass/(eta_j*4*pi*r^2) at random close packing eta_j=0.64.
  - where the local packing fraction eta > eta_c ("fluid"): two adjacent shells that approach within d_c collide
    ELASTICALLY: radial velocities are exchanged (equal masses -> exact KE conservation) and the pair is isotropized:
    each shell's (v_r, v_t) is rotated by a random angle, conserving its kinetic energy. Collisions transfer
    radial energy into tangential motion = pressure. No artificial viscosity, no jamming floor.
  - where eta < eta_c ("billiard balls"): shells cross freely and keep their j.
Kick: v_r = f * v_esc,local(r) * (1 + sigma*xi), xi ~ N(0,1) per shell; initial j from a small isotropic
thermal speed v_t0 = t0 * v_esc,local.
Usage: python3 rings_contact.py f sigma eta_c N tend seed [t0]
Output: unbound fraction (final self-consistent potential), energy drift, fraction of shells that ever collided,
density edges (>0.5 dex between adjacent 5%-mass bins), radial structure."""
import numpy as np, sys
f = float(sys.argv[1]); sig = float(sys.argv[2]); eta_c = float(sys.argv[3])
N = int(sys.argv[4]); tend = float(sys.argv[5]); seed = int(sys.argv[6])
t0 = float(sys.argv[7]) if len(sys.argv) > 7 else 0.05
rng = np.random.default_rng(seed)
eta0, eta_j, eps = 0.60, 0.64, 0.02
m = 1.0 / N
r = (np.arange(1, N + 1) / N) ** (1 / 3)
rho0 = 3 / (4 * np.pi); vol_per_mass = eta0 / rho0
Phi0 = -(3 - r**2) / 2.0
vesc = np.sqrt(2 * np.abs(Phi0))
v = f * vesc * (1 + sig * rng.standard_normal(N))
j = t0 * vesc * r
ever = np.zeros(N, bool)

def accel(r, j):
    order = np.argsort(r); Menc = np.empty(N); Menc[order] = m * np.arange(0, N)          # mass strictly inside (consistent with the pairwise potential)
    return -Menc * r / (r**2 + eps**2) ** 1.5 + j**2 * r / (r**2 + eps**2) ** 2, order, Menc

def pot_energy_pairwise(r):
    order = np.argsort(r); rs = r[order]; inv = 1 / np.sqrt(rs**2 + eps**2)
    # U = -sum_{i<k} m^2 / r_k  (outer shell k feels inner masses; softened)
    return -m * m * np.sum(np.arange(0, N) * inv)

def total_energy(r, v, j):
    return 0.5 * m * np.sum(v**2 + j**2 / (r**2 + eps**2)) + pot_energy_pairwise(r)

dt = 5e-5
a, order, Menc = accel(r, j)
E0 = total_energy(r, v, j)
t = 0.0; step = 0; ncoll = 0
while t < tend:
    v += 0.5 * dt * a
    r += dt * v
    neg = r < 0; r[neg] = -r[neg]; v[neg] = -v[neg]      # pass through the centre
    a, order, Menc = accel(r, j)
    v += 0.5 * dt * a
    # contact pass on sorted neighbours
    rs = r[order]; vs = v[order]; js = j[order]
    r_in = np.concatenate(([0.0], rs[:-1]))
    vol = 4 / 3 * np.pi * np.maximum(rs**3 - r_in**3, 1e-12)
    eta = m / vol * vol_per_mass
    d_c = m * vol_per_mass / (eta_j * 4 * np.pi * rs**2)
    for parity in (0, 1):
        i = np.arange(parity, N - 1, 2)
        close = (rs[i + 1] - rs[i] < d_c[i]) & (vs[i] > vs[i + 1]) & (np.maximum(eta[i], eta[i + 1]) > eta_c) & (eta_c < 1e8)
        k = i[close]
        if k.size:
            vs[k], vs[k + 1] = vs[k + 1].copy(), vs[k].copy()
            for kk in (k, k + 1):
                th = rng.uniform(0, 2 * np.pi, kk.size)
                vt = js[kk] / np.sqrt(rs[kk]**2 + eps**2); vr = vs[kk]
                vs[kk] = vr * np.cos(th) - vt * np.sin(th); js[kk] = np.abs(vr * np.sin(th) + vt * np.cos(th)) * np.sqrt(rs[kk]**2 + eps**2)
            ever[order[k]] = True; ever[order[k + 1]] = True; ncoll += k.size
    v[order] = vs; j[order] = js
    t += dt; step += 1

E1 = total_energy(r, v, j)
order = np.argsort(r); rs = r[order]; Menc = m * np.arange(0, N)
inv = 1 / np.sqrt(rs**2 + eps**2)
Phi = -Menc * inv - (np.cumsum((m * inv)[::-1])[::-1] - m * inv)
Etot = 0.5 * (v[order]**2 + j[order]**2 / (rs**2 + eps**2)) + Phi
unbound = m * np.sum(Etot > 0)
nb = 20; idx = (np.arange(N) * nb) // N
rb = np.array([rs[idx == k].max() for k in range(nb)]); rb_in = np.concatenate(([0.0], rb[:-1]))
eta_b = (m * np.bincount(idx) / (4 / 3 * np.pi * (rb**3 - rb_in**3))) * vol_per_mass
grad = np.abs(np.diff(np.log10(np.maximum(eta_b, 1e-30)))); edges = [int(i) for i in np.where(grad > 0.5)[0]]
Menc_b = np.cumsum(m * np.bincount(idx))
print(f"f={f} sigma={sig} eta_c={eta_c} N={N} seed={seed} t0={t0} t_end={t:.1f} steps={step} KE/|U|0={f**2*(1+sig**2)*2:.2f}")
print(f"unbound fraction = {unbound:.3f}   energy drift (E1-E0)/|E0| = {(E1-E0)/abs(E0):+.4f}   shells ever collided = {ever.mean():.2f}  collisions = {ncoll}")
print(f"edges (>0.5 dex between adjacent 5%-mass bins) = {len(edges)} at mass fractions {[round(float(Menc_b[i]),2) for i in edges]}")
for qf in (0.1, 0.3, 0.5, 0.7, 0.9, 1.0):
    i = min(int(qf * N) - 1, N - 1)
    print(f"  M<{qf:3.1f}: r={rs[i]:8.3f} eta_bin={eta_b[min(int(qf*nb)-1,nb-1)]:8.2e} E={Etot[i]:+7.3f} {'collided' if ever[order[i]] else 'ballistic'}")
