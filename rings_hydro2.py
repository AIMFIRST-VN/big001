"""Controlled version of the 1-D rings toy (replaces rings_hydro.py, whose Table-2 numbers were
shown to be the kick profile 1 - fbar^-3 and whose band edges were a median-threshold artefact).

Changes: (1) kick relative to the LOCAL escape speed, independent of radius:
            v(r) = f * sqrt(2|Phi(r)|),  Phi = -GM(3R^2 - r^2)/(2R^3)  (uniform sphere)
         so the analytic ballistic null is: all bound for f < 1, all unbound for f > 1.
         Any intermediate structure is then produced by pressure/waves, not by the kick.
         (2) thermal energy set by a pressure-to-gravity ratio beta = P0 / (G M^2 / R^4).
         (3) band edges defined by a fixed physical contrast: |d log10(eta)| > 0.5 dex over 5 shells.
         (4) resolution N and seed from the command line; ballistic null reported alongside.
Usage: python3 rings_hydro2.py f beta N tend
"""
import numpy as np, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
f = float(sys.argv[1]) if len(sys.argv) > 1 else 0.8
beta = float(sys.argv[2]) if len(sys.argv) > 2 else 0.3
N = int(sys.argv[3]) if len(sys.argv) > 3 else 100
tend = float(sys.argv[4]) if len(sys.argv) > 4 else 20.0
eta0 = 0.60
m = np.full(N, 1.0 / N); Menc = np.cumsum(m)
r = Menc ** (1 / 3)
rho0 = 3 / (4 * np.pi); vol_per_mass = eta0 / rho0
d_sphere = (eta0 / rho0) ** (1 / 3)


def cs_pressure(rho, e):
    eta = rho * vol_per_mass; eta_j = 0.64
    eta_c = np.clip(eta, 0, 0.55)
    Z = (1 + eta_c + eta_c**2 - eta_c**3) / (1 - eta_c) ** 3
    Z = np.where(eta > 0.55, Z * (eta_j - 0.55) / np.maximum(eta_j - eta, 1e-3), Z)
    return rho * (2.0 / 3.0) * e * Z, Z


# thermal energy from beta: P0 = rho0*(2/3)*e*Z(eta0) = beta * G M^2/R^4 (=beta in code units)
_, Z0 = cs_pressure(np.array([rho0]), np.array([1.0]))
e = np.full(N, beta / (rho0 * (2 / 3) * Z0[0]))
# kick relative to local escape speed
Phi = -(3 - r**2) / 2.0
v = f * np.sqrt(2 * np.abs(Phi))
t, dt = 0.0, 1e-4
step, passages, last_sign = 0, 0, 0
rec = []
while t < tend:
    r_in = np.concatenate(([0.0], r[:-1]))
    vol = 4 / 3 * np.pi * (r**3 - r_in**3); rho = m / vol
    P, Z = cs_pressure(rho, e)
    dv = v - np.concatenate(([0.0], v[:-1]))
    q = np.where(dv < 0, 2.0 * rho * dv**2, 0.0); Pt = P + q
    Pout = np.concatenate((Pt[1:], [0.0]))
    area = 4 * np.pi * r**2
    mass_face = 0.5 * (m + np.concatenate((m[1:], [m[-1]])))
    acc = -(Pout - Pt) * area / mass_face - Menc / r**2
    v = v + acc * dt
    r_new = np.maximum(r + v * dt, 1e-4)
    vmin3 = m * vol_per_mass / 0.639 * 3 / (4 * np.pi)
    r_prev = 0.0
    for k in range(N):
        r_k = max(r_new[k], (r_prev**3 + vmin3[k]) ** (1 / 3)); r_new[k] = r_k; r_prev = r_k
    v = np.where(r_new > r + v * dt + 1e-12, np.maximum(v, 0.0), v)
    r = r_new
    r_in = np.concatenate(([0.0], r[:-1])); vol_new = 4 / 3 * np.pi * (r**3 - r_in**3)
    e = np.maximum(e - Pt * (vol_new - vol) / m, 1e-8)
    t += dt
    cs = np.sqrt(np.maximum((5 / 3) * Pt / rho, 1e-12))
    dt = min(0.3 * np.min((r - r_in) / (cs + np.abs(v) + 1e-9)), 5e-3)
    step += 1
    if step % 200 == 0:
        s = np.sign(v[int(0.1 * N)])
        if s != 0 and s != last_sign and last_sign != 0: passages += 1
        last_sign = s

r_in = np.concatenate(([0.0], r[:-1])); vol = 4 / 3 * np.pi * (r**3 - r_in**3)
rho = m / vol; eta = rho * vol_per_mass
# self-consistent potential of the final configuration: Phi(r_i) = -Menc_i/r_i - sum_{k>i} m_k/r_k
Phi_f = -Menc / r - (np.cumsum((m / r)[::-1])[::-1] - m / r)
Etot = 0.5 * v**2 + e + Phi_f
unbound = m[Etot > 0].sum()
# ballistic null: a collective kick v = f*sqrt(2|Phi|) gives total KE = 2 f^2 |U| (since sum m|Phi| = 2|U|),
# so the whole system is unbound for f > 1/sqrt(2) = 0.707; the escaping fraction below that needs the
# pressureless run (beta -> 0) of this same code, which is reported as the null by running beta=1e-6.
null = 2 * f**2   # total kinetic / |U|; > 1 means E_total > 0
le = np.log10(np.maximum(eta, 1e-30))
grad = np.abs(le[5:] - le[:-5])          # contrast over 5 shells
edges = np.where(grad > 0.5)[0]
# merge adjacent edge indices into distinct edges
distinct = []
for i in edges:
    if not distinct or i - distinct[-1] > 5: distinct.append(int(i))
print(f"f_local={f} beta={beta} N={N} t_end={t:.1f} steps={step}")
print(f"unbound fraction = {unbound:.3f}   total KE/|U| = {null:.2f} (E_tot>0 if >1); pressureless run of this code is the ballistic null")
print(f"wave passages through core = {passages}")
print(f"physical edges (>0.5 dex over 5 shells) = {len(distinct)} at mass fractions {[round(float(Menc[i+2]),2) for i in distinct]}")
for qf in (0.1, 0.3, 0.5, 0.7, 0.9, 1.0):
    i = min(int(qf * N) - 1, N - 1)
    print(f"  M<{qf:3.1f}: r={r[i]:7.3f} eta={eta[i]:8.2e} E={Etot[i]:+8.3f}")
