"""Toy FLRW bounce model for the Big Pool Bangs paper.

Compares two mechanisms for halting collapse at the jamming density:
  (A) Classical GR + jamming (Salsburg-Wood) hard-sphere pressure:
      H^2 = (8piG/3) rho          -> pressure gravitates, NO bounce.
  (B) LQC-style modified Friedmann equation:
      H^2 = (8piG/3) rho (1 - rho/rho_c) -> bounce at rho = rho_c.

Planck units internally (G = hbar = c = 1), rho in units of the critical
(jamming) density rho_c = 0.64 * rho_P. Matter treated as the hard-sphere
fluid: rho evolves via continuity drho/dt = -3H(rho + P).
"""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ETA_RCP = 0.64          # random close packing
RHO_C   = ETA_RCP       # jamming density in Planck units (rho_P = 1)
V_TH2   = 1e-2          # thermal velocity^2 / c^2 (sets kinetic pressure scale)

def eta(rho):
    return np.clip(rho, 0, None) * ETA_RCP / RHO_C * ETA_RCP  # packing fraction ~ rho scaled
# simpler: eta proportional to rho, eta=0.64 at rho=RHO_C
def eta_of(rho):
    return ETA_RCP * rho / RHO_C

def P_jam(rho):
    """Salsburg-Wood free-volume EOS: P = n kT * d/(1 - (eta/eta_rcp)^(1/3)) approx.
    Diverges as eta -> eta_rcp. Kinetic (positive) pressure."""
    x = np.clip(eta_of(rho) / ETA_RCP, 0, 1 - 1e-9)
    return rho * V_TH2 * 3.0 / (1.0 - x ** (1.0 / 3.0))

def evolve(modified, t_end=80.0, dt=1e-4, rho0=1e-3 * RHO_C):
    """Integrate collapsing FLRW from low density; return t, a, rho."""
    rho, a, H = rho0, 1.0, None
    # start collapsing: H = -sqrt((8pi/3) rho) (classical branch valid at low rho)
    ts, as_, rhos = [], [], []
    t = 0.0
    sign = -1.0
    while t < t_end and a > 1e-6 and rho < 50 * RHO_C:
        P = P_jam(rho)
        f = rho * (1 - rho / RHO_C) if modified else rho
        if f < 0:
            f = 0.0
        H = sign * np.sqrt(8 * np.pi / 3 * f)
        # Raychaudhuri to detect turnaround in modified case:
        # in LQC, H passes through 0 at rho_c and expansion follows.
        drho = -3 * H * (rho + P)
        rho_new = rho + drho * dt
        if modified and rho_new >= RHO_C:
            rho_new = RHO_C * (1 - 1e-12)
            sign = +1.0  # bounce: contraction -> expansion
        a *= np.exp(H * dt)
        rho = rho_new
        ts.append(t); as_.append(a); rhos.append(rho)
        t += dt
    return np.array(ts), np.array(as_), np.array(rhos)

t1, a1, r1 = evolve(modified=False)
t2, a2, r2 = evolve(modified=True)

fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].plot(t1, a1, "r-", label="Classical GR + jamming pressure")
ax[0].plot(t2, a2, "b-", label="Modified Friedmann (LQC-type)")
ax[0].set_xlabel("t  [Planck times]"); ax[0].set_ylabel("scale factor a(t)")
ax[0].set_yscale("log"); ax[0].legend(fontsize=8); ax[0].set_title("Collapse vs. bounce")
ax[1].plot(t1, r1 / RHO_C, "r-"); ax[1].plot(t2, r2 / RHO_C, "b-")
ax[1].axhline(1.0, color="k", ls=":", lw=0.8, label=r"jamming density $\rho_c$")
ax[1].set_xlabel("t  [Planck times]"); ax[1].set_ylabel(r"$\rho/\rho_c$")
ax[1].set_yscale("log"); ax[1].legend(fontsize=8); ax[1].set_title("Density history")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "bounce.png"), dpi=150)

print("classical: final a = %.3e, max rho/rho_c = %.3e  (runaway collapse)" % (a1[-1], r1.max() / RHO_C))
print("modified : min  a = %.3e, max rho/rho_c = %.6f  (bounce at rho_c)" % (a2.min(), r2.max() / RHO_C))
