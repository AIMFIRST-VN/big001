"""Toy effective w(z) from the measured ejecta kinematics.

Construction (explicitly a toy; assumptions listed in the paper):
- After ejection, each particle free-streams: r(t) = v * (t - t_ej) ~ v*t.
  Velocity sorting => at cosmic time t the flow field is v(r) = r/t (exact
  for cold streams at late times), i.e. locally Hubble-like with H = 1/t.
- BUT the *density* profile is not uniform: n(r) at fixed t maps from the
  measured speed distribution f(v):  n(r) dr ~ f(v=r/t) dv.
- An interior observer fits FLRW to this. We compute the effective
  deceleration from the evolution of the flow they infer:
  For pure free streaming (no gravity), a(t) ~ t (coasting, q=0, w_eff=-1/3).
  Gravity of the interior mass M(<r) decelerates each shell:
  dv/dt = -G M(<r) / r^2. The evolving deceleration of the flow maps to
  q(z) and w_eff(z) = (2q(z) - 1)/3.
We integrate shells forward with self-gravity from the measured speed
distribution (spherical, mass-weighted), then compute q(t) of the flow at
the observer-weighted radius and convert to w_eff(z) with 1+z ~ a0/a.
Units: G = M_ej = 1, initial radius spread from t_ej spread (small) ignored.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# measured ejecta speeds (all 15 seeds, unbound, in units where v_esc=sqrt2)
spd = []
for s in range(1, 16):
    d = np.load(f"ejecta_full_{s}.npz")
    u = d["unbound"]
    spd.append(np.linalg.norm(d["v"][u], axis=1))
v0 = np.concatenate(spd)
n = len(v0)
m = 1.0 / n                      # equal masses, total ejecta mass = 1

# spherical shell integration: sort by speed => shells never cross (cold radial flow)
order = np.argsort(v0)
v = v0[order].copy()
Mcore = 15.0 / 85.0              # trapped furnace mass in ejecta-mass units
Menc = Mcore + m * np.arange(1, n + 1)  # mass enclosed per shell (sorted => fixed)
# shells bound in the spherical remap (E = v^2/2 - Menc/r0 < 0) rejoin the furnace
E = 0.5 * v ** 2 - Menc / 1.0
keep = E > 0
print(f"fallback (re-bound in spherical remap): {(~keep).mean()*100:.1f}% of ejecta")
v = v[keep]; Menc = Menc[keep]; n = keep.sum()
r = np.full(n, 1.0) + 1e-6 * np.arange(n)  # shells start at system radius

dt0, tend = 1e-3, 200.0
ts, qs, aa = [], [], []
t = 0.0
a_prev = None
while t < tend:
    dt = max(dt0, 0.01 * t)
    acc = -Menc / r ** 2
    v += acc * dt
    r += v * dt
    t += dt
    # observer-weighted flow: mass-weighted mean of (r,v); effective scale factor
    a_eff = np.average(r)                     # mean shell radius ~ a(t)
    H_eff = np.average(v / np.maximum(r, 1e-12))
    qdotterm = np.average(acc / np.maximum(r, 1e-12))
    q_eff = -qdotterm / H_eff ** 2 if H_eff > 0 else np.nan
    ts.append(t); qs.append(q_eff); aa.append(a_eff)
ts, qs, aa = map(np.array, (ts, qs, aa))
z = aa[-1] / aa - 1.0
w = (2 * qs - 1) / 3.0

sel = (z > 0.01) & (z < 3) & np.isfinite(w)
print("z     q_eff   w_eff")
for zz in (2.0, 1.0, 0.5, 0.2, 0.05):
    i = np.argmin(np.abs(z - zz))
    print(f"{z[i]:5.2f}  {qs[i]:+.3f}  {w[i]:+.3f}")
i0, i1 = np.argmin(np.abs(z - 1.0)), np.argmin(np.abs(z - 0.05))
print(f"\ndrift: w_eff(z=1)={w[i0]:+.3f} -> w_eff(z~0)={w[i1]:+.3f}",
      "(more negative today = 'strengthening dark energy', matches DES/DESI direction)" if w[i1] < w[i0] else "(WRONG direction vs DES/DESI)")

plt.figure(figsize=(7, 4.5))
plt.plot(z[sel], w[sel], "b-", label=r"toy $w_{\rm eff}(z)$ from measured ejecta shells")
plt.axhline(-1, color="k", ls=":", lw=0.8, label=r"$\Lambda$ ($w=-1$)")
plt.axhline(-1/3, color="0.5", ls="--", lw=0.8, label=r"coasting ($w=-1/3$)")
plt.errorbar([0.3], [-0.803], yerr=[0.054], fmt="rs", ms=5, label=r"DES+DESI $w_0$ (pivot, approx)")
plt.gca().invert_xaxis()
plt.xlabel("redshift z"); plt.ylabel(r"$w_{\rm eff}$")
plt.title("Effective equation of state from velocity-sorted ejecta shells")
plt.legend(fontsize=8); plt.tight_layout()
plt.savefig("wz_effective.png", dpi=150)
