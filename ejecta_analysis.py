"""Radial + tangential + coupling analysis of ejecta kinematics (15 seeds)."""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

runs = [np.load(f"ejecta_full_{s}.npz") for s in range(1, 16)]

# ---------- radial: ejection time & speed sorting ----------
tej, spd = [], []
for d in runs:
    u = d["unbound"] & ~np.isnan(d["t_ej"])
    tej.append(d["t_ej"][u]); spd.append(np.linalg.norm(d["v"][u], axis=1))
tej = np.concatenate(tej); spd = np.concatenate(spd) / np.sqrt(2)  # units of v_esc
print(f"n={len(tej)} ejecta with times; t_ej median={np.median(tej):.2f}, 10-90%={np.percentile(tej,10):.2f}-{np.percentile(tej,90):.2f}")
for lo, hi in [(0, 2), (2, 5), (5, 10), (10, 25)]:
    m = (tej >= lo) & (tej < hi)
    if m.sum() > 10:
        print(f"  t_ej {lo:>2}-{hi:<2}: {m.sum():5d} ejecta, v_med={np.median(spd[m]):.2f} v_esc")

# ---------- tangential: dipole + low-l power vs shot noise ----------
def dipole(dirs):
    return np.linalg.norm(dirs.mean(0))
dips, ns = [], []
for d in runs:
    u = d["unbound"]
    dirs = d["v"][u] / np.linalg.norm(d["v"][u], axis=1)[:, None]
    dips.append(dipole(dirs)); ns.append(u.sum())
dips, ns = np.array(dips), np.array(ns)
shot = np.sqrt(3 / ns)  # expected |mean| of n random unit vectors ~ sqrt(3/n)... use MC
rng = np.random.default_rng(0)
mc = np.array([[dipole(rng.normal(size=(n, 3)) / np.linalg.norm(rng.normal(size=(n, 3)), axis=1)[:, None]) for _ in range(200)] for n in ns[:1]])
print(f"\ndipole |<n>|: measured mean={dips.mean():.4f} +- {dips.std():.4f}  (15 skies)")
print(f"shot-noise MC (n={ns[0]}): {mc.mean():.4f} +- {mc.std():.4f}")
print(f"ratio measured/shot = {dips.mean()/mc.mean():.2f}")

# ---------- coupling: early vs late ejecta dipole directions ----------
cosangs = []
for d in runs:
    u = d["unbound"] & ~np.isnan(d["t_ej"])
    dirs = d["v"][u] / np.linalg.norm(d["v"][u], axis=1)[:, None]
    t = d["t_ej"][u]
    med = np.median(t)
    d1 = dirs[t <= med].mean(0); d2 = dirs[t > med].mean(0)
    cosangs.append(d1 @ d2 / np.linalg.norm(d1) / np.linalg.norm(d2))
cosangs = np.array(cosangs)
print(f"\nearly-vs-late dipole alignment: <cos angle> = {cosangs.mean():.3f} +- {cosangs.std():.3f}")
print("(+1 = same axis persists; 0 = uncorrelated shells)")

# ---------- plot ----------
fig, ax = plt.subplots(1, 2, figsize=(10, 4))
ax[0].hist(tej, bins=40, color="tab:blue", alpha=0.8)
ax[0].set_xlabel("ejection time [dyn. units]"); ax[0].set_ylabel("N ejecta"); ax[0].set_title("Ejection-time distribution")
h = ax[1].hist2d(tej, spd, bins=[30, 30], cmap="viridis")
ax[1].set_xlabel("ejection time"); ax[1].set_ylabel(r"final speed [$v_{esc}$]"); ax[1].set_title("Speed vs ejection time")
fig.colorbar(h[3], ax=ax[1])
fig.tight_layout(); fig.savefig("ejecta_kinematics.png", dpi=150)
