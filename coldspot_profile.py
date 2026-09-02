"""P1 test: azimuthally averaged radial temperature profile around the CMB Cold Spot.

Prediction (Big Pool Bangs, Sec. 9, P1): hot ring at 5-7 deg radius,
amplitude +20 to +35 microK, around the -70 microK core.

Uses Planck 2018 SMICA map, degraded to Nside 128. Error bars from the
same profile measured at 500 random sky positions (|b| > 30 deg, away
from the galactic plane).
"""
import numpy as np
import healpy as hp
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NSIDE = 128
CS_L, CS_B = 209.0, -57.0          # Cold Spot galactic coords
RMAX, NBIN = 15.0, 30              # profile out to 15 deg, 0.5 deg bins

print("loading map...")
m = hp.read_map("smica.fits", field=0)          # I (K_CMB)
m = hp.ud_grade(m, NSIDE) * 1e6                 # -> microK
np.save("smica_nside128.npy", m)                # keep small version, drop 2GB later

# vectorized angular distance from center for all pixels
px, py, pz = hp.pix2vec(NSIDE, np.arange(len(m)))
def radial_profile(l, b):
    vx, vy, vz = hp.ang2vec(l, b, lonlat=True)
    ang = np.degrees(np.arccos(np.clip(px*vx + py*vy + pz*vz, -1, 1)))
    bins = np.linspace(0, RMAX, NBIN + 1)
    idx = np.digitize(ang, bins) - 1
    prof = np.array([m[idx == i].mean() for i in range(NBIN)])
    return prof

r = (np.arange(NBIN) + 0.5) * RMAX / NBIN
cs = radial_profile(CS_L, CS_B)

print("random reference profiles...")
rng = np.random.default_rng(42)
ref = []
while len(ref) < 500:
    l = rng.uniform(0, 360); sb = rng.uniform(-1, 1); b = np.degrees(np.arcsin(sb))
    if abs(b) < 30: continue
    ref.append(radial_profile(l, b))
ref = np.array(ref)
sigma = ref.std(axis=0)

plt.figure(figsize=(8, 5))
plt.fill_between(r, -2*sigma, 2*sigma, color="0.85", label=r"random sky $\pm 2\sigma$")
plt.fill_between(r, -sigma, sigma, color="0.7", label=r"random sky $\pm 1\sigma$")
plt.plot(r, cs, "b-o", ms=3, label="Cold Spot profile")
plt.axvspan(5, 7, color="orange", alpha=0.2, label="predicted ring (P1)")
plt.axhspan(20, 35, color="red", alpha=0.12, label=r"predicted amplitude $+20$–$35\,\mu$K")
plt.axhline(0, color="k", lw=0.5)
plt.xlabel("radius from Cold Spot center [deg]")
plt.ylabel(r"$\langle \Delta T \rangle$  [$\mu$K]")
plt.title("Planck 2018 SMICA: radial profile around the Cold Spot")
plt.legend(fontsize=8); plt.tight_layout()
plt.savefig("coldspot_profile.png", dpi=150)

for ri, ci, si in zip(r, cs, sigma):
    flag = " <-- ring band" if 5 <= ri <= 7 else ""
    print(f"r={ri:5.2f}  dT={ci:+7.2f} uK   sigma={si:5.2f}   ({ci/si:+.1f} sigma){flag}")
