"""Item 3: model distance-redshift relation vs Pantheon+ SNe Ia.

Compares the Big Pool Bangs kinetic/coasting expansion (a ~ t, q=0, Milne)
against flat LambdaCDM (Om=0.334, Pantheon+ best fit) and matter-only EdS.
Absolute magnitude M marginalized analytically per model (chi2 with diagonal
errors — adequate for a first-pass discriminant).
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from numpy import trapz

d = np.genfromtxt("pantheon.dat", names=True, dtype=None, encoding=None)
mask = (d["zHD"] > 0.023) & (d["IS_CALIBRATOR"] == 0)   # Hubble-flow SNe only
z = d["zHD"][mask]; mu = d["MU_SH0ES"][mask]; err = d["MU_SH0ES_ERR_DIAG"][mask]
print(f"{mask.sum()} Hubble-flow SNe, z = {z.min():.3f} - {z.max():.3f}")

C, H0 = 299792.458, 70.0  # km/s, km/s/Mpc (absorbed into offset anyway)

def dl(zv, model):
    zs = np.linspace(0, zv, 200)
    if model == "lcdm":
        E = np.sqrt(0.334 * (1 + zs) ** 3 + 0.666)
    elif model == "eds":
        E = (1 + zs) ** 1.5
    elif model == "coast":
        E = 1 + zs                      # a ~ t  (Milne / kinetic)
    Dc = C / H0 * trapz(1 / E, zs)
    return (1 + zv) * Dc

def mu_model(model):
    return np.array([5 * np.log10(dl(zv, model)) + 25 for zv in z])

results = {}
for model in ("lcdm", "coast", "eds"):
    mm = mu_model(model)
    off = np.sum((mu - mm) / err**2) / np.sum(1 / err**2)  # marginalize offset
    chi2 = np.sum(((mu - mm - off) / err) ** 2)
    results[model] = (mm + off, chi2)
    print(f"{model:6s}  chi2 = {chi2:8.1f}  / {len(z)} SNe   chi2/dof = {chi2/len(z):.3f}")

# binned residual plot relative to coasting model
zb = np.logspace(np.log10(0.023), np.log10(2.3), 25)
ib = np.digitize(z, zb)
fig, ax = plt.subplots(figsize=(8, 5))
base = results["coast"][0]
for model, col, lab in (("lcdm", "tab:blue", r"$\Lambda$CDM ($\Omega_m=0.334$)"),
                        ("coast", "tab:red", "kinetic/coasting ($a \\propto t$, this work)"),
                        ("eds", "tab:green", "EdS (matter only)")):
    mm = results[model][0]
    zg = np.logspace(np.log10(0.023), np.log10(2.3), 100)
    mg = np.array([5 * np.log10(dl(zv, model)) + 25 for zv in zg])
    mgc = np.array([5 * np.log10(dl(zv, "coast")) + 25 for zv in zg])
    off = (results[model][0] - np.array([5*np.log10(dl(zv,model))+25 for zv in z])).mean()
    offc = (results["coast"][0] - np.array([5*np.log10(dl(zv,"coast"))+25 for zv in z])).mean()
    ax.plot(zg, (mg + off) - (mgc + offc), color=col, label=lab + f"  ($\\chi^2/dof={results[model][1]/len(z):.2f}$)")
res = mu - base
bz, bm, be = [], [], []
for i in range(1, len(zb)):
    s = ib == i
    if s.sum() > 3:
        bz.append(z[s].mean()); bm.append(np.average(res[s], weights=1/err[s]**2))
        be.append(1 / np.sqrt(np.sum(1/err[s]**2)))
ax.errorbar(bz, bm, yerr=be, fmt="ko", ms=4, capsize=2, label="Pantheon+ (binned)")
ax.set_xscale("log"); ax.axhline(0, color="0.6", lw=0.5)
ax.set_xlabel("redshift z"); ax.set_ylabel(r"$\Delta\mu$ relative to coasting model [mag]")
ax.set_title("Distance modulus residuals vs. kinetic (coasting) expansion")
ax.legend(fontsize=8); fig.tight_layout()
fig.savefig("sn_hubble.png", dpi=150)
