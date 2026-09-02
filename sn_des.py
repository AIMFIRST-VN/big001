"""DES-SN5YR (Dovekie recalibration) test of the coasting model, full covariance.

chi2 with offset marginalization: chi2 = r^T C^-1 r - (sum C^-1 r)^2 / (sum C^-1).
Cov = stat (diag MUERR^2) + SYS (packed lower/upper triangular from npz).
"""
import numpy as np

rows = [l.split()[1:] for l in open("DES-Dovekie_HD.csv") if l.startswith("SN:")]
arr = np.array([[r[2], r[4], r[5]] for r in rows], dtype=float)  # zHD MU MUERR
z, mu, muerr = arr[:, 0], arr[:, 1], arr[:, 2]
n = len(z)
p = np.load("STAT+SYS.npz")
nsn = int(p["nsn"][0]); flat = p["cov"]
print(f"{n} SNe, cov nsn={nsn}, packed len={len(flat)}")
assert nsn == n and len(flat) == n * (n + 1) // 2 or len(flat) == n * n

# npz stores the INVERSE covariance (stat+sys), upper triangular (per official likelihood)
Ci = np.zeros((n, n))
Ci[np.triu_indices(n)] = flat
il = np.tril_indices(n, -1)
Ci[il] = Ci.T[il]

CH, H0 = 299792.458, 70.0
def mu_model(model):
    out = np.empty(n)
    for i, zv in enumerate(z):
        zs = np.linspace(0, zv, 120)
        if model == "lcdm":  E = np.sqrt(0.334 * (1 + zs) ** 3 + 0.666)
        elif model == "coast": E = 1 + zs
        elif model == "eds":   E = (1 + zs) ** 1.5
        dl = (1 + zv) * CH / H0 * np.trapezoid(1 / E, zs)
        out[i] = 5 * np.log10(dl) + 25
    return out

one = np.ones(n)
for model in ("lcdm", "coast", "eds"):
    r = mu - mu_model(model)
    A = r @ Ci @ r; B = one @ Ci @ r; Cn = one @ Ci @ one
    chi2 = A - B * B / Cn
    print(f"{model:6s}  chi2 = {chi2:8.1f} / {n} SNe   chi2/dof = {chi2/n:.3f}")
