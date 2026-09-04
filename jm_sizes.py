"""Cell-size distribution of the Johnson-Mehl pool: seeds Poisson in space-time (rate I=1 -> L=1), growth at c=1.
Each cell = points reached first by that seed. Estimate every seed's cell volume with random sample points.
Reports the per-cell volume distribution (small cells = seeds born late in the gaps) and the volume-weighted one
(what a random observer sees). Also the age at birth vs final size."""
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
box, tmax = 4.0, 6.0        # transition complete long before tmax: untransformed fraction exp(-pi t^4/3)
n = rng.poisson((2 * box)**3 * tmax); x = rng.uniform(-box, box, (n, 3)); t = rng.uniform(0, tmax, n)
# keep only seeds that nucleate in still-untransformed space: seed i is real if no earlier seed j has t_j+|x_i-x_j| < t_i
arr = t[:, None] + np.linalg.norm(x[:, None, :] - x[None, :, :], axis=2)
cond = arr > t[None, :]; np.fill_diagonal(cond, True); real = np.all(cond, axis=0)
xr, tr = x[real], t[real]; m = len(xr)
P = rng.uniform(-box + 1.0, box - 1.0, (400000, 3))   # interior sample points (edge cells excluded)
win = np.concatenate([np.argmin(tr[None, :] + np.linalg.norm(Pc[:, None, :] - xr[None, :, :], axis=2), axis=1)
                      for Pc in np.array_split(P, 40)])   # chunked: ~120 MB peak instead of ~5 GB
inner = np.all(np.abs(xr) < box - 1.0, axis=1)
vol = np.bincount(win, minlength=m) * ((2 * box - 2.0)**3 / len(P))
v = vol[inner & (vol > 0)]; age = tr[inner & (vol > 0)]
q = lambda a: np.percentile(a, [5, 25, 50, 75, 95]).round(3)
print(f"real seeds {m} (of {n} attempted; the rest fell in already-transformed space); interior cells {len(v)}; mean cell volume {v.mean():.3f} L^3")
print(f"per-cell volume / L^3: 5/25/50/75/95% = {q(v)};  fraction of cells with V < 0.1 L^3: {np.mean(v<0.1):.2f}, < 0.01: {np.mean(v<0.01):.3f}")
w = v / v.sum(); order = np.argsort(v); cw = np.cumsum(w[order])
print(f"volume-weighted (random observer): P(own cell V < 0.1 L^3) = {cw[np.searchsorted(v[order],0.1)]:.3f}, < 0.5: {cw[np.searchsorted(v[order],0.5)]:.3f}")
for lo, hi in ((0, 0.3), (0.3, 0.6), (0.6, 1.0), (1.0, 2.0)):
    s = (age >= lo) & (age < hi); print(f"  seeds born at t in [{lo},{hi}) L/c: {s.sum():4d} cells, median volume {np.median(v[s]) if s.any() else float('nan'):.3f} L^3")
