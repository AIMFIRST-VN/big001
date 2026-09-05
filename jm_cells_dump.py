"""COPY of jm_cells.py that also saves per-realization arrays (wall, offs, nnb, dts) to symbolic_runs/jm_cells_samples_<seed>.npz."""
"""Johnson-Mehl (nucleation-and-growth) Monte Carlo for the pool of bangs.
Seeds appear as a Poisson process in space-time with rate I per unit volume per unit time; each grows at c = 1.
A point p belongs to the seed that reaches it first: arrival a_i(p) = t_i + |p - x_i|.
We put a random observer at the origin, find its cell (winner w), then along many directions find the first
distance s at which another seed beats the winner (the wall). Outputs, in units of the mean cell size
L = I^{-1/4} (space-time scaling): nearest-wall distance, number of distinct neighbours touched by rays,
age differences t_i - t_w across walls (negative = neighbour is OLDER), and the observer's offset from
its own seed (the 'Where is Nemo' distribution)."""
import numpy as np, sys, os
HERE = os.path.dirname(os.path.abspath(__file__))
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
nreal = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
ndir = 200
I = 1.0                       # rate -> L = 1
box = 3.0; tmax = 5.0         # seeds in [-box,box]^3 x [0,tmax]; expected count = I * (2box)^3 * tmax
dirs = rng.normal(size=(ndir, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
s_grid = np.linspace(0.005, 3.0, 400)
wall, nnb, older_frac, dts, offs, seedage = [], [], [], [], [], []
for _ in range(nreal):
    n = rng.poisson(I * (2 * box)**3 * tmax)
    x = rng.uniform(-box, box, size=(n, 3)); t = rng.uniform(0, tmax, size=n)
    arr0 = t + np.linalg.norm(x, axis=1)
    w = np.argmin(arr0)
    if arr0[w] > tmax - 3.0: continue      # avoid edge effects in time
    # points along rays
    P = s_grid[:, None, None] * dirs[None, :, :]                     # (S, D, 3)
    aw = t[w] + np.linalg.norm(P - x[w], axis=2)                     # (S, D)
    # arrival of all other seeds: (S, D, n) too big? n ~ 14000 -> S*D*n = 800*300*14000 = 3e9: too big.
    # restrict to seeds that could matter: those with t_i + |x_i| - |x_w| - t_w < 2*s_max
    cand = np.where((t + np.linalg.norm(x, axis=1)) < arr0[w] + 2 * s_grid[-1] + 0.5)[0]
    cand = cand[cand != w]
    best_s = np.full(ndir, np.inf); best_i = np.full(ndir, -1)
    for i in cand:
        ai = t[i] + np.linalg.norm(P - x[i], axis=2)                 # (S, D)
        beat = ai <= aw
        hit = beat.any(axis=0)
        if hit.any():
            first = np.argmax(beat, axis=0)
            s_hit = np.where(hit, s_grid[first], np.inf)
            better = s_hit < best_s
            best_s[better] = s_hit[better]; best_i[better] = i
    ok = np.isfinite(best_s)
    if not ok.all(): continue
    wall.append(best_s.min())
    nb = np.unique(best_i)
    nnb.append(len(nb))
    dt = t[nb] - t[w]
    dts.extend(dt.tolist()); older_frac.append(np.mean(dt < 0))
    offs.append(np.linalg.norm(x[w])); seedage.append(t[w])
wall, nnb, older_frac, dts, offs = map(np.array, (wall, nnb, older_frac, dts, offs))
seed_id = int(sys.argv[1]) if len(sys.argv) > 1 else 0
np.savez(os.path.join(HERE, "symbolic_runs", f"jm_cells_samples_{seed_id}.npz"), wall=wall, offs=offs, nnb=nnb, dts=dts, older_frac=older_frac, seedage=np.array(seedage))
q = lambda a: np.percentile(a, [10, 50, 90])
print(f"realizations used: {len(wall)}  (mean cell size L = 1, c = 1)")
print(f"nearest wall distance / L: 10/50/90% = {q(wall).round(3)}")
print(f"observer offset from own seed / L: 10/50/90% = {q(offs).round(3)}")
print(f"number of neighbours touched (distinct cells along {ndir} rays): 10/50/90% = {q(nnb).round(1)}")
print(f"fraction of a cell's neighbours that are OLDER: mean = {older_frac.mean():.2f}")
print(f"age difference across walls (t_neighbour - t_ours)/ (L/c): 10/50/90% = {q(dts).round(3)}; |dt| median = {np.median(np.abs(dts)):.3f}")
print(f"P(nearest wall < 0.25 L) = {np.mean(wall < 0.25):.2f}, < 0.5 L = {np.mean(wall < 0.5):.2f}, > 1 L = {np.mean(wall > 1.0):.2f}")
