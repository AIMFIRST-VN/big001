"""Johnson-Mehl (nucleation-and-growth) Monte Carlo for the pool of bangs.
Seeds appear as a Poisson process in space-time with rate I per unit volume per unit time; each grows at c = 1.
A point p belongs to the seed that reaches it first: arrival a_i(p) = t_i + |p - x_i|.
We put a random observer at the origin, find its cell (winner w), then along many directions find the first
distance s at which another seed beats the winner (the wall). Outputs, in units of the mean cell size
L = I^{-1/4} (space-time scaling): nearest-wall distance, number of distinct neighbours touched by rays,
age differences t_i - t_w across walls (negative = neighbour is OLDER), and the observer's offset from
its own seed (the 'Where is Nemo' distribution)."""
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
nreal = int(sys.argv[2]) if len(sys.argv) > 2 else 1000
ndir = 200
I = 1.0                       # rate -> L = 1
box = 3.0; tmax = 5.0         # seeds in [-box,box]^3 x [0,tmax]; expected count = I * (2box)^3 * tmax
dirs = rng.normal(size=(ndir, 3)); dirs /= np.linalg.norm(dirs, axis=1)[:, None]
s_grid = np.linspace(0.005, 3.0, 400)
wall, nnb, older_frac, dts, offs, seedage = [], [], [], [], [], []
D0s=[0.05,0.1,0.15,0.2,0.3]; res={d:[] for d in D0s}
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
    # per-neighbour nearest wall distance (min over rays hitting that neighbour)
    smin={}
    for i,sv in zip(best_i,best_s):
        smin[i]=min(smin.get(i,np.inf),sv)
    smins=np.array(list(smin.values()))
    for D0 in D0s:
        if smins.min() < D0: res[D0].append(None); continue        # a wall inside the horizon: excluded configuration
        n_band=int(np.sum((smins>=D0)&(smins<1.06*D0)))            # imprint band epsilon <= 0.06
        n_wide=int(np.sum((smins>=D0)&(smins<1.3*D0)))
        res[D0].append((n_band,n_wide))

print(f"realizations: {len(wall)}")
for D0 in D0s:
    r=res[D0]; ok=[x for x in r if x is not None]
    if not ok: print(f"D0/L={D0}: no configuration without a wall inside the horizon"); continue
    nb=np.array([x[0] for x in ok]); nw=np.array([x[1] for x in ok])
    print(f"D0/L={D0}: P(no wall inside horizon)={len(ok)/len(r):.2f}; neighbours with nearest wall within 6% beyond horizon: mean={nb.mean():.2f}, P(=8)={np.mean(nb==8):.3f}, P(6-10)={np.mean((nb>=6)&(nb<=10)):.3f}, P(>=1)={np.mean(nb>=1):.2f}; within 30%: mean={nw.mean():.2f}, P(=8)={np.mean(nw==8):.3f}, P(6-10)={np.mean((nw>=6)&(nw<=10)):.3f}")
