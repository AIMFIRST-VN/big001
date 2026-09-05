"""Mobile relic/anti-relic annihilation on a jammed packing, many reshuffles, vectorized (random priority matching
on a simple cubic lattice, z=6). Each round: random reshuffle of survivors' positions (full mixing between
recompressions), then annihilation of opposite-sign contact pairs. Reports the survivor fraction versus rounds and
fits the late-time decay law, to see how many reshuffles the required survivor fraction (1e-28..1e-13) needs.
Equal numbers of relics and anti-relics (exact CP), and a run with a 1e-3 excess (survivor floor = the excess)."""
import numpy as np, sys
rng = np.random.default_rng(int(sys.argv[1]) if len(sys.argv) > 1 else 0)
n = 64; N = n**3
def run(excess, rounds):
    signs = rng.choice([-1, 1], size=N); k = int(excess * N); signs[:k] = 1
    alive = np.ones(N, bool); hist = []
    for r in range(rounds):
        idx = np.where(alive)[0]
        if len(idx) < 2: hist.append(0.0); break
        # reshuffle: survivors randomly placed on the lattice
        pos = rng.choice(N, size=len(idx), replace=False)
        grid = np.zeros(N, np.int8); grid[pos] = signs[idx]
        g = grid.reshape(n, n, n); prio = rng.random(N).reshape(n, n, n)
        matched = np.zeros((n, n, n), bool)
        for axis in range(3):
            for sh in (1, -1):
                nb = np.roll(g, sh, axis); nbm = np.roll(matched, sh, axis); nbp = np.roll(prio, sh, axis)
                pair = (g != 0) & (nb == -g) & (~matched) & (~nbm) & (prio > nbp)
                # a site pairs with at most one neighbour: resolve by first-come in this loop order
                matched |= pair; matched |= np.roll(pair, -sh, axis)
        surv = (g != 0) & (~matched)
        keep = surv.reshape(-1)[pos]
        alive[idx[~keep]] = False
        hist.append(alive.sum() / N)
    return np.array(hist)
for excess in (0.0, 1e-3):
    h = run(excess, 400)
    print(f"excess={excess}: survivors after 1,3,10,30,100,200,400 rounds: " + " ".join(f"{h[i-1]:.2e}" for i in (1,3,10,30,100,200,400) if i <= len(h)))
    if excess == 0.0:
        rr = np.arange(1, len(h)+1); m = (rr > 50) & (h > 20 / N)
        if m.sum() > 5:
            p = np.polyfit(np.log(rr[m]), np.log(h[m]), 1); print(f"  late-time power law: survivors ~ rounds^({p[0]:.2f})")
            print(f"  rounds needed for 1e-13: {np.exp((np.log(1e-13)-p[1])/p[0]):.1e}; for 1e-28: {np.exp((np.log(1e-28)-p[1])/p[0]):.1e}")
