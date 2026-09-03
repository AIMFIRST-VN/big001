"""Contact annihilation of relic/anti-relic on a jammed packing, approximated by a simple cubic lattice (z=6)
with random signs. Immobile relics: only interface pairs annihilate; same-sign clusters survive (17.5% after one
round, unchanged thereafter). Mobile relics (random reshuffle between rounds): survivors fall slowly,
~1% after 30 rounds, versus the required 1e-28 (Planck-scale conversion). Result used to assess the
'annihilation on contact' channel: it does not remove the reshuffle/cross-section bottleneck of open problem 4(a)."""
import numpy as np
rng = np.random.default_rng(0); n = 40
s = rng.choice([-1, 1], size=(n, n, n)); alive = np.ones_like(s, bool)
def round_annihilate(s, alive):
    idx = np.argwhere(alive); rng.shuffle(idx)
    for (i, j, k) in idx:
        if not alive[i, j, k]: continue
        for d in ((1,0,0),(-1,0,0),(0,1,0),(0,-1,0),(0,0,1),(0,0,-1)):
            a, b, c = (i+d[0]) % n, (j+d[1]) % n, (k+d[2]) % n
            if alive[a, b, c] and s[a, b, c] == -s[i, j, k]:
                alive[i, j, k] = False; alive[a, b, c] = False; break
    return alive
N0 = alive.sum()
for r in range(1, 4):
    alive = round_annihilate(s, alive); print(f"immobile, round {r}: surviving fraction = {alive.sum()/N0:.3f}")
for r in range(1, 31):
    vals = s[alive]; flat = np.zeros(n**3, int); flat[:len(vals)] = vals; rng.shuffle(flat)
    s = flat.reshape(n, n, n); alive = s != 0
    alive = round_annihilate(s, alive)
    if r in (1, 3, 10, 30): print(f"mobile (reshuffle), round {r}: surviving fraction = {alive.sum()/N0:.4f}")
