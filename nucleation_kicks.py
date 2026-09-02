"""Nucleation-merger statistics -> kick spectrum prediction.

Bubbles nucleate Poisson-in-time, uniform in a unit ball; walls grow at w=1.
A sample point x is converted at t_c = min_i (t_i + |x-x_i|).  Its kick is
the vector sum of wall impulses from all bubbles arriving within a window
dt_w after first conversion, each directed away from that bubble's center,
weighted by the bubble's radius at arrival (wall momentum ~ R):
    v(x) = kappa * sum_i R_i * rhat_i
kappa is unknown (energy calibration) but cancels in the dimensionless
prediction  sigma_rel = <3D dispersion> / <mean radial drift>,
which the dark-matter fraction requires to be ~1.0.
Usage: python3 nucleation_kicks.py <gamma> (nucleation rate per vol per time)
"""
import sys
import numpy as np

def run(gamma, seed=0, NB_max=4000, NP=4000):
    rng = np.random.default_rng(seed)
    # nucleation events until the ball is (statistically) fully converted:
    # draw candidate sites/times; keep those not already inside another bubble
    # (thinning: a site inside converted phase cannot nucleate)
    T = 3.0
    n_cand = rng.poisson(gamma * (4/3*np.pi) * T)
    n_cand = min(n_cand, NB_max)
    tt = np.sort(rng.uniform(0, T, n_cand))
    xx = rng.normal(size=(n_cand, 3))
    xx /= np.linalg.norm(xx, axis=1)[:, None]
    xx *= rng.uniform(0, 1, n_cand)[:, None] ** (1/3)
    keep = np.ones(n_cand, bool)
    for i in range(n_cand):
        if not keep[i]: continue
        d = np.linalg.norm(xx[i+1:] - xx[i], axis=1)
        inside = d < (tt[i+1:] - tt[i])          # nucleates inside bubble i
        keep[i+1:] &= ~inside
    tt, xx = tt[keep], xx[keep]
    NBub = len(tt)
    # sample points in the ball
    p = rng.normal(size=(NP, 3)); p /= np.linalg.norm(p, axis=1)[:, None]
    p *= rng.uniform(0, 1, NP)[:, None] ** (1/3)
    # arrival times of every bubble wall at every point
    D = np.linalg.norm(p[:, None, :] - xx[None, :, :], axis=2)   # NP x NB
    arr = tt[None, :] + D
    t_first = arr.min(1)
    dt_w = np.median(np.diff(np.sort(tt))) if NBub > 3 else 0.05  # merger window
    v = np.zeros((NP, 3))
    for j in range(NBub):
        hit = arr[:, j] <= t_first + dt_w
        R = arr[hit, j] - tt[j]                                   # radius at arrival
        rhat = (p[hit] - xx[j]) / np.maximum(D[hit, j], 1e-9)[:, None]
        v[hit] += R[:, None] * rhat
    # decompose: radial (w.r.t. ball center) vs dispersion
    r = np.linalg.norm(p, axis=1)
    rhat_c = p / np.maximum(r, 1e-9)[:, None]
    vr = (v * rhat_c).sum(1)
    vres = v - vr[:, None] * rhat_c
    disp3d = np.sqrt(np.var(vr) + (vres ** 2).sum(1).mean())      # total 3D dispersion
    mean_r = vr.mean()
    return NBub, mean_r, disp3d, vr, v

if __name__ == "__main__":
    gamma = float(sys.argv[1]) if len(sys.argv) > 1 else 300.0
    print(f"gamma={gamma}")
    ratios = []
    for seed in range(5):
        NB, mr, dd, vr, v = run(gamma, seed)
        ratios.append(dd / mr if mr > 0 else np.nan)
        print(f"  seed {seed}: bubbles={NB:4d}  <v_r>={mr:.3f}  disp3D={dd:.3f}  sigma_rel={dd/mr:.2f}")
    ratios = np.array(ratios)
    print(f"  => sigma_rel = {np.nanmean(ratios):.2f} +- {np.nanstd(ratios):.2f}   (DM fraction requires ~1.0)")
    # shape check: normalized speed distribution vs Maxwellian
    sp = np.linalg.norm(v, axis=1)
    cv = sp.std() / sp.mean()
    print(f"  speed-distribution coefficient of variation = {cv:.2f} (3D Maxwellian: 0.42)")
