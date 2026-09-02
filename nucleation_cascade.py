"""Cascade nucleation: seed bubble at center; new bubbles are TRIGGERED near
existing walls (catalyzed), chain-reacting outward until a jamming stop
condition (converted fraction >= eta_stop of the ball). Kicks as before.
Compare sigma_rel = dispersion / <v_r> with Poisson result (~5-8) and the
DM-fraction requirement (~1.0).
Usage: python3 nucleation_cascade.py <d_trig> <lag>
"""
import sys
import numpy as np

def run(d_trig, lag, seed=0, NP=4000, eta_stop=0.64):
    rng = np.random.default_rng(seed)
    xs = [np.zeros(3)]; ts = [0.0]           # seed bubble at center
    t = 0.0
    for step in range(3000):
        # pick a random existing bubble weighted by current wall area (R^2)
        R = np.maximum(t - np.array(ts), 0)
        w = R ** 2 + 1e-12
        j = rng.choice(len(xs), p=w / w.sum())
        # trigger point: distance R_j + d_trig from bubble j, random direction
        d = rng.normal(size=3); d /= np.linalg.norm(d)
        x_new = xs[j] + (R[j] + d_trig) * d
        if np.linalg.norm(x_new) > 1.0:      # outside the patch
            t += lag; continue
        # skip if already converted
        if any(np.linalg.norm(x_new - x) < (t - tt) for x, tt in zip(xs, ts)):
            t += lag; continue
        t_new = t + lag
        xs.append(x_new); ts.append(t_new); t = t_new
        # stop condition: converted volume fraction (MC estimate every 25)
        if step % 25 == 0 and len(xs) > 10:
            pts = rng.normal(size=(400, 3)); pts /= np.linalg.norm(pts, axis=1)[:, None]
            pts *= rng.uniform(0, 1, 400)[:, None] ** (1/3)
            X = np.array(xs); T = np.array(ts)
            conv = ((np.linalg.norm(pts[:, None] - X[None], axis=2) < (t - T)[None]).any(1)).mean()
            if conv >= eta_stop: break
    X = np.array(xs); T = np.array(ts)
    t_end = t
    # kicks at sample points (same recipe as Poisson run)
    p = rng.normal(size=(NP, 3)); p /= np.linalg.norm(p, axis=1)[:, None]
    p *= rng.uniform(0, 1, NP)[:, None] ** (1/3)
    D = np.linalg.norm(p[:, None] - X[None], axis=2)
    arr = T[None, :] + D
    t_first = np.minimum(arr.min(1), t_end)
    dt_w = max(np.median(np.diff(np.sort(T))), 1e-3) if len(T) > 3 else 0.05
    v = np.zeros((NP, 3))
    for j in range(len(X)):
        hit = arr[:, j] <= t_first + dt_w
        Rj = arr[hit, j] - T[j]
        rhat = (p[hit] - X[j]) / np.maximum(D[hit, j], 1e-9)[:, None]
        v[hit] += Rj[:, None] * rhat
    r = np.linalg.norm(p, axis=1); rhat_c = p / np.maximum(r, 1e-9)[:, None]
    vr = (v * rhat_c).sum(1)
    vres = v - vr[:, None] * rhat_c
    disp = np.sqrt(np.var(vr) + (vres ** 2).sum(1).mean())
    sp = np.linalg.norm(v, axis=1)
    return len(X), vr.mean(), disp, sp.std() / sp.mean()

if __name__ == "__main__":
    d_trig, lag = float(sys.argv[1]), float(sys.argv[2])
    print(f"d_trig={d_trig} lag={lag}")
    rr = []
    for seed in range(5):
        nb, mr, dd, cv = run(d_trig, lag, seed)
        rr.append(dd / mr if mr > 0 else np.nan)
        print(f"  seed {seed}: bubbles={nb:4d}  <v_r>={mr:.3f}  disp={dd:.3f}  sigma_rel={dd/mr if mr>0 else float('nan'):.2f}  cv={cv:.2f}")
    rr = np.array(rr)
    print(f"  => sigma_rel = {np.nanmean(rr):.2f} +- {np.nanstd(rr):.2f}  (need ~1.0; Poisson gave 5-8)")
