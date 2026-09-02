"""Turbulent-kick N-body: does ~85% ejection arise naturally?

Kick model: radial rebound f_mean*v_esc plus isotropic turbulent component
drawn from a Maxwellian with 3D dispersion sigma_rel * f_mean * v_esc.
Sweep (f_mean, sigma_rel); report ejected fraction.
Usage: python3 nbody_turbulent.py <f_mean> <sigma_rel>
"""
import sys
import numpy as np
from nbody_ejection import N, EPS, DT, T_END, accel  # reuse integrator pieces

def run_turb(f_mean, sigma_rel, seed=1):
    rng = np.random.default_rng(seed)
    x = rng.normal(size=(N, 3)); x /= np.linalg.norm(x, axis=1)[:, None]
    x *= rng.uniform(0, 1, N)[:, None] ** (1 / 3)
    r = np.linalg.norm(x, axis=1)
    v_esc = np.sqrt(2.0)
    v = f_mean * v_esc * x / np.maximum(r, 1e-9)[:, None]        # radial rebound
    v += rng.normal(scale=sigma_rel * f_mean * v_esc / np.sqrt(3), size=(N, 3))  # turbulence
    a = accel(x)
    for _ in range(int(T_END / DT)):
        v += 0.5 * DT * a
        x += DT * v
        a = accel(x)
        v += 0.5 * DT * a
    d = x[:, None, :] - x[None, :, :]
    r2 = (d ** 2).sum(-1) + EPS ** 2
    inv = r2 ** -0.5
    np.fill_diagonal(inv, 0.0)
    pe = -inv.sum(1) / N
    ke = 0.5 * (v ** 2).sum(1)
    unbound = (ke + pe) > 0
    sp = np.linalg.norm(v[unbound], axis=1) if unbound.any() else np.array([0.0])
    return unbound.mean(), np.median(sp) / v_esc

if __name__ == "__main__":
    f, s = float(sys.argv[1]), float(sys.argv[2])
    frac, vmed = run_turb(f, s)
    print(f"{f:.2f} {s:.2f} {frac:.3f} {vmed:.2f}")

def run_turb_save(f_mean, sigma_rel, seed, outfile):
    """Same as run_turb but saves unit direction vectors of ejected particles."""
    import numpy as _np
    rng = _np.random.default_rng(seed)
    _np.random.seed(seed)
    # duplicate of run_turb with state capture
    x = rng.normal(size=(N, 3)); x /= _np.linalg.norm(x, axis=1)[:, None]
    x *= rng.uniform(0, 1, N)[:, None] ** (1 / 3)
    r = _np.linalg.norm(x, axis=1)
    v_esc = _np.sqrt(2.0)
    v = f_mean * v_esc * x / _np.maximum(r, 1e-9)[:, None]
    v += rng.normal(scale=sigma_rel * f_mean * v_esc / _np.sqrt(3), size=(N, 3))
    a = accel(x)
    for _ in range(int(T_END / DT)):
        v += 0.5 * DT * a
        x += DT * v
        a = accel(x)
        v += 0.5 * DT * a
    d = x[:, None, :] - x[None, :, :]
    r2 = (d ** 2).sum(-1) + EPS ** 2
    inv = r2 ** -0.5
    _np.fill_diagonal(inv, 0.0)
    pe = -inv.sum(1) / N
    ke = 0.5 * (v ** 2).sum(1)
    unbound = (ke + pe) > 0
    dirs = v[unbound] / _np.linalg.norm(v[unbound], axis=1)[:, None]
    _np.save(outfile, dirs)
    print(f"{seed} {unbound.mean():.3f} {unbound.sum()}")

def run_turb_full(f_mean, sigma_rel, seed, outfile):
    """Save full ejecta kinematics: velocities, positions, ejection times."""
    import numpy as _np
    rng = _np.random.default_rng(seed)
    x = rng.normal(size=(N, 3)); x /= _np.linalg.norm(x, axis=1)[:, None]
    x *= rng.uniform(0, 1, N)[:, None] ** (1 / 3)
    r = _np.linalg.norm(x, axis=1)
    v_esc = _np.sqrt(2.0)
    v = f_mean * v_esc * x / _np.maximum(r, 1e-9)[:, None]
    v += rng.normal(scale=sigma_rel * f_mean * v_esc / _np.sqrt(3), size=(N, 3))
    a = accel(x)
    t_ej = _np.full(N, _np.nan)
    for step in range(int(T_END / DT)):
        v += 0.5 * DT * a
        x += DT * v
        a = accel(x)
        v += 0.5 * DT * a
        if step % 50 == 0:
            d = x[:, None, :] - x[None, :, :]
            r2 = (d ** 2).sum(-1) + EPS ** 2
            inv = r2 ** -0.5
            _np.fill_diagonal(inv, 0.0)
            e = 0.5 * (v ** 2).sum(1) - inv.sum(1) / N
            t_ej[_np.isnan(t_ej) & (e > 0)] = step * DT
            t_ej[~_np.isnan(t_ej) & (e <= 0)] = _np.nan  # re-captured
    d = x[:, None, :] - x[None, :, :]
    r2 = (d ** 2).sum(-1) + EPS ** 2
    inv = r2 ** -0.5
    _np.fill_diagonal(inv, 0.0)
    e = 0.5 * (v ** 2).sum(1) - inv.sum(1) / N
    _np.savez(outfile, v=v, x=x, t_ej=t_ej, unbound=e > 0)
    print(f"{seed} {(e>0).mean():.3f} {int((e>0).sum())}")
