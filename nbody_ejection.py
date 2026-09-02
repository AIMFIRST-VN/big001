"""Toy N-body: ejecta fraction from a cold dense core with an outward radial kick.

Tests the paper's 85/15 ejected-halo / trapped-furnace split. Violent
relaxation is nearly scale-free, so N=2000 softened equal masses suffice
for a first estimate. Units G = M_tot = R_init = 1 (v_esc(surface) ~ sqrt(2)).

Sweep kick amplitude f = v_kick / v_esc; report bound vs unbound fraction
after the system relaxes.
"""
import numpy as np

N = 1000
EPS = 0.05          # softening
DT = 0.01
T_END = 25.0
rng = np.random.default_rng(1)

def accel(x):
    d = x[:, None, :] - x[None, :, :]
    r2 = (d ** 2).sum(-1) + EPS ** 2
    inv = r2 ** -1.5
    np.fill_diagonal(inv, 0.0)
    return -(d * inv[..., None]).sum(1) / N   # m_i = 1/N

def run(f_kick):
    # uniform sphere, radius 1
    x = rng.normal(size=(N, 3)); x /= np.linalg.norm(x, axis=1)[:, None]
    x *= rng.uniform(0, 1, N)[:, None] ** (1 / 3)
    r = np.linalg.norm(x, axis=1)
    v_esc = np.sqrt(2.0)                       # at surface, order unity inside
    v = f_kick * v_esc * x / np.maximum(r, 1e-9)[:, None]  # radial outward kick
    v += rng.normal(scale=0.02, size=(N, 3))   # tiny thermal noise
    a = accel(x)
    for _ in range(int(T_END / DT)):
        v += 0.5 * DT * a
        x += DT * v
        a = accel(x)
        v += 0.5 * DT * a
    # energies: KE + PE per particle
    d = x[:, None, :] - x[None, :, :]
    r2 = (d ** 2).sum(-1) + EPS ** 2
    inv = r2 ** -0.5
    np.fill_diagonal(inv, 0.0)
    pe = -inv.sum(1) / N
    ke = 0.5 * (v ** 2).sum(1)
    unbound = (ke + pe) > 0
    speeds = np.linalg.norm(v[unbound], axis=1)
    return unbound.mean(), (np.median(speeds) if unbound.any() else 0.0)

if __name__ == "__main__":
    print("f_kick   ejected_fraction   median_ejecta_speed/v_esc")
    for f in (0.6, 0.8, 1.0, 1.2, 1.5):
        frac, vmed = run(f)
        print(f"{f:5.2f}    {frac:6.3f}             {vmed / np.sqrt(2):5.2f}")
