"""Gemini idea 1 check: gravitational focusing (Sommerfeld-like) on Planck-relic mergers after the bounce.
sigma_eff = pi l_P^2 (1 + v_esc^2/v^2), v_esc^2 = 2 G m_P / l_P = 2 c^2 -> sigma v ~ 2 pi l_P^2 c^2 / v for v << c.
Kinetically decoupled relics: v = c (a_P/a). Radiation era: H = H_P (a_P/a)^2, n = f n_P (a_P/a)^3, n_P = 1/l_P^3 (jammed).
Then Gamma/H = n sigma v / H = 2 pi f (c/(l_P H_P)) = 2 pi f * k, k = O(1): constant in a except through f.
Depletion: df/dln a = -(Gamma/H) f = -2 pi k f^2  ->  f = 1 / (1/f0 + 2 pi k ln a).  Logarithmic, not exponential.
Also: force-chain collapse into N m_P black holes -> evaporation time t = N^3 t_P (order of magnitude)."""
import numpy as np
k=1.0; f0=1.0
for lna in (10, 30, 64, 100, 130):
    f=1/(1/f0+2*np.pi*k*lna); print(f"ln a = {lna:4d} (T ~ {1.2e28*np.exp(-lna):8.1e} eV):  f = {f:.2e}")
print("required f = T_eq/T_conv = 0.8 eV / 1.2e28 eV = 6.7e-29  -> gap remaining ~ 25 orders")
print("\nforce chains: N relics -> M = N m_P, evaporation time N^3 t_P (t_P = 5.4e-44 s)")
for N in (1e2,1e5,1e10,1e14,1e20):
    t=N**3*5.4e-44; print(f"  N={N:8.0e}  M={N*2.2e-8:8.1e} kg  t_evap={t:8.1e} s  ({'before BBN' if t<1 else 'after BBN' if t<3.8e5*3.15e7 else 'after recombination' if t<4.3e17 else 'not yet'})")
