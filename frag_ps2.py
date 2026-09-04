"""Rank-1/2 combined: horizon-crossing fragmentation of the re-jammed kernel with three spectra and the bookkeeping
that decides it: f_PBH = (permanent tail)/(unfragmented survivors), since the survivors are the dark matter (Sec 5),
and the energy injected after BBN, in the mu and y eras, and after recombination, relative to the radiation already
made by earlier explosions (PBH-dominated era: radiation at t is the mass exploded before t).
Spectra: 'kolmogorov' sigma = A (M/M_L)^(2/9); 'flat' sigma = A (M/M_L)^0.006; 'shot' sigma^2 = A^2 (M/M_L)^0.012 + m_P/M
(granular shot noise of discrete relics plus a flat floor). Planck-star clock t = (M/m_P)^2 t_P.
Usage: python3 frag_ps2.py"""
import numpy as np
from scipy.special import erfc
mP = 2.176e-8; tP = 5.39e-44; ML = 2.2e52 / mP
lnM = np.linspace(0, np.log(ML), 6000); M = np.exp(lnM); dl = lnM[1] - lnM[0]; kg = M * mP
t = (M) ** 2 * tP
eras = [('BBN 1-1e4 s (hadronic)', (t >= 1) & (t < 1e4), 1e-6), ('late BBN 1e4-1e6 s', (t >= 1e4) & (t < 1e6), 1e-4),
        ('mu 1e6-1e9 s', (t >= 1e6) & (t < 1e9), 1.2e-4), ('y 1e9-1e13 s', (t >= 1e9) & (t < 1e13), 6e-5),
        ('post-recomb', t >= 1e13, 1e-7)]
def sigma(spec, A):
    if spec == 'kolmogorov': return A * (M / ML) ** (2 / 9)
    if spec == 'flat': return A * (M / ML) ** 0.006
    return np.sqrt(A ** 2 * (M / ML) ** 0.012 + 1.0 / M)
print("spec        A     dc    f_surv    tail>6e22  f_PBH=tail/f  IMBH/f    " + " | ".join(e[0] for e in eras) + "   (limit)")
for spec in ('kolmogorov', 'flat', 'shot'):
    for A in ([0.1, 0.3, 1.0] if spec == 'kolmogorov' else [0.0, 0.1, 0.3, 0.5, 0.7, 1.0]):
        if spec != 'shot' and A == 0.0: continue
        for dc in (0.3, 0.45, 0.7, 1.0):
            beta = erfc(dc / (np.sqrt(2) * sigma(spec, A)))
            U = np.exp(-np.cumsum(beta) * dl); dfrac = U * beta * dl; f = U[-1]
            tail = dfrac[kg >= 6e22].sum(); imbh = dfrac[(kg >= 2e32) & (kg < 2e35)].sum()
            inj = []
            for name, m, lim in eras:
                before = dfrac[t < (t[m].min() if m.any() else np.inf)].sum()
                inj.append(dfrac[m].sum() / max(before, 1e-300))
            flag = ' KILL' if (f > 1e-13 or tail / max(f, 1e-300) > 1e-3 or any(i > l for i, (_, _, l) in zip(inj, eras))) else ' <-- alive'
            print(f"{spec:<10s} {A:<5.2f} {dc:<5.2f} {f:9.2e} {tail:9.2e} {tail/max(f,1e-300):11.2e} {imbh/max(f,1e-300):9.2e}   "
                  + " | ".join(f"{i:9.2e}" for i in inj) + flag)
print("limits per era (delta rho/rho):", [(e[0], e[2]) for e in eras], " rough; BBN hadronic from Kawasaki et al. 2018, mu/y from FIRAS")
