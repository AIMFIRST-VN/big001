"""Can any spectrum thread the needle? Blue power law sigma = (dc/2) (M/M_piv)^(-q): strong collapse below the pivot,
none far above. Requirements: survivors f in 1e-28..1e-13 (radiation era), permanent tail/f < 1e-3, injection limits per era.
Usage: python3 frag_blue.py"""
import numpy as np
from scipy.special import erfc
mP = 2.176e-8; tP = 5.39e-44; ML = 2.2e52 / mP
lnM = np.linspace(0, np.log(ML), 6000); M = np.exp(lnM); dl = lnM[1] - lnM[0]; kg = M * mP; t = M ** 2 * tP
eras = [('BBN', (t >= 1) & (t < 1e4), 1e-6), ('lateBBN', (t >= 1e4) & (t < 1e6), 1e-4), ('mu', (t >= 1e6) & (t < 1e9), 1.2e-4),
        ('y', (t >= 1e9) & (t < 1e13), 6e-5), ('post', t >= 1e13, 1e-7)]
print("q     M_piv[kg]  dc    f_surv    tail/f     " + " ".join(f"{e[0]:>9s}" for e in eras) + "   n=1+6q (blue)")
for q in (0.05, 0.1, 0.15, 0.2, 0.3, 0.5):
    for Mp in (1e8, 1e10, 1e12, 1e13):
        for dc in (0.45,):
            sig = (dc / 2) * (M * mP / Mp) ** (-q)
            beta = erfc(dc / (np.sqrt(2) * sig)); U = np.exp(-np.cumsum(beta) * dl); dfrac = U * beta * dl; f = U[-1]
            tail = dfrac[kg >= 6e22].sum() / f
            inj = [dfrac[m].sum() / max(dfrac[t < t[m].min()].sum(), 1e-300) for _, m, _ in eras]
            ok = (1e-28 < f < 1e-13) and tail < 1e-3 and all(i < l for i, (_, _, l) in zip(inj, eras))
            print(f"{q:<5.2f} {Mp:9.0e}  {dc:<5.2f} {f:9.2e} {tail:9.2e}  " + " ".join(f"{i:9.2e}" for i in inj) + f"   {1+6*q:5.2f}" + ("  ALIVE" if ok else ""))
