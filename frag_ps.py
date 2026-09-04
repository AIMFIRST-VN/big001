"""Rank-1 toy: fragment (PBH) mass function of the re-jammed kernel by horizon crossing (Press-Schechter).
At time t after re-jamming the largest coherent clump is the horizon mass M_H = m_P t/t_P; a region collapses when its
smoothed density contrast exceeds delta_c. Two spectra: Kolmogorov, sigma(M) = A (M/M_L)^(2/9) (delta ~ v^2/c^2 ~ k^-2/3),
and near-scale-invariant, sigma(M) = A (M/M_L)^((1-n_s)/6) with n_s = 0.965. Sequential depletion: the still-uncollapsed
fraction U obeys dU/dlnM = -U beta(M), beta = erfc(delta_c / (sqrt2 sigma)). Planck-star clock t = (M/m_P)^2 t_P sorts the
fragments into: explode before BBN (M<1e14 kg), before recombination (<3e20), before today (<6e22), permanent PBH.
Outputs the (A, delta_c) plane: unfragmented fraction f, mass fraction per window, and the 'window' where f is in the range
the paper needs (1e-28..1e-13, Sec 5) AND the permanent tail is >1e-32 (a solar mass per Hubble volume-ish) but <1e-3.
Usage: python3 frag_ps.py [spectrum: kolmogorov|flat]"""
import numpy as np, sys
from scipy.special import erfc
spec = sys.argv[1] if len(sys.argv) > 1 else 'kolmogorov'
mP = 2.176e-8; ML = 2.2e52 / mP                      # kernel (present Hubble volume) mass in m_P
lnM = np.linspace(0, np.log(ML), 4000); M = np.exp(lnM); dl = lnM[1] - lnM[0]
p = 2 / 9 if spec == 'kolmogorov' else (1 - 0.965) / 6
kg = M * mP
win = {'<1e14 kg (pre-BBN)': kg < 1e14, '1e14-3e20 (BBN-recomb)': (kg >= 1e14) & (kg < 3e20),
       '3e20-6e22 (recomb-today)': (kg >= 3e20) & (kg < 6e22), '>6e22 (permanent PBH)': kg >= 6e22,
       '>2e30 (>1 Msun)': kg >= 2e30, '2e32-2e35 (IMBH)': (kg >= 2e32) & (kg < 2e35)}
print(f"spectrum={spec}  sigma(M)=A (M/M_L)^{p:.3f}   M_L={ML:.2e} m_P")
print("A       dc    f_unfrag   " + "  ".join(f"{k:>24s}" for k in win))
good = []
for A in [0.01, 0.03, 0.1, 0.2, 0.3, 0.5, 0.7, 1.0]:
    for dc in [0.3, 0.45, 0.7, 1.0]:
        sig = A * (M / ML) ** p
        beta = erfc(dc / (np.sqrt(2) * sig))
        U = np.exp(-np.cumsum(beta) * dl)                 # uncollapsed fraction after crossing mass M
        dfrac = U * beta * dl                              # mass fraction collapsing into clumps of mass ~M
        f = U[-1]
        fr = {k: dfrac[m].sum() for k, m in win.items()}
        ok = (1e-28 < f < 1e-13) and (1e-32 < fr['>6e22 (permanent PBH)'] < 1e-3)
        good.append((A, dc, ok))
        print(f"{A:<7.2f} {dc:<5.2f} {f:9.2e}  " + "  ".join(f"{fr[k]:24.2e}" for k in win) + ("   <-- window" if ok else ""))
print("configurations in the window:", [(a, d) for a, d, o in good if o])
