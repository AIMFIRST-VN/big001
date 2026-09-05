"""Different contributions to the waves: decompose the spectrum into its components, each with its energy share and
the epoch (time spread) that made it. Smooth growing mode (the swell set at the bounce, all epochs); baryon acoustic
component (sound waves frozen at the drag epoch z~1060, a wiggle of a few percent on top of the smooth spectrum);
isocurvature (a second wave system: Planck bounds its fraction below a few percent at k~0.05/Mpc); tensor waves
(r < 0.036, BICEP/Keck+Planck: gravitational-wave power relative to scalar); and the neighbour tail of P9 (a fixed
l<=3 template carrying ~20% of the quadrupole power). Output: figures/wave_components.png"""
import numpy as np, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.signal import savgol_filter
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'pool_spectrum.py')).read()
g_ = {'np': np}; exec(src[src.index('h, Om, Ob, ns, Tcmb'):src.index('k = np.logspace(-4, 1.5, 1200)')], g_)
eh98, ns = g_['eh98'], g_['ns']
k = np.logspace(-3, 0.5, 2000); P = k**ns*eh98(k)**2; P /= P.max()
lnP = np.log(P); smooth = np.exp(savgol_filter(lnP, 301, 3))      # no-wiggle envelope by smoothing in log k
wig = P/smooth - 1
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
ax[0].loglog(k, P, 'k-', lw=2, label='total: growing mode (the swell, set at the bounce)')
ax[0].loglog(k, smooth, 'r--', lw=1, label='smooth part (no acoustic component)')
ax[0].loglog(k, smooth*0.03, 'c-', lw=1, label='isocurvature bound: < 3% of scalar power (Planck 2018 X)')
ax[0].loglog(k, smooth*0.036*0.5, 'm-', lw=1, label='tensor bound: r < 0.036 (BICEP/Keck 2021), shown as r/2 of scalar')
ax[0].set_xlabel('frequency k [h/Mpc]'); ax[0].set_ylabel('energy per unit frequency (peak = 1)'); ax[0].set_ylim(1e-5, 2)
ax[0].set_title('Wave systems in the matter spectrum: what is there and what is bounded'); ax[0].legend(fontsize=7, loc='lower left'); ax[0].grid(alpha=.3)
ax[1].semilogx(k, 100*wig, 'b-', lw=1.5, label='acoustic component: sound waves frozen at the drag epoch (z ≈ 1060)')
ax[1].axhline(0, color='k', lw=0.5); ax[1].set_xlim(0.01, 0.5); ax[1].set_ylim(-12, 12)
ax[1].set_xlabel('frequency k [h/Mpc]'); ax[1].set_ylabel('wiggle amplitude, % of the smooth spectrum'); ax[1].set_title('The fossil swell: BAO ringing at the sound-horizon frequency, 1/(147 Mpc)')
ax[1].legend(fontsize=8); ax[1].grid(alpha=.3)
plt.tight_layout(); plt.savefig(os.path.join(HERE, 'figures', 'wave_components.png'), dpi=130)
m = (k > 0.03) & (k < 0.3)
print(f"acoustic component: rms wiggle {100*np.sqrt(np.mean(wig[m]**2)):.1f}% of the smooth spectrum over k=0.03-0.3, first crest at k={k[m][np.argmax(wig[m])]:.3f} h/Mpc; energy share of the wiggles ~ {np.mean(wig[m]**2)*100:.2f}% of the power in that band")
print("component | energy share | epoch (time spread) | frequency spread")
print("growing mode (swell)        | ~100%       | set at the bounce, grows always     | all k, peak at k_eq")
print("acoustic (BAO)              | few % modulation, <0.1% energy | frozen at z~1060 (drag), width ~ Silk scale | k_BAO = 2pi/147 Mpc and harmonics")
print("isocurvature (2nd system)   | < 3% (bound) | would be set at the bounce           | scale-dependent; not detected")
print("tensor (gravitational waves)| < 1.8% (r<0.036)| bounce/horizon exit              | red, l<100; not detected")
print("vector (vorticity)          | decaying; omega/H < 1e-9 | storm only               | not detected")
print("neighbour tail (P9)         | ~20% of C_2, 8% of C_3 (template) | fossil, laid before the bounce ends | l<=3 only")
print("secondaries (lensing, SZ, ISW) | few % of CMB power at l>1000 | z<10                | high l")
