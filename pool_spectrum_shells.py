"""The frequency curve as we move out: angular power spectrum of the projected mass in thin shells at increasing
distance (Limber: C_l(z) = P(k=(l+1/2)/chi(z), z) / chi(z)^2 per unit comoving depth). Multipole l is angular
frequency: a feature of angular size theta sits at l ~ 180 deg/theta. Moving out does two things: the same physical
wavelength k appears at higher l (smaller angle) because chi grows, and the amplitude drops as the growth factor D(z)^2
because the shell is younger. The all-mass projection (Planck lensing) is the sum over shells weighted by the lensing
kernel, which is why it shows one blurred curve. Output: figures/pool_spectrum_shells.png"""
import numpy as np, os
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from scipy.integrate import quad
HERE = os.path.dirname(os.path.abspath(__file__))
src = open(os.path.join(HERE, 'pool_spectrum.py')).read()
g_ = {'np': np}; exec(src[src.index('h, Om, Ob, ns, Tcmb'):src.index('k = np.logspace(-4, 1.5, 1200)')], g_)
eh98, Om, ns, h = g_['eh98'], g_['Om'], g_['ns'], g_['h']
k = np.logspace(-4, 2, 3000); P0 = k**ns*eh98(k)**2
def sigmaR(P, R):
    x = k*R; W = 3*(np.sin(x)-x*np.cos(x))/x**3; return np.sqrt(np.trapz(k**2*P*W**2, k)/(2*np.pi**2))
P0 *= (0.81/sigmaR(P0, 8.0))**2
def growth(z):
    a = 1/(1+z); E = lambda a: np.sqrt(Om/a**3+(1-Om)); return 2.5*Om*E(a)*quad(lambda x: 1/(x*E(x))**3, 1e-6, a)[0]
D0 = growth(0)
def chi(z):  # comoving distance in Mpc/h
    return 2997.9*quad(lambda zz: 1/np.sqrt(Om*(1+zz)**3+(1-Om)), 0, z)[0]
l = np.arange(2, 3000)
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
shells = [(0.05, 'z=0.05, 150 Mpc/h (the local sea)'), (0.3, 'z=0.3'), (1.0, 'z=1'), (3.0, 'z=3'), (10.0, 'z=10 (21-cm era)'), (1100.0, 'z=1100 (last scattering, linear)')]
for z, lab in shells:
    c = chi(z); D = growth(z)/D0
    kk = (l+0.5)/c; Pz = np.interp(kk, k, P0, right=0)*D**2
    cl = Pz/c**2
    ax[0].loglog(l, cl, label=f'{lab}: $\\chi$={c:.0f} Mpc/h, D²={D**2:.2g}')
    ax[1].loglog(l, cl/np.max(cl), label=lab)
ax[0].set_xlabel('multipole $\\ell$ (angular frequency: $\\ell \\approx 180°/\\theta$)'); ax[0].set_ylabel('$C_\\ell$ per unit depth: energy per unit angular frequency (arbitrary)')
ax[0].set_title('Projected-mass spectrum of a thin shell at increasing radius'); ax[0].legend(fontsize=7); ax[0].grid(alpha=.3); ax[0].set_ylim(1e-14, 1e-2)
ax[1].set_xlabel('multipole $\\ell$'); ax[1].set_ylabel('same, each shell normalised to its peak'); ax[1].set_title('Shape: the peak walks to higher $\\ell$ as we move out')
ax[1].legend(fontsize=7); ax[1].grid(alpha=.3); ax[1].set_ylim(1e-3, 2)
for a in ax: a.axvline(180/1.0, color='k', ls=':', lw=1); a.text(185, a.get_ylim()[0]*3, '1°', fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(HERE, 'figures', 'pool_spectrum_shells.png'), dpi=130)
for z, lab in shells:
    c = chi(z); kk = (l+0.5)/c; Pz = np.interp(kk, k, P0, right=0)*(growth(z)/D0)**2; cl = Pz/c**2
    print(f"{lab:38s} chi={c:6.0f} Mpc/h  peak at l={l[np.argmax(cl)]:5d} (theta={180/l[np.argmax(cl)]:.1f} deg)  amplitude x{(growth(z)/D0)**2:.2g}")
