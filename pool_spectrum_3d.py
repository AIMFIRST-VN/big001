"""3-D view of the energy-frequency spectrum with time as the third axis.
Left: the universe -- linear matter power spectrum P(k, z) (Eisenstein-Hu 1998 transfer, growth factor D(z) for Planck-2018
LCDM), z from 1000 to 0, log energy vs log frequency vs log(1+z); the BAO ridge and the fixed peak at k_eq are visible, the
whole surface rises as D(z)^2 and the non-linear scale (where Delta^2 = 1, marked) walks to lower frequency with time.
Right: a sea developing under a steady wind -- JONSWAP spectrum as a function of fetch (Hasselmann et al. 1973):
f_p = 3.5 (g/U10) (g x / U10^2)^-0.33, alpha = 0.076 (g x/U10^2)^-0.22; the peak walks to lower frequency and the energy
grows with fetch, the ocean's version of growth with time. Output: figures/pool_spectrum_3d.png"""
import numpy as np, os, sys
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D  # noqa
HERE = os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0, HERE)
import importlib.util
spec = importlib.util.spec_from_file_location('ps', os.path.join(HERE, 'pool_spectrum.py'))
# reuse the EH98 transfer function without executing the plotting script: copy the function body
src = open(os.path.join(HERE, 'pool_spectrum.py')).read()
start = src.index('h, Om, Ob, ns, Tcmb'); end = src.index('k = np.logspace(-4, 1.5, 1200)')
g_ = {'np': np}; exec(src[start:end], g_); eh98 = g_['eh98']; Om = g_['Om']; ns = g_['ns']
from scipy.integrate import quad
def growth(z):
    a = 1/(1+z); E = lambda a: np.sqrt(Om/a**3 + (1-Om))
    I = quad(lambda x: 1/(x*E(x))**3, 1e-6, a)[0]; return 2.5*Om*E(a)*I
k = np.logspace(-3, 1.2, 300); zs = np.logspace(0, 3, 40)-1+1e-9   # 0..999
P0 = k**ns*eh98(k)**2
# normalise: sigma_8 = 0.81 today (top-hat) with h=0.6736
def sigmaR(P, R):
    x = k*R; W = 3*(np.sin(x)-x*np.cos(x))/x**3; return np.sqrt(np.trapz(k**2*P*W**2, k)/(2*np.pi**2))
P0 *= (0.81/sigmaR(P0, 8.0))**2
D0 = growth(0); Ds = np.array([growth(z)/D0 for z in zs])
Z = np.log10(np.outer(Ds**2, P0))           # (nz, nk)
KK, ZZ = np.meshgrid(np.log10(k), np.log10(1+zs))
# non-linear scale: Delta^2 = k^3 P/2pi^2 = 1
knl = []
for D in Ds:
    d2 = k**3*D**2*P0/(2*np.pi**2); i = np.where(d2 >= 1)[0]; knl.append(k[i[0]] if len(i) else np.nan)
knl = np.array(knl)
fig = plt.figure(figsize=(14, 6))
ax = fig.add_subplot(1, 2, 1, projection='3d')
ax.plot_surface(KK, ZZ, Z, cmap='viridis', rstride=1, cstride=4, linewidth=0, antialiased=True, alpha=0.95)
m = np.isfinite(knl); ax.plot(np.log10(knl[m]), np.log10(1+zs[m]), [np.interp(np.log10(kn), np.log10(k), Z[i]) for i, kn in zip(np.where(m)[0], knl[m])], 'r-', lw=2, label='non-linear scale (breaking line)')
ax.set_xlabel('log10 frequency  k [h/Mpc]'); ax.set_ylabel('log10(1+z)  (time runs toward 0)'); ax.set_zlabel('log10 energy per unit frequency  P(k,z)')
ax.set_title('The universe: spectrum growing with time (linear theory)'); ax.view_init(elev=28, azim=-125); ax.legend(loc='upper left', fontsize=8)
# ocean: JONSWAP vs fetch
g, U10 = 9.81, 15.0
f = np.logspace(-1.3, 0.3, 300); X = np.logspace(3, 6, 40)      # fetch 1 km .. 1000 km
S = []
for x in X:
    xt = g*x/U10**2; fp = 3.5*(g/U10)*xt**-0.33; alpha = 0.076*xt**-0.22
    sig = np.where(f < fp, 0.07, 0.09); r = np.exp(-(f-fp)**2/(2*sig**2*fp**2))
    S.append(alpha*g**2/(2*np.pi)**4*f**-5*np.exp(-1.25*(fp/f)**4)*3.3**r)
S = np.log10(np.maximum(np.array(S), 1e-8))
FF, XX = np.meshgrid(np.log10(f), np.log10(X/1e3))
ax2 = fig.add_subplot(1, 2, 2, projection='3d')
ax2.plot_surface(FF, XX, S, cmap='Blues_r', rstride=1, cstride=4, linewidth=0, antialiased=True, alpha=0.95)
fps = 3.5*(g/U10)*(g*X/U10**2)**-0.33
ax2.plot(np.log10(fps), np.log10(X/1e3), [S[i][np.argmin(np.abs(f-fp))] for i, fp in enumerate(fps)], 'r-', lw=2, label='peak (swell walks to lower frequency)')
ax2.set_xlabel('log10 frequency  f [Hz]'); ax2.set_ylabel('log10 fetch [km]  (time under the wind)'); ax2.set_zlabel('log10 energy per unit frequency  S(f) [m$^2$/Hz]')
ax2.set_title('The sea: JONSWAP spectrum developing with fetch (wind 15 m/s)'); ax2.view_init(elev=28, azim=-125); ax2.legend(loc='upper left', fontsize=8)
plt.tight_layout(); plt.savefig(os.path.join(HERE, 'figures', 'pool_spectrum_3d.png'), dpi=130)
print(f"non-linear scale today k_nl = {knl[0]:.3f} h/Mpc; at z=1: {np.interp(1, zs, np.nan_to_num(knl, nan=99)):.2f}; enters the grid (k<16) at z ~ {zs[np.isfinite(knl)].max():.0f}")
print(f"sea: peak frequency 1 km -> 1000 km fetch: {fps[0]:.2f} -> {fps[-1]:.3f} Hz; peak energy grows x{10**(S[-1].max()-S[0].max()):.0f}")
