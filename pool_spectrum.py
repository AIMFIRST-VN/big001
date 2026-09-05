"""Energy-versus-frequency spectrum of the universe's mass, projected onto a plane, compared with ocean-wave spectra.
Cosmic side: Eisenstein & Hu (1998) linear matter power spectrum (with baryon wiggles), Planck-2018 parameters, z=0;
  3-D 'energy per log frequency' Delta^2(k) = k^3 P(k)/2pi^2, and the 2-D projection of all mass onto a plane,
  P_2D(K) = int P(sqrt(K^2+kz^2)) dkz/2pi (slab projection), shown as K^2 P_2D/2pi.
  Measured projected mass: the Planck PR4 lensing convergence map (kappa, nside 256) with its mask, C_l^kk via anafast,
  shown as l(l+1)C_l/2pi (noise-dominated above l ~ 300; band-averaged).
Ocean side: Pierson-Moskowitz fully developed sea (no land, unlimited fetch), S(f) ~ f^-5 exp(-1.25 (f_p/f)^4);
  JONSWAP fetch-limited sea (finite pool, sharper peak, gamma = 3.3); and a sloping-seabed (shoaling) note:
  Green's law amplitude ~ h^-1/4 -> a depth slope makes an amplitude gradient across the pool, the analogue of the
  hemispherical power asymmetry / kernel gradient P12.
All spectra normalised to unit peak and plotted against frequency / peak frequency. Output: figures/pool_spectrum.png"""
import numpy as np, os, healpy as hp
import matplotlib; matplotlib.use('Agg'); import matplotlib.pyplot as plt
HERE = os.path.dirname(os.path.abspath(__file__))
# --- Eisenstein-Hu 1998 transfer function (with wiggles) ---
h, Om, Ob, ns, Tcmb = 0.6736, 0.3153, 0.0493, 0.9649, 2.7255
def eh98(k):   # k in h/Mpc
    om, ob = Om*h*h, Ob*h*h; th = Tcmb/2.7
    zeq = 2.5e4*om/th**4; keq = 7.46e-2*om/th**2
    b1 = 0.313*om**-0.419*(1+0.607*om**0.674); b2 = 0.238*om**0.223
    zd = 1291*om**0.251/(1+0.659*om**0.828)*(1+b1*ob**b2)
    Req = 31.5*ob/th**4/(zeq/1e3); Rd = 31.5*ob/th**4/(zd/1e3)
    s = 2/(3*keq)*np.sqrt(6/Req)*np.log((np.sqrt(1+Rd)+np.sqrt(Rd+Req))/(1+np.sqrt(Req)))
    ksilk = 1.6*ob**0.52*om**0.73*(1+(10.4*om)**-0.95)
    q = k*h/(13.41*keq)
    a1 = (46.9*om)**0.670*(1+(32.1*om)**-0.532); a2 = (12.0*om)**0.424*(1+(45.0*om)**-0.582)
    ac = a1**(-ob/om)*a2**(-(ob/om)**3)
    bb1 = 0.944/(1+(458*om)**-0.708); bb2 = (0.395*om)**-0.0266
    bc = 1/(1+bb1*((1-ob/om)**bb2-1))
    def T0(q, a, b):
        C = 14.2/a+386/(1+69.9*q**1.08); return np.log(np.e+1.8*b*q)/(np.log(np.e+1.8*b*q)+C*q*q)
    f = 1/(1+(k*h*s/5.4)**4)
    Tc = f*T0(q,1,bc)+(1-f)*T0(q,ac,bc)
    y = (1+zeq)/(1+zd); G = y*(-6*np.sqrt(1+y)+(2+3*y)*np.log((np.sqrt(1+y)+1)/(np.sqrt(1+y)-1)))
    ab = 2.07*keq*s*(1+Rd)**-0.75*G
    bnode = 8.41*om**0.435; st = s/(1+(bnode/(k*h*s))**3)**(1/3)
    bbb = 0.5+ob/om+(3-2*ob/om)*np.sqrt((17.2*om)**2+1)
    Tb = (T0(q,1,1)/(1+(k*h*s/5.2)**2)+ab/(1+(bbb/(k*h*s))**3)*np.exp(-(k*h/ksilk)**1.4))*np.sinc(k*h*st/np.pi)
    return ob/om*Tb+(1-ob/om)*Tc
k = np.logspace(-4, 1.5, 1200)
P = k**ns*eh98(k)**2; P *= 1/np.max(P)               # arbitrary normalisation
# ocean S(f) is variance per unit frequency; the cosmic counterpart per unit wavenumber is P(k) itself (peak at k_eq).
D2 = P.copy()
# 1-D cut through the field (what a single line of sight or a single buoy sees): P_1D(k) = (1/2pi) int_k^inf P(k') k' dk'
P1 = np.array([np.trapz((P*k)[k>=kk], k[k>=kk]) for kk in k])/(2*np.pi)
# 2-D slab projection
K = np.logspace(-3.5, 1, 300); P2 = []
for KK in K:
    kz = np.logspace(-5, 2, 4000); kk = np.sqrt(KK**2+kz**2)
    P2.append(2*np.trapz(np.interp(kk, k, P, right=0), kz)/(2*np.pi))
P2 = np.array(P2); D2_2d = P2.copy()
# --- measured projected mass: Planck PR4 kappa ---
kap = np.load(os.path.join(HERE, 'pr4_kappa_ns256.npy')); mask = hp.read_map(os.path.join(HERE, 'PR4_variations', 'mask.fits.gz'), verbose=False) if os.path.exists(os.path.join(HERE,'PR4_variations','mask.fits.gz')) else np.ones_like(kap)
mask = hp.ud_grade(mask, 256); fsky = mask.mean()
cl = hp.anafast(kap*mask, lmax=500)/fsky; l = np.arange(len(cl)); dl = cl.copy()
lb = np.array([8,16,32,64,128,256,400]); dlb = [dl[(l>=a)&(l<b)].mean() for a,b in zip(lb[:-1], lb[1:])]; lbc = np.sqrt(lb[:-1]*lb[1:])
# --- ocean spectra ---
f = np.logspace(-1, 1.5, 600)   # f/f_p
PM = f**-5*np.exp(-1.25*f**-4)
sig = np.where(f<1, 0.07, 0.09); r = np.exp(-(f-1)**2/(2*sig**2*1)); JS = PM*3.3**r
norm = lambda a: a/np.nanmax(a)
# peaks
kp = k[np.argmax(D2)]; Kp = K[np.argmax(D2_2d)]; lp = lbc[np.argmax(dlb)]; k1p = kp   # 1-D cut has no peak (monotone); plot it on the 3-D peak scale
fig, ax = plt.subplots(1, 2, figsize=(13, 5.2))
ax[0].loglog(f, norm(PM), 'b-', lw=2, label='ocean, no land (Pierson–Moskowitz): $f^{-5}$ tail')
ax[0].loglog(f, norm(JS), 'c--', lw=1.5, label='ocean, finite pool (JONSWAP, fetch-limited)')
ax[0].loglog(k/kp, norm(D2), 'k-', lw=2, label='universe, 3-D: $P(k)$, tail $\\propto k^{-3}\\ln^2 k$')
ax[0].loglog(K/Kp, norm(D2_2d), 'r-', lw=2, label='universe projected on a plane: $P_{2D}(K)$')
ax[0].loglog(k/k1p, norm(P1), 'g-', lw=1.5, label='universe, one line cut (a single buoy): $P_{1D}(k)$')
ax[0].axvspan(0.1, 1/(kp*1.4e4)*1e4 if False else 0.1, alpha=0); ax[0].axvline(7e-5/kp, color='k', ls=':', lw=1); ax[0].text(7e-5/kp*1.2, 1.2e-4, 'horizon: longer\nwaves unobservable', fontsize=7)
ax[0].axvline(1.0/kp, color='k', ls=':', lw=1); ax[0].text(1.05/kp, 1.2e-4, 'non-linear /\nresolution', fontsize=7)
ax[0].set_xlim(2e-3, 100); ax[0].set_ylim(1e-4, 2); ax[0].set_xlabel('frequency / peak frequency'); ax[0].set_ylabel('energy per unit frequency (peak = 1)')
ax[0].set_title('Energy–frequency spectra: ocean vs the mass of the universe'); ax[0].legend(fontsize=8, loc='lower left'); ax[0].grid(alpha=.3)
ax[1].loglog(l[2:], dl[2:], color='0.7', lw=0.8, label='Planck PR4 lensing $\\kappa$, raw $C_\\ell$ (noise above $\\ell\\sim300$)')
ax[1].loglog(lbc, dlb, 'ro-', label='band averages (measured projected mass)')
ax[1].loglog(K/Kp*lp, norm(D2_2d)*np.nanmax(dlb), 'r--', lw=1, label='projected linear theory, scaled to peak')
ax[1].set_xlabel('multipole $\\ell$ (angular frequency)'); ax[1].set_ylabel('$C_\\ell^{\\kappa\\kappa}$ (energy per unit angular frequency)'); ax[1].set_title('All mass projected on the sky: Planck lensing convergence')
ax[1].legend(fontsize=8, loc='lower left'); ax[1].grid(alpha=.3); ax[1].set_xlim(2, 500)
plt.tight_layout(); plt.savefig(os.path.join(HERE, 'figures', 'pool_spectrum.png'), dpi=130)
# slopes
def slope(x, y, lo, hi):
    m = (x>lo)&(x<hi)&(y>0); return np.polyfit(np.log(x[m]), np.log(y[m]), 1)[0] if m.sum()>2 else float('nan')
print(f"peak: k_p = {kp:.4f} h/Mpc (3-D), K_p = {Kp:.4f} h/Mpc (2-D), lensing band peak l ~ {lp:.0f}; fsky = {fsky:.2f}")
print(f"high-frequency slopes (f/f_p = 10..100): ocean PM {slope(f,PM,10,100):.2f}, universe 3-D {slope(k/kp,D2,10,100):.2f}, projected {slope(K/Kp,D2_2d,10,100):.2f}, 1-D cut {slope(k/k1p,P1,10,100):.2f}")
print(f"low-frequency slopes (f/f_p = 0.2..0.5): ocean PM {slope(f,PM,0.2,0.5):.2f} (exp cut-off), universe 3-D {slope(k/kp,D2,0.2,0.5):.2f} (k^n_s), projected {slope(K/Kp,D2_2d,0.2,0.5):.2f}, 1-D {slope(k/k1p,P1,0.2,0.5):.2f}")
print(f"horizon cut: k_H = 7e-5 h/Mpc = {7e-5/kp:.4f} of the peak frequency; everything to the left is the unobservable low-frequency tail")
print("shoaling (Green's law A ~ h^-1/4): a 7% amplitude asymmetry across the sky (Planck hemispherical asymmetry, l<60) corresponds to a 28% depth change across the pool")
