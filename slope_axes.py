"""Does a curvature slope along n (our offset from our own seed) line up with the sky?
A gradient-type multipole Y_l0 about n has zero m^2-dispersion about n and maximal about any axis
perpendicular to n (Wigner d^2_{20}(90deg)^2*2 = 3/4), so the de Oliveira-Costa 'axis' of the
quadrupole/octopole is predicted to lie on the great circle at 90 deg from n, azimuth unconstrained.
Tests, taking n = 2M++ void centroid: (a) Axis of Evil on the ring, (b) CF4 bulk flow on the ring? no -
the flow should be ALONG n (gradient), so test (b) is 'CF4 within angle of n'. Combined with Fisher's
method, with and without the void-centroid detection itself; dependencies stated."""
import numpy as np
from scipy.stats import chi2, norm
def u(l,b):
    l,b=np.radians(l),np.radians(b); return np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
n=u(305,-30)                      # 2M++ void centroid (this paper, Sec. 7)
aoe=u(260,60); cf4=u(298,-8); qso=u(238,29); dip=u(264,48)
def ang(a,b): return np.degrees(np.arccos(abs(a@b)))            # axis (sign-free)
def angd(a,b): return np.degrees(np.arccos(a@b))                 # direction
tests={}
# (a) AoE axis on the ring perpendicular to n: p = P(random axis within delta of the ring) = sin(delta)
d=90-ang(n,aoe); tests['Axis of Evil on ring (|90-sep|=%.0f deg)'%d]=np.sin(np.radians(d))
# (b) CF4 bulk flow along n: p = P(random direction within theta) = (1-cos theta)/2
th=angd(n,cf4); tests['CF4 bulk flow within %.0f deg of n'%th]=(1-np.cos(np.radians(th)))/2
# (c) quasar-dipole excess along n (only if isocurvature): direction test
hpa=u(227,-27); th=angd(n,hpa); tests['power-asymmetry axis within %.0f deg of n'%th]=(1-np.cos(np.radians(th)))/2
th=angd(n,qso); tests['CatWISE dipole within %.0f deg of n'%th]=(1-np.cos(np.radians(th)))/2
def fisher(ps):
    X=-2*np.sum(np.log(ps)); p=chi2.sf(X,2*len(ps)); return p, norm.isf(p)
print("individual tests (n = 2M++ void centroid):")
for k,v in tests.items(): print(f"  {k:45s} p={v:.3f}  ({norm.isf(v):.1f} sigma one-sided)")
keys=list(tests)
for label,sel in (("AoE + CF4",keys[:2]),("AoE + CF4 + power asymmetry",keys[:3]),("AoE + CF4 + PA + quasar",keys),("AoE + CF4 + void centroid (p=0.05, pre look-elsewhere)",keys[:2]+['void']),
                  ("AoE + CF4 + void centroid (p=0.3, post look-elsewhere)",keys[:2]+['void2'])):
    ps=[tests[k] if k in tests else (0.05 if k=='void' else 0.3) for k in sel]
    p,s=fisher(ps); print(f"Fisher {label:58s}: p={p:.3f}  {s:.1f} sigma")
print("caveats: CF4 flow and 2M++ voids trace the same local structure (not independent); the quasar test needs an\n"
      "isocurvature slope; the ring test has no azimuth prediction and does not explain the l=2/l=3 mutual alignment.")

# --- post-hoc cross-product NOTE (not counted in the paper) (added after the 87-deg quasar result was seen) ---
cs=u(203,-56)
x=np.cross(n,cs); x/=np.linalg.norm(x); d=ang(x,qso)
p_raw=1-np.cos(np.radians(d)); p_le=min(1,9*p_raw)   # trials: 3 cross products x 3 target axes tried/possible
print(f"\n(n x ColdSpot) axis vs CatWISE dipole axis: {d:.0f} deg; p_raw={p_raw:.4f}, x9 trials p={p_le:.4f}")
for label,ps in (("AoE + CF4 + cross (raw)",[tests[keys[0]],tests[keys[1]],p_raw]),("AoE + CF4 + cross (trials x9)",[tests[keys[0]],tests[keys[1]],p_le])):
    p,s=fisher(ps); print(f"Fisher {label:40s}: p={p:.4f}  {s:.1f} sigma")
print("caveat: a perpendicular dipole needs a parity-odd (rotation x gradient) term; global rotation is bounded at omega/H < 1e-9,\n"
      "so the geometry is matched but no mechanism supplies the amplitude.")
