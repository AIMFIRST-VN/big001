"""Fit the slope direction n instead of assuming it: for every n on a HEALPix grid, combine (Fisher) the p-values of
the pre-stated tests -- Axis of Evil on the ring perpendicular to n; CF4 bulk flow, hemispherical power asymmetry
and the cluster H0 dipole along n -- and find the best n. Look-elsewhere: repeat with the four data axes rotated
randomly (jointly, preserving their mutual geometry? no -- independently, which is the null of 'no common axis')."""
import numpy as np, healpy as hp
from scipy.stats import chi2, norm
def u(l,b):
    l,b=np.radians(l),np.radians(b); return np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
aoe=u(260,60); along={'CF4':u(298,-8),'HPA':u(227,-27),'clusterH0':u(280,-15)}
NS=32; npix=hp.nside2npix(NS); V=np.array(hp.pix2vec(NS,np.arange(npix))).T
def score(aoe,along):
    d=np.degrees(np.arccos(np.abs(V@aoe))); p_ring=np.sin(np.radians(np.abs(90-d)))
    X=-2*np.log(np.clip(p_ring,1e-12,1))
    for v in along.values():
        th=np.arccos(np.clip(V@v,-1,1)); X+=-2*np.log(np.clip((1-np.cos(th))/2,1e-12,1))
    return chi2.sf(X,2*(1+len(along)))
p=score(aoe,along); i=np.argmin(p); l,b=hp.pix2ang(NS,i,lonlat=True)
print(f"best-fit slope direction: (l,b)=({l:.0f},{b:.0f}), joint p={p[i]:.2e} ({norm.isf(p[i]):.1f} sigma before look-elsewhere)")
nv=u(305,-30); print(f"  angle to void centroid (305,-30): {np.degrees(np.arccos(V[i]@nv)):.0f} deg; joint p at void centroid = {score(aoe,along)[hp.vec2pix(NS,*nv)]:.3f}")
# 68% region size
print(f"  sky fraction with p < 3*p_best: {np.mean(p<3*p[i]):.3f}")
rng=np.random.default_rng(0); best=[]
for _ in range(400):
    R=lambda: (lambda v: v/np.linalg.norm(v))(rng.normal(size=3))
    best.append(score(R(),{k:R() for k in along}).min())
best=np.array(best); ple=np.mean(best<=p[i])
print(f"look-elsewhere (400 nulls with four random axes): P(best p <= observed) = {ple:.3f}  -> {norm.isf(max(ple,1e-4)):.1f} sigma global")
