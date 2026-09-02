"""Falsifiable predictions of the 'surrounding neighbours' picture (H4/H7).
Each spot is modelled as the Sachs-Wolfe imprint of a point mass just beyond the
last-scattering sphere at comoving D=(1+eps)d.  Profile: f(t)=1/sqrt(eps^2+2(1+eps)(1-cos t)).
Half-max angle th -> eps = th/sqrt(3).  The same mass produces a 1/r tail over the
whole sky -> predicted l=1..3 pattern, and a tidal pull -> predicted bulk-flow direction.
Compare with the observed SMICA low multipoles and the CF4 bulk flow."""
import os
HERE = os.path.dirname(os.path.abspath(__file__))
import numpy as np, healpy as hp

spots = [  # l, b, dT peak (uK), half-max (deg)
 (157.1,-70.5,+170,2.5),(79.5,-33.2,-170,3.5),(203.2,-56.3,-168,2.5),(155.0,-29.3,-163,5.5),
 (304.5,-29.0,+153,2.5),(210.6,-35.0,-153,2.5),(170.3,-46.6,+162,4.5),(184.2,-54.3,+155,4.5)]
NS=64; npix=hp.nside2npix(NS); vec=np.array(hp.pix2vec(NS,np.arange(npix))).T
def unit(l,b):
    l,b=np.radians(l),np.radians(b); return np.array([np.cos(b)*np.cos(l),np.cos(b)*np.sin(l),np.sin(b)])
def tail(l,b,dT,th):
    eps=np.radians(th)/np.sqrt(3); c=vec@unit(l,b)
    f=1/np.sqrt(eps**2+2*(1+eps)*(1-c)); return dT*eps*f, eps
tot=np.zeros(npix); g=np.zeros(3)
print("spot     eps    D/d    tail@90deg(uK)  tail@antipode(uK)")
for l,b,dT,th in spots:
    m,eps=tail(l,b,dT,th); tot+=m
    print("(%5.1f,%5.1f) %.4f %.4f  %6.2f  %6.2f"%(l,b,eps,1+eps,dT*eps/np.sqrt(2*(1+eps)),dT*eps/np.sqrt(eps**2+4*(1+eps))))
    # tidal pull toward mass: M ∝ dT*eps*d (relative), pull ∝ M/D^2 ~ dT*eps ; sign: cold=attractive, hot=repulsive(negative mass)
    g+= (dT*eps/(1+eps)**2)*unit(l,b)*(-1)   # cold (dT<0) -> +toward
alm=hp.map2alm(tot,lmax=3)
cl=hp.alm2cl(alm)
print("predicted C1,C2,C3 (uK^2): %.2f %.2f %.2f  -> rms l=2: %.2f uK, l=3: %.2f uK"%(cl[1],cl[2],cl[3],np.sqrt(5*cl[2]/(4*np.pi)),np.sqrt(7*cl[3]/(4*np.pi))))
gl=np.degrees(np.arctan2(g[1],g[0]))%360; gb=np.degrees(np.arcsin(g[2]/np.linalg.norm(g)))
print("predicted bulk-flow (tidal pull) direction: l=%.0f b=%.0f"%(gl,gb))
cf4=unit(297,-6); print("angle to CF4 bulk flow (297,-6): %.0f deg"%np.degrees(np.arccos(g@cf4/np.linalg.norm(g))))
# axis of predicted l=2,3 vs observed SMICA l=2,3
def axis(alm,l):
    # maximize angular-momentum dispersion: brute force over directions
    best=None
    for i in range(hp.nside2npix(16)):
        n=hp.pix2vec(16,i); th,ph=hp.vec2ang(np.array(n))
        r=hp.Rotator(rot=[np.degrees(ph[0]),90-np.degrees(th[0]),0],deg=True,inv=True)
        a=r.rotate_alm(alm.copy())
        s=sum(m*m*abs(a[hp.Alm.getidx(3,l,m)])**2*(2 if m>0 else 1) for m in range(l+1))
        if best is None or s>best[0]: best=(s,i)
    return hp.pix2ang(16,best[1],lonlat=True)
sm=np.load(os.path.join(HERE, "smica_nside128.npy"))
if np.std(sm)<1e-2: sm*=1e6
almo=hp.map2alm(sm-sm.mean(),lmax=3); clo=hp.alm2cl(almo)
print("observed (unmasked SMICA) C2,C3: %.0f %.0f uK^2"%(clo[2],clo[3]))
for l in (2,3):
    ap=axis(alm,l); ao=axis(almo,l)
    ang=np.degrees(np.arccos(abs(unit(*ap)@unit(*ao))))
    print("l=%d axis: predicted (%.0f,%.0f)  observed (%.0f,%.0f)  separation %.0f deg (axes, <=90)"%(l,ap[0],ap[1],ao[0],ao[1],ang))
# cross-correlation of predicted tail template with observed l=2,3 maps
for l in (2,3):
    fp=np.zeros(4); fp[l]=1
    mp=hp.alm2map(hp.almxfl(alm,fp),NS); mo=hp.alm2map(hp.almxfl(almo,fp),NS)
    print("l=%d template/observed correlation r=%.2f"%(l,np.corrcoef(mp,mo)[0,1]))
# flip sign convention (hot attractive) for the bulk-flow prediction
g2=-g; print("sign-flipped bulk-flow direction: l=%.0f b=%.0f, angle to CF4 %.0f deg"%(np.degrees(np.arctan2(g2[1],g2[0]))%360,np.degrees(np.arcsin(g2[2]/np.linalg.norm(g2))),np.degrees(np.arccos(g2@cf4/np.linalg.norm(g2)))))
