"""'Smudging' test for H5 (Where Is Nemo): smooth the 2M++ density field (Carrick et al. 2015) with a Gaussian
of radius R_s and follow the mass-centroid offset (dipole of the smoothed field) within R < 150 Mpc/h as R_s grows.
A genuine large-scale radial gradient keeps its direction and amplitude as R_s grows; local structure's dipole
falls and its direction wanders. Null: Gaussian random fields with the data's own power spectrum, same mask."""
import numpy as np, os
from scipy.ndimage import gaussian_filter
HERE = os.path.dirname(os.path.abspath(__file__)); rng = np.random.default_rng(1)
den = np.load(os.path.join(HERE, 'twompp_density.npy')); ax = np.linspace(-200, 200, 257); h = ax[1] - ax[0]
X, Y, Z = np.meshgrid(ax, ax, ax, indexing='ij'); R = np.sqrt(X**2 + Y**2 + Z**2)
RMAX = 150.0; sel = (R < RMAX) & (R > 5)
def lb(u): return np.degrees(np.arctan2(u[1], u[0])) % 360, np.degrees(np.arcsin(u[2] / np.linalg.norm(u)))
def centroid(field, Rs):
    fs = gaussian_filter(field, Rs / h, mode='constant') if Rs > 0 else field
    w = 1 + fs[sel]; c = np.array([np.sum(w * X[sel]), np.sum(w * Y[sel]), np.sum(w * Z[sel])]) / np.sum(w)
    return c
# power spectrum of the data (for the null)
F = np.fft.rfftn(den); k = np.sqrt(sum(np.meshgrid(*[np.fft.fftfreq(257, h)]*2, np.fft.rfftfreq(257, h), indexing='ij'))[0]**2 for _ in [0]) if False else None
kx = np.fft.fftfreq(257, h); kz = np.fft.rfftfreq(257, h)
KX, KY, KZ = np.meshgrid(kx, kx, kz, indexing='ij'); K = np.sqrt(KX**2 + KY**2 + KZ**2)
P2 = np.abs(F)**2; kb = np.logspace(np.log10(K[K > 0].min()), np.log10(K.max()), 40)
idx = np.clip(np.digitize(K, kb), 0, len(kb) - 1); Pk = np.bincount(idx.ravel(), P2.ravel()) / np.maximum(np.bincount(idx.ravel()), 1)
amp_k = np.sqrt(Pk[idx])
def mock():
    ph = rng.standard_normal(F.shape) + 1j * rng.standard_normal(F.shape)
    f = np.fft.irfftn(amp_k * ph / np.sqrt(2), s=den.shape); f *= den[sel].std() / f[sel].std()
    return np.maximum(f, -1.0)
NM = 25; Rs_list = (0, 10, 20, 40, 60, 80)
mocks = [mock() for _ in range(NM)]
print(f"2M++ mass-centroid offset within R<{RMAX:.0f} Mpc/h vs Gaussian smoothing radius (null: {NM} Gaussian mocks with the data power spectrum)")
prev = None
for Rs in Rs_list:
    c = centroid(den, Rs); a = np.linalg.norm(c); l, b = lb(c)
    am = np.array([np.linalg.norm(centroid(m, Rs)) for m in mocks])
    ang = '' if prev is None else f"  turn from previous R_s: {np.degrees(np.arccos(np.clip(c @ prev / a / np.linalg.norm(prev), -1, 1))):.0f} deg"
    print(f"  R_s={Rs:3d}: offset {a:5.1f} Mpc/h toward ({l:.0f},{b:.0f});  mocks {am.mean():5.1f}+-{am.std():4.1f}, p(>=data)={np.mean(am >= a):.2f}{ang}")
    prev = c
