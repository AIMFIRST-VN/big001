"""Sign of dR/dt through the LQC-form dust bounce (rho = rho_c/(1+6 pi G rho_c t^2), Planck units)
and across a cascade wall (false vacuum w=-1 -> dust), for gravitational baryogenesis (Davoudiasl et al. 2004,
coupling dR/dt * J_B). Result: dR/dt > 0 during contraction, 0 at the bounce, < 0 for all t > 0; across the wall
R drops from 32 pi G rho to 8 pi G rho, so dR/dt < 0 there as well. Without a contracting branch every curvature
change in the model has the same sign -> a universal matter/antimatter sign (loaded coin)."""
import numpy as np
rc = 0.41; k = 6 * np.pi * rc
t = np.linspace(-6, 6, 240001); dt = t[1] - t[0]
a = (1 + k * t**2) ** (1 / 3)
H = np.gradient(np.log(a), dt); Hd = np.gradient(H, dt)
R = 6 * (Hd + 2 * H**2); Rd = np.gradient(R, dt)
for tt in (-1, -0.2, 0, 0.2, 1, 4):
    i = np.argmin(abs(t - tt)); print(f"t={tt:5.1f} t_P  R={R[i]:8.3f}  dR/dt={Rd[i]:8.3f}")
print("R(vacuum) = 32 pi rho_c = %.2f ; R(dust at rho_c) = 8 pi rho_c = %.2f -> dR/dt < 0 across the wall" % (32 * np.pi * rc, 8 * np.pi * rc))
