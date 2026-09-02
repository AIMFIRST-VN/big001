"""Bounce-to-equality energy budget for Planck-relic ejecta.
Components at the bounce (per unit total relic rest energy):
  ejected relics: fraction f_ej, Lorentz factor g_ej  (momentum p ∝ 1/a)
  trapped relics: fraction 1-f_ej, Lorentz factor g_t; a fraction eps of their
                  total energy is converted to radiation by grinding in Kernel 0,0,0 (at the bounce, best case)
  extra radiation R0 from the crystallization itself (latent heat), not from relics
Question: is there a history with a radiation era (rho_rad >> rho_m down to T~MeV),
relics cold at matter-radiation equality (p/mc << 1), and surviving relic fraction ~0.84?"""
import numpy as np

def history(f_ej, g_ej, g_t=1.0, eps=1.0, R0=0.0):
    # per-particle momentum in units of m c at the bounce
    p_ej0 = np.sqrt(g_ej**2 - 1); p_t0 = np.sqrt(g_t**2 - 1)
    f_t = (1 - f_ej) * (1 - eps)            # trapped relics that survive
    rad0 = R0 + eps * (1 - f_ej) * g_t      # Kernel 0,0,0 conversion radiation (rest+kinetic of converted relics)
    a = np.logspace(0, 40, 4001)
    E_ej = f_ej * np.sqrt(1 + (p_ej0 / a)**2)      # energy density * a^3 (ejected)
    E_t = f_t * np.sqrt(1 + (p_t0 / a)**2)
    rad = rad0 / a
    tot = E_ej + E_t + rad
    # equality: radiation-like part (photons + relativistic kinetic) equals rest mass part
    rest = f_ej + f_t
    radlike = tot - rest
    i = np.argmax(radlike < rest) if np.any(radlike < rest) else -1
    a_eq = a[i]
    p_eq = p_ej0 / a_eq
    return a_eq, p_eq, rest, radlike[0]

print("case 1: radiation only from relics (kinetic + Kernel 0,0,0 rest-mass conversion)")
for g in (2, 10, 1e3, 1e10, 1e20):
    for eps in (0.0, 1.0):
        a_eq, p_eq, fs, A = history(0.88, g, g_t=1.0, eps=eps)
        print(f"  g_ej={g:>7.0e} eps={eps}: rad/rest at bounce={A:9.3g}  a_eq/a_b={a_eq:9.3g}  p/mc at eq={p_eq:6.3f}  surviving relic fraction={fs:.2f}")
print("\ncase 2: trapped relics hotter than ejecta (violates velocity sorting) g_t = 100 g_ej")
a_eq, p_eq, fs, A = history(0.88, 10, g_t=1000, eps=1.0)
print(f"  p/mc at eq={p_eq:.3f} surviving={fs:.2f}  (needs Kernel 0,0,0 energy >> ejecta energy)")
print("\ncase 3: radiation from crystallization latent heat, relics a trace")
for R0 in (1e3, 1e6, 1e28):
    a_eq, p_eq, fs, A = history(0.88, 10, eps=0.0, R0=R0)
    print(f"  R0={R0:8.0e}: a_eq/a_b={a_eq:9.3g}  p/mc at eq={p_eq:8.2e}  surviving={fs:.2f}  relic packing fraction at bounce ~ {0.64/(1+R0):.1e}")

print("\ncase 4 (H8 Precious Rings): core-only accounting. Kernel 0,0,0 converts fraction (1-f) of core relics to light at T_conv;")
print("survivors f are the local cold dark matter. Consistency with T_eq = 0.8 eV requires f = T_eq/T_conv.")
T_eq=0.8  # eV
for name,T_conv in (("Planck",1.2e28),("GUT 1e16 GeV",1e25),("1e10 GeV",1e19),("TeV",1e12),("MeV",1e6)):
    f=T_eq/T_conv
    for g in (2,1e3):
        R=(1-f)/f; p_eq=np.sqrt(g**2-1)/R
        print(f"  T_conv={name:>13}: survivor fraction f={f:8.1e}  R={R:8.1e}  p/mc at eq (gamma={g:g})={p_eq:8.1e}  cold={'yes' if p_eq<1e-3 else 'no'}")
print("  The 88% ejection fraction concerns the outer rings (shell) and does not enter the core budget.")
