"""Decay law for the kernel's turbulence. Free decaying turbulence loses energy as E ~ t^-10/7 (Kolmogorov/Batchelor-
Townsend; Saffman spectra give t^-6/5), measured in eddy turnovers. In the expanding kernel the number of turnovers
between the bounce and recombination is N = ln(a_rec/a_bounce)/nu_t with nu_t the turnover-to-expansion ratio (the
paper fits nu_t = 0.105 from n_s). Vortical velocity does NOT redshift away in the radiation era (v ~ a^0 for w=1/3), so
only the cascade's own decay removes it. Requirement at last scattering: residual turbulent v/c < 1e-5 (Doppler).
Outputs: residual v/c for v0 = c and nu_t = 0.105; the nu_t that would suffice; the v0 that would suffice at nu_t = 0.105."""
import numpy as np
lnA = np.log(1e29)                    # bounce a ~ 1e-32 to recombination a ~ 1e-3
for p in (10/7, 6/5):
    for nu in (0.105, 0.01, 1e-3):
        N = lnA/nu; v = N**(-p/2)
        print(f"E~t^-{p:.2f}: nu_t={nu:<6} turnovers N={N:9.0f}  residual v/c (v0=c) = {v:.1e}")
    Nreq = (1e-5)**(-2/p); print(f"   need N >= {Nreq:.1e} turnovers -> nu_t <= {lnA/Nreq:.1e}; or at nu_t=0.105, v0/c <= {1e-5/(lnA/0.105)**(-p/2):.1e}")
