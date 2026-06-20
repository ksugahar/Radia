# -*- coding: utf-8 -*-
"""Sufficient surface mesh vs Kelvin-region surface resolution:
the resolvable multipole band (eigenvalue rel_err < 0.5%) as a function of the
truncation-sphere surface DOF, vs the raw angular-Nyquist ceiling sqrt(ndof)."""
import math
from radia_mcp.radia_ngsolve.bem_integral import dtn_spectrum_vs_mesh

study = dtn_spectrum_vs_mesh(R=1.0, maxh_list=(0.5, 0.4, 0.32, 0.28),
                             order=1, intorder=10, nmax=5, band_tol=0.005)
nd = study["ndof_list"]
band = study["accurate_band"]
et = study["error_table"]
# normalise error_table keys to int
et = {int(k): v for k, v in et.items()}
degrees = sorted(et)

print(f"{'ndof':>6} {'n_acc(<0.5%)':>12} {'n_Nyquist≈√ndof-1':>18}   per-degree rel_err")
for j, (n_, b_) in enumerate(zip(nd, band)):
    nyq = math.sqrt(n_) - 1
    errs = "  ".join(f"n{n}:{et[n][j]:.1e}" for n in degrees)
    print(f"{n_:>6d} {('n<='+str(b_)):>12} {nyq:>18.1f}   {errs}")

print("\ncoarse_low_mode_max_err (n<=2, coarsest mesh):",
      f"{study['coarse_low_mode_max_err']:.2e}")
print("accurate_band per mesh:", band)
print("DONE")
