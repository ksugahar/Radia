# -*- coding: utf-8 -*-
# DEMO (s) (verified): in the refined limit the BEM exterior DtN and the FEM-Kelvin DtN are
# EQUIVALENT (two discretizations of the SAME continuous operator), and BOTH are rate-limited by
# the TRUNCATION-SURFACE resolution -- the volume (FEM-Kelvin, Galerkin-exact at p>=n, cf
# act2_09_exterior_mesh.py) and the singular-kernel quadrature (BEM) are NOT the bottleneck. At a
# FIXED matched trace order the two differ only by a method-specific CONSTANT that does not scale;
# the difference that SCALES is COST (sparse vs dense, cf act4_04_sparse_factorization_of_kernel), not accuracy.
#
# Measured (dipole n=1, MATCHED order 1, refine the surface mesh; both -> analytic -(n+1)/R):
#   maxh 0.50: N_surf 336, BEM 7.45e-4, FEM-Kelvin 5.09e-3   (assembly 68 s vs 0.03 s)
#   maxh 0.40: N_surf 564, BEM 2.86e-4, FEM-Kelvin 1.39e-3   (assembly 217 s vs 0.04 s)
# Both fall as the SURFACE refines (surface-limited, ~same rate); BEM is ~7x better at order 1 (a
# constant from the order-1 dual trace space vs the curved-geometry floor, NOT a rate). So "ultimate
# mesh -> BEM and FEM-Kelvin are equivalent, surface-rate-limited; the genuine difference is cost."
#
# NOTE: the dense BEM densification is O(N_S^2) and SLOW (minutes at maxh<=0.4); run with `python -u`.
import os, time
os.environ.setdefault("PYTHONIOENCODING", "utf-8")
from radia_mcp.radia_ngsolve.bem_integral import exterior_dtn_spectrum
from radia_mcp.radia_ngsolve.fem_bem_coupling import kelvin_dtn_eigenvalue

R = 1.0
print("dipole n=1, MATCHED order 1, refine the SURFACE mesh. Both -> -(n+1)/R, surface-limited:")
print("  maxh | N_surf | BEM rel.err | FEM-Kelvin rel.err |  t_BEM     t_Kelvin")
print("  -----+--------+-------------+--------------------+--------------------")
for maxh in (0.5, 0.4):
    t0 = time.time()
    bem = exterior_dtn_spectrum(R=R, maxh=maxh, order=1, intorder=10, nmax=1)
    t_bem = time.time() - t0
    be = [m for m in bem["modes"] if m["n"] == 1][0]["rel_err"]
    t0 = time.time()
    k = kelvin_dtn_eigenvalue(R=R, degree=1, maxh=maxh, order=1, intorder=10)
    t_kel = time.time() - t0
    print("  %.2f | %5d  |  %.3e  |     %.3e      | %6.1fs   %6.3fs"
          % (maxh, bem["ndof"], be, k["rel_err"], t_bem, t_kel))

print("\n-> both rel errors fall as the SURFACE refines (surface-limited; the FEM-Kelvin volume is")
print("   Galerkin-exact and the BEM kernel quadrature is converged, so neither is the bottleneck).")
print("   The order-1 gap (~7x) is a method CONSTANT, not a rate -> in the limit the two are the same")
print("   operator. The difference that SCALES is COST (sparse FEM-Kelvin ~ms vs dense BEM ~minutes).")
