"""radia.moment_galerkin -- the SYMMETRIC moment-Galerkin MMMM demag solver.

MMMM (Multipole Magnetic Moment Method) made SYMMETRIC by a Bubnov-Galerkin (moment) test functional instead
of point collocation.  Collocation MMMM (radia's `RadHACApKMomentSystem`, the moment-yano route) is
non-symmetric (||A - A^T||/||A|| ~ 0.07-1.6) because the test is the field at the centroid != the trial moment
distribution -> loop modes that ill-condition BiCGSTAB.  The moment-Galerkin coupling is the moment-moment
mutual-energy tensor N = B^T G B, SYMMETRIC by Green-kernel reciprocity (validated 1e-16 on the loop-heavy
C-yoke) -> the loop modes are field-null by construction -> mu_r-independent / loop-free convergence.

This is the de-risk-validated dipole-level solver (3 DOF/hex, constant M per element), built on the EXISTING
C++ charge-Gram H-matrix (`_ChargeGramHMatrix`) + mass-Riesz CG (`solve_linear_material_mass_riesz`) -- the
same kernels HDiv-VIM (radia.vim) ships.  No new C++ kernel: the heavy compute (Gram + solve) is C++, the
moment-basis assembly (B, M_mass) is this thin Python layer.

Public API:
  moment_galerkin_demag_solve(hexes, mu_r=/chi=, H_ext=..) -> dict(M, iters, demag_factor)
      Linear isotropic soft-iron demag on a hexahedral body (the validated production entry).
  assemble_moment_system(hexes, ...) -> dict(G, B, M_mass, vols, ...)   (the sparse pieces + C++ Gram)
  solve_assembled(sys, H_ext, chi)   -> (M, iters)                       (solve a pre-assembled system)
  demag_factor(sys, kdir=2)          -> float                            (operator demag factor, cube -> 1/3)

Scope (this dipole-level increment): a uniform M per hex (the lowest moment order = the standard demag).
The higher moment modes (the 2 quad residual-eigenmodes per hex, the '6-DOF' set) are a separate increment.
Loop-EXCITING sources (azimuthal / transformer drive) route through HDiv-VIM (radia.vim) per the loop-free
architecture decision.
"""
from ._assemble import assemble_moment_system, HEX_FACES  # noqa: F401
from ._solve import (  # noqa: F401
    moment_galerkin_demag_solve,
    solve_assembled,
    demag_factor,
)

__all__ = [
    "moment_galerkin_demag_solve",
    "assemble_moment_system",
    "solve_assembled",
    "demag_factor",
    "HEX_FACES",
]
