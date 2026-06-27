"""radia.moment_galerkin -- the SYMMETRIC moment-Galerkin MMMM demag solver.

MMMM (Multipole Magnetic Moment Method) made SYMMETRIC by a Bubnov-Galerkin (moment) test functional instead
of point collocation.  Collocation MMMM (radia's `RadHACApKMomentSystem`, the moment-yano route) is
non-symmetric (||A - A^T||/||A|| ~ 0.07-1.6) because the test is the field at the centroid != the trial moment
distribution -> loop modes that ill-condition BiCGSTAB.  The moment-Galerkin coupling is the moment-moment
mutual-energy tensor N = B^T G B, SYMMETRIC by Green-kernel reciprocity (validated 1e-16 on the loop-heavy
C-yoke) -> the loop modes are field-null by construction -> mu_r-independent / loop-free convergence.

Built on the EXISTING C++ charge-Gram H-matrix (`_ChargeGramHMatrix`) + Krylov solve (mass-Riesz CG for the
dipole, diagonal-Jacobi auto_prec for the quad) -- the same kernels HDiv-VIM (radia.vim) ships.  No new C++
kernel: the heavy compute (Gram + solve) is C++, the moment-basis assembly (B, M_mass) is this thin Python
layer.

Moment orders (`quad=`):
  quad=False (default) -- 3 DOF/hex: constant magnetization per hex (the dipole / standard demag).
  quad=True            -- 5 DOF/hex: + 2 quad residual-eigenmode amplitudes (higher per-element order; the
      quad modes improve skew / gradient loads -- the iron-yoke case -- and are ~null for axial/symmetric
      loads, never worse).

Public API:
  moment_galerkin_demag_solve(hexes, mu_r=/chi=, H_ext=.., quad=False)
      -> dict(M (n,3) dipole magnetization, m raw amplitudes, quad_amps, iters, demag_factor, n_hex, sys).
  reconstruct_field(sys, m, probes) -> external B (Tesla) from the solved amplitudes (field-from-sigma,
      anchored to rad.Fld ~1e-13).
  assemble_moment_system(hexes, quad=False, ...) -> dict(G, B, M_mass, vols, ndof_per, all_tris, ...).
  solve_assembled(sys, H_ext, chi) -> (M_dipole (n,3), iters);  solve_raw -> (m, iters).
  demag_factor(sys, kdir=2) -> float (operator demag factor, cube -> 1/3).

Loop-EXCITING sources (azimuthal / transformer drive) route through HDiv-VIM (radia.vim) per the loop-free
architecture decision.
"""
from ._assemble import assemble_moment_system, HEX_FACES  # noqa: F401
from ._solve import (  # noqa: F401
    moment_galerkin_demag_solve,
    solve_assembled,
    solve_raw,
    reconstruct_field,
    demag_factor,
)

__all__ = [
    "moment_galerkin_demag_solve",
    "assemble_moment_system",
    "solve_assembled",
    "solve_raw",
    "reconstruct_field",
    "demag_factor",
    "HEX_FACES",
]
