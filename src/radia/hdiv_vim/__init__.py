"""radia.hdiv_vim -- HDiv-type VIM demag operator.

The FEEC H(div) RT alternative to the collocation MMM/MSC kernel, and the candidate replacement for the
yano-type distortion elements: a SYMMETRIC demag operator N = B^T G B whose loop modes are field-null by
construction (de Rham), giving mu_r-independent convergence with no hand-crafted loop-star.

This is the PRODUCTION home (productionization milestone M1): the validated core was promoted here from
examples/feec_vim.  Canonical docs: docs/hdiv_vim/README.md; roadmap: docs/hdiv_vim/PRODUCTIONIZATION.md.

Public API (validated solve primitives):
  build_demag(mesh, nsub=4, wilton_surface=False, analytic_gram=False)
      -> dict(N, M_mass, B, loops, ...): the linear demag operator N = B^T G B.
         wilton_surface=True : exact analytic SURFACE Gram (uniform-M linear demag, demag factor 1/3).
         analytic_gram=True  : full analytic volume Gram (REQUIRED for non-uniform / nonlinear, div M!=0).
  demag_factor(d) -> the demag factor (Rayleigh quotient) from a build_demag result.
  solve_nonlinear_newton(mesh, chi0, Msat, H0, analytic_gram=..., bh_table=..., require_convergence=True)
      -> (M_avg, n_iter, D): damped Newton (fail-loud on non-convergence).
  Gram building blocks: tri_potential, phi_tet, wilton_surface_block, analytic_charge_gram.

  ngsolve.bem-STYLE API (._vim): DemagOperator(fes, intorder=, eps=) -- construct from an HDiv FESpace
      (the order comes from the fes); `.mat` is the H-matrix-backed NGSolve BaseMatrix N = B^T G B, which
      composes with NGSolve's solvers / BlockMatrix exactly like ngsolve.bem's SingleLayerPotentialOperator.
      order=0 (RT0) and order=p go through ONE call.  `.DemagFactor(M_cf)` -> the demag factor (~1/3).

NOTE: importing this package imports `radia` (the C++ core).  The NGSolve-side HDiv-VIM solve itself does
not require the C++ core, but the production home is the radia package.
"""
from . import _core, _nonlinear  # noqa: F401
from ._core import (  # noqa: F401
    build_demag,
    demag_factor,
    tri_potential,
    phi_tet,
    wilton_surface_block,
    analytic_charge_gram,
    build_near_correction,
    C_TRI,
)
from ._nonlinear import (  # noqa: F401
    solve_nonlinear_newton,
    solve_nonlinear_newton_scalable,
    solve_nonlinear,
)
from ._vim import DemagOperator, build_charge_gram  # noqa: F401  (ngsolve.bem-style operator + .mat)
from ._field import reconstruct_field, reconstruct_field_polynomial  # noqa: F401  (field-at-points from solved M; NOT M_mass^-1 N m)

__all__ = [
    "build_demag", "demag_factor", "tri_potential", "phi_tet", "wilton_surface_block",
    "analytic_charge_gram", "build_near_correction", "C_TRI",
    "solve_nonlinear_newton", "solve_nonlinear_newton_scalable", "solve_nonlinear",
    "DemagOperator", "build_charge_gram", "reconstruct_field", "reconstruct_field_polynomial",
    "_core", "_nonlinear", "_vim", "_field",
]
