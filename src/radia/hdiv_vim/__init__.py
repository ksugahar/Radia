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
from ._field import (  # noqa: F401  (field-at-points from solved M; NOT M_mass^-1 N m)
    reconstruct_field,
    reconstruct_field_polynomial,  # Step 1: EXTERNAL polynomial-charge field (tet + hex)
    reconstruct_field_internal,    # Step 2: INTERNAL/near field (self-volume spherical + analytic surface)
    flat_triangle_charge_field,    # Step-2 building block: exact uniform-triangle field (surface near-field)
    tet_self_volume_field,         # Step-2 building block: tet self volume-charge field (spherical ray-trace)
    triangle_potential_const,      # degree-1 building block: INT_T 1/R dS' (Wilton)
    triangle_potential_moment,     # degree-1 building block: INT_T r'/R dS' (first moment)
    tet_newtonian_potential,       # degree-1 building block: INT_V 1/R dV' (PhiTet, pure-Python)
    tet_volume_field_linear,       # EXACT closed-form LINEAR volume-charge field (order-2 -div M term)
    linear_triangle_charge_field,  # EXACT closed-form LINEAR surface-charge field (order-2 M.n term)
    triangle_potential_moment2,    # degree-2 building block: INT_T r'(x)r'/R dS' (second moment)
    tet_newtonian_moment,          # degree-2 building block: INT_V r'/R dV' (volume first moment)
    tet_volume_field_quadratic,    # EXACT closed-form QUADRATIC volume-charge field
    quadratic_triangle_charge_field,  # EXACT closed-form QUADRATIC surface-charge field
    triangle_inplane_moments,         # general surface moment dicts A_k (1/R), B_k (1/R^3), any degree
    polynomial_triangle_charge_field,  # ARBITRARY-degree surface-charge field (general assembler)
    tet_volume_field_polynomial,      # ARBITRARY-degree volume-charge field (general assembler)
    tet_boundary_triangles,           # flat-faced polytope helper: tet -> 4 (tri, outward n)
    hex_boundary_triangles,           # flat-faced polytope helper: hex -> 12 (tri, outward n)
    polytope_newtonian_potential,     # INT_V 1/R over any flat-faced polytope
    polytope_newtonian_moment,        # INT_V r'/R over any flat-faced polytope
    polytope_volume_field_quadratic,  # quadratic volume-charge field over any flat-faced polytope
    hex_volume_field_linear,          # linear volume-charge field over an (affine) hex
    hex_volume_field_quadratic,       # quadratic volume-charge field over an (affine) hex
    make_t6_surface_map,              # build a T6 (quadratic) curved-triangle parametrization
    curved_triangle_charge_field,     # CURVED-face surface-charge field (singularity subtraction + Duffy)
)

__all__ = [
    "build_demag", "demag_factor", "tri_potential", "phi_tet", "wilton_surface_block",
    "analytic_charge_gram", "build_near_correction", "C_TRI",
    "solve_nonlinear_newton", "solve_nonlinear_newton_scalable", "solve_nonlinear",
    "DemagOperator", "build_charge_gram", "reconstruct_field", "reconstruct_field_polynomial",
    "reconstruct_field_internal", "flat_triangle_charge_field", "tet_self_volume_field",
    "triangle_potential_const", "triangle_potential_moment", "tet_newtonian_potential",
    "tet_volume_field_linear", "linear_triangle_charge_field",
    "triangle_potential_moment2", "tet_newtonian_moment", "tet_volume_field_quadratic",
    "quadratic_triangle_charge_field", "triangle_inplane_moments",
    "polynomial_triangle_charge_field", "tet_volume_field_polynomial",
    "tet_boundary_triangles", "hex_boundary_triangles", "polytope_newtonian_potential",
    "polytope_newtonian_moment", "polytope_volume_field_quadratic",
    "hex_volume_field_linear", "hex_volume_field_quadratic",
    "make_t6_surface_map", "curved_triangle_charge_field",
    "_core", "_nonlinear", "_vim", "_field",
]
