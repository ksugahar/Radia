"""radia.vim -- HDiv-type VIM demag operator.

The FEEC H(div) RT alternative/complement to the canonical multipole-moment MMM MSC kernel: a SYMMETRIC demag
operator N = B^T G B whose loop modes are field-null by construction (de Rham), giving mu_r-independent
convergence with no hand-crafted loop-star.

This is the PRODUCTION home (productionization milestone M1): the validated core was promoted here from
the retired Python prototype inventory.  Canonical docs: docs/hdiv_vim/README.md; roadmap:
docs/hdiv_vim/PRODUCTIONIZATION.md; retirement ledger: docs/hdiv_vim/vim_examples_retirement.ipynb.

Public API (NGSolve-aligned validated solve primitives):
  Solve(mesh, mu_r=/bh_table=, H_ext=..)
      -> the production one-call HDiv-VIM demag solve.  This is the preferred method-level entry,
         matching the NGSolve convention that solve/operator objects use short CamelCase names inside
         their owning namespace.
  DemagOperator(HDiv(mesh, order=1), intorder=, eps=, gram_backend=)
      -> an ngsolve.bem-style operator.  `.mat` is the H-matrix-backed NGSolve BaseMatrix
         N = B^T G B, which composes with NGSolve solvers / BlockMatrix exactly like ngsolve.bem's
         SingleLayerPotentialOperator.  `.DemagFactor(M_cf)` -> the demag factor (~1/3).
  ChargeGram(HDiv(mesh, order=1), ...)
      -> (B, G, M_mass), the charge map, charge-Gram H-matrix, and HDiv mass used by DemagOperator.
  MeshSoftIron(mesh, mu_r=/bh_table=) / VolSoftIron(path, mu_r=/bh_table=)
      -> method-layer constructors for mesh-backed Radia soft iron.  For ordinary user code prefer the
         user-intent API `rad.SoftIron(geometry, mu_r=...).solve(...)`.
  PlanarSolve(...) and PlanarDemagBody(...)
      -> the 2D planar tri/quad layer.

Lower-level diagnostics:
  build_demag(mesh, nsub=4)
      -> dict(M_mass, B, cell_verts, face_verts, poly, ...): the SPARSE pieces + the C++ charge-Gram
         geometry.  There is NO dense N: the C++ `_ChargeGramHMatrix` kernel IS the demag operator
         (N v = B^T (H.matvec(B v))) and folds IMA in via image_masks/image_signs (see Solve).
         The dense Python charge-Gram path (the O(N^2) Gram, the SVD loop basis, the analytic / image /
         Wilton dense Gram builders) was REMOVED.
  solve_nonlinear_newton(mesh, chi0, Msat, H0, bh_table=..., require_convergence=True)
      -> (M_avg, n_iter, D): damped Newton on the C++ scalable charge-Gram operator (fail-loud).
  Gram building block: tri_potential (the exact flat-triangle potential).

NOTE: importing this package imports `radia` (the C++ core).  The NGSolve-side HDiv-VIM solve itself does
not require the C++ core, but the production home is the radia package.
"""
import inspect as _inspect

from . import _core, _nonlinear  # noqa: F401
from ._core import (  # noqa: F401
    build_demag,
    tri_potential,
    build_near_correction,
    C_TRI,
)
from ._nonlinear import (  # noqa: F401
    solve_nonlinear_newton,
    solve_nonlinear_newton_scalable,
)
from ._vim import (  # noqa: F401  (ngsolve.bem-style operator + .mat)
    DemagOperator,
    build_charge_gram as _build_charge_gram,
    build_charge_gauss as _build_charge_gauss,
)
from ._solve import hdiv_demag_solve as _hdiv_demag_solve  # noqa: F401  (production demag solve)
from ._vim2d import (  # noqa: F401  (2D planar motor-cross-section layer; hdiv_demag_solve dispatches here)
    PlanarDemagBody,
    maxwell_torque_circle,
    solve_planar_demag as _solve_planar_demag,
)
from ._radsolve import (  # noqa: F401  (.vol/mesh -> both-backend iron)
    soft_iron_from_mesh as _soft_iron_from_mesh,
    soft_iron_from_vol as _soft_iron_from_vol,
)
from ._shapes import soft_iron_box, soft_iron_hex, magnet_box, magnet_hex  # noqa: F401  (mesh-less-SHAPE intent constructors: soft iron -> HDiv-VIM; PM -> analytic)
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
    make_t10_tet_map,                 # build a T10 (quadratic) curved-tet parametrization
    curved_tet_volume_field,          # CURVED-tet VOLUME-charge field (1->8 subdivision + flat closed form)
    assemble_demag_field,             # Step 2: fast charge-coefficient assembly (analytic C++ kernels)
    solve_demag_picard,               # Step 3: field-based nonlinear demag solve (under-relaxed Picard)
    solve_demag_newton,               # (4 step 2): matrix-free Newton-Krylov on the analytic field (robust stiff)
    saturating_tangent,               # (M_of_H, dM_of_H) for the saturating isotropic law (rank-1 tangent)
)


def Solve(*args, **kwargs):
    """NGSolve-style production HDiv-VIM one-call solve.
    """
    return _hdiv_demag_solve(*args, **kwargs)


def ChargeGram(*args, **kwargs):
    """NGSolve-style charge-Gram builder for an HDiv finite element space.

    Returns ``(B, G, M_mass)``.
    """
    return _build_charge_gram(*args, **kwargs)


def ChargeGramGauss(*args, **kwargs):
    """Retired Gauss point-operator charge Gram entry.

    Kept only so diagnostics fail through the same fail-loud path as the
    internal Gauss builder.
    """
    return _build_charge_gauss(*args, **kwargs)


def MeshSoftIron(*args, **kwargs):
    """Build a mesh-backed Radia soft iron container from an NGSolve mesh.

    User-facing workflows should usually use ``rad.SoftIron(mesh, ...)``.
    """
    return _soft_iron_from_mesh(*args, **kwargs)


def VolSoftIron(*args, **kwargs):
    """Build a mesh-backed Radia soft iron container from a Netgen ``.vol`` file.

    User-facing workflows should usually use ``rad.SoftIron(path, ...)``.
    """
    return _soft_iron_from_vol(*args, **kwargs)


def PlanarSolve(*args, **kwargs):
    """NGSolve-style alias for the 2D planar HDiv-VIM solve."""
    return _solve_planar_demag(*args, **kwargs)


def SolveNonlinearNewton(*args, **kwargs):
    """NGSolve-style alias for the nonlinear Newton diagnostic solve."""
    return solve_nonlinear_newton(*args, **kwargs)


def SolveNonlinearNewtonScalable(*args, **kwargs):
    """NGSolve-style alias for the scalable nonlinear Newton diagnostic solve."""
    return solve_nonlinear_newton_scalable(*args, **kwargs)


for _new, _old in [
    (Solve, _hdiv_demag_solve),
    (ChargeGram, _build_charge_gram),
    (ChargeGramGauss, _build_charge_gauss),
    (MeshSoftIron, _soft_iron_from_mesh),
    (VolSoftIron, _soft_iron_from_vol),
    (PlanarSolve, _solve_planar_demag),
    (SolveNonlinearNewton, solve_nonlinear_newton),
    (SolveNonlinearNewtonScalable, solve_nonlinear_newton_scalable),
]:
    _new.__signature__ = _inspect.signature(_old)

__all__ = [
    "build_demag", "tri_potential", "build_near_correction", "C_TRI",
    "Solve", "DemagOperator", "ChargeGram", "ChargeGramGauss",
    "MeshSoftIron", "VolSoftIron", "PlanarSolve",
    "SolveNonlinearNewton", "SolveNonlinearNewtonScalable",
    "solve_nonlinear_newton", "solve_nonlinear_newton_scalable",
    "PlanarDemagBody", "maxwell_torque_circle",
    "soft_iron_box", "soft_iron_hex",
    "magnet_box", "magnet_hex",
    "reconstruct_field", "reconstruct_field_polynomial",
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
    "make_t10_tet_map", "curved_tet_volume_field", "assemble_demag_field",
    "solve_demag_picard", "solve_demag_newton", "saturating_tangent",
    "_core", "_nonlinear", "_vim", "_field", "_solve", "_radsolve",
]
