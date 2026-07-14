"""radia.vim -- HDiv-type VIM demag operator.

The FEEC H(div) RT production demag path: a SYMMETRIC demag
operator N = B^T G B whose loop modes are field-null by construction (de Rham), giving mu_r-independent
convergence with no hand-crafted loop-star.

This is the PRODUCTION home for the validated HDiv-VIM core.  Canonical docs:
docs/hdiv_vim/README.md; roadmap:
docs/hdiv_vim/PRODUCTIONIZATION.md.

Public API (NGSolve-aligned validated solve primitives):
  Solve(mesh, mu_r=/bh_table=, H_ext=..)
      -> the production one-call HDiv-VIM demag solve.  This is the preferred method-level entry,
         matching the NGSolve convention that solve/operator objects use short CamelCase names inside
         their owning namespace.
  SolveHysteresis(mesh, h_steps, play=(K, eta, f_k_tables))
      -> quasi-static B-input hysteresis stepping on the SAME RT1 charge Gram: the chi-free
         H-matrix is built ONCE and reused by every step / nonlinear iteration; per-element
         committed play states advance step-by-step (save -> restore-before-eval -> commit).
         Duck-typed material protocol (state0/forward/commit/nu_B0) so lab-local B-input models
         plug in without touching radia.
  DemagOperator(HDiv(mesh, order=1), intorder=, eps=)
      -> an ngsolve.bem-style operator.  `.mat` is the H-matrix-backed NGSolve BaseMatrix
         N = B^T G B, which composes with NGSolve solvers / BlockMatrix exactly like ngsolve.bem's
         SingleLayerPotentialOperator.  `.DemagFactor(M_cf)` -> the demag factor (~1/3).
  ChargeGram(HDiv(mesh, order=1), ...)
      -> (B, G, M_mass), the charge map, charge-Gram H-matrix, and HDiv mass used by DemagOperator.
  MeshSoftIron(mesh, mu_r=/bh_table=) / VolSoftIron(path, mu_r=/bh_table=)
      -> method-layer constructors for mesh-backed Radia soft iron.  For ordinary user code prefer the
         user-intent API `rad.SoftIron(geometry, mu_r=...).solve(...)`.  When `rad.Solve(..., image=...)`
          is used on a MeshSoftIron, `rad.Fld(iron, ...)` evaluates the solved RT field and its reflected
          IMA contributions directly; `M_avg_reduced` is the
         reduced-domain diagnostic and `M_avg` is the physical full-domain average.  Unconstrained explicit
         full-solve `rad.Fld` parity is a separate 10-eps validation target, not a percent-level tolerance.
  PlanarSolve(...) and PlanarDemagBody(...)
      -> the 2D planar tri/quad layer.

NOTE: importing this package imports `radia` (the C++ core).  The NGSolve-side HDiv-VIM solve itself does
not require the C++ core, but the production home is the radia package.
"""
import inspect as _inspect

from . import _nonlinear  # noqa: F401
from ._vim import (  # noqa: F401  (ngsolve.bem-style operator + .mat)
    DemagOperator,
    build_charge_gram as _charge_gram_impl,
)
from ._solve import hdiv_demag_solve as _solve_impl  # noqa: F401  (production demag solve)
from ._vim2d import (  # noqa: F401  (2D planar motor-cross-section layer; vim.Solve dispatches here)
    PlanarDemagBody,
    maxwell_torque_circle,
    solve_planar_demag as _solve_planar_demag,
)
from ._radsolve import (  # noqa: F401  (.vol/mesh -> both-backend iron)
    soft_iron_from_mesh as _mesh_soft_iron_impl,
    soft_iron_from_vol as _vol_soft_iron_impl,
)
from ._hysteresis import (  # noqa: F401  (B-input hysteresis stepping: ONE Gram build, per-step W-CG)
    PlayHysteresisMaterial,
    SolveHysteresis,
)
from ._field_batch import (  # noqa: F401  (batch exterior field of the ORDER-1 solution, no constant-M collapse)
    field_from_solution as _field_from_solution_impl,
)
from ._shapes import soft_iron_box, soft_iron_hex, magnet_box, magnet_hex  # noqa: F401  (mesh-less-SHAPE intent constructors: soft iron -> HDiv-VIM; PM -> analytic)
def Solve(*args, **kwargs):
    """NGSolve-style production HDiv-VIM one-call solve.
    """
    return _solve_impl(*args, **kwargs)


def ChargeGram(*args, **kwargs):
    """NGSolve-style charge-Gram builder for an HDiv finite element space.

    Returns ``(B, G, M_mass)``.
    """
    return _charge_gram_impl(*args, **kwargs)


def MeshSoftIron(*args, **kwargs):
    """Build a mesh-backed Radia soft iron container from an NGSolve mesh.

    User-facing workflows should usually use ``rad.SoftIron(mesh, ...)``.
    """
    return _mesh_soft_iron_impl(*args, **kwargs)


def VolSoftIron(*args, **kwargs):
    """Build a mesh-backed Radia soft iron container from a Netgen ``.vol`` file.

    User-facing workflows should usually use ``rad.SoftIron(path, ...)``.
    """
    return _vol_soft_iron_impl(*args, **kwargs)


def PlanarSolve(*args, **kwargs):
    """NGSolve-style alias for the 2D planar HDiv-VIM solve."""
    return _solve_planar_demag(*args, **kwargs)


def FieldFromSolution(*args, **kwargs):
    """Batch demagnetizing H (A/m) at points from the ORDER-1 HDiv solution directly
    (no per-element constant-M collapse -- none of the near-surface ripple of
    ``rad.Fld`` on the write-back elements).  Pass ``vim.Solve``'s result dict.
    """
    return _field_from_solution_impl(*args, **kwargs)


for _new, _old in [
    (Solve, _solve_impl),
    (ChargeGram, _charge_gram_impl),
    (MeshSoftIron, _mesh_soft_iron_impl),
    (VolSoftIron, _vol_soft_iron_impl),
    (PlanarSolve, _solve_planar_demag),
    (FieldFromSolution, _field_from_solution_impl),
]:
    _new.__signature__ = _inspect.signature(_old)

__all__ = [
    "Solve", "DemagOperator", "ChargeGram",
    "MeshSoftIron", "VolSoftIron", "PlanarSolve",
    "PlanarDemagBody", "maxwell_torque_circle",
    "soft_iron_box", "soft_iron_hex",
    "magnet_box", "magnet_hex",
    "SolveHysteresis", "PlayHysteresisMaterial", "FieldFromSolution",
    "_nonlinear", "_vim", "_solve", "_radsolve", "_hysteresis",
]
