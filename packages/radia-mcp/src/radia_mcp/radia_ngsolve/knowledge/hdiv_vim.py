"""HDiv-type VIM knowledge for radia_mcp.radia_ngsolve.

This module describes the current Radia soft-iron direction: mesh-backed
H(div) volume-integral demagnetization using NGSolve meshes and Radia's C++
charge-Gram H-matrix.  It intentionally avoids retired solver history so MCP
answers steer agents toward the live implementation.
"""

_OVERVIEW = r"""
# HDiv-type VIM demag operator

Radia's soft-iron demag path is the HDiv-type volume integral method (VIM).
The unknown magnetization is represented in an H(div) space on an NGSolve mesh;
the demag operator is

    N = B^T G B

where B maps magnetization to volume/surface magnetic charge, and G is the
Laplace single-layer/Coulomb Gram.  The important engineering properties are:

- loop modes are field-null by construction because loops live in ker(B);
- the operator is symmetric and compatible with MINRES/CG-style solvers;
- material state, fields, and source terms live in NGSolve Mesh,
  GridFunction, CoefficientFunction, and BilinearForm vocabulary;
- this makes reduced FEM coupling cleaner than sampling a separate object
  field back into FEM.

Use this route for mesh-backed TET, HEX, and WEDGE soft iron, nonlinear BH
curves, and planar motor cross sections.  Mesh-less soft-iron solves are not a
supported Radia production path; create an NGSolve mesh and call the HDiv API.
"""

_IMPLEMENTATION = r"""
# Implementation

Primary Python entry points:

- `radia.vim.MeshSoftIron(mesh, mu_r=... | bh_table=...)`
- `radia.Solve(model, prec, maxiter, method, demag_backend="hdiv")`
- `radia.vim.Solve(mesh, mu_r=... | bh_table=..., H_ext=..., image=...)`
- `radia.vim.build_demag(mesh, order=..., image=...)` for diagnostics

Core pieces:

- `src/radia/vim/` handles mesh ingestion, dispatch, material setup, image
  symmetry contracts, and field reconstruction.
- `src/core/rad_hdiv_vim.*` contains structured and unstructured HDiv assembly
  helpers.
- `src/core/rad_hacapk_hdiv.*` contains `_ChargeGramHMatrix`, the C++ H-matrix
  backend for the Coulomb Gram.
- `src/radia/planar_geometry.py`, `planar_materials.py`, `planar_charges.py`,
  `planar_hysteresis.py`, and `planar_aniso.py` provide 2D planar shared
  geometry/material helpers.

TaskManager is assumed.  NGSolve assembly should run under
`with ngsolve.TaskManager():`, and the C++ kernels use parallel loops for
charge gather, dot products, preconditioner/vector updates, and sparse scatters.
"""

_SCALING = r"""
# Scaling

The costly object is the charge Gram G.  Radia builds it through HACApK as a
charge H-matrix, then applies the material operator as B^T G B without
materializing a dense N for production runs.

Record these quantities in validation and benchmark artifacts:

- number of magnetic elements;
- H(div) unknown count and charge count;
- H-matrix build time and compression;
- solve iterations and residual;
- peak memory when available;
- machine label (`LAB` smoke vs `mdx` validation).

Small problems are allowed to be simply "interactive".  The scaling question
matters at engineering size, where charge count and matrix build dominate.
Timing claims should be taken on mdx when it is idle.
"""

_VERIFICATION = r"""
# Verification

Fast tests should cover API contracts and small deterministic checks:

- backend selection rejects unsupported mesh-less soft iron;
- pure TET/HEX/WEDGE mesh-backed soft iron dispatches to HDiv;
- `rad.Fld` after `rad.Solve(..., image=...)` matches an explicitly mirrored
  full model for truly symmetric meshes to near roundoff;
- 2D planar helpers preserve material labels and PM source regions;
- public solver names and config keys match the current API.

Validation-class tests live under `validation_test/feec/` and should cover:

- sphere/cube demag factors and convergence trends;
- nonlinear BH curves with convergence metadata;
- IMA/image symmetry for TET/HEX/WEDGE;
- curved and high-order geometry where analytic demag truth exists;
- reduced-FEM handoff through NGSolve fields/CoefficientFunctions.
"""

_NONLINEAR = r"""
# Nonlinear Material Solve

For BH curves, the HDiv route updates a material state on the mesh and solves
the demag equation with a robust nonlinear iteration.  Engineering defaults:

- use tolerances that match observable mesh/discretization error;
- record convergence status, iteration count, and max update;
- fail loudly on non-convergence;
- keep the solve under TaskManager.

Deep saturation can require safeguarded nonlinear steps.  Do not judge the
method from a single scalar residual; inspect the material-state update and the
field observable used by the application.
"""

_CURVED = r"""
# Curved And High-Order Geometry

Curved geometry is one of the main reasons Radia keeps soft iron in the HDiv /
NGSolve lane.  `mesh.Curve(p)`, Piola mappings, curved boundary normals, and
high-order integration are shared with the reduced FEM side.

Good validation cases:

- sphere and spheroid demag factors against analytic values;
- curved sphere external field against the exact dipole;
- curved high-order element convergence compared with flat low-order faceting;
- GMSH/Netgen export checks when the mesh originates from Cubit.
"""

_SYMMETRY = r"""
# Image Symmetry

Image symmetry is part of the HDiv field contract.  A reduced model and an
explicit full model should agree to near roundoff when the mesh is geometrically
and topologically symmetric.  Percent-level agreement is a warning sign for
asymmetric mesh cuts, incorrect image signs, wrong materialization of images, or
quadrature/charge-basis mismatch.

For reduced models, record:

- image string and reflected axes;
- real charge count and image count;
- whether `rad.Fld` was evaluated through the reduced model or a materialized
  full model;
- max/mean field difference at probes.
"""

_CROSS_METHOD = r"""
# Cross-Method Checks

Prefer analytic truth first: ellipsoid demag factors, cuboid permanent-magnet
fields, dipole limits, and closed-form thin/axisymmetric cases.  When analytic
truth is unavailable, use independent formulations:

- HDiv VIM on the NGSolve mesh;
- volume FEM A/phi or reduced-potential FEM where appropriate;
- boundary-element single-layer checks for surface-charge problems;
- direct full-model image materialization for symmetry tests.

Do not put local third-party provenance into public artifacts.  Public docs
should state the analytic convention and the reproduced number, not internal
comparison file names.
"""

_REFERENCE_AUDIT = r"""
# Reference Audit Ladder

When a disagreement appears:

1. Inspect mesh materials, boundaries, and finite-element spaces before solving.
2. Verify image signs and whether the mesh cut is exactly symmetric.
3. Compare charge maps and field evaluation before nonlinear iteration.
4. Check the same observable through two evaluators (`M_avg`, `rad.Fld`, probe
   grid, energy) before changing solver tolerances.
5. Move heavy sweeps to mdx and label the result as validation, not LAB smoke.
"""

_STATUS = r"""
# Status

Current direction:

- Radia soft iron: HDiv-VIM.
- Planar 2D support: HDiv/planar shared geometry and material helpers.
- Public docs: result-bearing HDiv notebooks plus synchronized JSON.
- MCP: teach the live HDiv API and reduced-FEM coupling path.

Open work:

- extend 2D and 3D validation coverage around `rad.Fld`;
- harden image-symmetry roundoff contracts;
- continue charge-Gram H-matrix performance checks on mdx;
- keep Cubit/GMSH mesh-export artifacts aligned with the HDiv API.
"""

_SECTIONS = {
    "overview": _OVERVIEW,
    "implementation": _IMPLEMENTATION,
    "scaling": _SCALING,
    "verification": _VERIFICATION,
    "nonlinear": _NONLINEAR,
    "curved": _CURVED,
    "symmetry": _SYMMETRY,
    "cross_method": _CROSS_METHOD,
    "reference_audit": _REFERENCE_AUDIT,
    "status": _STATUS,
}


def get_hdiv_vim_documentation(topic: str = "overview") -> str:
    """Return HDiv-VIM knowledge for a topic; use 'all' for every section."""
    t = (topic or "overview").strip().lower()
    if t == "all":
        return "\n\n".join(_SECTIONS[key] for key in _SECTIONS)
    if t in _SECTIONS:
        return _SECTIONS[t]
    return (
        f"Unknown topic '{topic}'. Options: "
        + ", ".join(_SECTIONS.keys())
        + ", all.\n\n"
        + _OVERVIEW
    )
