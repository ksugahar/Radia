"""Application-neutral electromagnetic-force knowledge router.

The detailed documents predate the standalone Force server and remain in
their original modules for import compatibility.  This module is now the
canonical public routing layer.  Application servers must import these
functions instead of reaching into another application's knowledge module.
"""

from __future__ import annotations

from radia_mcp.differential_forms.em_force_extras_knowledge import (
    get_em_force_extras,
)
from radia_mcp.differential_forms.em_force_ngsolve_recipe_knowledge import (
    get_em_force_ngsolve_recipe,
)
from radia_mcp.differential_forms.forces_knowledge import (
    get_forces_documentation,
)
from radia_mcp.radia_ngsolve.knowledge.force_validation import (
    get_force_validation_documentation,
)

TOPICS: dict[str, str] = {
    "overview": "Shared Force-domain boundary, workflow, and Motor/MagLev ownership",
    "methods": "Unified electromagnetic-force theory and seven-method catalog",
    "recipe": "Practical method-selection and high-order NGSolve recipe",
    "extras": "Permanent-magnet, energy/coenergy, shape-derivative, Lorentz, and Meissner methods",
    "validation": "Method map, eggshell guidance, cross-checks, and reference-result contract",
    "motor": "How Motor consumes common Force methods while retaining motor-specific torque gates",
    "maglev": "How MagLev consumes common Force methods while retaining levitation dynamics",
    "source_quadrature": "Quadrature of an analytic source field on elements near the source: distance classes, measured errors, near-element subdivision",
    "all": "Concatenation of all common Force knowledge",
}


OVERVIEW = r"""
# radia-mcp.force — common electromagnetic-force layer

`radia_mcp.force` is the shared front door for electromagnetic force and
torque post-processing.  Motor and MagLev consume this layer; neither owns the
general Maxwell-stress, Lorentz-force, virtual-work, or validation recipe.

## Boundary

| Layer | Owner |
|---|---|
| Field solution, mesh, material law, quadrature samples | Radia, NGSolve, BEM/PEEC, or an application solver |
| Solver-independent sample integration | `radia.force` |
| Method choice, common formulas, pitfalls, cross-validation | `radia_mcp.force` |
| Motor geometry, periodic torque, cogging, rotating-frame gates | `radia_mcp.motor` |
| Lift dynamics, settling, stability, TEAM 28 motion gates | `radia_mcp.maglev` |

## Minimal workflow

1. Solve the electromagnetic field with the application-appropriate solver.
2. Select the force method with `force_recipe("method_choice")` and normalize
   solver-owned resultants with `force_result`.
3. For a unit-permeability conductor, integrate `J x B` with
   `force_lorentz`.  For a body enclosed by a closed air surface, integrate
   Maxwell traction with `force_maxwell_surface`. Supplying sample positions
   and a pivot returns force and torque together.
4. For harmonic fields, use `force_time_average_lorentz` or
   `force_time_average_maxwell_surface` and declare peak or RMS phasors.
5. For magnetic material, prefer weighted stress or constant-current
   coenergy virtual work. Use `force_virtual_work`, `force_coenergy_torque`,
   and `force_air_gap_torque` for solver-independent tables and estimates.
6. Check method suitability, independent-method agreement, Newton's third
   law, lift/weight equilibrium, path/surface independence, and mesh
   refinement with the four `force_*_gate` tools and
   `force_validation_guide`.

Static tools accept SI data.  Vectors use Cartesian components; current
density is A/m^2, flux density is T, volume weights are m^3, surface weights
are m^2, returned force is N, and torque is N m. Phasor tools require an
explicit peak/RMS convention and implement the conjugated time-average
identity; static tools reject complex samples.
"""


MOTOR_GUIDANCE = r"""
# Motor handoff to the common Force layer

Use `force_recipe`, `force_extras`, and `force_validation_guide` for general
force/torque extraction.  `motor_em_force_recipe` and
`motor_em_force_extras` are compatibility aliases that forward here.

Keep motor-only concerns in `radia_mcp.motor`: rotor/stator selection,
periodicity, cogging torque, rotating-frame covariance, air-gap sampling, and
the motor result-artifact contract. The compatibility tool
`motor_force_torque_method_agreement_gate` delegates to the common Force gate.
`HDivReducedMotor.sweep` emits normalized Maxwell-surface,
magnetization-volume, and virtual-work torque records for every angle, so those
independent routes can be checked directly by the common agreement gate.
"""


MAGLEV_GUIDANCE = r"""
# MagLev handoff to the common Force layer

Use `force_recipe`, `force_extras`, and `force_validation_guide` for Maxwell
stress, Lorentz force, and virtual-work extraction.  The MagLev
`force_computation` topic forwards to this common layer and then adds
levitation-specific notes.

Keep MagLev-only concerns in `radia_mcp.maglev`: lift/weight equilibrium,
stability, periodic eddy-current settling, motion coupling, and TEAM 28
cycle-averaged dynamics. Its force-method agreement and lift/weight tools are
application aliases of the common Force gates.
For quantitative 3-D plate forces, build the divergence-free HCurl-VIM current
model and use `compute_lorentz_force_result_via_hcurl_vim`; it emits the
conductor/source reaction pair in the common peak-phasor schema.  The older
`compute_lorentz_force_result_via_foster` is only a reduced solve of a scalar
local-reaction approximation.  Independent 3-D evidence rejects that scalar
model quantitatively; `compute_lorentz_force_via_foster_verified` can control
its high-frequency truncation but cannot repair its physical-model error.
`PositionForceCurve.force_result_at` converts a sampled position-force
interpolation into the same common result schema.
"""


SOURCE_QUADRATURE = r"""
# Source-field quadrature near the source

Force and load integrals often evaluate an analytic source field at the
quadrature points of meshed elements: the reduced-potential load
`int mu H_s . grad v dx`, the Lorentz force `F = int J x B_ext dV` on a
conductor, and drive terms such as `int A_ext . psi`.  On elements close to
the source this integrand varies quickly or is not smooth.  It is the
reciprocal of evaluating a layer potential near its source surface, where the
source element is split at the projection of the target point.

## Classify elements by distance

Use the signed distance `d` from the element to the source conductor
(negative inside) over the longest element edge `h`:

| class | integrand | convergence in the rule order |
|---|---|---|
| straddles the conductor surface, or inside it | `grad H_s` jumps on the surface | algebraic; raising the order is inefficient |
| `0 < d/h < 0.5` | smooth, near-singular | fast |
| `d/h >= 1` | smooth | small error in the recorded case; no general distance-only error bound |

A thin or filament source (a `1/r` line singularity) is worse than the solid
conductor measured below and should be checked on its own.

## Measured case: C-type magnet, coarse Kelvin mesh

Element-wise source integrals against an order-20 rule, production
tetrahedral order 8, solid rectangular coil
(`validation_test/c_type_three_engine/source_load_quadrature_20260927/`):

- Air (the coil passes through unmeshed air elements): 534 elements straddle
  the conductor and 58 lie inside it; together they carry about 97% of the
  squared error.  Median relative element error 6.5e-4, maximum 4.5e-2.
  Elements with `d/h >= 1` are at or below 1e-11.  The global relative error
  falls only from 5.3e-3 to 1.3e-3 between orders 6 and 16.
- Iron (no element straddles the coil): global relative error 1.9e-6 at
  order 8, all of it in the 104 elements with `d/h < 0.5`, falling to 4.4e-8
  at order 16.

What matters is where the observable sits relative to those elements.  The
gap-core B of that magnet moved by less than 1e-6 under every load-only
change (air rule order 20, face-flux load, whole-solver bonus, exact Kelvin
exterior), because the erroneous elements sit at the coil and not at the
gap.  A force integral whose largest weights are on the elements nearest the
source, such as `J x B_ext` on a plate close to a coil, is the opposite case.

## Remedy: subdivide only the near elements

Keep the production rule on far elements and replace it on the near set
(straddling, inside, `d/h < 1`) with a composite rule: split the reference
tetrahedron into `8^L` Bey children and put the production rule on each.
Assemble the two parts with `dx(definedonelements=...)` and check that
near + far with the unchanged rule reproduces the whole-mesh load to
round-off before comparing.  On a kinked test integrand each level reduced
the error about tenfold (order 8: 1.2e-4, 1.3e-5, 1.3e-6 for L = 0, 1, 2); an
order-8 tetrahedral rule has 125 points, so L = 2 costs 8000 points per near
element and L = 3 exceeds NGSolve's default local heap.

Measured effect on the C-type magnet (1,325 near air elements, reference
L = 2): B on lines through the coil leg moved from 3-5e-4 relative RMS
(production order 8) to 4-7e-5 with a single split level, the same accuracy
as order 20 on every element. The reported near-only L1 assembly was about
one thirteenth of whole-mesh order-20 assembly, excluding the far-element
assembly; this is not a total assembly or solver speedup.  The gap-core field changed by less than 1e-6 either way.

## Invariance to check first

In the mixed total/reduced Omega route the iron total-Hodge projection of
`H_s` cannot change B: a change of `Phi_s` inside the order-p space is
absorbed by `phi_total` through the interface jump.  Its quadrature bonus
moved the harmonic norm but changed gap B by 1e-14.  Vary a load that B
actually depends on before drawing conclusions.
"""


def get_force_methods(topic: str = "all") -> str:
    """Return unified electromagnetic-force theory by legacy subtopic."""

    return get_forces_documentation(topic)


def get_force_recipe(topic: str = "method_choice") -> str:
    """Return practical method-selection and NGSolve implementation guidance."""

    return get_em_force_ngsolve_recipe(topic)


def get_force_extras(topic: str = "all") -> str:
    """Return specialized electromagnetic-force formulations."""

    return get_em_force_extras(topic)


def get_force_validation(topic: str = "all") -> str:
    """Return the common force validation and evidence contract."""

    return get_force_validation_documentation(topic)


def get_force_knowledge(topic: str = "overview") -> str:
    """Dispatch the top-level common Force topics."""

    key = (topic or "overview").strip().lower().replace("-", "_")
    if key in {"overview", "intro", ""}:
        return OVERVIEW
    if key in {"methods", "theory", "force_methods", "differential_forms"}:
        return get_force_methods("all")
    if key in {"recipe", "method_choice", "implementation"}:
        return get_force_recipe("method_choice")
    if key in {"extras", "specialized", "specialised"}:
        return get_force_extras("all")
    if key in {"validation", "cross_validation", "evidence"}:
        return get_force_validation("all")
    if key in {"motor", "torque"}:
        return MOTOR_GUIDANCE
    if key in {"maglev", "levitation", "lift"}:
        return MAGLEV_GUIDANCE
    if key in {"source_quadrature", "near_source", "near_field_quadrature",
               "source_load_quadrature"}:
        return SOURCE_QUADRATURE
    if key == "all":
        return "\n\n".join(
            [
                OVERVIEW,
                get_force_methods("all"),
                get_force_recipe("all"),
                get_force_extras("all"),
                get_force_validation("all"),
                MOTOR_GUIDANCE,
                MAGLEV_GUIDANCE,
                SOURCE_QUADRATURE,
            ]
        )
    return f"Unknown topic '{topic}'. Available: {', '.join(TOPICS)}."
