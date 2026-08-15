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
`compute_lorentz_force_result_via_foster` emits conductor/source reaction
records and their residual, while `PositionForceCurve.force_result_at` converts
a CLN position-force interpolation into the same common result schema.
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
            ]
        )
    return f"Unknown topic '{topic}'. Available: {', '.join(TOPICS)}."
