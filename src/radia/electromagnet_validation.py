"""Shared formulation contracts for static soft-iron electromagnets.

The HDiv-VIM acceptance route is intentionally formulation-independent.  A
static electromagnet is accepted only when the identical physical model has
been evaluated by HDiv-MMM, HCurl reduced-A, and the H1 TOSCA-style mixed
total/reduced Omega formulation.
"""
from __future__ import annotations

from collections.abc import Mapping


STATIC_ELECTROMAGNET_THREE_ENGINES = (
    "hdiv_mmm",
    "reduced_a",
    "mixed_total_reduced_omega",
)

STATIC_ELECTROMAGNET_FORMULATIONS = {
    "hdiv_mmm": "HDiv-MMM",
    "reduced_a": "HCurl reduced-A",
    "mixed_total_reduced_omega": "H1 TOSCA mixed total/reduced Omega",
}

# The seven source cases are the shipped ESRF Radia notebooks.  They are a
# traceable physical corpus, not a convenient collection of element fixtures.
# A fixed permanent magnet is a source case; a hybrid/yoked magnet has an
# additional nonlinear material lane.  The H1 route is mandatory whenever a
# source drives material response, but a source-only example must not be
# dressed up as a synthetic three-PDE comparison.
ESRF_RADIA_SEVEN_CASES = (
    {
        "number": 1,
        "slug": "pm_cube",
        "source_notebook": "Example#1.nb",
        "source_kind": "fixed_magnetization",
        "nonlinear_material": False,
        "requires_three_engine_comparison": False,
    },
    {
        "number": 2,
        "slug": "racetrack_coils",
        "source_notebook": "Example#2.nb",
        "source_kind": "coilbuilder_current",
        "nonlinear_material": False,
        "requires_three_engine_comparison": False,
    },
    {
        "number": 3,
        "slug": "hybrid_undulator",
        "source_notebook": "Example#3.nb",
        "source_kind": "hybrid_fixed_magnetization",
        "nonlinear_material": True,
        "requires_three_engine_comparison": True,
    },
    {
        "number": 4,
        "slug": "magnetized_sphere",
        "source_notebook": "Example#4.nb",
        "source_kind": "fixed_magnetization",
        "nonlinear_material": False,
        "requires_three_engine_comparison": False,
    },
    {
        "number": 5,
        "slug": "c_dipole",
        "source_notebook": "Example#5.nb",
        "source_kind": "coilbuilder_current",
        "nonlinear_material": True,
        "requires_three_engine_comparison": True,
    },
    {
        "number": 6,
        "slug": "quadrupole",
        "source_notebook": "Example#6.nb",
        "source_kind": "coilbuilder_current",
        "nonlinear_material": True,
        "requires_three_engine_comparison": True,
    },
    {
        "number": 7,
        "slug": "esrf_storage_ring_quadrupole",
        "source_notebook": "Example#7.nb",
        "source_kind": "coilbuilder_current",
        "nonlinear_material": True,
        "requires_three_engine_comparison": True,
    },
)


def static_electromagnet_three_engine_contract() -> dict[str, object]:
    """Return the only accepted three-formulation set.

    A historical global reduced-Omega calculation remains useful as diagnostic
    archaeology, but is never a substitute for the mixed H1 formulation.
    """
    return {
        "schema": "radia.static-electromagnet-three-engine/v1",
        "engines": list(STATIC_ELECTROMAGNET_THREE_ENGINES),
        "formulations": dict(STATIC_ELECTROMAGNET_FORMULATIONS),
        "h1_acceptance_route": "mixed_total_reduced_omega",
        "global_reduced_omega_acceptance_forbidden": True,
    }


def require_static_electromagnet_three_engine_contract(
    diagnostics: Mapping[str, Mapping[str, object]],
) -> dict[str, object]:
    """Reject a result that omits, relabels, or replaces an engine."""
    if not isinstance(diagnostics, Mapping):
        raise TypeError("static electromagnet engine diagnostics must be a mapping")
    expected = set(STATIC_ELECTROMAGNET_THREE_ENGINES)
    observed = set(diagnostics)
    if observed != expected:
        raise ValueError(
            "static electromagnet acceptance needs exactly "
            f"{list(STATIC_ELECTROMAGNET_THREE_ENGINES)}; "
            f"missing={sorted(expected - observed)}, extra={sorted(observed - expected)}"
        )
    for name, formulation in STATIC_ELECTROMAGNET_FORMULATIONS.items():
        diagnostic = diagnostics[name]
        if not isinstance(diagnostic, Mapping):
            raise TypeError(f"{name} diagnostics must be a mapping")
        if diagnostic.get("formulation") != formulation:
            raise ValueError(
                f"{name} must report formulation {formulation!r}; got "
                f"{diagnostic.get('formulation')!r}"
            )
    return static_electromagnet_three_engine_contract()


def esrf_radia_seven_case_contract() -> dict[str, object]:
    """Return completion requirements for the seven original Radia examples."""
    return {
        "schema": "radia.hdiv-mmm-esrf-seven-case-suite/v1",
        "cases": list(ESRF_RADIA_SEVEN_CASES),
        "three_engine_contract": static_electromagnet_three_engine_contract(),
        "required_evidence": {
            "shared_physical_input": True,
            "fixed_magnetization_source_projection_when_present": True,
            "linear_mesh_convergence": True,
            "nonlinear_mesh_convergence_when_material_is_nonlinear": True,
            "result_bearing_documentation": True,
        },
        "formulation_scope": {
            "three_engine_comparison": (
                "required only when a permanent-magnet or CoilBuilder source "
                "drives a material response"
            ),
            "source_only_examples": (
                "require the native Radia source oracle and source-field "
                "fidelity; a synthetic three-PDE comparison is forbidden"
            ),
        },
    }


__all__ = [
    "STATIC_ELECTROMAGNET_FORMULATIONS",
    "STATIC_ELECTROMAGNET_THREE_ENGINES",
    "ESRF_RADIA_SEVEN_CASES",
    "esrf_radia_seven_case_contract",
    "require_static_electromagnet_three_engine_contract",
    "static_electromagnet_three_engine_contract",
]
