"""Regression gates for the static-electromagnet three-way contract."""

from pathlib import Path

import pytest

from radia.electromagnet_validation import (
    ESRF_RADIA_SEVEN_CASES,
    STATIC_ELECTROMAGNET_FORMULATIONS,
    esrf_radia_seven_case_contract,
    require_static_electromagnet_three_engine_contract,
)
from radia.esrf_examples import build_esrf_coils, list_esrf_example_specs
from radia.static_electromagnet import StaticElectromagnetMixedDomain


ROOT = Path(__file__).resolve().parents[1]


def _diagnostics():
    return {
        name: {"formulation": formulation}
        for name, formulation in STATIC_ELECTROMAGNET_FORMULATIONS.items()
    }


def test_static_electromagnets_require_all_three_formulations():
    contract = require_static_electromagnet_three_engine_contract(_diagnostics())
    assert contract["h1_acceptance_route"] == "mixed_total_reduced_omega"
    assert contract["global_reduced_omega_acceptance_forbidden"] is True


def test_esrf_seven_case_corpus_is_source_traceable_and_physical():
    contract = esrf_radia_seven_case_contract()
    cases = contract["cases"]
    assert cases == list(ESRF_RADIA_SEVEN_CASES)
    assert [case["number"] for case in cases] == list(range(1, 8))
    assert [case["source_notebook"] for case in cases] == [
        f"Example#{number}.nb" for number in range(1, 8)
    ]
    assert cases[0]["source_kind"] == "fixed_magnetization"
    assert cases[0]["requires_three_engine_comparison"] is False
    assert cases[2]["source_kind"] == "hybrid_fixed_magnetization"
    assert cases[2]["requires_three_engine_comparison"] is True
    assert cases[4]["slug"] == "c_dipole"
    assert cases[6]["slug"] == "esrf_storage_ring_quadrupole"
    assert [case["number"] for case in cases if case["requires_three_engine_comparison"]] == [
        3,
        5,
        6,
        7,
    ]
    assert contract["required_evidence"][
        "nonlinear_mesh_convergence_when_material_is_nonlinear"
    ] is True
    assert contract["required_evidence"][
        "fixed_magnetization_source_projection_when_present"
    ] is True


def test_esrf_validation_registry_cannot_drift_from_source_definitions():
    registry = esrf_radia_seven_case_contract()["cases"]
    sources = list_esrf_example_specs()
    assert [(case["number"], case["slug"], case["source_notebook"])
            for case in registry] == [
        (source.number, source.slug, source.source_notebook)
        for source in sources
    ]


def test_esrf_source_kinds_match_the_executable_source_builders():
    cases = {
        case["number"]: case
        for case in esrf_radia_seven_case_contract()["cases"]
    }
    assert {number for number in cases if build_esrf_coils(number)} == {2, 5, 6, 7}
    assert cases[1]["source_kind"] == "fixed_magnetization"
    assert cases[3]["source_kind"] == "hybrid_fixed_magnetization"
    assert cases[4]["source_kind"] == "fixed_magnetization"


def test_global_reduced_omega_cannot_be_relabelled_as_mixed_h1():
    diagnostics = _diagnostics()
    diagnostics["mixed_total_reduced_omega"]["formulation"] = (
        "H1 global reduced Omega"
    )
    with pytest.raises(ValueError, match="TOSCA mixed"):
        require_static_electromagnet_three_engine_contract(diagnostics)


def test_static_electromagnet_contract_rejects_an_incomplete_result():
    diagnostics = _diagnostics()
    diagnostics.pop("reduced_a")
    with pytest.raises(ValueError, match="exactly"):
        require_static_electromagnet_three_engine_contract(diagnostics)


def test_mixed_h1_domain_requires_interfaces_and_kelvin_vertex_gauge():
    domain = StaticElectromagnetMixedDomain(
        reduced_materials=("air",),
        total_materials=("iron", "kelvin"),
        nonlinear_materials=("iron",),
    )
    domain.validate_mesh_labels(
        ("air", "iron", "kelvin"),
        ("iron_air_interface", "kelvin_int"),
        ("GND",),
    )
    assert domain.as_dict()["reduced_total_interface"] == "iron_air_interface"
    with pytest.raises(ValueError, match="missing_boundaries"):
        domain.validate_mesh_labels(("air", "iron", "kelvin"), (), ("GND",))
    with pytest.raises(ValueError, match="missing_bbboundaries"):
        domain.validate_mesh_labels(
            ("air", "iron", "kelvin"),
            ("iron_air_interface", "kelvin_int"),
        )
    with pytest.raises(ValueError, match="exhaustive declared material partition"):
        domain.validate_mesh_labels(
            ("air", "iron", "kelvin", "coil"),
            ("iron_air_interface", "kelvin_int"),
            ("GND",),
        )


def test_mixed_h1_uses_kelvin_vertex_gauge_not_surface_dirichlet_label():
    source = (ROOT / "src" / "radia" / "kelvin_solver.py").read_text(
        encoding="utf-8"
    )
    mixed = source.split(
        "def solve_magnetostatic_mixed_total_reduced_omega_kelvin(", 1
    )[1].split(
        "def solve_magnetostatic_mixed_total_reduced_omega_picard_kelvin(", 1
    )[0]
    assert 'dirichlet_bbbnd="GND"' in mixed
    assert "dirichlet_bbbnd=dirichlet_bbbnd" in mixed
    assert "dirichlet_bbnd=" not in mixed
