# -*- coding: utf-8 -*-
"""Public AC/DC learning catalog should be complete and provenance-scrubbed."""

from radia_mcp.radia_ngsolve.acdc_cross_learning import (
    acdc_family_counts,
    acdc_problem_catalog_manifest_gate,
    public_acdc_problem_catalog,
)
from radia_mcp.radia_ngsolve.knowledge.ngsolve import get_ngsolve_documentation
from radia_mcp.radia_ngsolve.server import ngsolve_usage


EXPECTED_COUNTS = {
    "bem_and_open_boundary": 4,
    "bioelectric_and_em_tissue": 5,
    "coil_excitation_inductance": 12,
    "eddy_current_and_braking": 2,
    "electric_currents_conduction": 5,
    "electrostatics_capacitance": 3,
    "field_exposure_and_sensors": 3,
    "general_acdc_multiphysics": 9,
    "general_interface_selection": 1,
    "induction_heating_nonlinear": 7,
    "lumped_circuit_extraction": 4,
    "magnet_signature_geophysics": 3,
    "material_models_micromagnetics": 5,
    "mesh_geometry_postprocess": 2,
    "multiphysics_coupling": 9,
    "optimization_design": 1,
    "particle_charged_device": 3,
    "power_line_and_cable_fields": 6,
    "rotating_machine_motor": 8,
    "transformer_and_magnetic_circuit": 8,
}


def _first_case(catalog, family):
    return next(case for case in catalog["cases"] if case["family"] == family)


def test_public_acdc_problem_catalog_has_100_scrubbed_cases():
    catalog = public_acdc_problem_catalog(created_at_utc="2026-07-05T00:00:00Z")
    gate = acdc_problem_catalog_manifest_gate(catalog)

    assert catalog["count"] == 100
    assert len(catalog["cases"]) == 100
    assert gate["status"] == "ok"
    assert gate["checks"]["public_provenance_scrubbed"] is True
    assert acdc_family_counts(catalog) == EXPECTED_COUNTS


def test_public_acdc_problem_catalog_does_not_publish_source_native_provenance():
    catalog = public_acdc_problem_catalog(created_at_utc="2026-07-05T00:00:00Z")
    text = repr(catalog).lower()

    assert ("co" + "msol") not in text
    assert ("ht" + "tp") not in text
    assert "www." not in text
    assert all(case["source_provenance"] == "private_crossval_only" for case in catalog["cases"])


def test_acdc_catalog_records_solver_ready_validation_gates_for_key_families():
    catalog = public_acdc_problem_catalog(created_at_utc="2026-07-05T00:00:00Z")

    coil = _first_case(catalog, "coil_excitation_inductance")
    assert "inductance energy/flux linkage reciprocity" in coil["validation_gate"]
    assert coil["radia_ngsolve_lane"] == "magnetostatic_hcurl_or_reduced_coil"

    eddy = _first_case(catalog, "eddy_current_and_braking")
    assert "skin-depth/loss monotonicity" in eddy["validation_gate"]
    assert eddy["result_observable"] == "Joule loss, phase lag, and braking force direction"

    transformer = _first_case(catalog, "transformer_and_magnetic_circuit")
    assert "open/short circuit equivalent parameters" in transformer["validation_gate"]

    motor = _first_case(catalog, "rotating_machine_motor")
    assert "torque/coenergy identity" in motor["validation_gate"]
    assert "AGE/HDiv" in motor["validation_gate"]

    electrostatic = _first_case(catalog, "electrostatics_capacitance")
    assert "capacitance matrix symmetry/positive semidefinite" in electrostatic["validation_gate"]

    open_boundary = _first_case(catalog, "bem_and_open_boundary")
    assert "high-order impedance boundary" in open_boundary["validation_gate"]


def test_acdc_cross_learning_topic_is_exposed_to_ngsolve_usage():
    doc = get_ngsolve_documentation("acdc_problem_catalog")
    tool_doc = ngsolve_usage("acdc_cross_learning")

    assert "public-safe 100-case learning queue" in doc
    assert "public_acdc_problem_catalog()" in doc
    assert "AGE/HDiv lane agreement" in tool_doc
    assert "Unknown topic" not in tool_doc[:120]
