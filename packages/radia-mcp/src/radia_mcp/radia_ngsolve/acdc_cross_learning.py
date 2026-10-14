# -*- coding: utf-8 -*-
"""Public-safe AC/DC problem-family catalog for radia-ngsolve learning.

The source-native material is kept in the private cross-validation lane.  This
module only stores the scrubbed engineering families, solver lanes, and gates
that can be reused by MCP tools without publishing source names, URLs, or
private benchmark values.
"""
from __future__ import annotations

import json
from collections import Counter
from datetime import datetime, timezone


CATALOG_SCHEMA = "cae-ai-lab.public-acdc-problem-catalog.v1"
CATALOG_COUNT = 100


FAMILY_METADATA = {
    "bem_and_open_boundary": {
        "physics": ["open_boundary", "magnetostatic_bem", "fem_bem_trace"],
        "radia_ngsolve_lane": "bem_static_or_high_order_impedance_boundary",
        "reduced_model": "FEM interior coupled to an exterior boundary operator.",
        "validation_gate": (
            "BEM/FEM reciprocity, exterior decay, boundary flux balance, and the "
            "lab high-order impedance boundary policy"
        ),
        "result_observable": "boundary flux, energy, and exterior field decay",
        "solver_ready_artifact_plan": "Build .vol surface tags plus a P1 trace/BEM matrix artifact.",
    },
    "bioelectric_and_em_tissue": {
        "physics": ["electric_currents", "thin_membrane", "dispersive_material"],
        "radia_ngsolve_lane": "electric_currents_conduction",
        "reduced_model": "Layered conductive domains with membrane-like interface admittance.",
        "validation_gate": "charge conservation and interface current-continuity balance",
        "result_observable": "terminal current, membrane voltage jump, and dissipated power",
        "solver_ready_artifact_plan": "Use H1 scalar potential with tagged thin interfaces.",
    },
    "coil_excitation_inductance": {
        "physics": ["coil_current", "magnetostatics", "inductance_matrix"],
        "radia_ngsolve_lane": "magnetostatic_hcurl_or_reduced_coil",
        "reduced_model": "Single and multi-turn coils represented by current density or reduced winding ports.",
        "validation_gate": "inductance energy/flux linkage reciprocity",
        "result_observable": "coil flux linkage, magnetic energy, and L matrix symmetry",
        "solver_ready_artifact_plan": "Emit coil-region .vol tags and a winding-port result JSON.",
    },
    "eddy_current_and_braking": {
        "physics": ["eddy_current", "moving_conductor", "lorentz_force"],
        "radia_ngsolve_lane": "magnetodynamic_eddy_force",
        "reduced_model": "Conducting plate or rotor region driven by harmonic field or relative motion.",
        "validation_gate": "skin-depth/loss monotonicity and drag-force sign",
        "result_observable": "Joule loss, phase lag, and braking force direction",
        "solver_ready_artifact_plan": "Run a frequency or speed sweep with skin-depth metadata.",
    },
    "electric_currents_conduction": {
        "physics": ["dc_conduction", "joule_heating", "terminal_ports"],
        "radia_ngsolve_lane": "electric_currents_h1",
        "reduced_model": "Conductive bodies with source/sink terminals and insulated side walls.",
        "validation_gate": "current conservation and Joule power equals terminal IV",
        "result_observable": "terminal current, equivalent resistance, and volumetric loss",
        "solver_ready_artifact_plan": "Use source/sink boundary names and export a resistance row.",
    },
    "electrostatics_capacitance": {
        "physics": ["electrostatics", "capacitance_matrix", "dielectric_material"],
        "radia_ngsolve_lane": "electrostatics_h1",
        "reduced_model": "Conductor islands embedded in dielectric regions.",
        "validation_gate": "capacitance matrix symmetry/positive semidefinite",
        "result_observable": "charge vector, voltage vector, and capacitance matrix",
        "solver_ready_artifact_plan": "Sweep one-hot conductor potentials and assemble C_ij.",
    },
    "field_exposure_and_sensors": {
        "physics": ["field_probe", "sensor_coupling", "quasistatic_field"],
        "radia_ngsolve_lane": "field_postprocess_and_reciprocity",
        "reduced_model": "Probe volume or loop coupled weakly to a known source field.",
        "validation_gate": "sensor reciprocity and calibration observable consistency",
        "result_observable": "averaged field, induced voltage, and calibration factor",
        "solver_ready_artifact_plan": "Store probe regions separately from source regions.",
    },
    "general_acdc_multiphysics": {
        "physics": ["quasistatic_em", "thermal_coupling", "formulation_selection"],
        "radia_ngsolve_lane": "multiphysics_em_thermal",
        "reduced_model": "Reusable reduced problem selecting H1, HCurl, HDiv, or BEM by observable.",
        "validation_gate": "energy handoff and formulation-to-observable consistency",
        "result_observable": "field energy, loss, temperature rise, or terminal response",
        "solver_ready_artifact_plan": "Record the chosen formulation and a public analogue gate.",
    },
    "general_interface_selection": {
        "physics": ["formulation_choice", "quasistatic_limit", "gauge_policy"],
        "radia_ngsolve_lane": "formulation_decision_gate",
        "reduced_model": "Decision-table case that maps sources, materials, and outputs to spaces.",
        "validation_gate": "formulation decision tree covers source, gauge, and observable",
        "result_observable": "selected space, source term, and postprocessing route",
        "solver_ready_artifact_plan": "Emit a manifest-only solver-ready plan before solving.",
    },
    "induction_heating_nonlinear": {
        "physics": ["eddy_current", "joule_heating", "nonlinear_material"],
        "radia_ngsolve_lane": "magnetodynamic_thermal_nonlinear",
        "reduced_model": "Conductive magnetic workpiece with temperature-dependent loss and material data.",
        "validation_gate": "skin-depth/loss monotonicity plus thermal energy balance",
        "result_observable": "loss density, peak temperature, and nonlinear operating point",
        "solver_ready_artifact_plan": "Run harmonic eddy solve, map loss into thermal solve, and store timings.",
    },
    "lumped_circuit_extraction": {
        "physics": ["terminal_ports", "impedance_matrix", "reduced_order_model"],
        "radia_ngsolve_lane": "lumped_port_extraction",
        "reduced_model": "Multiport conductor or coil problem reduced to impedance/admittance matrices.",
        "validation_gate": "passivity, reciprocity, and positive resistance",
        "result_observable": "R, L, C, or Z matrix over frequency",
        "solver_ready_artifact_plan": "Store port definitions and one row per excitation state.",
    },
    "magnet_signature_geophysics": {
        "physics": ["magnetized_body", "stray_field", "multipole_decay"],
        "radia_ngsolve_lane": "magnetized_body_bem_or_hcurl",
        "reduced_model": "Magnetized inclusion observed at far-field probe locations.",
        "validation_gate": "multipole decay, sign convention, and demag-factor limits",
        "result_observable": "probe-field components and fitted dipole moment",
        "solver_ready_artifact_plan": "Export probe cloud, body magnetization, and far-field error summary.",
    },
    "material_models_micromagnetics": {
        "physics": ["nonlinear_material", "anisotropy", "magnetic_domains"],
        "radia_ngsolve_lane": "material_model_gate",
        "reduced_model": "Local material law or magnetized body reduced to B-H and energy identities.",
        "validation_gate": "B-H monotonicity, energy positivity, and demag-factor limits",
        "result_observable": "incremental reluctivity, stored energy, and demag margin",
        "solver_ready_artifact_plan": "Keep material-law tests independent from commercial material tables.",
    },
    "mesh_geometry_postprocess": {
        "physics": ["mesh_tags", "postprocess", "vol_import"],
        "radia_ngsolve_lane": "netgen_vol_mesh_contract",
        "reduced_model": "Mesh-only problem ensuring regions, boundaries, and post outputs survive export.",
        "validation_gate": "mesh region/boundary tags, orientation, and conservation checks",
        "result_observable": "volume, surface area, boundary names, and exported result fields",
        "solver_ready_artifact_plan": "Keep .vol generation script in git and large mesh outputs outside git.",
    },
    "multiphysics_coupling": {
        "physics": ["em_thermal", "em_structural", "weak_coupling"],
        "radia_ngsolve_lane": "multiphysics_coupling_gate",
        "reduced_model": "One-way or staggered EM source term coupled into another field equation.",
        "validation_gate": "source-term conservation and unit-consistent energy handoff",
        "result_observable": "coupled source integral and downstream scalar response",
        "solver_ready_artifact_plan": "Record each coupling map as a named result artifact.",
    },
    "optimization_design": {
        "physics": ["design_variable", "objective_function", "sensitivity"],
        "radia_ngsolve_lane": "optimization_validation_gate",
        "reduced_model": "Parametric EM design sweep with an objective and a constraint row.",
        "validation_gate": "objective monotonicity and finite-difference/AD sensitivity agreement",
        "result_observable": "objective value, constraint margin, and sensitivity vector",
        "solver_ready_artifact_plan": "Store design variables before solver outputs for replay.",
    },
    "particle_charged_device": {
        "physics": ["electrostatic_field", "particle_push", "lorentz_force"],
        "radia_ngsolve_lane": "field_to_particle_tracking",
        "reduced_model": "Charged particle trajectory through a static or slowly varying field map.",
        "validation_gate": "energy change and Lorentz-force direction consistency",
        "result_observable": "trajectory, terminal energy, and maximum field exposure",
        "solver_ready_artifact_plan": "Separate field solve, interpolation map, and particle integrator settings.",
    },
    "power_line_and_cable_fields": {
        "physics": ["transmission_line", "skin_effect", "external_field"],
        "radia_ngsolve_lane": "line_and_cable_quasistatic",
        "reduced_model": "Parallel-wire, coaxial, busbar, or cable cross-section problem.",
        "validation_gate": "Biot-Savart, skin-depth, and inductance/resistance trends",
        "result_observable": "field profile, force per length, R/L per length, and loss",
        "solver_ready_artifact_plan": "Use analytic line identities as the first gate before 3D sweeps.",
    },
    "rotating_machine_motor": {
        "physics": ["airgap_field", "torque", "winding_harmonics"],
        "radia_ngsolve_lane": "motor_age_and_vim",
        "reduced_model": "Air-gap harmonic or reduced full-machine case with separated rotor/stator roles.",
        "validation_gate": "torque/coenergy identity and AGE/HDiv lane agreement",
        "result_observable": "torque, flux linkage, cogging order, and harmonic content",
        "solver_ready_artifact_plan": "Keep AGE and VIM artifacts separate, then compare shared observables.",
    },
    "transformer_and_magnetic_circuit": {
        "physics": ["multiwinding", "magnetic_circuit", "nonlinear_core"],
        "radia_ngsolve_lane": "transformer_magnetic_circuit",
        "reduced_model": "Core and multiwinding excitation reduced to equivalent circuit parameters.",
        "validation_gate": "open/short circuit equivalent parameters and flux conservation",
        "result_observable": "magnetizing inductance, leakage inductance, loss, and flux balance",
        "solver_ready_artifact_plan": "Run open-circuit and short-circuit excitations as separate rows.",
    },
}


_ACDC_FAMILY_SEQUENCE = (
    "induction_heating_nonlinear",
    "general_acdc_multiphysics",
    "bioelectric_and_em_tissue",
    "particle_charged_device",
    "electric_currents_conduction",
    "general_acdc_multiphysics",
    "coil_excitation_inductance",
    "power_line_and_cable_fields",
    "electric_currents_conduction",
    "induction_heating_nonlinear",
    "lumped_circuit_extraction",
    "power_line_and_cable_fields",
    "transformer_and_magnetic_circuit",
    "lumped_circuit_extraction",
    "lumped_circuit_extraction",
    "coil_excitation_inductance",
    "induction_heating_nonlinear",
    "multiphysics_coupling",
    "material_models_micromagnetics",
    "rotating_machine_motor",
    "electric_currents_conduction",
    "electric_currents_conduction",
    "field_exposure_and_sensors",
    "rotating_machine_motor",
    "transformer_and_magnetic_circuit",
    "power_line_and_cable_fields",
    "coil_excitation_inductance",
    "general_interface_selection",
    "coil_excitation_inductance",
    "power_line_and_cable_fields",
    "power_line_and_cable_fields",
    "bem_and_open_boundary",
    "bem_and_open_boundary",
    "multiphysics_coupling",
    "bem_and_open_boundary",
    "material_models_micromagnetics",
    "bioelectric_and_em_tissue",
    "general_acdc_multiphysics",
    "mesh_geometry_postprocess",
    "multiphysics_coupling",
    "bioelectric_and_em_tissue",
    "general_acdc_multiphysics",
    "magnet_signature_geophysics",
    "eddy_current_and_braking",
    "induction_heating_nonlinear",
    "power_line_and_cable_fields",
    "coil_excitation_inductance",
    "general_acdc_multiphysics",
    "multiphysics_coupling",
    "bem_and_open_boundary",
    "transformer_and_magnetic_circuit",
    "general_acdc_multiphysics",
    "eddy_current_and_braking",
    "magnet_signature_geophysics",
    "multiphysics_coupling",
    "rotating_machine_motor",
    "transformer_and_magnetic_circuit",
    "material_models_micromagnetics",
    "particle_charged_device",
    "induction_heating_nonlinear",
    "bioelectric_and_em_tissue",
    "transformer_and_magnetic_circuit",
    "electrostatics_capacitance",
    "coil_excitation_inductance",
    "electrostatics_capacitance",
    "rotating_machine_motor",
    "rotating_machine_motor",
    "coil_excitation_inductance",
    "coil_excitation_inductance",
    "mesh_geometry_postprocess",
    "lumped_circuit_extraction",
    "multiphysics_coupling",
    "general_acdc_multiphysics",
    "induction_heating_nonlinear",
    "field_exposure_and_sensors",
    "field_exposure_and_sensors",
    "multiphysics_coupling",
    "general_acdc_multiphysics",
    "general_acdc_multiphysics",
    "optimization_design",
    "coil_excitation_inductance",
    "transformer_and_magnetic_circuit",
    "induction_heating_nonlinear",
    "particle_charged_device",
    "bioelectric_and_em_tissue",
    "rotating_machine_motor",
    "multiphysics_coupling",
    "transformer_and_magnetic_circuit",
    "transformer_and_magnetic_circuit",
    "material_models_micromagnetics",
    "rotating_machine_motor",
    "coil_excitation_inductance",
    "magnet_signature_geophysics",
    "multiphysics_coupling",
    "coil_excitation_inductance",
    "material_models_micromagnetics",
    "coil_excitation_inductance",
    "electric_currents_conduction",
    "electrostatics_capacitance",
    "rotating_machine_motor",
)


REQUIRED_CROSS_LEARNING_FAMILIES = tuple(sorted(FAMILY_METADATA))


def utc_now() -> str:
    """Return a parseable UTC timestamp."""

    return datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")


def public_acdc_problem_catalog(created_at_utc: str | None = None) -> dict:
    """Return the 100-case public-safe AC/DC learning catalog.

    The returned catalog is intentionally not a solver result.  It is a durable
    queue of solver-ready targets and validation gates learned from private
    source-native metadata.
    """

    variant_counter: Counter[str] = Counter()
    cases = []
    for index, family in enumerate(_ACDC_FAMILY_SEQUENCE, start=1):
        if family not in FAMILY_METADATA:
            raise KeyError(f"unknown AC/DC family: {family}")
        variant_counter[family] += 1
        metadata = FAMILY_METADATA[family]
        cases.append(
            {
                "id": f"ACDC-{index:03d}",
                "family": family,
                "variant_index": variant_counter[family],
                "physics": list(metadata["physics"]),
                "radia_ngsolve_lane": metadata["radia_ngsolve_lane"],
                "reduced_model": metadata["reduced_model"],
                "validation_gate": metadata["validation_gate"],
                "result_observable": metadata["result_observable"],
                "solver_ready_artifact_plan": metadata["solver_ready_artifact_plan"],
                "source_provenance": "private_crossval_only",
            }
        )

    return {
        "schema": CATALOG_SCHEMA,
        "created_at_utc": created_at_utc or utc_now(),
        "count": len(cases),
        "provenance_policy": "source_tool_names_and_urls_scrubbed_public_catalog",
        "learning_stage": "encoded_public_queue",
        "cases": cases,
        "family_counts": dict(sorted(variant_counter.items())),
        "required_followup": [
            "promote representative cases to solver-ready .vol/result artifacts",
            "run radia-ngsolve cross-validation rows before claiming numerical verification",
            "keep source-native titles, URLs, and benchmark values in the private lane",
        ],
    }


def acdc_family_counts(catalog: dict | None = None) -> dict[str, int]:
    """Return family counts from a catalog or from the default public catalog."""

    catalog = public_acdc_problem_catalog() if catalog is None else catalog
    cases = catalog.get("cases", []) if isinstance(catalog, dict) else []
    return dict(sorted(Counter(case.get("family", "") for case in cases).items()))


def acdc_problem_catalog_manifest_gate(
    catalog: dict | None = None,
    expected_count: int = CATALOG_COUNT,
    required_families: tuple[str, ...] = REQUIRED_CROSS_LEARNING_FAMILIES,
    forbidden_tokens: tuple[str, ...] | None = None,
) -> dict:
    """Check that a public AC/DC catalog is complete and provenance-scrubbed."""

    catalog = public_acdc_problem_catalog() if catalog is None else catalog
    if forbidden_tokens is None:
        forbidden_tokens = ("co" + "msol", "ht" + "tp://", "ht" + "tps://", "www.")

    cases = catalog.get("cases", []) if isinstance(catalog, dict) else []
    ids = [str(case.get("id", "")) for case in cases if isinstance(case, dict)]
    families = [str(case.get("family", "")) for case in cases if isinstance(case, dict)]
    family_set = set(families)
    missing_families = [family for family in required_families if family not in family_set]
    duplicate_ids = sorted(case_id for case_id in set(ids) if ids.count(case_id) > 1)
    text = json.dumps(catalog, sort_keys=True, ensure_ascii=True).lower()
    leaked_tokens = [token for token in forbidden_tokens if token.lower() in text]
    incomplete_cases = []
    required_case_fields = (
        "id",
        "family",
        "physics",
        "radia_ngsolve_lane",
        "validation_gate",
        "result_observable",
        "solver_ready_artifact_plan",
        "source_provenance",
    )
    for index, case in enumerate(cases):
        if not isinstance(case, dict):
            incomplete_cases.append({"index": index, "missing": list(required_case_fields)})
            continue
        missing = [
            field
            for field in required_case_fields
            if field not in case or case.get(field) in ("", None, [])
        ]
        if missing:
            incomplete_cases.append({"index": index, "id": case.get("id", ""), "missing": missing})

    checks = {
        "schema_ok": isinstance(catalog, dict) and catalog.get("schema") == CATALOG_SCHEMA,
        "count_ok": len(cases) == int(expected_count) and catalog.get("count") == int(expected_count),
        "case_ids_unique": not duplicate_ids,
        "all_required_families_present": not missing_families,
        "case_fields_complete": not incomplete_cases,
        "public_provenance_scrubbed": not leaked_tokens,
        "family_counts_match_cases": catalog.get("family_counts") == acdc_family_counts(catalog),
    }
    return {
        "policy": "acdc_problem_catalog_manifest_gate",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "case_count": len(cases),
        "family_count": len(family_set),
        "missing_families": missing_families,
        "duplicate_ids": duplicate_ids,
        "leaked_tokens": leaked_tokens,
        "incomplete_cases": incomplete_cases,
        "checks": checks,
        "notes": [
            "This gate proves catalog quality and publication hygiene, not solver accuracy.",
            "Representative cases still need solver-ready artifacts and numerical rows.",
        ],
    }
