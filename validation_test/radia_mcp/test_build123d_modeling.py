# -*- coding: utf-8 -*-
r"""Parametric modeling ops for build123d (radia_mcp.build123d.modeling) -- geometry-gated.

Each op is checked against a closed-form volume / count where possible, for OCCT validity, for the
region label, and -- the point of "CAE-safe" -- that the result MESHES cleanly in Netgen (the
build123d -> Netgen -> Radia/NGSolve tet pipeline).
"""
import math
import os
import sys
import json

import pytest

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)
from radia_mcp.build123d.modeling import (annular_segment, tube, racetrack_coil, polar_array,
                                          linear_array, mirrored, assembly,
                                          shape_envelope_row, enclosing_box,
                                          enclosure_clearance_row, enclosure_difference_region,
                                          shape_measurement_row, shape_measurement_rows,
                                          box_through_cylinder_reference_row,
                                          mounting_plate_boss_reference_row,
                                          keyed_terminal_plate_reference_row,
                                          flanged_sleeve_reference_row,
                                          coax_annular_sleeve_reference_row,
                                          ribbed_busbar_heat_sink_reference_row,
                                          three_phase_busbar_snubber_plate_reference_row,
                                          rcd_snubber_heat_spreader_reference_row,
                                          rcd_snubber_capacitance_sweep_rows,
                                          thermal_robin_cooling_plate_reference_row,
                                          motor_housing_radial_fin_reference_row,
                                          v_type_ipm_rotor_coupon_reference_row,
                                          box_face_vector_area_rows,
                                          box_face_pressure_force_rows,
                                          box_face_pressure_moment_rows,
                                          box_face_pressure_resultant_summary,
                                          box_face_traction_moment_rows,
                                          compare_boundary_vector_area_rows,
                                          compare_shape_measurement_rows,
                                          compare_shape_volume_rows,
                                          shape_name_identity_gate,
                                          shape_role_metadata_gate,
                                          shape_transition_role_metadata_gate,
                                          shape_cubit_meshing_scheme_intent_gate,
                                          cst_cad_volume_export_manifest_gate,
                                          shape_volume_crosscheck_summary,
                                          shape_perforated_prism_roundtrip_gate,
                                          shape_volume_crosscheck_source_coverage_gate,
                                          shape_volume_crosscheck_source_identity_gate,
                                          shape_external_cad_volume_evidence_package_gate,
                                          shape_cad_route_source_contract_gate,
                                          shape_mass_property_crosscheck_summary,
                                          shape_cubit_export_package_handoff_gate,
                                          shape_cubit_quality_package_handoff_gate,
                                          shape_cubit_quality_ledger_handoff_gate,
                                          shape_cubit_solver_route_handoff_gate,
                                          shape_cad_handoff_manifest_gate,
                                          shape_submodel_cad_handoff_gate,
                                          shape_curvilinear_mesh_intent_gate,
                                          shape_mesh_environment_handoff_gate,
                                          shape_measurement_comparison_summary,
                                          shape_measurement_inventory_summary,
                                          worst_shape_measurement_comparison_rows,
                                          shape_measurement_health_summary,
                                          shape_bbox_pair_clearance_summary,
                                          shape_parameter_sweep_summary)
from radia_mcp.build123d.build123d_knowledge import get_build123d_documentation
from radia_mcp.cubit.vol_inventory import (
    cubit_export_package_identity_gate,
    cubit_headless_installation_route_gate,
    cubit_mass_property_sidecar_gate,
    cubit_meshing_scheme_trace_gate,
    cubit_mixed_solver_route_manifest_gate,
)
from build123d import Box, Compound, Pos
from radia_mcp.build123d.server import build123d_motor_housing_thermal_reference


def test_motor_housing_radial_fin_reference_closes_live_cubit_roundtrip():
    row = motor_housing_radial_fin_reference_row(42.0, 50.0, 90.0, 8, 12.0, 3.0)
    assert row["policy"] == "analytic_motor_housing_radial_fin_multibody_reference"
    assert row["body_count"] == 9
    assert row["volume"] == pytest.approx(234019.097373788)
    assert row["area"] == pytest.approx(78825.19872953116)

    measured = {"cubit": [{"name": row["name"], "volume": 233997.78602906506}]}
    summary = shape_volume_crosscheck_summary(
        [row], measured, rtol=row["roundtrip_tolerances"]["curved_step_volume_rtol"]
    )
    assert summary["status"] == "ok"
    assert summary["max_volume_rel_error"] == pytest.approx(9.10666905484026e-5)

    payload = json.loads(build123d_motor_housing_thermal_reference())
    assert payload["body_count"] == 9
    assert payload["roundtrip_tolerances"]["body_count_exact"] is True


def test_motor_housing_radial_fin_reference_rejects_overlap_and_invalid_radii():
    with pytest.raises(ValueError, match="outer_radius"):
        motor_housing_radial_fin_reference_row(50.0, 42.0, 90.0, 8, 12.0, 3.0)
    with pytest.raises(ValueError, match="separated radial fins"):
        motor_housing_radial_fin_reference_row(42.0, 50.0, 90.0, 8, 12.0, 40.0)


def _vol(obj):
    return sum(s.volume for s in obj.solids()) if isinstance(obj, Compound) else obj.volume


def _valid(obj):
    return all(s.is_valid for s in obj.solids()) if isinstance(obj, Compound) else obj.is_valid


def test_build123d_lab_policy_routes_tet_to_netgen_and_mixed_to_cubit():
    doc = get_build123d_documentation("lab_policy")

    assert "tet-only path" in doc
    assert "build123d → Netgen (tet) → Radia" in doc
    assert "hex path + mixed-mesh CAD fallback" in doc
    assert "pyramid transition elements" in doc
    assert "cubit_vol_inventory" in doc
    assert "box_through_cylinder_reference_row" in doc
    assert "volume as the common currency" in doc
    assert "build123d_volume_crosscheck" in doc
    assert "shape_volume_crosscheck_source_coverage_gate" in doc
    assert "(\"cubit\", \"cst_import\")" in doc
    assert "Slot331 extends the volume crosscheck with source-identity metadata" in doc
    assert "shape_volume_crosscheck_source_identity_gate" in doc
    assert "build123d_volume_crosscheck_source_identity_gate" in doc
    assert "measurement_method" in doc
    assert "source_artifact_id" in doc
    assert "shape_external_cad_volume_evidence_package_gate" in doc
    assert "build123d_external_cad_volume_evidence_package" in doc
    assert "Slot339 bundles coverage, source identity, and route contract" in doc
    assert "mesh-volume-after-import" in doc
    assert "Slot391 extends the same volume crosscheck" in doc
    assert "parameter_set_artifact_id" in doc
    assert "parameter_set_digest" in doc
    assert "parameter_set_path" in doc
    assert "objective_observable_id" in doc
    assert "objective_observable_family" in doc
    assert "stale sizing parameters" in doc
    assert "Slot398 extends the downstream mesh-environment handoff" in doc
    assert "shape_mesh_environment_handoff_gate" in doc
    assert "sanitized `version_probe_summary`" in doc
    assert "stale summary version line" in doc
    assert "Slot441 keeps Cubit headless teardown warnings separate" in doc
    assert "do not treat the exit code alone as\na build123d volume mismatch" in doc
    assert "source`, `sink`, `sibc`, `coil`, `workpiece`, and\n`air" in doc
    assert "Slot449 adds the one-loop learning closure rule for build123d" in doc
    assert "learning_lanes.source_tool=verified" in doc
    assert "CST are\ncross-kernel volume witnesses" in doc
    assert "shape_mass_property_crosscheck_summary" in doc
    assert "build123d_mass_property_crosscheck" in doc
    assert "Slot211 tightens that bridge" in doc
    assert "volume_unit" in doc
    assert "mm^3" in doc
    assert "shape_cubit_export_package_handoff_gate" in doc
    assert "shape_cubit_quality_package_handoff_gate" in doc
    assert "headless mesh-quality package" in doc
    assert "Slot370 adds the stricter Cubit quality-ledger handoff" in doc
    assert "shape_cubit_quality_ledger_handoff_gate" in doc
    assert "cubit_quality_ledger_handoff" in doc
    assert "cubit_quality_ledger_json" in doc
    assert "mesh_digest" in doc
    assert "Slot384 carries the Coreform mixed-route reader contract" in doc
    assert "require_solver_contract_artifact=True" in doc
    assert "solver_contract_digest" in doc
    assert "mixed-element reader contract digest is stale" in doc
    assert "Slot419 extends the same mixed-route handoff" in doc
    assert "solver_route_convention_schema_id" in doc
    assert "require_solver_route_convention_schema=True" in doc
    assert "value-only or missing route\nconvention" in doc
    assert "shape_cad_handoff_manifest_gate" in doc
    assert "Slot347 adds the Cubit solver-route handoff" in doc
    assert "shape_cubit_solver_route_handoff_gate" in doc
    assert "cubit_solver_route_handoff" in doc
    assert "cubit_mixed_solver_route_manifest_gate" in doc
    assert "no_implicit_tetization=true" in doc
    assert "tet_only_owner=netgen_tri_tet_path" in doc
    assert "Slot251 adds the local submodel CAD preflight" in doc
    assert "shape_submodel_cad_handoff_gate" in doc
    assert "crop_box" in doc
    assert "boundary-handoff gate" in doc
    assert "Slot259 extends this handoff" in doc
    assert "transition_handoff" in doc
    assert "volume_kind_counts` containing a `pyramid" in doc
    assert "shape_curvilinear_mesh_intent_gate" in doc
    assert "cubit_curvilinear_handoff" in doc
    assert "projection_error_within_tolerance" in doc
    assert "negative_jacobian_count_zero" in doc
    assert "Slot323 also lets this gate consume" in doc
    assert "cubit_mixed_order_series_inventory_gate" in doc
    assert "non-curved first-order inventory" in doc
    assert "shape_mesh_environment_handoff_gate" in doc
    assert "require_export_inventory=True" in doc
    assert "export_inventory" in doc
    assert "cubit_headless_installation_route_gate" in doc
    assert "release note is only a watchlist" in doc
    assert "Slot227 extends this bridge" in doc
    assert 'license_status="ValidStudent"' in doc
    assert "coreform_cubit.com -version" in doc
    assert "not CAD evidence" in doc
    assert "Slot355 carries the Slot354 console-binary rule" in doc
    assert "binary_path_is_console_com=True" in doc
    assert "version_probe_uses_recorded_binary=True" in doc
    assert "coreform_cubit.exe" in doc
    assert "external volume summary JSON" in doc
    assert "final build123d-side preflight" in doc
    assert "stale STEP or stale quality JSON" in doc
    assert "Slot315 extends the same CAD handoff" in doc
    assert "cad_measurement_convention" in doc
    assert "occt_closed_solid_mass_properties" in doc
    assert "mesh_volume_after_import" in doc
    assert "Slot426 extends the same CAD handoff" in doc
    assert "cad_measurement_postprocess_row_convention_schema_id" in doc
    assert "build123d_occt_mass_property_row_convention_v1" in doc
    assert "require_measurement_postprocess_row_convention_schema=True" in doc
    assert "Slot433 extends the same CAD handoff" in doc
    assert "cad_measurement_component_basis_schema_id" in doc
    assert "build123d_occt_volume_area_bbox_component_basis_v1" in doc
    assert "require_measurement_component_basis_schema=True" in doc
    assert "stale scalar-volume\ncomponent basis" in doc
    assert "Slot235 extends this handoff" in doc
    assert "tri/tet-only" in doc
    assert "explicit `geometry_id`" in doc
    assert "stale sidecars and wrong-geometry" in doc
    assert "Slot405 extends that CAD-to-Cubit handoff" in doc
    assert "require_sidecar_inventory_counts=True" in doc
    assert "vol_sidecar_element_count_matches_inventory" in doc
    assert "vol_sidecar_order_matches_expected" in doc
    assert "Slot412 extends the same handoff" in doc
    assert "vol_sidecar_schema_id" in doc
    assert "require_sidecar_schema=True" in doc
    assert "coreform_cubit.com -nographics -batch" in doc
    assert '"box_hole", "volume": 22.994690350851265' in doc
    assert '"l_bracket_two_holes", "volume": 2.8982123980236905' in doc
    assert "boolean union/overlap accounting" in doc
    assert "mounting_plate_boss_reference_row" in doc
    assert "three_phase_busbar_snubber_plate_reference_row" in doc
    assert '"mounting_plate_boss_five_holes", "volume": 12.786811880091562' in doc
    assert "min_edge_over_characteristic" in doc
    assert "max_volume_rel_error = 0.0" in doc
    assert "stepped cylindrical spacer" in doc
    assert "rel_error = 5.854366888206992e-6" in doc
    assert "shape_volume_crosscheck_summary(..., rtol=1e-5)" in doc
    assert "shape_name_identity_gate" in doc
    assert "shape_role_metadata_gate" in doc
    assert "shape_transition_role_metadata_gate" in doc
    assert "hex_region" in doc
    assert "transition_kind` set to `\"pyramid\"" in doc
    assert "Slot267 extends the same bridge with surface-family intent" in doc
    assert 'required_surface_kinds=("quad", "triangle")' in doc
    assert 'expected_surface_kinds=["quad", "triangle"]' in doc
    assert "Slot275 extends the same bridge with material/block label intent" in doc
    assert "downstream_material_name" in doc
    assert "allowed_zero_downstream_material_names" in doc
    assert "Slot291 extends the bridge with Cubit meshing-scheme intent" in doc
    assert "shape_cubit_meshing_scheme_intent_gate" in doc
    assert "cubit_meshing_scheme_trace_gate" in doc
    assert "hex_region -> map" in doc
    assert "export netgen" in doc
    assert "Slot363 carries the Slot362 export-artifact rule into build123d" in doc
    assert "require_downstream_export_output_artifact=True" in doc
    assert "export_output_artifact_id" in doc
    assert "export_output_digest" in doc
    assert "export_output_path" in doc
    assert "cubit_meshing_scheme_handoff" in doc
    assert "old Cubit `.vol` digest" in doc
    assert "missing, extra, duplicate, or unnamed imported solids" in doc
    assert "keyed terminal plate" in doc
    assert '"keyed_terminal_plate_two_bosses", "volume": 9.364087557965556' in doc
    assert "flanged sleeve" in doc
    assert '"flanged_sleeve_four_bolt_holes", "volume": 5.085955762823052' in doc
    assert "coax_annular_sleeve_reference_row" in doc
    assert '"coax_annular_sleeve", "volume": 38.453094079939064' in doc
    assert "max_volume_rel_error = 3.693299710979559e-5" in doc
    assert "`1e-4` external-CAD volume gate" in doc
    assert "Ribbed busbar / heat-sink" in doc
    assert "ribbed_busbar_heat_sink_reference_row" in doc
    assert '"ribbed_busbar_heat_sink_four_holes", "volume": 13.617497357233166' in doc
    assert "RCD snubber heat-spreader" in doc
    assert "rcd_snubber_heat_spreader_reference_row" in doc
    assert '"rcd_snubber_heat_spreader", "volume": 15.52205629192717' in doc
    assert "RCD capacitance sweep design-table gate" in doc
    assert "rcd_snubber_capacitance_sweep_rows" in doc
    assert "parameter_key=\"capacitance_uF\"" in doc
    assert "a `0.114` volume mismatch" in doc
    assert "Thermal Robin cooling plate" in doc
    assert "thermal_robin_cooling_plate_reference_row" in doc
    assert '"thermal_robin_cooling_plate", "volume": 15.11290974078757' in doc
    assert "multi-line dict literals" in doc
    assert "V-type IPM rotor-coupon" in doc
    assert "v_type_ipm_rotor_coupon_reference_row" in doc
    assert '"v_type_ipm_rotor_coupon", "volume": 13.238339620676824' in doc
    assert "Cubit mass-property\nsidecar" in doc
    assert "cubit_mass_property_sidecar_gate" in doc
    assert "surface area\nand bbox dimensions" in doc


def test_annular_segment_volume_and_validity():
    seg = annular_segment(40, 55, 20, 0, 30, label="seg0")
    ana = (30/360)*math.pi*(55**2-40**2)*20
    assert abs(_vol(seg)-ana)/ana < 1e-6, f"wedge volume {_vol(seg)} vs {ana}"
    assert seg.is_valid and seg.label == "seg0"


def test_tube_volume():
    t = tube(40, 55, 20)
    assert abs(_vol(t) - math.pi*(55**2-40**2)*20)/_vol(t) < 1e-9 and t.is_valid


def test_racetrack_coil_builds():
    c = racetrack_coil(120, 80, 8, 15, 20, label="rc")
    assert c.is_valid and c.volume > 0 and c.label == "rc"


def test_polar_array_full_ring_sums_to_annulus():
    seg = annular_segment(40, 55, 20, 0, 30, label="seg")
    ring = polar_array(seg, 12, 360.0, label="halbach")
    assert len(ring.solids()) == 12 and _valid(ring)
    assert abs(_vol(ring) - math.pi*(55**2-40**2)*20)/_vol(ring) < 1e-6, "12 x 30deg = full ring"
    # build123d keeps region labels on the Compound's CHILDREN (.solids() flattens & drops them)
    assert [c.label for c in ring.children] == [f"halbach_{k:02d}" for k in range(12)]


def test_polar_array_partial_fan():
    fan = polar_array(Box(10, 4, 4).translate((30, 0, 0)), 4, total_angle=90.0)
    assert len(fan.solids()) == 4 and _valid(fan)


def test_linear_array_count_and_spacing():
    arr = linear_array(Box(10, 10, 10), 5, 15.0, (1, 0, 0), label="row")
    assert len(arr.solids()) == 5 and _valid(arr)
    assert abs(_vol(arr) - 5*1000.0) < 1e-6
    xs = sorted(s.center().X for s in arr.solids())
    assert abs((xs[-1]-xs[0]) - 4*15.0) < 1e-6, "5 copies span 4 spacings"


def test_mirrored_doubles_volume():
    half = Box(20, 10, 10).translate((0, 20, 0))
    both = mirrored(half, label="sym")
    assert len(both.solids()) == 2 and abs(_vol(both) - 2*2000.0) < 1e-6


def test_assembly_keeps_regions_separate():
    a = annular_segment(40, 55, 20, 0, 90, label="iron")
    b = tube(10, 20, 20, label="coil")
    asm = assembly(a, b, label="machine")
    assert len(asm.solids()) == 2 and asm.label == "machine"
    labels = {c.label for c in asm.children}
    assert labels == {"iron", "coil"}, "regions stay separate and labelled on the children (not fused)"


def test_assembly_accepts_raw_primitives():
    """A raw build123d Part (Box / Cylinder / extrude result) is a Compound with its solid in .solids()
    but EMPTY .children -- assembly must still keep it (fall back to .solids()), not silently drop it,
    and a label set on the raw Part must survive onto its child solid."""
    from build123d import Box, Cylinder, Pos
    core = Cylinder(0.5, 2); core.label = "core"
    asm = assembly(Box(1, 1, 1), Pos(3, 0, 0) * core, label="mixed")
    assert len(asm.solids()) == 2 and len(asm.children) == 2, "raw primitives are kept, not dropped"
    assert abs(asm.volume - (1.0 + math.pi * 0.5 ** 2 * 2)) < 1e-6
    assert "core" in {c.label for c in asm.children}, "a label on a raw Part carries onto its solid"


def test_shape_measurement_row_matches_box_geometry():
    box = Box(2, 3, 4).solid()
    box.label = "box"
    row = shape_measurement_row(box)

    assert row["name"] == "box"
    assert row["is_valid"]
    assert row["volume"] == pytest.approx(24.0)
    assert row["area"] == pytest.approx(52.0)
    assert row["faces"] == 6
    assert row["edges"] == 12
    assert row["vertices"] == 8
    assert row["solids"] == 1
    assert row["bounding_box"]["min"] == pytest.approx([-1.0, -1.5, -2.0])
    assert row["bounding_box"]["max"] == pytest.approx([1.0, 1.5, 2.0])
    assert row["bounding_box"]["center"] == pytest.approx([0.0, 0.0, 0.0])
    assert row["bounding_box"]["size"] == pytest.approx([2.0, 3.0, 4.0])
    assert row["bounding_box"]["diagonal"] == pytest.approx(math.sqrt(29.0))
    assert row["characteristic_length"] == pytest.approx(4.0)


def test_build123d_measurement_row_can_feed_cubit_mass_property_sidecar_gate():
    box = Box(1.5, 2.0, 0.75).solid()
    box.label = "hex_brick"
    row = shape_measurement_row(box)

    gate = cubit_mass_property_sidecar_gate(
        [row],
        expected_total_volume=2.25,
        expected_total_area=11.25,
        expected_bbox_size=[1.5, 2.0, 0.75],
        rel_tol=1.0e-12,
        abs_tol=1.0e-12,
    )

    assert gate["status"] == "ok"
    assert gate["row_names"] == ["hex_brick"]
    assert gate["total_volume"] == pytest.approx(2.25)
    assert gate["total_area"] == pytest.approx(11.25)
    assert gate["bbox_size"] == pytest.approx([1.5, 2.0, 0.75])


def test_build123d_cubit_export_package_handoff_matches_geometry_identity():
    box = Box(1.5, 2.0, 0.75).solid()
    box.label = "hex_brick"
    row = shape_measurement_row(box)
    row["geometry_id"] = "hex_brick_v1"
    vol_path = "artifacts/slot147_hex_brick_o3.vol"
    package = cubit_export_package_identity_gate(
        [
            {
                "kind": "vol",
                "path": vol_path,
                "export_id": "slot147_hex_brick_o3",
                "geometry_id": "hex_brick_v1",
                "order": 3,
            },
            {
                "kind": "vol_sidecar",
                "path": vol_path + ".json",
                "export_id": "slot147_hex_brick_o3",
                "geometry_id": "hex_brick_v1",
                "order": 3,
                "n_elements": 12,
                "n_points": 13,
                "vol_sidecar_schema_id": "coreform_netgen_vol_sidecar_inventory_v1",
            },
            {
                "kind": "raw_result",
                "path": "artifacts/slot147_hex_brick_raw.json",
                "export_id": "slot147_hex_brick_o3",
                "geometry_id": "hex_brick_v1",
            },
        ],
        expected_export_id="slot147_hex_brick_o3",
        expected_geometry_id="hex_brick_v1",
        expected_order=3,
        expected_routing_hint="cubit_hex_or_mixed_path",
        expected_vol_sidecar_schema_id="coreform_netgen_vol_sidecar_inventory_v1",
        require_vol_sidecar_schema=True,
        require_vol_sidecar_inventory_counts=True,
        inventory={
            "source": vol_path,
            "routing_hint": "cubit_hex_or_mixed_path",
            "volume_elements": 12,
            "points": 13,
        },
    )

    handoff = shape_cubit_export_package_handoff_gate(
        [row],
        package,
        expected_export_id="slot147_hex_brick_o3",
        require_sidecar_schema=True,
        require_sidecar_inventory_counts=True,
    )

    assert handoff["policy"] == "build123d_cubit_export_package_handoff_gate"
    assert handoff["status"] == "ok"
    assert handoff["checks"]["geometry_id_matches_package"] is True
    assert handoff["checks"]["package_vol_sidecar_pairs_vol"] is True
    assert handoff["checks"]["package_sidecar_schema_recorded"] is True
    assert handoff["checks"]["package_sidecar_schema_matches_expected"] is True
    assert handoff["checks"]["package_sidecar_inventory_counts_recorded"] is True
    assert handoff["checks"]["package_sidecar_element_count_matches_inventory"] is True
    assert handoff["checks"]["package_sidecar_point_count_matches_inventory"] is True
    assert handoff["checks"]["package_sidecar_order_matches_expected"] is True
    assert handoff["checks"]["shape_rows_have_volume_area_bbox"] is True

    wrong_row = dict(row)
    wrong_row["geometry_id"] = "hex_brick_old"
    bad_geometry = shape_cubit_export_package_handoff_gate([wrong_row], package)
    assert bad_geometry["status"] == "needs_attention"
    assert bad_geometry["checks"]["geometry_id_matches_package"] is False

    stale_package = dict(package)
    stale_package["status"] = "needs_attention"
    stale_package["checks"] = {**package["checks"], "vol_sidecar_pairs_vol": False}
    bad_package = shape_cubit_export_package_handoff_gate([row], stale_package)
    assert bad_package["status"] == "needs_attention"
    assert bad_package["checks"]["package_gate_ok"] is False
    assert bad_package["checks"]["package_vol_sidecar_pairs_vol"] is False

    stale_sidecar_package = dict(package)
    stale_sidecar_package["checks"] = {
        **package["checks"],
        "vol_sidecar_element_count_matches_inventory": False,
    }
    bad_sidecar_counts = shape_cubit_export_package_handoff_gate(
        [row],
        stale_sidecar_package,
        require_sidecar_inventory_counts=True,
    )
    assert bad_sidecar_counts["status"] == "needs_attention"
    assert bad_sidecar_counts["checks"]["package_sidecar_element_count_matches_inventory"] is False

    stale_schema_package = dict(package)
    stale_schema_package["checks"] = {
        **package["checks"],
        "expected_vol_sidecar_schema_id_matches": False,
    }
    bad_sidecar_schema = shape_cubit_export_package_handoff_gate(
        [row],
        stale_schema_package,
        require_sidecar_schema=True,
    )
    assert bad_sidecar_schema["status"] == "needs_attention"
    assert bad_sidecar_schema["checks"]["package_sidecar_schema_matches_expected"] is False


def test_build123d_cubit_quality_package_handoff_matches_geometry_identity():
    box = Box(1.0, 1.0, 1.0).solid()
    box.label = "unit_brick"
    row = shape_measurement_row(box)
    row["geometry_id"] = "unit_brick_mapped_hex_v1"
    package = {
        "policy": "cubit_headless_batch_quality_package_gate",
        "status": "ok",
        "export_id": "slot154_headless_hex_quality_A",
        "geometry_id": "unit_brick_mapped_hex_v1",
        "quality_count": 64,
        "export_inventory_source": r"artifacts/cubit/slot154_hex.vol",
        "export_inventory_volume_kind_counts": {"hex": 64},
        "export_inventory_routing_hint": "cubit_hex_or_mixed_path",
        "export_inventory_is_tri_tet_only": False,
        "checks": {
            "headless_command_recorded": True,
            "quality_count_positive": True,
            "export_inventory_recorded": True,
            "export_inventory_volume_elements_positive": True,
            "export_inventory_routing_hint_matches_expected": True,
            "export_inventory_count_matches_quality": True,
            "export_inventory_contains_quality_element": True,
            "export_inventory_not_tri_tet_only_for_cubit_hex_route": True,
        },
    }

    handoff = shape_cubit_quality_package_handoff_gate(
        [row],
        package,
        expected_export_id="slot154_headless_hex_quality_A",
        require_export_inventory=True,
    )

    assert handoff["policy"] == "build123d_cubit_quality_package_handoff_gate"
    assert handoff["status"] == "ok"
    assert handoff["checks"]["geometry_id_matches_quality_package"] is True
    assert handoff["checks"]["quality_package_headless"] is True
    assert handoff["checks"]["quality_package_count_positive"] is True
    assert handoff["checks"]["quality_package_export_inventory_present"] is True
    assert handoff["checks"]["quality_package_export_inventory_count_matches"] is True
    assert handoff["checks"]["quality_package_export_inventory_contains_quality_element"] is True
    assert handoff["checks"]["quality_package_export_inventory_not_tri_tet_only_for_cubit_route"] is True
    assert handoff["quality_package_export_inventory_volume_kind_counts"] == {"hex": 64}
    assert handoff["quality_package_export_inventory_is_tri_tet_only"] is False

    wrong_row = dict(row)
    wrong_row["geometry_id"] = "unit_brick_old"
    bad_geometry = shape_cubit_quality_package_handoff_gate([wrong_row], package)
    assert bad_geometry["status"] == "needs_attention"
    assert bad_geometry["checks"]["geometry_id_matches_quality_package"] is False

    gui_package = dict(package)
    gui_package["checks"] = {**package["checks"], "headless_command_recorded": False}
    bad_gui = shape_cubit_quality_package_handoff_gate([row], gui_package)
    assert bad_gui["status"] == "needs_attention"
    assert bad_gui["checks"]["quality_package_headless"] is False

    zero_count = dict(package)
    zero_count["quality_count"] = 0
    zero_count["checks"] = {**package["checks"], "quality_count_positive": False}
    bad_count = shape_cubit_quality_package_handoff_gate([row], zero_count)
    assert bad_count["status"] == "needs_attention"
    assert bad_count["checks"]["quality_package_count_positive"] is False

    missing_inventory = dict(package)
    missing_inventory["checks"] = {
        key: value for key, value in package["checks"].items()
        if not key.startswith("export_inventory")
    }
    bad_inventory = shape_cubit_quality_package_handoff_gate(
        [row],
        missing_inventory,
        require_export_inventory=True,
    )
    assert bad_inventory["status"] == "needs_attention"
    assert bad_inventory["checks"]["quality_package_export_inventory_present"] is False

    tri_tet_inventory = dict(package)
    tri_tet_inventory["export_inventory_volume_kind_counts"] = {"tet": 64}
    tri_tet_inventory["export_inventory_routing_hint"] = "netgen_tri_tet_path"
    tri_tet_inventory["export_inventory_is_tri_tet_only"] = True
    tri_tet_inventory["checks"] = {
        **package["checks"],
        "export_inventory_routing_hint_matches_expected": False,
        "export_inventory_contains_quality_element": False,
        "export_inventory_not_tri_tet_only_for_cubit_hex_route": False,
    }
    bad_route = shape_cubit_quality_package_handoff_gate(
        [row],
        tri_tet_inventory,
        expected_export_id="slot154_headless_hex_quality_A",
        require_export_inventory=True,
    )
    assert bad_route["status"] == "needs_attention"
    assert bad_route["checks"]["quality_package_export_inventory_routing_ok"] is False
    assert bad_route["checks"]["quality_package_export_inventory_contains_quality_element"] is False
    assert bad_route["checks"]["quality_package_export_inventory_not_tri_tet_only_for_cubit_route"] is False
    assert bad_route["quality_package_export_inventory_is_tri_tet_only"] is True


def test_build123d_cubit_quality_ledger_handoff_binds_mesh_digest_and_route():
    box = Box(1.0, 1.0, 1.0).solid()
    box.label = "unit_brick_quality_ledger"
    row = shape_measurement_row(box)
    row["geometry_id"] = "unit_brick_quality_ledger_v1"
    row["mesh_route"] = "cubit_hex_or_mixed_path"
    ledger = {
        "policy": "cubit_mesh_quality_ledger_identity_gate",
        "status": "ok",
        "quality_artifact_id": "slot369_hex_quality_ledger_json_v1",
        "quality_digest": "sha256:slot369-quality-ledger",
        "metric_set_id": "cubit_scaled_jacobian_hex_v1",
        "export_id": "slot369_hex_quality_ledger",
        "geometry_id": "unit_brick_quality_ledger_v1",
        "mesh_artifact_id": "slot369_hex_quality_ledger_vol_v1",
        "mesh_digest": "sha256:slot369-hex-quality-vol",
        "routing_hint": "cubit_hex_or_mixed_path",
        "min_scaled_jacobian": 1.0,
        "negative_jacobian_count": 0,
        "element_type_counts": {"hex": 64},
        "inventory_is_tri_tet_only": False,
        "checks": {
            "quality_artifact_id_recorded": True,
            "quality_digest_recorded": True,
            "metric_set_id_recorded": True,
            "mesh_artifact_id_recorded": True,
            "mesh_digest_recorded": True,
            "min_scaled_jacobian_above_threshold": True,
            "negative_jacobian_count_zero": True,
            "hex_or_mixed_volume_family_present": True,
            "not_tri_tet_only_for_cubit_quality_ledger": True,
            "created_at_utc_recorded_when_required": True,
            "created_at_utc_parseable_when_present": True,
            "version_recorded_when_required": True,
            "elapsed_s_recorded_when_required": True,
            "elapsed_s_finite_nonnegative_when_present": True,
            "timing_breakdown_recorded_when_required": True,
            "timing_breakdown_has_required_stage_count": True,
            "timing_breakdown_values_finite_nonnegative": True,
            "timing_breakdown_total_within_elapsed_when_present": True,
        },
    }

    handoff = shape_cubit_quality_ledger_handoff_gate(
        [row],
        ledger,
        expected_quality_artifact_id="slot369_hex_quality_ledger_json_v1",
        expected_quality_digest="sha256:slot369-quality-ledger",
        expected_metric_set_id="cubit_scaled_jacobian_hex_v1",
        expected_export_id="slot369_hex_quality_ledger",
        expected_mesh_artifact_id="slot369_hex_quality_ledger_vol_v1",
        expected_mesh_digest="sha256:slot369-hex-quality-vol",
        require_quality_execution_metadata=True,
    )

    assert handoff["policy"] == "build123d_cubit_quality_ledger_handoff_gate"
    assert handoff["status"] == "ok"
    assert handoff["checks"]["geometry_id_matches_quality_ledger"] is True
    assert handoff["checks"]["expected_quality_digest_matches"] is True
    assert handoff["checks"]["expected_mesh_digest_matches"] is True
    assert handoff["checks"]["mesh_route_matches_expected"] is True
    assert handoff["checks"]["routing_hint_matches_expected"] is True
    assert handoff["checks"]["hex_or_mixed_volume_family_present"] is True
    assert handoff["checks"]["not_tri_tet_only_for_cubit_quality_ledger"] is True
    assert handoff["checks"]["quality_ledger_execution_metadata_ok"] is True
    assert handoff["require_quality_execution_metadata"] is True
    assert handoff["quality_digest"] == "sha256:slot369-quality-ledger"
    assert handoff["mesh_digest"] == "sha256:slot369-hex-quality-vol"

    wrong_row = dict(row)
    wrong_row["geometry_id"] = "unit_brick_old"
    bad_geometry = shape_cubit_quality_ledger_handoff_gate([wrong_row], ledger)
    assert bad_geometry["status"] == "needs_attention"
    assert bad_geometry["checks"]["geometry_id_matches_quality_ledger"] is False

    stale_ledger = dict(ledger)
    stale_ledger["status"] = "needs_attention"
    stale_ledger["quality_digest"] = "sha256:old-quality-ledger"
    stale_ledger["checks"] = {**ledger["checks"], "expected_quality_digest_matches": False}
    bad_digest = shape_cubit_quality_ledger_handoff_gate(
        [row],
        stale_ledger,
        expected_quality_digest="sha256:slot369-quality-ledger",
    )
    assert bad_digest["status"] == "needs_attention"
    assert bad_digest["checks"]["quality_ledger_gate_ok"] is False
    assert bad_digest["checks"]["expected_quality_digest_matches"] is False

    tri_tet_ledger = {
        **ledger,
        "routing_hint": "netgen_tri_tet_path",
        "element_type_counts": {"tet": 64},
        "inventory_is_tri_tet_only": True,
        "checks": {
            **ledger["checks"],
            "hex_or_mixed_volume_family_present": False,
            "not_tri_tet_only_for_cubit_quality_ledger": False,
        },
    }
    bad_route = shape_cubit_quality_ledger_handoff_gate([row], tri_tet_ledger)
    assert bad_route["status"] == "needs_attention"
    assert bad_route["checks"]["routing_hint_matches_expected"] is False
    assert bad_route["checks"]["hex_or_mixed_volume_family_present"] is False
    assert bad_route["checks"]["not_tri_tet_only_for_cubit_quality_ledger"] is False

    weak_mesh = {**ledger, "min_scaled_jacobian": 0.05, "negative_jacobian_count": 2}
    bad_quality = shape_cubit_quality_ledger_handoff_gate([row], weak_mesh)
    assert bad_quality["status"] == "needs_attention"
    assert bad_quality["checks"]["min_scaled_jacobian_above_threshold"] is False
    assert bad_quality["checks"]["negative_jacobian_count_zero"] is False

    missing_execution = {
        **ledger,
        "checks": {
            **ledger["checks"],
            "timing_breakdown_has_required_stage_count": False,
        },
    }
    bad_execution = shape_cubit_quality_ledger_handoff_gate(
        [row],
        missing_execution,
        require_quality_execution_metadata=True,
    )
    assert bad_execution["status"] == "needs_attention"
    assert bad_execution["checks"]["quality_ledger_execution_metadata_ok"] is False


def test_build123d_cad_handoff_manifest_bundles_volume_files_and_quality():
    box = Box(1.0, 2.0, 3.0).solid()
    box.label = "terminal_block"
    row = shape_measurement_row(box)
    row["geometry_id"] = "terminal_block_v1"
    row["mesh_route"] = "cubit_hex_or_mixed_path"
    row["role"] = "hex_region"
    row["expected_cubit_scheme"] = "map"
    row["downstream_meshing_trace_id"] = "slot363_terminal_block_scheme_trace"
    row["expected_cubit_command_fragments"] = ["create brick", "volume 1 scheme map", "export netgen"]
    row["expected_cubit_export_order"] = 2
    row["units"] = {"length": "mm", "area": "mm^2", "volume": "mm^3"}
    row["cad_measurement_convention"] = "occt_closed_solid_mass_properties"
    row["cad_measurement_postprocess_row_convention_schema_id"] = (
        "build123d_occt_mass_property_row_convention_v1"
    )
    row["cad_measurement_component_basis_schema_id"] = (
        "build123d_occt_volume_area_bbox_component_basis_v1"
    )
    external_volume = shape_volume_crosscheck_summary(
        [row],
        {
            "cubit": [{"name": "terminal_block", "volume": row["volume"]}],
            "cst": [{"name": "terminal_block", "volume": row["volume"]}],
        },
        rtol=1.0e-12,
    )
    quality_package = {
        "policy": "cubit_headless_batch_quality_package_gate",
        "status": "ok",
        "export_id": "slot163_terminal_block_hex_quality_A",
        "geometry_id": "terminal_block_v1",
        "quality_count": 27,
        "checks": {
            "headless_command_recorded": True,
            "quality_count_positive": True,
        },
    }
    quality_handoff = shape_cubit_quality_package_handoff_gate(
        [row],
        quality_package,
        expected_export_id="slot163_terminal_block_hex_quality_A",
    )
    quality_ledger_gate = {
        "policy": "cubit_mesh_quality_ledger_identity_gate",
        "status": "ok",
        "quality_artifact_id": "slot369_terminal_block_quality_ledger_v1",
        "quality_digest": "sha256:slot369-terminal-block-quality-ledger",
        "metric_set_id": "cubit_scaled_jacobian_hex_v1",
        "export_id": "slot163_terminal_block_hex_quality_A",
        "geometry_id": "terminal_block_v1",
        "mesh_artifact_id": "slot363_terminal_block_vol_v1",
        "mesh_digest": "sha256:slot363-terminal-block-vol",
        "routing_hint": "cubit_hex_or_mixed_path",
        "min_scaled_jacobian": 1.0,
        "negative_jacobian_count": 0,
        "element_type_counts": {"hex": 27},
        "inventory_is_tri_tet_only": False,
        "checks": {
            "quality_artifact_id_recorded": True,
            "quality_digest_recorded": True,
            "metric_set_id_recorded": True,
            "mesh_artifact_id_recorded": True,
            "mesh_digest_recorded": True,
            "min_scaled_jacobian_above_threshold": True,
            "negative_jacobian_count_zero": True,
            "hex_or_mixed_volume_family_present": True,
            "not_tri_tet_only_for_cubit_quality_ledger": True,
        },
    }
    quality_ledger_handoff = shape_cubit_quality_ledger_handoff_gate(
        [row],
        quality_ledger_gate,
        expected_quality_artifact_id="slot369_terminal_block_quality_ledger_v1",
        expected_quality_digest="sha256:slot369-terminal-block-quality-ledger",
        expected_metric_set_id="cubit_scaled_jacobian_hex_v1",
        expected_export_id="slot163_terminal_block_hex_quality_A",
        expected_mesh_artifact_id="slot363_terminal_block_vol_v1",
        expected_mesh_digest="sha256:slot363-terminal-block-vol",
    )
    solver_route_gate = cubit_mixed_solver_route_manifest_gate(
        {
            "volume_kind_counts": {"hex": 1, "pyramid": 1, "tet": 1},
            "surface_kind_counts": {"quad": 1, "triangle": 1},
            "routing_hint": "cubit_hex_or_mixed_path",
        },
        {
            "solver_route_package_id": "slot347_terminal_block_solver_route_v1",
            "routing_hint": "cubit_hex_or_mixed_path",
            "route_policy": "hex_primary_pyramid_transition_tet_compatibility",
            "downstream_solver": "NGSolve/radia-ngsolve",
            "downstream_solver_contract_artifact_id": "slot384_terminal_block_reader_contract_v1",
            "downstream_solver_contract_digest": "sha256:slot384-terminal-block-reader-contract",
            "downstream_solver_contract_path": r"artifacts/cubit/slot384_terminal_block_reader_contract.json",
            "solver_route_convention_schema_id": "coreform_mixed_hex_pyramid_tet_route_convention_v1",
            "tet_only_owner": "netgen_tri_tet_path",
            "no_implicit_tetization": True,
            "volume_routes": [
                {"volume_kind": "hex", "solver_role": "primary_volume_fem"},
                {"volume_kind": "pyramid", "solver_role": "transition_bridge", "not_primary_region": True},
                {"volume_kind": "tet", "solver_role": "compatibility_subregion_volume_fem"},
            ],
            "surface_routes": [
                {"surface_kind": "quad", "solver_role": "hex_boundary_trace"},
                {"surface_kind": "triangle", "solver_role": "tet_boundary_trace"},
            ],
        },
        expected_package_id="slot347_terminal_block_solver_route_v1",
        expected_solver_contract_artifact_id="slot384_terminal_block_reader_contract_v1",
        expected_solver_contract_digest="sha256:slot384-terminal-block-reader-contract",
        expected_solver_contract_path=r"artifacts/cubit/slot384_terminal_block_reader_contract.json",
        expected_solver_route_convention_schema_id="coreform_mixed_hex_pyramid_tet_route_convention_v1",
        require_solver_contract_artifact=True,
        require_solver_route_convention_schema=True,
    )
    solver_route_handoff = shape_cubit_solver_route_handoff_gate(
        [row],
        solver_route_gate,
        expected_solver_route_package_id="slot347_terminal_block_solver_route_v1",
        expected_solver_contract_artifact_id="slot384_terminal_block_reader_contract_v1",
        expected_solver_contract_digest="sha256:slot384-terminal-block-reader-contract",
        expected_solver_contract_path=r"artifacts/cubit/slot384_terminal_block_reader_contract.json",
        expected_solver_route_convention_schema_id="coreform_mixed_hex_pyramid_tet_route_convention_v1",
        require_solver_contract_artifact=True,
        require_solver_route_convention_schema=True,
    )
    scheme_trace = cubit_meshing_scheme_trace_gate(
        {
            "trace_id": "slot363_terminal_block_scheme_trace",
            "command_digest": "sha256:slot363-terminal-block-map-export",
            "commands": [
                "create brick x 1 y 2 z 3",
                "volume 1 scheme map",
                'export netgen "artifacts/cubit/slot363_terminal_block.vol" order 2 overwrite',
            ],
            "volume_schemes": {"1": "map"},
            "export_order": 2,
            "export_output_artifact_id": "slot363_terminal_block_vol_v1",
            "export_output_digest": "sha256:slot363-terminal-block-vol",
            "export_output_path": r"artifacts/cubit/slot363_terminal_block.vol",
        },
        expected_trace_id="slot363_terminal_block_scheme_trace",
        expected_command_digest="sha256:slot363-terminal-block-map-export",
        expected_volume_schemes={"1": "map"},
        required_command_fragments=("create brick", "volume 1 scheme map", "export netgen"),
        expected_export_order=2,
        expected_export_output_artifact_id="slot363_terminal_block_vol_v1",
        expected_export_output_digest="sha256:slot363-terminal-block-vol",
        expected_export_output_path=r"artifacts/cubit/slot363_terminal_block.vol",
        require_export_output_artifact=True,
    )
    scheme_handoff = shape_cubit_meshing_scheme_intent_gate(
        [row],
        scheme_trace_gate=scheme_trace,
        required_roles=("hex_region",),
        expected_scheme_by_role={"hex_region": "map"},
        required_command_fragments=("create brick", "volume 1 scheme map", "export netgen"),
        expected_trace_id="slot363_terminal_block_scheme_trace",
        expected_export_order=2,
        expected_export_output_artifact_id="slot363_terminal_block_vol_v1",
        expected_export_output_digest="sha256:slot363-terminal-block-vol",
        expected_export_output_path=r"artifacts/cubit/slot363_terminal_block.vol",
        require_downstream_export_output_artifact=True,
    )
    manifest = [
        {
            "kind": "step",
            "path": r"artifacts/build123d/terminal_block.step",
            "cad_output_artifact_id": "terminal_block_cad_output_step_v1",
            "cad_output_digest": "sha256:terminal_block_cad_output_step_v1",
            "cad_output_path": r"artifacts/build123d/terminal_block.step",
        },
        {
            "kind": "build123d_measurement_json",
            "path": r"artifacts/build123d/terminal_block_measure.json",
            "cad_output_artifact_id": "terminal_block_cad_output_step_v1",
            "cad_output_digest": "sha256:terminal_block_cad_output_step_v1",
            "cad_output_path": r"artifacts/build123d/terminal_block.step",
        },
        {
            "kind": "external_volume_summary_json",
            "path": r"artifacts/build123d/terminal_block_volume.json",
            "cad_output_artifact_id": "terminal_block_cad_output_step_v1",
            "cad_output_digest": "sha256:terminal_block_cad_output_step_v1",
            "cad_output_path": r"artifacts/build123d/terminal_block.step",
        },
        {
            "kind": "cubit_quality_json",
            "path": r"artifacts/cubit/terminal_block_quality.json",
            "cad_output_artifact_id": "terminal_block_cad_output_step_v1",
            "cad_output_digest": "sha256:terminal_block_cad_output_step_v1",
            "cad_output_path": r"artifacts/build123d/terminal_block.step",
        },
        {
            "kind": "cubit_quality_ledger_json",
            "path": r"artifacts/cubit/terminal_block_quality_ledger.json",
            "cad_output_artifact_id": "terminal_block_cad_output_step_v1",
            "cad_output_digest": "sha256:terminal_block_cad_output_step_v1",
            "cad_output_path": r"artifacts/build123d/terminal_block.step",
        },
        {
            "kind": "cubit_solver_route_json",
            "path": r"artifacts/cubit/terminal_block_solver_route.json",
            "cad_output_artifact_id": "terminal_block_cad_output_step_v1",
            "cad_output_digest": "sha256:terminal_block_cad_output_step_v1",
            "cad_output_path": r"artifacts/build123d/terminal_block.step",
        },
    ]
    for item in manifest:
        item["cad_observable_id"] = "terminal_block_cad_handoff_measurements_v1"
        item["cad_observable_family"] = "cad_mass_properties_handoff"
        item["length_unit"] = "mm"
        item["area_unit"] = "mm^2"
        item["volume_unit"] = "mm^3"
        item["cad_measurement_convention"] = "occt_closed_solid_mass_properties"
        item["cad_measurement_postprocess_row_convention_schema_id"] = (
            "build123d_occt_mass_property_row_convention_v1"
        )
        item["cad_measurement_component_basis_schema_id"] = (
            "build123d_occt_volume_area_bbox_component_basis_v1"
        )

    gate = shape_cad_handoff_manifest_gate(
        [row],
        file_manifest=manifest,
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        cubit_quality_ledger_handoff=quality_ledger_handoff,
        cubit_solver_route_handoff=solver_route_handoff,
        cubit_meshing_scheme_handoff=scheme_handoff,
        required_file_kinds=(
            "step",
            "build123d_measurement_json",
            "external_volume_summary_json",
            "cubit_quality_json",
            "cubit_quality_ledger_json",
        ),
        expected_geometry_ids=("terminal_block_v1",),
        expected_cad_output_artifact_id="terminal_block_cad_output_step_v1",
        expected_cad_output_digest="sha256:terminal_block_cad_output_step_v1",
        require_cad_output_artifact=True,
        expected_cad_observable_id="terminal_block_cad_handoff_measurements_v1",
        expected_cad_observable_family="cad_mass_properties_handoff",
        require_cad_observable=True,
        expected_length_unit="mm",
        expected_area_unit="mm^2",
        expected_volume_unit="mm^3",
        expected_measurement_convention="occt_closed_solid_mass_properties",
        expected_measurement_postprocess_row_convention_schema_id=(
            "build123d_occt_mass_property_row_convention_v1"
        ),
        expected_measurement_component_basis_schema_id=(
            "build123d_occt_volume_area_bbox_component_basis_v1"
        ),
        require_measurement_postprocess_row_convention_schema=True,
        require_measurement_component_basis_schema=True,
    )

    assert gate["policy"] == "build123d_cad_handoff_manifest_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["external_volume_summary_ok"] is True
    assert gate["checks"]["cubit_quality_handoff_ok"] is True
    assert gate["checks"]["cubit_quality_ledger_handoff_ok"] is True
    assert gate["checks"]["cubit_solver_route_handoff_ok"] is True
    assert gate["checks"]["cubit_meshing_scheme_handoff_ok"] is True
    assert solver_route_handoff["checks"]["solver_contract_artifact_id_recorded_when_required"] is True
    assert solver_route_handoff["checks"]["solver_contract_digest_recorded_when_required"] is True
    assert solver_route_handoff["checks"]["solver_contract_path_recorded_when_required"] is True
    assert solver_route_handoff["checks"]["expected_solver_contract_artifact_id_matches"] is True
    assert solver_route_handoff["checks"]["expected_solver_contract_digest_matches"] is True
    assert solver_route_handoff["checks"]["expected_solver_contract_path_matches"] is True
    assert solver_route_handoff["checks"]["solver_route_convention_schema_id_recorded_when_required"] is True
    assert solver_route_handoff["checks"]["expected_solver_route_convention_schema_id_matches"] is True
    assert solver_route_handoff["solver_contract_artifact_id"] == "slot384_terminal_block_reader_contract_v1"
    assert solver_route_handoff["solver_contract_digest"] == "sha256:slot384-terminal-block-reader-contract"
    assert solver_route_handoff["solver_contract_path"] == r"artifacts/cubit/slot384_terminal_block_reader_contract.json"
    assert solver_route_handoff["solver_route_convention_schema_id"] == "coreform_mixed_hex_pyramid_tet_route_convention_v1"
    assert gate["checks"]["required_file_kinds_present"] is True
    assert gate["checks"]["cad_output_artifact_id_recorded_when_required"] is True
    assert gate["checks"]["cad_output_digest_recorded_when_required"] is True
    assert gate["checks"]["cad_output_path_recorded_when_required"] is True
    assert gate["checks"]["expected_cad_output_artifact_id_matches"] is True
    assert gate["checks"]["expected_cad_output_digest_matches"] is True
    assert gate["checks"]["cad_observable_id_recorded_when_required"] is True
    assert gate["checks"]["expected_cad_observable_id_matches"] is True
    assert gate["checks"]["expected_cad_observable_family_matches"] is True
    assert gate["checks"]["cad_unit_metadata_consistent_when_present"] is True
    assert gate["checks"]["expected_cad_length_unit_matches"] is True
    assert gate["checks"]["expected_cad_area_unit_matches"] is True
    assert gate["checks"]["expected_cad_volume_unit_matches"] is True
    assert gate["checks"]["cad_measurement_convention_consistent_when_present"] is True
    assert gate["checks"]["expected_cad_measurement_convention_matches"] is True
    assert gate["checks"]["cad_measurement_postprocess_row_convention_schema_id_consistent_when_present"] is True
    assert gate["checks"]["cad_measurement_postprocess_row_convention_schema_id_recorded_when_required"] is True
    assert gate["checks"]["cad_measurement_postprocess_row_convention_schema_id_recorded_when_expected"] is True
    assert (
        gate["checks"][
            "expected_cad_measurement_postprocess_row_convention_schema_id_matches"
        ]
        is True
    )
    assert gate["checks"]["cad_measurement_component_basis_schema_id_consistent_when_present"] is True
    assert gate["checks"]["cad_measurement_component_basis_schema_id_recorded_when_required"] is True
    assert gate["checks"]["cad_measurement_component_basis_schema_id_recorded_when_expected"] is True
    assert (
        gate["checks"]["expected_cad_measurement_component_basis_schema_id_matches"]
        is True
    )
    assert gate["cad_output_artifact_id"] == "terminal_block_cad_output_step_v1"
    assert gate["cad_output_digest"] == "sha256:terminal_block_cad_output_step_v1"
    assert gate["cad_output_path"] == r"artifacts/build123d/terminal_block.step"
    assert gate["cad_observable_id"] == "terminal_block_cad_handoff_measurements_v1"
    assert gate["cad_observable_family"] == "cad_mass_properties_handoff"
    assert gate["units"] == {"length": ["mm"], "area": ["mm^2"], "volume": ["mm^3"]}
    assert gate["cad_measurement_convention"] == "occt_closed_solid_mass_properties"
    assert gate["cad_measurement_postprocess_row_convention_schema_id"] == (
        "build123d_occt_mass_property_row_convention_v1"
    )
    assert gate["cad_measurement_postprocess_row_convention_schema_ids"] == [
        "build123d_occt_mass_property_row_convention_v1"
    ]
    assert gate["require_measurement_postprocess_row_convention_schema"] is True
    assert gate["cad_measurement_component_basis_schema_id"] == (
        "build123d_occt_volume_area_bbox_component_basis_v1"
    )
    assert gate["cad_measurement_component_basis_schema_ids"] == [
        "build123d_occt_volume_area_bbox_component_basis_v1"
    ]
    assert gate["require_measurement_component_basis_schema"] is True
    assert gate["external_volume_sources"] == ["cubit", "cst"]
    assert gate["cubit_quality_ledger_handoff_policy"] == "build123d_cubit_quality_ledger_handoff_gate"
    assert gate["cubit_solver_route_handoff_policy"] == "build123d_cubit_solver_route_handoff_gate"
    assert gate["cubit_meshing_scheme_handoff_policy"] == "build123d_cubit_meshing_scheme_intent_gate"
    assert scheme_handoff["checks"]["downstream_export_output_artifact_id_recorded_when_required"] is True
    assert scheme_handoff["checks"]["downstream_export_output_digest_matches"] is True
    assert scheme_handoff["checks"]["downstream_export_output_path_matches"] is True
    assert solver_route_handoff["checks"]["solver_route_pyramid_transition_role_recorded"] is True
    assert solver_route_handoff["checks"]["solver_route_no_implicit_tetization"] is True

    stale_volume = {**external_volume, "status": "needs_attention"}
    bad_volume = shape_cad_handoff_manifest_gate(
        [row],
        file_manifest=manifest,
        external_volume_summary=stale_volume,
        cubit_quality_handoff=quality_handoff,
        required_file_kinds=("step", "build123d_measurement_json", "external_volume_summary_json", "cubit_quality_json"),
    )
    assert bad_volume["status"] == "needs_attention"
    assert bad_volume["checks"]["external_volume_summary_ok"] is False

    missing_file = shape_cad_handoff_manifest_gate(
        [row],
        file_manifest=manifest[:2],
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        required_file_kinds=("step", "build123d_measurement_json", "external_volume_summary_json", "cubit_quality_json"),
    )
    assert missing_file["status"] == "needs_attention"
    assert missing_file["checks"]["required_file_kinds_present"] is False

    stale_output = [dict(item) for item in manifest]
    stale_output[1]["cad_output_artifact_id"] = "terminal_block_old_step"
    stale_output_gate = shape_cad_handoff_manifest_gate(
        [row],
        file_manifest=stale_output,
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        required_file_kinds=("step", "build123d_measurement_json", "external_volume_summary_json", "cubit_quality_json"),
        expected_cad_output_artifact_id="terminal_block_cad_output_step_v1",
        require_cad_output_artifact=True,
    )
    assert stale_output_gate["status"] == "needs_attention"
    assert stale_output_gate["checks"]["cad_output_artifact_id_consistent_when_present"] is False
    assert stale_output_gate["checks"]["expected_cad_output_artifact_id_matches"] is True

    stale_digest = [dict(item) for item in manifest]
    stale_digest[1]["cad_output_digest"] = "sha256:terminal_block_old_step"
    stale_digest_gate = shape_cad_handoff_manifest_gate(
        [row],
        file_manifest=stale_digest,
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        required_file_kinds=("step", "build123d_measurement_json", "external_volume_summary_json", "cubit_quality_json"),
        expected_cad_output_digest="sha256:terminal_block_cad_output_step_v1",
        require_cad_output_artifact=True,
    )
    assert stale_digest_gate["status"] == "needs_attention"
    assert stale_digest_gate["checks"]["cad_output_digest_consistent_when_present"] is False
    assert stale_digest_gate["checks"]["expected_cad_output_digest_matches"] is True

    missing_output_path = [dict(item) for item in manifest]
    for item in missing_output_path:
        item.pop("cad_output_path")
    missing_output_path_gate = shape_cad_handoff_manifest_gate(
        [row],
        file_manifest=missing_output_path,
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        required_file_kinds=("step", "build123d_measurement_json", "external_volume_summary_json", "cubit_quality_json"),
        require_cad_output_artifact=True,
    )
    assert missing_output_path_gate["status"] == "needs_attention"
    assert missing_output_path_gate["checks"]["cad_output_path_recorded_when_required"] is False

    stale_observable = [dict(item) for item in manifest]
    stale_observable[1]["cad_observable_id"] = "aa_mesh_inventory_v1"
    stale_observable_gate = shape_cad_handoff_manifest_gate(
        [row],
        file_manifest=stale_observable,
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        required_file_kinds=("step", "build123d_measurement_json", "external_volume_summary_json", "cubit_quality_json"),
        expected_cad_observable_id="terminal_block_cad_handoff_measurements_v1",
        require_cad_observable=True,
    )
    assert stale_observable_gate["status"] == "needs_attention"
    assert stale_observable_gate["checks"]["cad_observable_id_consistent_when_present"] is False
    assert stale_observable_gate["checks"]["expected_cad_observable_id_matches"] is False

    stale_observable_family = [dict(item) for item in manifest]
    stale_observable_family[1]["cad_observable_family"] = "aa_mesh_inventory"
    stale_observable_family_gate = shape_cad_handoff_manifest_gate(
        [row],
        file_manifest=stale_observable_family,
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        required_file_kinds=("step", "build123d_measurement_json", "external_volume_summary_json", "cubit_quality_json"),
        expected_cad_observable_family="cad_mass_properties_handoff",
        require_cad_observable=True,
    )
    assert stale_observable_family_gate["status"] == "needs_attention"
    assert stale_observable_family_gate["checks"]["cad_observable_family_consistent_when_present"] is False
    assert stale_observable_family_gate["checks"]["expected_cad_observable_family_matches"] is False

    stale_volume_unit = [dict(item) for item in manifest]
    stale_volume_unit[1]["volume_unit"] = "m^3"
    stale_volume_unit_gate = shape_cad_handoff_manifest_gate(
        [row],
        file_manifest=stale_volume_unit,
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        required_file_kinds=("step", "build123d_measurement_json", "external_volume_summary_json", "cubit_quality_json"),
        expected_volume_unit="mm^3",
    )
    assert stale_volume_unit_gate["status"] == "needs_attention"
    assert stale_volume_unit_gate["checks"]["cad_unit_metadata_consistent_when_present"] is False
    assert stale_volume_unit_gate["checks"]["expected_cad_volume_unit_matches"] is False

    wrong_convention = [dict(item) for item in manifest]
    wrong_convention[2]["cad_measurement_convention"] = "mesh_volume_after_import"
    wrong_convention_gate = shape_cad_handoff_manifest_gate(
        [row],
        file_manifest=wrong_convention,
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        required_file_kinds=("step", "build123d_measurement_json", "external_volume_summary_json", "cubit_quality_json"),
        expected_measurement_convention="occt_closed_solid_mass_properties",
    )
    assert wrong_convention_gate["status"] == "needs_attention"
    assert wrong_convention_gate["checks"]["cad_measurement_convention_consistent_when_present"] is False
    assert wrong_convention_gate["checks"]["expected_cad_measurement_convention_matches"] is False

    stale_measurement_row_convention = [dict(item) for item in manifest]
    stale_measurement_row_convention[2]["cad_measurement_postprocess_row_convention_schema_id"] = (
        "build123d_scalar_volume_row_v0"
    )
    stale_measurement_row_convention_gate = shape_cad_handoff_manifest_gate(
        [row],
        file_manifest=stale_measurement_row_convention,
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        required_file_kinds=("step", "build123d_measurement_json", "external_volume_summary_json", "cubit_quality_json"),
        expected_measurement_convention="occt_closed_solid_mass_properties",
        expected_measurement_postprocess_row_convention_schema_id=(
            "build123d_occt_mass_property_row_convention_v1"
        ),
        require_measurement_postprocess_row_convention_schema=True,
    )
    assert stale_measurement_row_convention_gate["status"] == "needs_attention"
    assert stale_measurement_row_convention_gate["checks"]["expected_cad_measurement_convention_matches"] is True
    assert (
        stale_measurement_row_convention_gate["checks"][
            "cad_measurement_postprocess_row_convention_schema_id_consistent_when_present"
        ]
        is False
    )
    assert (
        stale_measurement_row_convention_gate["checks"][
            "expected_cad_measurement_postprocess_row_convention_schema_id_matches"
        ]
        is False
    )

    missing_measurement_row_convention_row = {
        key: value
        for key, value in row.items()
        if key != "cad_measurement_postprocess_row_convention_schema_id"
    }
    missing_measurement_row_convention_manifest = [
        {
            key: value
            for key, value in item.items()
            if key != "cad_measurement_postprocess_row_convention_schema_id"
        }
        for item in manifest
    ]
    missing_measurement_row_convention_gate = shape_cad_handoff_manifest_gate(
        [missing_measurement_row_convention_row],
        file_manifest=missing_measurement_row_convention_manifest,
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        required_file_kinds=("step", "build123d_measurement_json", "external_volume_summary_json", "cubit_quality_json"),
        expected_measurement_postprocess_row_convention_schema_id=(
            "build123d_occt_mass_property_row_convention_v1"
        ),
        require_measurement_postprocess_row_convention_schema=True,
    )
    assert missing_measurement_row_convention_gate["status"] == "needs_attention"
    assert (
        missing_measurement_row_convention_gate["checks"][
            "cad_measurement_postprocess_row_convention_schema_id_recorded_when_required"
        ]
        is False
    )
    assert (
        missing_measurement_row_convention_gate["checks"][
            "cad_measurement_postprocess_row_convention_schema_id_recorded_when_expected"
        ]
        is False
    )

    stale_measurement_component_basis = [dict(item) for item in manifest]
    stale_measurement_component_basis[2]["cad_measurement_component_basis_schema_id"] = (
        "build123d_scalar_volume_component_basis_v0"
    )
    stale_measurement_component_basis_gate = shape_cad_handoff_manifest_gate(
        [row],
        file_manifest=stale_measurement_component_basis,
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        required_file_kinds=("step", "build123d_measurement_json", "external_volume_summary_json", "cubit_quality_json"),
        expected_measurement_convention="occt_closed_solid_mass_properties",
        expected_measurement_postprocess_row_convention_schema_id=(
            "build123d_occt_mass_property_row_convention_v1"
        ),
        expected_measurement_component_basis_schema_id=(
            "build123d_occt_volume_area_bbox_component_basis_v1"
        ),
        require_measurement_component_basis_schema=True,
    )
    assert stale_measurement_component_basis_gate["status"] == "needs_attention"
    assert stale_measurement_component_basis_gate["checks"]["expected_cad_measurement_convention_matches"] is True
    assert (
        stale_measurement_component_basis_gate["checks"][
            "expected_cad_measurement_postprocess_row_convention_schema_id_matches"
        ]
        is True
    )
    assert (
        stale_measurement_component_basis_gate["checks"][
            "cad_measurement_component_basis_schema_id_consistent_when_present"
        ]
        is False
    )
    assert (
        stale_measurement_component_basis_gate["checks"][
            "expected_cad_measurement_component_basis_schema_id_matches"
        ]
        is False
    )

    missing_measurement_component_basis_row = {
        key: value
        for key, value in row.items()
        if key != "cad_measurement_component_basis_schema_id"
    }
    missing_measurement_component_basis_manifest = [
        {
            key: value
            for key, value in item.items()
            if key != "cad_measurement_component_basis_schema_id"
        }
        for item in manifest
    ]
    missing_measurement_component_basis_gate = shape_cad_handoff_manifest_gate(
        [missing_measurement_component_basis_row],
        file_manifest=missing_measurement_component_basis_manifest,
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        required_file_kinds=("step", "build123d_measurement_json", "external_volume_summary_json", "cubit_quality_json"),
        expected_measurement_component_basis_schema_id=(
            "build123d_occt_volume_area_bbox_component_basis_v1"
        ),
        require_measurement_component_basis_schema=True,
    )
    assert missing_measurement_component_basis_gate["status"] == "needs_attention"
    assert (
        missing_measurement_component_basis_gate["checks"][
            "cad_measurement_component_basis_schema_id_recorded_when_required"
        ]
        is False
    )
    assert (
        missing_measurement_component_basis_gate["checks"][
            "cad_measurement_component_basis_schema_id_recorded_when_expected"
        ]
        is False
    )

    stale_solver_contract_handoff = shape_cubit_solver_route_handoff_gate(
        [row],
        {
            **solver_route_gate,
            "solver_contract_digest": "sha256:stale-reader-contract",
            "checks": {
                **solver_route_gate["checks"],
                "expected_solver_contract_digest_matches": False,
            },
        },
        expected_solver_route_package_id="slot347_terminal_block_solver_route_v1",
        expected_solver_contract_digest="sha256:slot384-terminal-block-reader-contract",
        require_solver_contract_artifact=True,
    )
    assert stale_solver_contract_handoff["status"] == "needs_attention"
    assert stale_solver_contract_handoff["checks"]["expected_solver_contract_digest_matches"] is False

    stale_route_convention_handoff = shape_cubit_solver_route_handoff_gate(
        [row],
        {
            **solver_route_gate,
            "solver_route_convention_schema_id": "coreform_value_only_mixed_route_v0",
            "checks": {
                **solver_route_gate["checks"],
                "expected_solver_route_convention_schema_id_matches": False,
            },
        },
        expected_solver_route_package_id="slot347_terminal_block_solver_route_v1",
        expected_solver_route_convention_schema_id="coreform_mixed_hex_pyramid_tet_route_convention_v1",
        require_solver_route_convention_schema=True,
    )
    assert stale_route_convention_handoff["status"] == "needs_attention"
    assert (
        stale_route_convention_handoff["checks"]["expected_solver_route_convention_schema_id_matches"]
        is False
    )

    missing_route_convention_handoff = shape_cubit_solver_route_handoff_gate(
        [row],
        {
            **solver_route_gate,
            "solver_route_convention_schema_id": "",
            "checks": {
                **solver_route_gate["checks"],
                "solver_route_convention_schema_id_recorded_when_required": False,
            },
        },
        expected_solver_route_package_id="slot347_terminal_block_solver_route_v1",
        require_solver_route_convention_schema=True,
    )
    assert missing_route_convention_handoff["status"] == "needs_attention"
    assert (
        missing_route_convention_handoff["checks"]["solver_route_convention_schema_id_recorded_when_required"]
        is False
    )

    stale_solver_route_handoff = shape_cubit_solver_route_handoff_gate(
        [row],
        {
            **solver_route_gate,
            "status": "needs_attention",
            "checks": {
                **solver_route_gate["checks"],
                "pyramid_transition_role_recorded": False,
                "no_implicit_tetization_recorded": False,
            },
        },
        expected_solver_route_package_id="slot347_terminal_block_solver_route_v1",
    )
    bad_solver_route_gate = shape_cad_handoff_manifest_gate(
        [row],
        file_manifest=manifest,
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        cubit_solver_route_handoff=stale_solver_route_handoff,
        required_file_kinds=("step", "build123d_measurement_json", "external_volume_summary_json", "cubit_quality_json"),
    )
    assert stale_solver_route_handoff["status"] == "needs_attention"
    assert stale_solver_route_handoff["checks"]["solver_route_gate_ok"] is False
    assert stale_solver_route_handoff["checks"]["solver_route_pyramid_transition_role_recorded"] is False
    assert stale_solver_route_handoff["checks"]["solver_route_no_implicit_tetization"] is False
    assert bad_solver_route_gate["status"] == "needs_attention"
    assert bad_solver_route_gate["checks"]["cubit_solver_route_handoff_ok"] is False

    stale_quality_ledger = {
        **quality_ledger_handoff,
        "status": "needs_attention",
        "checks": {
            **quality_ledger_handoff["checks"],
            "expected_mesh_digest_matches": False,
        },
    }
    bad_quality_ledger_gate = shape_cad_handoff_manifest_gate(
        [row],
        file_manifest=manifest,
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        cubit_quality_ledger_handoff=stale_quality_ledger,
        required_file_kinds=(
            "step",
            "build123d_measurement_json",
            "external_volume_summary_json",
            "cubit_quality_json",
            "cubit_quality_ledger_json",
        ),
    )
    assert stale_quality_ledger["status"] == "needs_attention"
    assert bad_quality_ledger_gate["status"] == "needs_attention"
    assert bad_quality_ledger_gate["checks"]["cubit_quality_ledger_handoff_ok"] is False

    stale_scheme_trace = dict(scheme_trace)
    stale_scheme_trace["export_output_digest"] = "sha256:old-slot362-vol"
    stale_scheme_handoff = shape_cubit_meshing_scheme_intent_gate(
        [row],
        scheme_trace_gate=stale_scheme_trace,
        required_roles=("hex_region",),
        expected_scheme_by_role={"hex_region": "map"},
        required_command_fragments=("create brick", "volume 1 scheme map", "export netgen"),
        expected_trace_id="slot363_terminal_block_scheme_trace",
        expected_export_order=2,
        expected_export_output_artifact_id="slot363_terminal_block_vol_v1",
        expected_export_output_digest="sha256:slot363-terminal-block-vol",
        expected_export_output_path=r"artifacts/cubit/slot363_terminal_block.vol",
        require_downstream_export_output_artifact=True,
    )
    stale_scheme_gate = shape_cad_handoff_manifest_gate(
        [row],
        file_manifest=manifest,
        external_volume_summary=external_volume,
        cubit_quality_handoff=quality_handoff,
        cubit_meshing_scheme_handoff=stale_scheme_handoff,
        required_file_kinds=("step", "build123d_measurement_json", "external_volume_summary_json", "cubit_quality_json"),
    )
    assert stale_scheme_handoff["status"] == "needs_attention"
    assert stale_scheme_handoff["checks"]["downstream_export_output_digest_matches"] is False
    assert stale_scheme_gate["status"] == "needs_attention"
    assert stale_scheme_gate["checks"]["cubit_meshing_scheme_handoff_ok"] is False


def test_build123d_submodel_cad_handoff_gate_binds_crop_recipe_and_boundary_contract():
    local = Box(1.0, 1.0, 1.0).solid()
    local.label = "slot251_local_tip_crop"
    row = shape_measurement_row(local)
    row["geometry_id"] = "slot251_local_tip_crop_v1"
    boundary_handoff = {
        "policy": "cubit_submodel_boundary_handoff_mesh_package_gate",
        "status": "ok",
        "submodel_region_id": "slot251_zoom_region_tip_01",
        "zoom_boundary_id": "zoom_boundary_outer",
        "boundary_transfer_error_estimate": 0.018,
        "checks": {
            "boundary_transfer_error_estimate_recorded": True,
        },
    }
    files = [
        {"kind": "step", "path": r"artifacts/build123d/slot251_local_tip_crop.step"},
        {"kind": "build123d_measurement_json", "path": r"artifacts/build123d/slot251_local_tip_crop_measure.json"},
    ]

    gate = shape_submodel_cad_handoff_gate(
        [row],
        recipe_id="slot251_local_box_crop_recipe",
        parent_model_id="slot249_global_plate_bending_coarse_v1",
        submodel_region_id="slot251_zoom_region_tip_01",
        crop_box={"min": [-1.0, -1.0, -1.0], "max": [1.0, 1.0, 1.0]},
        export_id="slot251_local_tip_crop_step_v1",
        unit="mm",
        file_manifest=files,
        boundary_handoff=boundary_handoff,
        expected_geometry_ids=("slot251_local_tip_crop_v1",),
    )

    assert gate["policy"] == "build123d_submodel_cad_handoff_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["recipe_id_recorded"] is True
    assert gate["checks"]["shape_bboxes_inside_crop"] is True
    assert gate["checks"]["step_file_present"] is True
    assert gate["checks"]["measurement_json_present"] is True
    assert gate["checks"]["boundary_handoff_ok"] is True
    assert gate["checks"]["boundary_handoff_submodel_matches"] is True
    assert gate["checks"]["boundary_handoff_error_recorded"] is True

    crop_too_small = shape_submodel_cad_handoff_gate(
        [row],
        recipe_id="slot251_local_box_crop_recipe",
        parent_model_id="slot249_global_plate_bending_coarse_v1",
        submodel_region_id="slot251_zoom_region_tip_01",
        crop_box={"min": [-0.25, -0.25, -0.25], "max": [0.25, 0.25, 0.25]},
        export_id="slot251_local_tip_crop_step_v1",
        unit="mm",
        file_manifest=files,
        boundary_handoff=boundary_handoff,
    )
    assert crop_too_small["status"] == "needs_attention"
    assert crop_too_small["checks"]["shape_bboxes_inside_crop"] is False

    wrong_boundary = shape_submodel_cad_handoff_gate(
        [row],
        recipe_id="slot251_local_box_crop_recipe",
        parent_model_id="slot249_global_plate_bending_coarse_v1",
        submodel_region_id="slot251_zoom_region_tip_01",
        crop_box={"min": [-1.0, -1.0, -1.0], "max": [1.0, 1.0, 1.0]},
        export_id="slot251_local_tip_crop_step_v1",
        unit="mm",
        file_manifest=files,
        boundary_handoff={**boundary_handoff, "submodel_region_id": "other_region"},
    )
    assert wrong_boundary["status"] == "needs_attention"
    assert wrong_boundary["checks"]["boundary_handoff_submodel_matches"] is False

    missing_files = shape_submodel_cad_handoff_gate(
        [row],
        recipe_id="slot251_local_box_crop_recipe",
        parent_model_id="slot249_global_plate_bending_coarse_v1",
        submodel_region_id="slot251_zoom_region_tip_01",
        crop_box={"min": [-1.0, -1.0, -1.0], "max": [1.0, 1.0, 1.0]},
        export_id="slot251_local_tip_crop_step_v1",
        unit="mm",
        file_manifest=[],
        boundary_handoff=boundary_handoff,
    )
    assert missing_files["status"] == "needs_attention"
    assert missing_files["checks"]["file_manifest_present"] is False


def test_build123d_submodel_cad_handoff_gate_pairs_transition_intent_with_pyramid_boundary():
    local = Box(1.0, 1.0, 1.0).solid()
    local.label = "slot259_local_mixed_crop"
    row = shape_measurement_row(local)
    row["geometry_id"] = "slot259_local_mixed_crop_v1"
    boundary_handoff = {
        "policy": "cubit_submodel_boundary_handoff_mesh_package_gate",
        "status": "ok",
        "submodel_region_id": "slot259_zoom_region_tip_01",
        "zoom_boundary_id": "zoom_boundary_outer",
        "boundary_transfer_error_estimate": 0.012,
        "volume_kind_counts": {"hex": 1, "pyramid": 1, "tet": 1},
        "surface_kind_counts": {"quad": 6, "triangle": 4},
        "transition_policy": "keep pyramid bridge as an explicit conformal hex-to-tet transition",
        "checks": {
            "boundary_transfer_error_estimate_recorded": True,
            "transition_policy_recorded_when_present": True,
        },
    }
    transition_rows = [
        {
            "name": "slot259_hex_core",
            "role": "hex_region",
            "material": "core_steel",
            "volume": 1.0,
        },
        {
            "name": "slot259_pyramid_transition_envelope",
            "role": "mesh_transition",
            "material": "transition_air",
            "transition_kind": "pyramid",
            "connects_roles": ["hex_region", "tet_region"],
            "expected_surface_kinds": ["quad", "triangle"],
            "volume": 0.25,
        },
        {
            "name": "slot259_tet_region",
            "role": "tet_region",
            "material": "air",
            "volume": 1.0,
        },
    ]
    transition_handoff = shape_transition_role_metadata_gate(
        transition_rows,
        required_surface_kinds=("quad", "triangle"),
        source_label="slot259_build123d",
    )
    files = [
        {"kind": "step", "path": r"artifacts/build123d/slot259_local_mixed_crop.step"},
        {"kind": "build123d_measurement_json", "path": r"artifacts/build123d/slot259_local_mixed_crop_measure.json"},
    ]

    gate = shape_submodel_cad_handoff_gate(
        [row],
        recipe_id="slot259_local_mixed_crop_recipe",
        parent_model_id="slot249_global_plate_bending_coarse_v1",
        submodel_region_id="slot259_zoom_region_tip_01",
        crop_box={"min": [-1.0, -1.0, -1.0], "max": [1.0, 1.0, 1.0]},
        export_id="slot259_local_mixed_crop_step_v1",
        unit="mm",
        file_manifest=files,
        boundary_handoff=boundary_handoff,
        transition_handoff=transition_handoff,
        expected_geometry_ids=("slot259_local_mixed_crop_v1",),
    )

    assert gate["status"] == "ok"
    assert gate["checks"]["transition_handoff_ok"] is True
    assert gate["checks"]["transition_handoff_present_for_pyramid_boundary"] is True
    assert gate["checks"]["transition_handoff_kind_matches_boundary"] is True
    assert gate["checks"]["transition_handoff_connects_hex_tet"] is True
    assert gate["checks"]["transition_handoff_surface_kinds_recorded"] is True
    assert gate["checks"]["transition_handoff_surface_kinds_match_boundary"] is True
    assert gate["boundary_volume_kind_counts"]["pyramid"] == 1
    assert gate["boundary_surface_kind_counts"] == {"quad": 6, "triangle": 4}
    assert gate["transition_handoff_kinds"] == ["pyramid"]
    assert gate["transition_handoff_surface_kinds"] == ["quad", "triangle"]

    missing_transition = shape_submodel_cad_handoff_gate(
        [row],
        recipe_id="slot259_local_mixed_crop_recipe",
        parent_model_id="slot249_global_plate_bending_coarse_v1",
        submodel_region_id="slot259_zoom_region_tip_01",
        crop_box={"min": [-1.0, -1.0, -1.0], "max": [1.0, 1.0, 1.0]},
        export_id="slot259_local_mixed_crop_step_v1",
        unit="mm",
        file_manifest=files,
        boundary_handoff=boundary_handoff,
    )
    assert missing_transition["status"] == "needs_attention"
    assert missing_transition["checks"]["transition_handoff_present_for_pyramid_boundary"] is False

    wrong_kind = {
        **transition_handoff,
        "transition_kinds": ["wedge"],
        "status": "ok",
    }
    wrong_kind_gate = shape_submodel_cad_handoff_gate(
        [row],
        recipe_id="slot259_local_mixed_crop_recipe",
        parent_model_id="slot249_global_plate_bending_coarse_v1",
        submodel_region_id="slot259_zoom_region_tip_01",
        crop_box={"min": [-1.0, -1.0, -1.0], "max": [1.0, 1.0, 1.0]},
        export_id="slot259_local_mixed_crop_step_v1",
        unit="mm",
        file_manifest=files,
        boundary_handoff=boundary_handoff,
        transition_handoff=wrong_kind,
    )
    assert wrong_kind_gate["status"] == "needs_attention"
    assert wrong_kind_gate["checks"]["transition_handoff_kind_matches_boundary"] is False

    missing_quad = {
        **transition_handoff,
        "surface_kinds": ["triangle"],
        "required_surface_kinds": ["quad", "triangle"],
        "status": "ok",
    }
    missing_quad_gate = shape_submodel_cad_handoff_gate(
        [row],
        recipe_id="slot259_local_mixed_crop_recipe",
        parent_model_id="slot249_global_plate_bending_coarse_v1",
        submodel_region_id="slot259_zoom_region_tip_01",
        crop_box={"min": [-1.0, -1.0, -1.0], "max": [1.0, 1.0, 1.0]},
        export_id="slot259_local_mixed_crop_step_v1",
        unit="mm",
        file_manifest=files,
        boundary_handoff=boundary_handoff,
        transition_handoff=missing_quad,
    )
    assert missing_quad_gate["status"] == "needs_attention"
    assert missing_quad_gate["checks"]["transition_handoff_surface_kinds_match_boundary"] is False


def test_build123d_submodel_cad_handoff_gate_pairs_material_intent_with_cubit_sidecar_labels():
    local = Box(1.0, 1.0, 1.0).solid()
    local.label = "slot275_local_mixed_crop"
    row = shape_measurement_row(local)
    row["geometry_id"] = "slot275_local_mixed_crop_v1"
    boundary_handoff = {
        "policy": "cubit_submodel_boundary_handoff_mesh_package_gate",
        "status": "ok",
        "submodel_region_id": "slot275_zoom_region_tip_01",
        "zoom_boundary_id": "zoom_boundary_outer",
        "boundary_transfer_error_estimate": 0.01,
        "volume_kind_counts": {"hex": 1, "pyramid": 1, "tet": 1},
        "surface_kind_counts": {"quad": 6, "triangle": 4},
        "material_names": ["hex_core", "pyramid_transition", "tet_region"],
        "allowed_zero_measurement_names": ["pyramid_transition"],
        "roles_present": ["hex_to_transition", "transition_to_tet"],
        "checks": {"boundary_transfer_error_estimate_recorded": True},
    }
    transition_rows = [
        {
            "name": "slot275_hex_body",
            "role": "hex_region",
            "material": "core_steel",
            "downstream_material_name": "hex_core",
            "volume": 1.0,
        },
        {
            "name": "slot275_transition_envelope",
            "role": "mesh_transition",
            "material": "transition_air",
            "transition_kind": "pyramid",
            "connects_roles": ["hex_region", "tet_region"],
            "expected_surface_kinds": ["quad", "triangle"],
            "expected_interface_roles": ["hex_to_transition", "transition_to_tet"],
            "downstream_material_name": "pyramid_transition",
            "volume": 0.25,
        },
        {
            "name": "slot275_tet_body",
            "role": "tet_region",
            "material": "air",
            "downstream_material_name": "tet_region",
            "volume": 1.0,
        },
    ]
    transition_handoff = shape_transition_role_metadata_gate(
        transition_rows,
        required_surface_kinds=("quad", "triangle"),
        required_interface_roles=("hex_to_transition", "transition_to_tet"),
        expected_downstream_material_names=("hex_core", "pyramid_transition", "tet_region"),
        allowed_zero_downstream_material_names=("pyramid_transition",),
        source_label="slot275_build123d",
    )
    files = [
        {"kind": "step", "path": r"artifacts/build123d/slot275_local_mixed_crop.step"},
        {"kind": "build123d_measurement_json", "path": r"artifacts/build123d/slot275_local_mixed_crop_measure.json"},
    ]

    gate = shape_submodel_cad_handoff_gate(
        [row],
        recipe_id="slot275_local_mixed_crop_recipe",
        parent_model_id="slot274_coreform_material_sidecar_v1",
        submodel_region_id="slot275_zoom_region_tip_01",
        crop_box={"min": [-1.0, -1.0, -1.0], "max": [1.0, 1.0, 1.0]},
        export_id="slot275_local_mixed_crop_step_v1",
        unit="mm",
        file_manifest=files,
        boundary_handoff=boundary_handoff,
        transition_handoff=transition_handoff,
        expected_geometry_ids=("slot275_local_mixed_crop_v1",),
    )

    assert gate["status"] == "ok"
    assert gate["checks"]["boundary_material_names_recorded"] is True
    assert gate["checks"]["transition_handoff_material_names_recorded"] is True
    assert gate["checks"]["transition_handoff_material_names_match_boundary"] is True
    assert gate["checks"]["transition_handoff_zero_material_names_match_boundary"] is True
    assert gate["checks"]["boundary_interface_roles_recorded"] is True
    assert gate["checks"]["transition_handoff_interface_roles_recorded"] is True
    assert gate["checks"]["transition_handoff_interface_roles_match_boundary"] is True
    assert gate["boundary_material_names"] == ["hex_core", "pyramid_transition", "tet_region"]
    assert gate["transition_handoff_material_names"] == ["hex_core", "pyramid_transition", "tet_region"]
    assert gate["boundary_interface_roles"] == ["hex_to_transition", "transition_to_tet"]
    assert gate["transition_handoff_interface_roles"] == ["hex_to_transition", "transition_to_tet"]

    stale_materials = {
        **transition_handoff,
        "downstream_material_names": ["hex_core", "tet_region"],
        "status": "ok",
    }
    stale_gate = shape_submodel_cad_handoff_gate(
        [row],
        recipe_id="slot275_local_mixed_crop_recipe",
        parent_model_id="slot274_coreform_material_sidecar_v1",
        submodel_region_id="slot275_zoom_region_tip_01",
        crop_box={"min": [-1.0, -1.0, -1.0], "max": [1.0, 1.0, 1.0]},
        export_id="slot275_local_mixed_crop_step_v1",
        unit="mm",
        file_manifest=files,
        boundary_handoff=boundary_handoff,
        transition_handoff=stale_materials,
        expected_geometry_ids=("slot275_local_mixed_crop_v1",),
    )
    assert stale_gate["status"] == "needs_attention"
    assert stale_gate["checks"]["transition_handoff_material_names_match_boundary"] is False

    stale_zero_contract = {
        **transition_handoff,
        "allowed_zero_downstream_material_names": [],
        "status": "ok",
    }
    stale_zero_gate = shape_submodel_cad_handoff_gate(
        [row],
        recipe_id="slot275_local_mixed_crop_recipe",
        parent_model_id="slot274_coreform_material_sidecar_v1",
        submodel_region_id="slot275_zoom_region_tip_01",
        crop_box={"min": [-1.0, -1.0, -1.0], "max": [1.0, 1.0, 1.0]},
        export_id="slot275_local_mixed_crop_step_v1",
        unit="mm",
        file_manifest=files,
        boundary_handoff=boundary_handoff,
        transition_handoff=stale_zero_contract,
        expected_geometry_ids=("slot275_local_mixed_crop_v1",),
    )
    assert stale_zero_gate["status"] == "needs_attention"
    assert stale_zero_gate["checks"]["transition_handoff_zero_material_names_match_boundary"] is False

    stale_interface_roles = {
        **transition_handoff,
        "interface_roles": ["hex_to_transition"],
        "required_interface_roles": ["hex_to_transition"],
        "status": "ok",
    }
    stale_interface_gate = shape_submodel_cad_handoff_gate(
        [row],
        recipe_id="slot275_local_mixed_crop_recipe",
        parent_model_id="slot274_coreform_material_sidecar_v1",
        submodel_region_id="slot275_zoom_region_tip_01",
        crop_box={"min": [-1.0, -1.0, -1.0], "max": [1.0, 1.0, 1.0]},
        export_id="slot275_local_mixed_crop_step_v1",
        unit="mm",
        file_manifest=files,
        boundary_handoff=boundary_handoff,
        transition_handoff=stale_interface_roles,
        expected_geometry_ids=("slot275_local_mixed_crop_v1",),
    )
    assert stale_interface_gate["status"] == "needs_attention"
    assert stale_interface_gate["checks"]["transition_handoff_interface_roles_match_boundary"] is False


def test_build123d_curvilinear_mesh_intent_gate_keeps_cubit_route_explicit():
    box = Box(2.0, 1.0, 1.0).solid()
    box.label = "curved_hex_block"
    row = shape_measurement_row(box)
    row["geometry_id"] = "slot171_curved_hex_block_v1"
    row["role"] = "hex_region"
    row["mesh_route"] = "cubit_hex_or_mixed_path"
    row["downstream_handoff"] = "cubit_curvilinear_handoff"
    manifest_gate = {
        "policy": "cubit_curvilinear_handoff_manifest_gate",
        "status": "ok",
        "checks": {
            "projection_error_recorded": True,
            "projection_error_within_tolerance": True,
            "negative_jacobian_count_recorded": True,
            "negative_jacobian_count_zero": True,
        },
    }
    order_series_gate = {
        "policy": "cubit_mixed_order_series_inventory_gate",
        "status": "ok",
        "checks": {
            "volume_kind_counts_invariant": True,
            "surface_kind_counts_invariant": True,
            "routing_hint_is_cubit_mixed": True,
            "first_order_inventory_present": True,
            "first_order_inventory_not_curved": True,
        },
    }

    gate = shape_curvilinear_mesh_intent_gate(
        [row],
        downstream_manifest_gate=manifest_gate,
        downstream_order_series_gate=order_series_gate,
    )

    assert gate["policy"] == "build123d_curvilinear_mesh_intent_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["required_roles_present"] is True
    assert gate["checks"]["route_is_cubit_hex_or_mixed"] is True
    assert gate["checks"]["handoff_label_recorded"] is True
    assert gate["checks"]["downstream_manifest_ok"] is True
    assert gate["checks"]["downstream_projection_error_ok"] is True
    assert gate["checks"]["downstream_negative_jacobian_zero"] is True
    assert gate["checks"]["downstream_order_series_ok"] is True
    assert gate["checks"]["downstream_order_series_policy_known"] is True
    assert gate["checks"]["downstream_order_series_topology_invariant"] is True
    assert gate["checks"]["downstream_order_series_route_matches"] is True
    assert gate["checks"]["downstream_order_series_first_order_inventory"] is True
    assert gate["roles"] == ["hex_region"]
    assert gate["routes"] == ["cubit_hex_or_mixed_path"]
    assert gate["handoffs"] == ["cubit_curvilinear_handoff"]
    assert gate["downstream_order_series_policy"] == "cubit_mixed_order_series_inventory_gate"

    tet_row = dict(row)
    tet_row["mesh_route"] = "netgen_tri_tet_path"
    bad_route = shape_curvilinear_mesh_intent_gate([tet_row])
    assert bad_route["status"] == "needs_attention"
    assert bad_route["checks"]["route_is_cubit_hex_or_mixed"] is False
    assert bad_route["checks"]["not_tet_only_route"] is False

    missing_role = dict(row)
    missing_role["role"] = "tet_only"
    bad_role = shape_curvilinear_mesh_intent_gate([missing_role])
    assert bad_role["status"] == "needs_attention"
    assert bad_role["checks"]["required_roles_present"] is False
    assert bad_role["checks"]["not_tet_only_route"] is False

    stale_manifest = {**manifest_gate, "status": "needs_attention"}
    bad_manifest = shape_curvilinear_mesh_intent_gate(
        [row],
        downstream_manifest_gate=stale_manifest,
    )
    assert bad_manifest["status"] == "needs_attention"
    assert bad_manifest["checks"]["downstream_manifest_ok"] is False

    poor_projection_manifest = {
        **manifest_gate,
        "checks": {
            **manifest_gate["checks"],
            "projection_error_within_tolerance": False,
        },
    }
    poor_projection = shape_curvilinear_mesh_intent_gate(
        [row],
        downstream_manifest_gate=poor_projection_manifest,
    )
    assert poor_projection["status"] == "needs_attention"
    assert poor_projection["checks"]["downstream_projection_error_ok"] is False
    assert "projection error within tolerance" in " ".join(poor_projection["issues"])

    inverted_manifest = {
        **manifest_gate,
        "checks": {
            **manifest_gate["checks"],
            "negative_jacobian_count_zero": False,
        },
    }
    inverted = shape_curvilinear_mesh_intent_gate(
        [row],
        downstream_manifest_gate=inverted_manifest,
    )
    assert inverted["status"] == "needs_attention"
    assert inverted["checks"]["downstream_negative_jacobian_zero"] is False
    assert "zero negative-Jacobian" in " ".join(inverted["issues"])

    stale_order_series = {
        **order_series_gate,
        "checks": {
            **order_series_gate["checks"],
            "volume_kind_counts_invariant": False,
        },
    }
    stale_order_gate = shape_curvilinear_mesh_intent_gate(
        [row],
        downstream_manifest_gate=manifest_gate,
        downstream_order_series_gate=stale_order_series,
    )
    assert stale_order_gate["status"] == "needs_attention"
    assert stale_order_gate["checks"]["downstream_order_series_topology_invariant"] is False
    assert "topology invariant" in " ".join(stale_order_gate["issues"])

    wrong_order_policy = {**order_series_gate, "policy": "cubit_mixed_transition_metadata_gate"}
    wrong_order_gate = shape_curvilinear_mesh_intent_gate(
        [row],
        downstream_order_series_gate=wrong_order_policy,
    )
    assert wrong_order_gate["status"] == "needs_attention"
    assert wrong_order_gate["checks"]["downstream_order_series_policy_known"] is False

    missing_first_order = {
        **order_series_gate,
        "checks": {
            **order_series_gate["checks"],
            "first_order_inventory_present": False,
        },
    }
    missing_first_gate = shape_curvilinear_mesh_intent_gate(
        [row],
        downstream_order_series_gate=missing_first_order,
    )
    assert missing_first_gate["status"] == "needs_attention"
    assert missing_first_gate["checks"]["downstream_order_series_first_order_inventory"] is False


def test_build123d_mesh_environment_handoff_gate_binds_cad_rows_to_installed_cubit_evidence():
    box = Box(2.0, 1.0, 0.5).solid()
    box.label = "mesh_ready_block"
    row = shape_measurement_row(box)
    row["geometry_id"] = "slot187_mesh_ready_block_v1"
    row["mesh_route"] = "cubit_hex_or_mixed_path"
    environment = cubit_headless_installation_route_gate(
        {
            "installed_version": "2025.12",
            "binary_path": "C:/Program Files/Coreform Cubit 2025.12/bin/coreform_cubit.com",
            "binary_exists": True,
            "headless_flags": ["-nographics", "-batch"],
            "gui_policy": "headless_no_gui_daemon_by_default",
            "allow_gui_daemon": False,
            "release_note_version": "2026.6",
            "release_note_status": "watchlist",
            "live_claimed_release_version": "2025.12",
            "license_status": "ValidStudent",
            "version_probe_command": "C:/Program Files/Coreform Cubit 2025.12/bin/coreform_cubit.com -version",
            "version_probe_summary": {
                "license_status": "ValidStudent",
                "version_line": "Coreform Cubit Version 2025.12 Build 3d8d3af7",
                "binary_kind": "coreform_cubit.com synchronous console probe",
            },
        }
    )

    gate = shape_mesh_environment_handoff_gate([row], environment)

    assert gate["policy"] == "build123d_mesh_environment_handoff_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["shape_rows_have_volume_area_bbox"] is True
    assert gate["checks"]["mesh_route_matches_expected"] is True
    assert gate["checks"]["mesh_environment_gate_ok"] is True
    assert gate["checks"]["headless_flags_present"] is True
    assert gate["checks"]["live_claim_matches_installed"] is True
    assert gate["checks"]["release_note_watchlist_not_live_claim"] is True
    assert gate["checks"]["license_status_allows_headless_probe"] is True
    assert gate["checks"]["version_probe_is_synchronous_console"] is True
    assert gate["checks"]["binary_path_is_console_com"] is True
    assert gate["checks"]["version_probe_uses_recorded_binary"] is True
    assert gate["checks"]["version_probe_summary_records_installed_version"] is True
    assert gate["checks"]["version_probe_summary_records_license_status"] is True
    assert gate["binary_path"].endswith("coreform_cubit.com")
    assert gate["license_status"] == "ValidStudent"
    assert "coreform_cubit.com -version" in gate["version_probe_command"]
    assert gate["version_probe_summary"]["version_line"].startswith("Coreform Cubit Version 2025.12")

    gui_binary_environment = cubit_headless_installation_route_gate(
        {
            "installed_version": "2025.12",
            "binary_path": "C:/Program Files/Coreform Cubit 2025.12/bin/coreform_cubit.exe",
            "binary_exists": True,
            "headless_flags": ["-nographics", "-batch"],
            "gui_policy": "headless_no_gui_daemon_by_default",
            "allow_gui_daemon": False,
            "release_note_version": "2026.6",
            "release_note_status": "watchlist",
            "live_claimed_release_version": "2025.12",
            "license_status": "ValidStudent",
            "version_probe_command": "C:/Program Files/Coreform Cubit 2025.12/bin/coreform_cubit.com -version",
        }
    )
    gui_binary = shape_mesh_environment_handoff_gate([row], gui_binary_environment)
    assert gui_binary["status"] == "needs_attention"
    assert gui_binary["checks"]["binary_path_is_console_com"] is False
    assert gui_binary["checks"]["version_probe_uses_recorded_binary"] is False
    assert any("coreform_cubit.com" in issue for issue in gui_binary["issues"])

    overclaim_environment = cubit_headless_installation_route_gate(
        {
            "installed_version": "2025.12",
            "binary_path": "C:/Program Files/Coreform Cubit 2025.12/bin/coreform_cubit.exe",
            "binary_exists": True,
            "headless_flags": ["-nographics", "-batch"],
            "gui_policy": "headless_no_gui_daemon_by_default",
            "release_note_version": "2026.6",
            "release_note_status": "installed",
            "live_claimed_release_version": "2026.6",
        }
    )
    overclaim = shape_mesh_environment_handoff_gate([row], overclaim_environment)
    assert overclaim["status"] == "needs_attention"
    assert overclaim["checks"]["mesh_environment_gate_ok"] is False
    assert overclaim["checks"]["live_claim_matches_installed"] is False

    stale_probe_summary_environment = cubit_headless_installation_route_gate(
        {
            "installed_version": "2025.12",
            "binary_path": "C:/Program Files/Coreform Cubit 2025.12/bin/coreform_cubit.com",
            "binary_exists": True,
            "headless_flags": ["-nographics", "-batch"],
            "gui_policy": "headless_no_gui_daemon_by_default",
            "allow_gui_daemon": False,
            "live_claimed_release_version": "2025.12",
            "license_status": "ValidStudent",
            "version_probe_summary": {
                "license_status": "ValidStudent",
                "version_line": "Coreform Cubit Version 2026.6 Build future",
            },
        }
    )
    stale_probe_summary = shape_mesh_environment_handoff_gate(
        [row], stale_probe_summary_environment
    )
    assert stale_probe_summary["status"] == "needs_attention"
    assert stale_probe_summary["checks"]["mesh_environment_gate_ok"] is False
    assert stale_probe_summary["checks"]["version_probe_summary_records_installed_version"] is False
    assert stale_probe_summary["checks"]["version_probe_summary_records_license_status"] is True
    assert any("version-probe summary" in issue for issue in stale_probe_summary["issues"])

    tet_row = dict(row)
    tet_row["mesh_route"] = "netgen_tri_tet_path"
    bad_route = shape_mesh_environment_handoff_gate([tet_row], environment)
    assert bad_route["status"] == "needs_attention"
    assert bad_route["checks"]["mesh_route_matches_expected"] is False


def test_box_through_cylinder_reference_matches_build123d_mass_properties():
    from build123d import Align, Cylinder

    x, y, z, radius = 4.0, 3.0, 2.0, 0.4
    part = Box(x, y, z) - Cylinder(
        radius=radius,
        height=1.5 * z,
        align=(Align.CENTER, Align.CENTER, Align.CENTER),
    )
    part = part.solid()
    part.label = "box_hole"

    reference = [box_through_cylinder_reference_row(x, y, z, radius, axis="z", label="box_hole")]
    measured = [shape_measurement_row(part)]
    health = shape_measurement_health_summary(reference, measured, rtol=1.0e-12, bbox_atol=1.0e-12)

    assert health["status"] == "ok"
    assert health["comparison_summary"]["max_volume_rel_error"] < 1.0e-14
    assert health["comparison_summary"]["max_area_rel_error"] < 1.0e-14
    assert health["comparison_summary"]["max_bbox_abs_error"] == pytest.approx(0.0)
    assert reference[0]["policy"] == "analytic_box_through_cylinder_mass_property_reference"


def test_shape_measurement_rows_follow_assembly_children():
    left = Box(1, 2, 3).solid()
    left.label = "left"
    right = Box(2, 2, 1).solid()
    right.label = "right"
    asm = assembly(left, right, label="two_region")
    rows = shape_measurement_rows(asm)

    assert [row["name"] for row in rows] == ["left", "right"]
    assert [row["volume"] for row in rows] == pytest.approx([6.0, 4.0])
    assert [row["area"] for row in rows] == pytest.approx([22.0, 16.0])


def test_shape_measurement_inventory_summary_reports_assembly_fractions():
    left = (Pos(-2, 0, 0) * Box(2, 2, 2)).solid()
    left.label = "left"
    right = (Pos(2, 0, 0) * Box(1, 2, 2)).solid()
    right.label = "right"
    rows = shape_measurement_rows(assembly(left, right, label="two_region"))

    summary = shape_measurement_inventory_summary(rows)

    assert summary["n_shapes"] == 2
    assert summary["n_valid"] == 2
    assert summary["total_volume"] == pytest.approx(12.0)
    assert summary["total_area"] == pytest.approx(40.0)
    assert summary["bounding_box"]["min"] == pytest.approx([-3.0, -1.0, -1.0])
    assert summary["bounding_box"]["max"] == pytest.approx([2.5, 1.0, 1.0])
    assert summary["bbox_volume"] == pytest.approx(22.0)
    assert summary["bbox_fill_fraction"] == pytest.approx(12.0 / 22.0)
    assert summary["largest_volume_name"] == "left"
    assert summary["smallest_volume_name"] == "right"
    fractions = {row["name"]: row for row in summary["volume_fraction_rows"]}
    assert fractions["left"]["volume_fraction"] == pytest.approx(8.0 / 12.0)
    assert fractions["right"]["volume_fraction"] == pytest.approx(4.0 / 12.0)


def test_shape_bbox_pair_clearance_summary_flags_overlap_for_precise_check():
    left = Box(1, 1, 1).solid()
    left.label = "left"
    right = (Pos(2.0, 0, 0) * Box(1, 1, 1)).solid()
    right.label = "right"
    overlap = (Pos(0.4, 0, 0) * Box(1, 1, 1)).solid()
    overlap.label = "overlap"
    rows = shape_measurement_rows(assembly(left, right, overlap, label="three_region"))

    summary = shape_bbox_pair_clearance_summary(rows)
    pairs = {row["pair"]: row for row in summary["pair_rows"]}

    assert summary["status"] == "needs_attention"
    assert summary["n_pairs"] == 3
    assert summary["separated_pair_count"] == 2
    assert summary["bbox_overlap_pair_count"] == 1
    assert summary["touching_pair_count"] == 0
    assert summary["min_positive_gap"] == pytest.approx(0.6)

    assert pairs["left::right"]["status"] == "separated"
    assert pairs["left::right"]["axis_gaps"]["x"] == pytest.approx(1.0)
    assert pairs["left::overlap"]["status"] == "bbox_overlap_needs_precise_check"
    assert pairs["left::overlap"]["bbox_intersection_volume"] == pytest.approx(0.6)
    assert pairs["left::overlap"]["axis_overlaps"]["x"] == pytest.approx(0.6)


def test_shape_parameter_sweep_summary_tracks_monotonic_metrics_and_limits():
    rows = []
    for height in [1.0, 2.0, 3.0, 4.0]:
        box = Box(2.0, 3.0, height).solid()
        row = shape_measurement_row(box, name=f"h_{height:g}")
        row["height"] = height
        rows.append(row)

    summary = shape_parameter_sweep_summary(
        rows,
        "height",
        metric_keys=("volume", "area"),
        limits_by_metric={"volume": {"min": 12.0, "max": 24.0}},
    )
    metrics = {row["metric"]: row for row in summary["metric_rows"]}

    assert summary["parameter_values"] == [1.0, 2.0, 3.0, 4.0]
    assert summary["parameter_strictly_increasing"] is True
    assert summary["status"] == "needs_attention"
    assert summary["constraint_violation_count"] == 1
    assert summary["constraint_violations"][0]["kind"] == "below_min"
    assert metrics["volume"]["min"] == pytest.approx(6.0)
    assert metrics["volume"]["max"] == pytest.approx(24.0)
    assert metrics["volume"]["monotonic_non_decreasing"] is True
    assert metrics["volume"]["min_step_delta"] == pytest.approx(6.0)
    assert metrics["area"]["first"] == pytest.approx(22.0)
    assert metrics["area"]["last"] == pytest.approx(52.0)
    assert metrics["area"]["monotonic_non_decreasing"] is True

    clean = shape_parameter_sweep_summary(rows, "height", metric_keys=("volume",))
    assert clean["status"] == "ok"
    assert clean["ok_for_design_table"] is True

    duplicate = shape_parameter_sweep_summary(rows + [dict(rows[-1])], "height", metric_keys=("volume",))
    assert duplicate["duplicate_parameter_values"] is True
    assert duplicate["status"] == "needs_attention"


def test_box_face_vector_area_rows_match_centered_box():
    box = Box(2, 3, 5).solid()
    measurement = shape_measurement_row(box)
    rows = box_face_vector_area_rows(measurement["bounding_box"]["size"])
    by_name = {row["name"]: row for row in rows}

    assert sum(row["surface_area"] for row in rows) == pytest.approx(measurement["area"])
    assert by_name["xmin"]["vector_area"] == pytest.approx((-15.0, 0.0, 0.0))
    assert by_name["xmax"]["vector_area"] == pytest.approx((15.0, 0.0, 0.0))
    assert by_name["ymin"]["vector_area"] == pytest.approx((0.0, -10.0, 0.0))
    assert by_name["ymax"]["vector_area"] == pytest.approx((0.0, 10.0, 0.0))
    assert by_name["zmin"]["vector_area"] == pytest.approx((0.0, 0.0, -6.0))
    assert by_name["zmax"]["vector_area"] == pytest.approx((0.0, 0.0, 6.0))
    assert all(row["vector_area_norm_over_area"] == pytest.approx(1.0) for row in rows)


def test_compare_boundary_vector_area_rows_marks_direction_mismatch():
    reference = box_face_vector_area_rows((2, 3, 5))
    measured = [dict(row) for row in reference]
    measured[0]["vector_area"] = (15.0, 0.0, 0.0)
    measured[0]["unit_normal"] = (1.0, 0.0, 0.0)

    rows = compare_boundary_vector_area_rows(reference, measured, vector_atol=1.0e-12)
    by_name = {row["name"]: row for row in rows}

    assert not by_name["xmin"]["passed"]
    assert by_name["xmin"]["vector_abs_error"] == pytest.approx(30.0)
    assert by_name["xmin"]["unit_normal_abs_error"] == pytest.approx(2.0)
    assert by_name["xmax"]["passed"]


def test_box_face_pressure_force_rows_integrate_pressure_loads():
    uniform = box_face_pressure_force_rows(
        (2, 3, 5),
        {"xmin": 2.0, "xmax": 2.0, "ymin": 2.0, "ymax": 2.0, "zmin": 2.0, "zmax": 2.0},
    )
    total = [sum(row["force_N"][axis] for row in uniform) for axis in range(3)]
    assert total == pytest.approx([0.0, 0.0, 0.0])

    zmax = box_face_pressure_force_rows((2, 3, 5), {"zmax": 2.0}, default_pressure=0.0)
    by_name = {row["name"]: row for row in zmax}
    assert by_name["zmax"]["force_N"] == pytest.approx((0.0, 0.0, 12.0))
    assert by_name["zmax"]["force_magnitude_N"] == pytest.approx(12.0)
    assert by_name["zmax"]["pressure_source"] == "name"
    assert by_name["xmin"]["pressure_source"] == "default"

    by_index = box_face_pressure_force_rows((2, 3, 5), {6: 3.0}, default_pressure=0.0)
    assert {row["name"]: row for row in by_index}["zmax"]["force_N"] == pytest.approx((0.0, 0.0, 18.0))

    with pytest.raises(KeyError):
        box_face_pressure_force_rows((2, 3, 5), {"zmax": 2.0}, default_pressure=None)


def test_box_face_pressure_moment_rows_integrate_pivot_moments():
    uniform = box_face_pressure_moment_rows(
        (2, 3, 5),
        {"xmin": 2.0, "xmax": 2.0, "ymin": 2.0, "ymax": 2.0, "zmin": 2.0, "zmax": 2.0},
        center=(1.0, 1.5, 2.5),
    )
    total_force = [sum(row["force_N"][axis] for row in uniform) for axis in range(3)]
    total_moment = [sum(row["moment_about_pivot_Nm"][axis] for row in uniform) for axis in range(3)]
    assert total_force == pytest.approx((0.0, 0.0, 0.0))
    assert total_moment == pytest.approx((0.0, 0.0, 0.0))

    zmax = box_face_pressure_moment_rows(
        (2, 3, 5),
        {"zmax": 2.0},
        center=(1.0, 1.5, 2.5),
        default_pressure=0.0,
    )
    by_name = {row["name"]: row for row in zmax}
    assert by_name["zmax"]["face_center"] == pytest.approx((1.0, 1.5, 5.0))
    assert by_name["zmax"]["force_N"] == pytest.approx((0.0, 0.0, 12.0))
    assert by_name["zmax"]["moment_about_pivot_Nm"] == pytest.approx((18.0, -12.0, 0.0))

    shifted = box_face_pressure_moment_rows(
        (2, 3, 5),
        {"zmax": 2.0},
        center=(1.0, 1.5, 2.5),
        default_pressure=0.0,
        pivot_m=(1.0, 1.5, 0.0),
    )
    assert {row["name"]: row for row in shifted}["zmax"]["moment_about_pivot_Nm"] == pytest.approx((0.0, 0.0, 0.0))

    with pytest.raises(KeyError):
        box_face_pressure_moment_rows((2, 3, 5), {"zmax": 2.0}, default_pressure=None)


def test_box_face_pressure_resultant_summary_matches_closed_box_balance():
    box = (Pos(1.0, 1.5, 2.5) * Box(2, 3, 5)).solid()
    measurement = shape_measurement_row(box)
    size = measurement["bounding_box"]["size"]
    center = measurement["bounding_box"]["center"]

    uniform = box_face_pressure_resultant_summary(size, {}, center=center, default_pressure=2.0)
    assert uniform["boundary_count"] == 6
    assert uniform["box_size"] == pytest.approx((2.0, 3.0, 5.0))
    assert uniform["box_center"] == pytest.approx((1.0, 1.5, 2.5))
    assert uniform["total_force_N"] == pytest.approx((0.0, 0.0, 0.0))
    assert uniform["total_moment_about_pivot_Nm"] == pytest.approx((0.0, 0.0, 0.0))
    assert uniform["absolute_force_sum_N"] == pytest.approx(124.0)
    assert uniform["force_balance_ratio"] == pytest.approx(0.0)
    assert uniform["surface_vector_area"] == pytest.approx((0.0, 0.0, 0.0))

    zmax = box_face_pressure_resultant_summary(
        size,
        {"zmax": 2.0},
        center=center,
        default_pressure=0.0,
    )
    assert zmax["total_force_N"] == pytest.approx((0.0, 0.0, 12.0))
    assert zmax["total_moment_about_pivot_Nm"] == pytest.approx((18.0, -12.0, 0.0))
    assert zmax["force_balance_ratio"] == pytest.approx(1.0)

    shifted = box_face_pressure_resultant_summary(
        size,
        {"zmax": 2.0},
        center=center,
        default_pressure=0.0,
        pivot_m=(1.0, 1.5, 0.0),
    )
    assert shifted["total_moment_about_pivot_Nm"] == pytest.approx((0.0, 0.0, 0.0))

    with pytest.raises(KeyError):
        box_face_pressure_resultant_summary(size, {"zmax": 2.0}, default_pressure=None)


def test_box_face_traction_moment_rows_integrate_vector_tractions():
    rows = box_face_traction_moment_rows(
        (2, 3, 5),
        {"zmax": (1.0, -2.0, 3.0)},
        center=(1.0, 1.5, 2.5),
        default_traction=(0.0, 0.0, 0.0),
    )
    by_name = {row["name"]: row for row in rows}
    assert by_name["zmax"]["face_center"] == pytest.approx((1.0, 1.5, 5.0))
    assert by_name["zmax"]["traction_N_per_m2"] == pytest.approx((1.0, -2.0, 3.0))
    assert by_name["zmax"]["force_N"] == pytest.approx((6.0, -12.0, 18.0))
    assert by_name["zmax"]["moment_about_pivot_Nm"] == pytest.approx((87.0, 12.0, -21.0))
    assert by_name["zmax"]["traction_source"] == "name"
    assert by_name["xmin"]["traction_source"] == "default"

    shifted = box_face_traction_moment_rows(
        (2, 3, 5),
        {"zmax": (1.0, -2.0, 3.0)},
        center=(1.0, 1.5, 2.5),
        default_traction=(0.0, 0.0, 0.0),
        pivot_m=(1.0, 1.5, 5.0),
    )
    assert {row["name"]: row for row in shifted}["zmax"]["moment_about_pivot_Nm"] == pytest.approx((0.0, 0.0, 0.0))

    by_index = box_face_traction_moment_rows((2, 3, 5), {6: (0.0, 0.0, 2.0)}, default_traction=(0.0, 0.0, 0.0))
    assert {row["name"]: row for row in by_index}["zmax"]["force_N"] == pytest.approx((0.0, 0.0, 12.0))

    with pytest.raises(KeyError):
        box_face_traction_moment_rows((2, 3, 5), {"zmax": (1.0, -2.0, 3.0)}, default_traction=None)
    with pytest.raises(ValueError):
        box_face_traction_moment_rows((2, 3, 5), {"zmax": (1.0, -2.0)}, default_traction=(0.0, 0.0, 0.0))


def test_shape_envelope_row_and_enclosing_box_use_union_bbox_margin():
    left = (Pos(-2, 0, 0) * Box(2, 4, 6)).solid()
    right = (Pos(3, 0, 1) * Box(2, 2, 2)).solid()

    row = shape_envelope_row([left, right], margin=(1, 2, 3), name="outer")
    assert row["n_shapes"] == 2
    assert row["min"] == pytest.approx([-4.0, -4.0, -6.0])
    assert row["max"] == pytest.approx([5.0, 4.0, 6.0])
    assert row["size"] == pytest.approx([9.0, 8.0, 12.0])
    assert row["center"] == pytest.approx([0.5, 0.0, 0.0])
    assert row["volume"] == pytest.approx(864.0)

    outer = enclosing_box([left, right], margin=(1, 2, 3), label="outer")
    assert outer.label == "outer"
    assert outer.volume == pytest.approx(864.0)
    size = outer.bounding_box().size
    assert [size.X, size.Y, size.Z] == pytest.approx([9.0, 8.0, 12.0])


def test_enclosure_difference_region_volume_area_and_clearance():
    inner = Box(2, 2, 2).solid()
    inner.label = "solid"
    outer = enclosing_box([inner], margin=1.0, label="outer")

    clear = enclosure_clearance_row(outer, [inner])
    assert clear["contained_by_bbox"]
    assert clear["min_clearance"] == pytest.approx(1.0)
    assert clear["clearances"] == pytest.approx({
        "xmin": 1.0,
        "xmax": 1.0,
        "ymin": 1.0,
        "ymax": 1.0,
        "zmin": 1.0,
        "zmax": 1.0,
    })
    assert clear["enclosure_volume"] == pytest.approx(64.0)
    assert clear["inner_volume_sum"] == pytest.approx(8.0)
    assert clear["nominal_void_volume"] == pytest.approx(56.0)
    assert clear["inner_volume_fraction"] == pytest.approx(0.125)

    void = enclosure_difference_region(outer, [inner], label="void")
    assert void.label == "void"
    assert void.is_valid
    assert void.volume == pytest.approx(56.0)
    assert void.area == pytest.approx(96.0 + 24.0)


def test_compare_shape_measurement_rows_marks_pass_fail_and_missing():
    reference = [
        {"name": "ok", "volume": 10.0, "area": 20.0},
        {"name": "bad", "volume": 10.0, "area": 20.0},
        {"name": "missing", "volume": 1.0, "area": 2.0},
    ]
    measured = [
        {"name": "ok", "volume": 10.00001, "area": 20.00001},
        {"name": "bad", "volume": 10.2, "area": 20.0},
    ]
    rows = compare_shape_measurement_rows(reference, measured, rtol=1.0e-4, measured_label="external")
    by_name = {row["name"]: row for row in rows}

    assert by_name["ok"]["passed"]
    assert by_name["ok"]["measured_label"] == "external"
    assert by_name["ok"]["volume_rel_error"] < 1.0e-4
    assert not by_name["bad"]["passed"]
    assert by_name["bad"]["reason"] == "outside tolerance"
    assert not by_name["missing"]["passed"]
    assert by_name["missing"]["reason"] == "missing measured row"


def test_shape_volume_crosscheck_summary_accepts_cubit_and_external_cad_rows():
    reference = [
        {"name": "iron", "volume": 100.0, "area": 220.0},
        {"name": "coil", "volume": 25.0, "area": 60.0},
    ]
    measured = {
        "cubit": [
            {"name": "iron", "volume": 100.0000001},
            {"name": "coil", "volume": 24.9999999},
        ],
        "external_cad": {
            "rows": [
                {"name": "iron", "volume": 100.0},
                {"name": "coil", "volume": 25.0},
            ],
        },
    }

    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-8)

    assert summary["policy"] == "build123d_external_cad_volume_crosscheck"
    assert summary["status"] == "ok"
    assert summary["ok_for_cad_roundtrip_volume"] is True
    assert summary["sources"] == ["cubit", "external_cad"]
    assert summary["max_volume_rel_error"] < 1.0e-8

    rows = compare_shape_volume_rows(reference, [{"name": "iron", "volume": 100.0}], measured_label="cubit")
    by_name = {row["name"]: row for row in rows}
    assert by_name["iron"]["passed"]
    assert not by_name["coil"]["passed"]
    assert by_name["coil"]["reason"] == "missing measured row"


def test_shape_volume_crosscheck_source_coverage_gate_requires_cubit_and_cst_rows():
    reference = [
        {"name": "stator", "volume": 100.0},
        {"name": "coil", "volume": 25.0},
    ]
    summary = shape_volume_crosscheck_summary(
        reference,
        {
            "cubit": [
                {"name": "stator", "volume": 100.0},
                {"name": "coil", "volume": 25.0},
            ],
            "cst_import": [
                {"name": "stator", "volume": 100.00000001},
                {"name": "coil", "volume": 24.99999999},
            ],
        },
        rtol=1.0e-8,
    )

    gate = shape_volume_crosscheck_source_coverage_gate(
        summary,
        required_sources=("cubit", "cst_import"),
        max_allowed_volume_rel_error=1.0e-8,
    )

    assert gate["policy"] == "build123d_volume_crosscheck_source_coverage_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["required_sources_present"] is True
    assert gate["checks"]["volume_error_within_limit"] is True

    missing_cst = shape_volume_crosscheck_summary(
        reference,
        {"cubit": [{"name": "stator", "volume": 100.0}, {"name": "coil", "volume": 25.0}]},
        rtol=1.0e-8,
    )
    bad = shape_volume_crosscheck_source_coverage_gate(missing_cst)
    assert bad["status"] == "needs_attention"
    assert bad["missing_sources"] == ["cst_import"]
    assert bad["checks"]["required_sources_present"] is False


def test_shape_volume_crosscheck_source_identity_gate_binds_methods_body_keys_and_artifacts():
    reference = [
        {"name": "stator", "volume": 100.0},
        {"name": "coil", "volume": 25.0},
    ]
    summary = shape_volume_crosscheck_summary(
        reference,
        [
            {
                "source": "cubit",
                "rows": [
                    {"name": "stator", "volume": 100.0},
                    {"name": "coil", "volume": 25.0},
                ],
                "measurement_method": "coreform_cubit_volume_command",
                "body_identity_key": "name",
                "source_artifact_id": "slot331_cubit_volume_rows_v1",
                "parameter_set_artifact_id": "slot391_cad_parameter_set_v1",
                "parameter_set_digest": "sha256:slot391-cad-parameter-set",
                "parameter_set_path": r"artifacts/build123d/slot391_cad_parameter_set.json",
                "objective_observable_id": "slot391_volume_quality_objective_v1",
                "objective_observable_family": "cad_volume_crosscheck",
            },
            {
                "source": "cst_import",
                "rows": [
                    {"name": "stator", "volume": 100.00000001},
                    {"name": "coil", "volume": 24.99999999},
                ],
                "measurement_method": "cst_modeler_solid_volume_export",
                "body_identity_key": "name",
                "source_artifact_id": "slot331_cst_volume_rows_v1",
                "parameter_set_artifact_id": "slot391_cad_parameter_set_v1",
                "parameter_set_digest": "sha256:slot391-cad-parameter-set",
                "parameter_set_path": r"artifacts/build123d/slot391_cad_parameter_set.json",
                "objective_observable_id": "slot391_volume_quality_objective_v1",
                "objective_observable_family": "cad_volume_crosscheck",
            },
        ],
        rtol=1.0e-8,
    )

    assert summary["status"] == "ok"
    by_source = {item["source"]: item for item in summary["comparison_sets"]}
    assert by_source["cubit"]["measurement_method"] == "coreform_cubit_volume_command"
    assert by_source["cst_import"]["source_artifact_id"] == "slot331_cst_volume_rows_v1"

    gate = shape_volume_crosscheck_source_identity_gate(
        summary,
        expected_measurement_methods={
            "cubit": "coreform_cubit_volume_command",
            "cst_import": "cst_modeler_solid_volume_export",
        },
        expected_body_identity_keys={"cubit": "name", "cst_import": "name"},
        expected_source_artifact_ids={
            "cubit": "slot331_cubit_volume_rows_v1",
            "cst_import": "slot331_cst_volume_rows_v1",
        },
        expected_parameter_set_artifact_ids={
            "cubit": "slot391_cad_parameter_set_v1",
            "cst_import": "slot391_cad_parameter_set_v1",
        },
        expected_parameter_set_digests={
            "cubit": "sha256:slot391-cad-parameter-set",
            "cst_import": "sha256:slot391-cad-parameter-set",
        },
        expected_parameter_set_paths={
            "cubit": r"artifacts/build123d/slot391_cad_parameter_set.json",
            "cst_import": r"artifacts/build123d/slot391_cad_parameter_set.json",
        },
        expected_objective_observable_ids={
            "cubit": "slot391_volume_quality_objective_v1",
            "cst_import": "slot391_volume_quality_objective_v1",
        },
        expected_objective_observable_families={
            "cubit": "cad_volume_crosscheck",
            "cst_import": "cad_volume_crosscheck",
        },
    )

    assert gate["policy"] == "build123d_volume_crosscheck_source_identity_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["expected_measurement_methods_match"] is True
    assert gate["checks"]["expected_body_identity_keys_match"] is True
    assert gate["checks"]["expected_source_artifact_ids_match"] is True
    assert gate["checks"]["expected_parameter_set_artifact_ids_match"] is True
    assert gate["checks"]["expected_parameter_set_digests_match"] is True
    assert gate["checks"]["expected_parameter_set_paths_match"] is True
    assert gate["checks"]["expected_objective_observable_ids_match"] is True
    assert gate["checks"]["expected_objective_observable_families_match"] is True

    stale_method = shape_volume_crosscheck_source_identity_gate(
        summary,
        expected_measurement_methods={"cubit": "mesh_volume_after_import"},
        expected_body_identity_keys={"cubit": "name"},
        expected_source_artifact_ids={"cubit": "slot331_cubit_volume_rows_v1"},
    )
    assert stale_method["status"] == "needs_attention"
    assert stale_method["checks"]["expected_measurement_methods_match"] is False

    missing_artifact_summary = {
        **summary,
        "comparison_sets": [
            {key: value for key, value in item.items() if key != "source_artifact_id"}
            for item in summary["comparison_sets"]
        ],
    }
    missing_artifact = shape_volume_crosscheck_source_identity_gate(
        missing_artifact_summary,
        expected_source_artifact_ids={"cubit": "slot331_cubit_volume_rows_v1"},
    )
    assert missing_artifact["status"] == "needs_attention"
    assert missing_artifact["checks"]["source_artifact_ids_recorded_when_expected"] is False

    stale_parameter_digest = shape_volume_crosscheck_source_identity_gate(
        summary,
        expected_parameter_set_digests={"cubit": "sha256:stale-cad-parameter-set"},
    )
    assert stale_parameter_digest["status"] == "needs_attention"
    assert stale_parameter_digest["checks"]["expected_parameter_set_digests_match"] is False

    missing_parameter_path_summary = {
        **summary,
        "comparison_sets": [
            {key: value for key, value in item.items() if key != "parameter_set_path"}
            for item in summary["comparison_sets"]
        ],
    }
    missing_parameter_path = shape_volume_crosscheck_source_identity_gate(
        missing_parameter_path_summary,
        expected_parameter_set_paths={
            "cubit": r"artifacts/build123d/slot391_cad_parameter_set.json",
        },
    )
    assert missing_parameter_path["status"] == "needs_attention"
    assert missing_parameter_path["checks"]["parameter_set_paths_recorded_when_expected"] is False

    wrong_objective_family = shape_volume_crosscheck_source_identity_gate(
        summary,
        expected_objective_observable_families={"cubit": "mesh_quality"},
    )
    assert wrong_objective_family["status"] == "needs_attention"
    assert wrong_objective_family["checks"]["expected_objective_observable_families_match"] is False


def test_shape_external_cad_volume_evidence_package_binds_coverage_identity_and_route():
    box = Box(2, 3, 4).solid()
    row = shape_measurement_row(box, name="cad_block")
    row.update({
        "geometry_id": "cad_block_revB",
        "recipe_id": "slot339_box_recipe",
        "cad_kernel": "occt",
        "cad_kernel_version": "7.8-test",
        "script_path": r"artifacts/build123d/slot339_build123d.py",
        "export_id": "slot339_box_step",
        "authoring_source": "build123d_occt",
        "mesh_route": "cubit_hex_or_mixed_path",
        "length_unit": "m",
        "area_unit": "m^2",
        "volume_unit": "m^3",
    })
    measured = [
        {
            "source": "cubit",
            "rows": [{"name": "cad_block", "volume": 24.0}],
            "measurement_method": "coreform_cubit_volume_command",
            "body_identity_key": "name",
            "source_artifact_id": "slot339_cubit_volume_rows_v1",
            "parameter_set_artifact_id": "slot391_cad_parameter_set_v1",
            "parameter_set_digest": "sha256:slot391-cad-parameter-set",
            "parameter_set_path": r"artifacts/build123d/slot391_cad_parameter_set.json",
            "objective_observable_id": "slot391_volume_quality_objective_v1",
            "objective_observable_family": "cad_volume_crosscheck",
        },
        {
            "source": "cst_import",
            "rows": [{"name": "cad_block", "volume": 24.0}],
            "measurement_method": "cst_modeler_solid_volume_export",
            "body_identity_key": "name",
            "source_artifact_id": "slot339_cst_volume_rows_v1",
            "parameter_set_artifact_id": "slot391_cad_parameter_set_v1",
            "parameter_set_digest": "sha256:slot391-cad-parameter-set",
            "parameter_set_path": r"artifacts/build123d/slot391_cad_parameter_set.json",
            "objective_observable_id": "slot391_volume_quality_objective_v1",
            "objective_observable_family": "cad_volume_crosscheck",
        },
    ]
    summary = shape_volume_crosscheck_summary([row], measured, rtol=1.0e-12)

    package = shape_external_cad_volume_evidence_package_gate(
        [row],
        summary,
        required_sources=("cubit", "cst_import"),
        max_allowed_volume_rel_error=1.0e-12,
        expected_measurement_methods={
            "cubit": "coreform_cubit_volume_command",
            "cst_import": "cst_modeler_solid_volume_export",
        },
        expected_body_identity_keys={"cubit": "name", "cst_import": "name"},
        expected_source_artifact_ids={
            "cubit": "slot339_cubit_volume_rows_v1",
            "cst_import": "slot339_cst_volume_rows_v1",
        },
        expected_parameter_set_artifact_ids={
            "cubit": "slot391_cad_parameter_set_v1",
            "cst_import": "slot391_cad_parameter_set_v1",
        },
        expected_parameter_set_digests={
            "cubit": "sha256:slot391-cad-parameter-set",
            "cst_import": "sha256:slot391-cad-parameter-set",
        },
        expected_parameter_set_paths={
            "cubit": r"artifacts/build123d/slot391_cad_parameter_set.json",
            "cst_import": r"artifacts/build123d/slot391_cad_parameter_set.json",
        },
        expected_objective_observable_ids={
            "cubit": "slot391_volume_quality_objective_v1",
            "cst_import": "slot391_volume_quality_objective_v1",
        },
        expected_objective_observable_families={
            "cubit": "cad_volume_crosscheck",
            "cst_import": "cad_volume_crosscheck",
        },
        expected_route="cubit_hex_or_mixed_path",
        expected_length_unit="m",
        expected_area_unit="m^2",
        expected_volume_unit="m^3",
        required_metadata_fields=("recipe_id", "cad_kernel", "cad_kernel_version", "script_path", "export_id"),
    )

    assert package["policy"] == "build123d_external_cad_volume_evidence_package_gate"
    assert package["status"] == "ok"
    assert package["coverage_gate_status"] == "ok"
    assert package["source_identity_gate_status"] == "ok"
    assert package["route_contract_gate_status"] == "ok"
    assert package["checks"]["source_identity_metadata_complete"] is True
    assert package["checks"]["gate_source_sets_consistent"] is True
    assert package["source_identity_gate"]["checks"]["expected_parameter_set_digests_match"] is True
    assert package["source_identity_gate"]["checks"]["expected_objective_observable_families_match"] is True

    missing_artifact_summary = {
        **summary,
        "comparison_sets": [
            {key: value for key, value in item.items() if not (item["source"] == "cst_import" and key == "source_artifact_id")}
            for item in summary["comparison_sets"]
        ],
    }
    missing_artifact = shape_external_cad_volume_evidence_package_gate([row], missing_artifact_summary)
    assert missing_artifact["status"] == "needs_attention"
    assert missing_artifact["checks"]["source_identity_metadata_complete"] is False

    tet_row = dict(row)
    tet_row["mesh_route"] = "netgen_tri_tet_path"
    wrong_route = shape_external_cad_volume_evidence_package_gate([tet_row], summary)
    assert wrong_route["status"] == "needs_attention"
    assert wrong_route["checks"]["route_contract_gate_ok"] is False


def test_shape_cad_route_source_contract_gate_requires_route_and_source_groups():
    box = Box(2, 3, 4).solid()
    row = shape_measurement_row(box, name="hex_block")
    row.update({
        "geometry_id": "hex_block_revA",
        "recipe_id": "slot243_box_recipe",
        "cad_kernel": "occt",
        "cad_kernel_version": "7.8-test",
        "script_path": r"artifacts/build123d/slot243_build123d.py",
        "export_id": "slot243_box_step",
        "authoring_source": "build123d_occt",
        "mesh_route": "cubit_hex_or_mixed_path",
        "length_unit": "m",
        "area_unit": "m^2",
        "volume_unit": "m^3",
    })
    summary = shape_mass_property_crosscheck_summary(
        [row],
        {
            "coreform_cubit": [dict(row)],
            "external_cad": [dict(row)],
        },
        rtol=1.0e-12,
        bbox_atol=1.0e-12,
    )

    gate = shape_cad_route_source_contract_gate(
        [row],
        summary,
        expected_length_unit="m",
        expected_area_unit="m^2",
        expected_volume_unit="m^3",
        required_metadata_fields=(
            "recipe_id",
            "cad_kernel",
            "cad_kernel_version",
            "script_path",
            "export_id",
        ),
    )

    assert gate["policy"] == "build123d_cad_route_source_contract_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["authoring_source_is_build123d"] is True
    assert gate["checks"]["mesh_route_matches_expected"] is True
    assert gate["checks"]["required_source_groups_present"] is True
    assert gate["checks"]["shape_length_unit_expected_ok"] is True
    assert gate["checks"]["shape_area_unit_expected_ok"] is True
    assert gate["checks"]["shape_volume_unit_expected_ok"] is True
    assert gate["checks"]["required_shape_metadata_present"] is True
    assert gate["required_metadata_fields"] == [
        "recipe_id",
        "cad_kernel",
        "cad_kernel_version",
        "script_path",
        "export_id",
    ]
    assert gate["missing_required_metadata_by_shape"] == {}
    assert gate["units"] == {"length": ["m"], "area": ["m^2"], "volume": ["m^3"]}
    assert [item["matched"] for item in gate["matched_source_groups"]] == [
        "coreform_cubit",
        "external_cad",
    ]

    tet_row = dict(row)
    tet_row["mesh_route"] = "netgen_tri_tet_path"
    bad_route = shape_cad_route_source_contract_gate([tet_row], summary)
    assert bad_route["status"] == "needs_attention"
    assert bad_route["checks"]["mesh_route_matches_expected"] is False
    assert bad_route["checks"]["disallowed_routes_absent"] is False

    missing_external = shape_mass_property_crosscheck_summary(
        [row],
        {"coreform_cubit": [dict(row)]},
        rtol=1.0e-12,
        bbox_atol=1.0e-12,
    )
    bad_source = shape_cad_route_source_contract_gate([row], missing_external)
    assert bad_source["status"] == "needs_attention"
    assert bad_source["checks"]["required_source_groups_present"] is False
    assert bad_source["missing_source_groups"] == [["external_cad", "cst_import"]]

    wrong_unit = dict(row)
    wrong_unit["volume_unit"] = "mm^3"
    bad_unit = shape_cad_route_source_contract_gate(
        [wrong_unit],
        summary,
        expected_volume_unit="m^3",
    )
    assert bad_unit["status"] == "needs_attention"
    assert bad_unit["checks"]["shape_volume_unit_expected_ok"] is False
    assert "expected CAD units" in " ".join(bad_unit["issues"])

    missing_metadata = dict(row)
    missing_metadata.pop("recipe_id")
    bad_metadata = shape_cad_route_source_contract_gate(
        [missing_metadata],
        summary,
        required_metadata_fields=("recipe_id", "cad_kernel"),
    )
    assert bad_metadata["status"] == "needs_attention"
    assert bad_metadata["checks"]["required_shape_metadata_present"] is False
    assert bad_metadata["missing_required_metadata_by_shape"] == {"hex_block": ["recipe_id"]}
    assert "required CAD provenance metadata" in " ".join(bad_metadata["issues"])


def test_cst_cad_volume_export_manifest_gate_normalizes_rows_for_crosscheck():
    manifest = {
        "source_tool": "CST Studio Suite",
        "project_id": "cad_volume_widget",
        "run_id": "run_geom_001",
        "export_id": "cst_volume_export_D",
        "geometry_id": "widget_revA",
        "volume_unit": "mm^3",
        "volume_rows": [
            {"solid_name": "stator", "volume": 100.0e9, "status": "exported"},
            {"solid_name": "coil", "volume": 25.0e9, "status": "exported"},
        ],
    }
    gate = cst_cad_volume_export_manifest_gate(
        manifest,
        expected_geometry_id="widget_revA",
        expected_export_id="cst_volume_export_D",
        required_shape_names=("stator", "coil"),
    )

    assert gate["policy"] == "cst_cad_volume_export_manifest_gate"
    assert gate["status"] == "ok"
    assert gate["normalized_rows"][0]["volume_m3"] == pytest.approx(100.0)
    assert gate["normalized_rows"][1]["volume_m3"] == pytest.approx(25.0)
    assert all(gate["checks"].values())

    summary = shape_volume_crosscheck_summary(
        [{"name": "stator", "volume": 100.0}, {"name": "coil", "volume": 25.0}],
        {
            "cst_import": gate["normalized_rows"],
            "cubit": [{"name": "stator", "volume": 100.0}, {"name": "coil", "volume": 25.0}],
        },
        rtol=1.0e-12,
    )
    coverage = shape_volume_crosscheck_source_coverage_gate(summary)
    assert summary["status"] == "ok"
    assert coverage["status"] == "ok"

    duplicate = dict(manifest)
    duplicate["volume_rows"] = [dict(row) for row in manifest["volume_rows"]]
    duplicate["volume_rows"][1]["solid_name"] = "stator"
    duplicate_gate = cst_cad_volume_export_manifest_gate(duplicate)
    assert duplicate_gate["status"] == "needs_attention"
    assert duplicate_gate["checks"]["shape_names_unique"] is False

    unknown_unit = dict(manifest)
    unknown_unit["volume_unit"] = "litre"
    unit_gate = cst_cad_volume_export_manifest_gate(unknown_unit)
    assert unit_gate["status"] == "needs_attention"
    assert unit_gate["checks"]["volume_units_known"] is False


def test_shape_name_identity_gate_rejects_extra_missing_and_duplicate_shapes():
    reference = [
        {"name": "stator", "volume": 100.0, "area": 220.0},
        {"name": "coil", "volume": 25.0, "area": 60.0},
    ]
    measured = [
        {"name": "stator", "volume": 100.0, "area": 220.0},
        {"name": "coil", "volume": 25.0, "area": 60.0},
    ]

    ok = shape_name_identity_gate(reference, measured, measured_label="cubit")
    assert ok["status"] == "ok"
    assert ok["policy"] == "build123d_cad_roundtrip_named_shape_identity_gate"
    assert ok["checks"]["same_name_multiset"] is True
    assert ok["version_note"].startswith("Use this before volume/area/bbox gates")

    bad = shape_name_identity_gate(
        reference,
        [
            {"name": "stator", "volume": 100.0, "area": 220.0},
            {"name": "stator", "volume": 100.0, "area": 220.0},
            {"name": "bolt", "volume": 1.0, "area": 2.0},
        ],
        measured_label="cst_import",
    )
    assert bad["status"] == "needs_attention"
    assert bad["measured_label"] == "cst_import"
    assert bad["missing_names"] == ["coil"]
    assert bad["extra_names"] == ["bolt", "stator"]
    assert bad["duplicate_measured_names"] == ["stator"]
    assert bad["checks"]["measured_names_unique"] is False

    summary = shape_mass_property_crosscheck_summary(
        reference,
        {"cst_import": [
            {"name": "stator", "volume": 100.0, "area": 220.0},
            {"name": "coil", "volume": 25.0, "area": 60.0},
            {"name": "bolt", "volume": 1.0, "area": 2.0},
        ]},
        rtol=1.0e-12,
    )
    assert summary["status"] == "needs_attention"
    assert summary["checks"]["all_sources_present_and_within_tolerance"] is True
    assert summary["checks"]["all_sources_preserve_named_shape_identity"] is False
    assert "missing, extra, duplicate, or unnamed shapes" in summary["issues"][0]
    assert summary["comparison_sets"][0]["name_identity_gate"]["extra_names"] == ["bolt"]


def test_shape_role_metadata_gate_requires_solver_handoff_semantics():
    rows = [
        {"name": "iron_core", "role": "magnetic_core", "material": "electrical_steel"},
        {"name": "phase_a_coil", "role": "conductor", "material_name": "copper"},
        {"name": "air_box", "solver_role": "air_region", "mat": "air"},
    ]

    ok = shape_role_metadata_gate(
        rows,
        required_names=["iron_core", "phase_a_coil"],
        required_roles=["magnetic_core", "conductor"],
        required_materials=["electrical_steel", "copper"],
        source_label="slot123_build123d",
    )

    assert ok["policy"] == "build123d_solver_handoff_role_material_metadata_gate"
    assert ok["status"] == "ok"
    assert ok["source_label"] == "slot123_build123d"
    assert ok["checks"]["required_roles_present"] is True
    assert ok["checks"]["required_materials_present"] is True
    assert ok["materials"] == ["air", "copper", "electrical_steel"]
    assert ok["version_note"].startswith("Run this after shape_name_identity_gate")

    bad = shape_role_metadata_gate(
        [
            {"name": "iron_core", "role": "magnetic_core", "material": "electrical_steel"},
            {"name": "iron_core", "role": "fixture", "material": "steel"},
            {"name": "phase_a_coil", "role": "conductor"},
            {"name": "air_box", "material": "air"},
        ],
        required_names=["iron_core", "phase_a_coil", "air_box", "shaft"],
        required_roles=["magnetic_core", "conductor", "air_region"],
        required_materials=["electrical_steel", "copper", "air"],
    )

    assert bad["status"] == "needs_attention"
    assert bad["duplicate_names"] == ["iron_core"]
    assert bad["rows_missing_material"] == ["phase_a_coil"]
    assert bad["rows_missing_role"] == ["air_box"]
    assert bad["missing_required_names"] == ["shaft"]
    assert bad["missing_required_roles"] == ["air_region"]
    assert bad["missing_required_materials"] == ["copper"]


def test_shape_transition_role_metadata_gate_preserves_hex_tet_handoff_intent():
    hex_body = (Pos(-1.25, 0, 0) * Box(1.0, 1.0, 1.0)).solid()
    hex_body.label = "hex_core"
    transition = Box(0.5, 1.0, 1.0).solid()
    transition.label = "pyramid_transition_envelope"
    tet_body = (Pos(1.25, 0, 0) * Box(1.0, 1.0, 1.0)).solid()
    tet_body.label = "tet_region"

    rows = shape_measurement_rows(assembly(hex_body, transition, tet_body, label="handoff"))
    by_name = {row["name"]: row for row in rows}
    by_name["hex_core"].update({"role": "hex_region", "material": "core_steel"})
    by_name["pyramid_transition_envelope"].update({
        "role": "mesh_transition",
        "material": "transition_air",
        "transition_kind": "pyramid",
        "connects_roles": ["hex_region", "tet_region"],
        "expected_surface_kinds": ["quad", "triangle"],
        "expected_interface_roles": ["hex_to_transition", "transition_to_tet"],
    })
    by_name["tet_region"].update({"role": "tet_region", "material": "air"})

    ok = shape_transition_role_metadata_gate(
        rows,
        required_surface_kinds=("quad", "triangle"),
        required_interface_roles=("hex_to_transition", "transition_to_tet"),
        source_label="slot131_build123d",
    )

    assert ok["policy"] == "build123d_hex_tet_transition_role_metadata_gate"
    assert ok["status"] == "ok"
    assert ok["source_label"] == "slot131_build123d"
    assert ok["checks"]["required_roles_present"] is True
    assert ok["checks"]["transition_kind_matches"] is True
    assert ok["checks"]["transition_connects_required_roles"] is True
    assert ok["roles"] == ["hex_region", "mesh_transition", "tet_region"]
    assert ok["transition_kinds"] == ["pyramid"]
    assert ok["connected_roles"] == ["hex_region", "tet_region"]
    assert ok["checks"]["required_surface_kinds_present"] is True
    assert ok["surface_kinds"] == ["quad", "triangle"]
    assert ok["checks"]["interface_roles_recorded"] is True
    assert ok["checks"]["required_interface_roles_present"] is True
    assert ok["interface_roles"] == ["hex_to_transition", "transition_to_tet"]

    wrong_kind = [dict(row) for row in rows]
    for row in wrong_kind:
        if row["name"] == "pyramid_transition_envelope":
            row["transition_kind"] = "unknown"
    bad_kind = shape_transition_role_metadata_gate(wrong_kind)
    assert bad_kind["status"] == "needs_attention"
    assert bad_kind["checks"]["transition_kind_matches"] is False

    missing_connection = [dict(row) for row in rows]
    for row in missing_connection:
        if row["name"] == "pyramid_transition_envelope":
            row["connects_roles"] = ["hex_region"]
    bad_connection = shape_transition_role_metadata_gate(missing_connection)
    assert bad_connection["status"] == "needs_attention"
    assert bad_connection["checks"]["transition_connects_required_roles"] is False

    missing_surface_family = [dict(row) for row in rows]
    for row in missing_surface_family:
        if row["name"] == "pyramid_transition_envelope":
            row["expected_surface_kinds"] = ["triangle"]
    bad_surface_family = shape_transition_role_metadata_gate(
        missing_surface_family,
        required_surface_kinds=("quad", "triangle"),
    )
    assert bad_surface_family["status"] == "needs_attention"
    assert bad_surface_family["checks"]["required_surface_kinds_present"] is False

    missing_interface_role = [dict(row) for row in rows]
    for row in missing_interface_role:
        if row["name"] == "pyramid_transition_envelope":
            row["expected_interface_roles"] = ["hex_to_transition"]
    bad_interface_role = shape_transition_role_metadata_gate(
        missing_interface_role,
        required_surface_kinds=("quad", "triangle"),
        required_interface_roles=("hex_to_transition", "transition_to_tet"),
    )
    assert bad_interface_role["status"] == "needs_attention"
    assert bad_interface_role["checks"]["required_interface_roles_present"] is False
    assert bad_interface_role["missing_required_interface_roles"] == ["transition_to_tet"]


def test_shape_transition_role_metadata_gate_records_downstream_material_labels():
    rows = [
        {
            "name": "slot275_hex_body",
            "role": "hex_region",
            "material": "core_steel",
            "downstream_material_name": "hex_core",
            "volume": 1.0,
        },
        {
            "name": "slot275_transition_envelope",
            "role": "mesh_transition",
            "material": "transition_air",
            "transition_kind": "pyramid",
            "connects_roles": ["hex_region", "tet_region"],
            "expected_surface_kinds": ["quad", "triangle"],
            "expected_interface_roles": ["hex_to_transition", "transition_to_tet"],
            "downstream_material_name": "pyramid_transition",
            "volume": 0.25,
        },
        {
            "name": "slot275_tet_body",
            "role": "tet_region",
            "material": "air",
            "downstream_material_name": "tet_region",
            "volume": 1.0,
        },
    ]

    gate = shape_transition_role_metadata_gate(
        rows,
        required_surface_kinds=("quad", "triangle"),
        required_interface_roles=("hex_to_transition", "transition_to_tet"),
        expected_downstream_material_names=("hex_core", "pyramid_transition", "tet_region"),
        allowed_zero_downstream_material_names=("pyramid_transition",),
        source_label="slot275_build123d",
    )

    assert gate["status"] == "ok"
    assert gate["checks"]["downstream_material_names_recorded"] is True
    assert gate["checks"]["expected_downstream_material_names_present"] is True
    assert gate["checks"]["allowed_zero_downstream_material_names_declared"] is True
    assert gate["downstream_material_names"] == ["hex_core", "pyramid_transition", "tet_region"]
    assert gate["allowed_zero_downstream_material_names"] == ["pyramid_transition"]

    missing_downstream_name = [dict(row) for row in rows]
    missing_downstream_name[1].pop("downstream_material_name")
    bad = shape_transition_role_metadata_gate(
        missing_downstream_name,
        required_surface_kinds=("quad", "triangle"),
        required_interface_roles=("hex_to_transition", "transition_to_tet"),
        expected_downstream_material_names=("hex_core", "pyramid_transition", "tet_region"),
        allowed_zero_downstream_material_names=("pyramid_transition",),
    )

    assert bad["status"] == "needs_attention"
    assert bad["checks"]["downstream_material_names_recorded"] is False
    assert bad["checks"]["expected_downstream_material_names_present"] is False
    assert bad["missing_expected_downstream_material_names"] == ["pyramid_transition"]
    assert bad["rows_missing_downstream_material_name"] == ["slot275_transition_envelope"]


def test_shape_cubit_meshing_scheme_intent_gate_binds_cad_roles_to_downstream_trace():
    rows = [
        {
            "name": "slot291_hex_body",
            "role": "hex_region",
            "material": "core_steel",
            "mesh_route": "cubit_hex_or_mixed_path",
            "expected_cubit_scheme": "map",
            "downstream_meshing_trace_id": "slot291_mixed_scheme_trace",
            "expected_cubit_command_fragments": ["imprint all", "merge all", "export netgen"],
            "expected_cubit_export_order": 2,
            "volume": 1.0,
        },
        {
            "name": "slot291_transition_envelope",
            "role": "mesh_transition",
            "material": "transition_air",
            "mesh_route": "cubit_hex_or_mixed_path",
            "expected_cubit_scheme": "tetmesh",
            "downstream_meshing_trace_id": "slot291_mixed_scheme_trace",
            "expected_cubit_command_fragments": ["imprint all", "merge all", "export netgen"],
            "expected_cubit_export_order": 2,
            "volume": 0.25,
        },
        {
            "name": "slot291_tet_body",
            "role": "tet_region",
            "material": "air",
            "mesh_route": "cubit_hex_or_mixed_path",
            "expected_cubit_scheme": "tetmesh",
            "downstream_meshing_trace_id": "slot291_mixed_scheme_trace",
            "expected_cubit_command_fragments": ["imprint all", "merge all", "export netgen"],
            "expected_cubit_export_order": 2,
            "volume": 1.0,
        },
    ]
    scheme_trace = cubit_meshing_scheme_trace_gate(
        {
            "trace_id": "slot291_mixed_scheme_trace",
            "command_digest": "sha256:slot291-imprint-merge-map-tet-export",
            "commands": [
                "imprint all",
                "merge all",
                "volume 1 scheme map",
                "volume 2 scheme tetmesh",
                "volume 3 scheme tetmesh",
                "export netgen \"slot291_mixed.vol\" order 2 overwrite",
            ],
            "volume_schemes": {"1": "map", "2": "tetmesh", "3": "tetmesh"},
            "export_order": 2,
            "export_output_artifact_id": "slot291_mixed_vol_v1",
            "export_output_digest": "sha256:slot291-mixed-vol",
            "export_output_path": r"artifacts/cubit/slot291_mixed.vol",
        },
        expected_trace_id="slot291_mixed_scheme_trace",
        expected_command_digest="sha256:slot291-imprint-merge-map-tet-export",
        expected_volume_schemes={"1": "map", "2": "tetmesh", "3": "tetmesh"},
        expected_export_order=2,
        expected_export_output_artifact_id="slot291_mixed_vol_v1",
        expected_export_output_digest="sha256:slot291-mixed-vol",
        expected_export_output_path=r"artifacts/cubit/slot291_mixed.vol",
        require_export_output_artifact=True,
    )

    gate = shape_cubit_meshing_scheme_intent_gate(
        rows,
        scheme_trace_gate=scheme_trace,
        expected_scheme_by_role={
            "hex_region": "map",
            "mesh_transition": "tetmesh",
            "tet_region": "tetmesh",
        },
        expected_trace_id="slot291_mixed_scheme_trace",
        expected_export_order=2,
        expected_export_output_artifact_id="slot291_mixed_vol_v1",
        expected_export_output_digest="sha256:slot291-mixed-vol",
        expected_export_output_path=r"artifacts/cubit/slot291_mixed.vol",
        require_downstream_export_output_artifact=True,
        source_label="slot291_build123d",
    )

    assert gate["policy"] == "build123d_cubit_meshing_scheme_intent_gate"
    assert gate["status"] == "ok"
    assert gate["role_scheme_intent"] == {
        "hex_region": "map",
        "mesh_transition": "tetmesh",
        "tet_region": "tetmesh",
    }
    assert gate["checks"]["downstream_scheme_trace_gate_ok"] is True
    assert gate["checks"]["downstream_trace_id_matches"] is True
    assert gate["checks"]["required_command_fragments_present"] is True
    assert gate["checks"]["downstream_export_output_artifact_id_matches"] is True
    assert gate["checks"]["downstream_export_output_digest_matches"] is True
    assert gate["checks"]["downstream_export_output_path_matches"] is True

    stale_scheme = [dict(row) for row in rows]
    stale_scheme[0]["expected_cubit_scheme"] = "tetmesh"
    bad_scheme = shape_cubit_meshing_scheme_intent_gate(
        stale_scheme,
        scheme_trace_gate=scheme_trace,
        expected_scheme_by_role={
            "hex_region": "map",
            "mesh_transition": "tetmesh",
            "tet_region": "tetmesh",
        },
        expected_trace_id="slot291_mixed_scheme_trace",
        expected_export_order=2,
    )
    assert bad_scheme["status"] == "needs_attention"
    assert bad_scheme["checks"]["expected_scheme_by_role_matches"] is False
    assert bad_scheme["missing_expected_scheme_roles"] == ["hex_region"]

    missing_export_fragment = [dict(row) for row in rows]
    for row in missing_export_fragment:
        row["expected_cubit_command_fragments"] = ["imprint all", "merge all"]
    bad_fragment = shape_cubit_meshing_scheme_intent_gate(
        missing_export_fragment,
        scheme_trace_gate=scheme_trace,
        expected_trace_id="slot291_mixed_scheme_trace",
        expected_export_order=2,
    )
    assert bad_fragment["status"] == "needs_attention"
    assert bad_fragment["checks"]["required_command_fragments_present"] is False

    stale_trace = dict(scheme_trace)
    stale_trace["trace_id"] = "slot290_old_scheme_trace"
    bad_trace = shape_cubit_meshing_scheme_intent_gate(
        rows,
        scheme_trace_gate=stale_trace,
        expected_trace_id="slot291_mixed_scheme_trace",
        expected_export_order=2,
    )
    assert bad_trace["status"] == "needs_attention"
    assert bad_trace["checks"]["downstream_trace_id_matches"] is False

    stale_digest_trace = dict(scheme_trace)
    stale_digest_trace["export_output_digest"] = "sha256:old-mixed-vol"
    bad_output = shape_cubit_meshing_scheme_intent_gate(
        rows,
        scheme_trace_gate=stale_digest_trace,
        expected_trace_id="slot291_mixed_scheme_trace",
        expected_export_order=2,
        expected_export_output_digest="sha256:slot291-mixed-vol",
        require_downstream_export_output_artifact=True,
    )
    assert bad_output["status"] == "needs_attention"
    assert bad_output["checks"]["downstream_export_output_digest_matches"] is False


def test_build123d_l_bracket_slot_volume_crosscheck_accepts_cubit_roundtrip():
    reference = [{"name": "l_bracket_two_holes", "volume": 2.898212398023691}]
    measured = {"cubit": [{"name": "l_bracket_two_holes", "volume": 2.8982123980236905}]}

    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert summary["status"] == "ok"
    assert summary["ok_for_cad_roundtrip_volume"] is True
    assert summary["max_volume_rel_error"] == pytest.approx(1.5322866265870984e-16)


def test_mounting_plate_boss_reference_matches_slot19_volume_gate():
    reference = [mounting_plate_boss_reference_row(6.0, 4.0, 0.5, 0.8, 0.6, 0.25, 0.18, 2.2, 1.2)]
    measured = {
        "build123d": [{"name": "mounting_plate_boss_five_holes", "volume": 12.78681188009156}],
        "cubit": [{"name": "mounting_plate_boss_five_holes", "volume": 12.786811880091562}],
    }

    row = reference[0]
    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert row["volume"] == pytest.approx(12.786811880091562)
    assert row["terms"]["base"] == pytest.approx(12.0)
    assert row["terms"]["central_hole"] < 0.0
    assert row["policy"] == "analytic_mounting_plate_boss_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([6.0, 4.0, 1.1])
    assert summary["status"] == "ok"
    assert summary["sources"] == ["build123d", "cubit"]
    assert summary["max_volume_rel_error"] < 2.0e-16


def test_stepped_spacer_slot27_uses_external_kernel_roundtrip_tolerance():
    reference = [{"name": "stepped_spacer_four_bolt_holes", "volume": 4.8536349860900865}]
    measured = {
        "build123d": [{"name": "stepped_spacer_four_bolt_holes", "volume": 4.8536349860900865}],
        "cubit": [{"name": "stepped_spacer_four_bolt_holes", "volume": 4.853663401216389}],
    }

    strict_summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-6)
    roundtrip_summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-5)

    assert strict_summary["status"] == "needs_attention"
    assert strict_summary["comparison_sets"][0]["status"] == "ok"
    assert strict_summary["comparison_sets"][1]["status"] == "needs_attention"
    assert roundtrip_summary["status"] == "ok"
    assert roundtrip_summary["ok_for_cad_roundtrip_volume"] is True
    assert roundtrip_summary["max_volume_rel_error"] == pytest.approx(5.854366888206992e-6)


def test_keyed_terminal_plate_slot35_matches_volume_roundtrip_gate():
    reference = [keyed_terminal_plate_reference_row(
        8.0, 3.0, 0.4,
        0.55, 0.5, 2.5, 0.22,
        1.8, 0.8,
        0.7, 1.0,
        0.15, 3.4, -0.9,
    )]
    measured = {
        "build123d": [{"name": "keyed_terminal_plate_two_bosses", "volume": 9.364087557965554}],
        "cubit": [{"name": "keyed_terminal_plate_two_bosses", "volume": 9.364087557965554}],
    }

    row = reference[0]
    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert row["volume"] == pytest.approx(9.364087557965556)
    assert row["terms"]["base"] == pytest.approx(9.6)
    assert row["terms"]["rectangular_window"] < 0.0
    assert row["terms"]["edge_key_slot"] < 0.0
    assert row["policy"] == "analytic_keyed_terminal_plate_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([8.0, 3.0, 0.9])
    assert summary["status"] == "ok"
    assert summary["sources"] == ["build123d", "cubit"]
    assert summary["max_volume_rel_error"] < 2.0e-16


def test_flanged_sleeve_slot43_matches_volume_roundtrip_gate():
    reference = [flanged_sleeve_reference_row(
        1.8, 0.35,
        0.75, 1.1,
        0.28,
        1.25, 0.12,
        bolt_count=4,
    )]
    measured = {
        "build123d": [{"name": "flanged_sleeve_four_bolt_holes", "volume": 5.085955762823052}],
        "cubit": [{"name": "flanged_sleeve_four_bolt_holes", "volume": 5.085958320184421}],
    }

    row = reference[0]
    strict_summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-7)
    roundtrip_summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-6)

    assert row["volume"] == pytest.approx(5.085955762823052)
    assert row["terms"]["flange_annulus"] == pytest.approx(3.476360766756322)
    assert row["terms"]["hub_annulus"] == pytest.approx(1.6729295039631007)
    assert row["terms"]["bolt_holes"] < 0.0
    assert row["bolt_count"] == 4
    assert row["policy"] == "analytic_flanged_sleeve_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([3.6, 3.6, 1.45])
    assert strict_summary["status"] == "needs_attention"
    assert roundtrip_summary["status"] == "ok"
    assert roundtrip_summary["sources"] == ["build123d", "cubit"]
    assert roundtrip_summary["max_volume_rel_error"] == pytest.approx(5.028278267134254e-7)


def test_coax_annular_sleeve_slot91_uses_physical_rc_geometry_volume_gate():
    reference = [coax_annular_sleeve_reference_row(0.9, 2.1, 3.4)]
    sleeve = tube(0.9, 2.1, 3.4, label="coax_annular_sleeve")
    measured = {
        "build123d": [{"name": "coax_annular_sleeve", "volume": _vol(sleeve)}],
        "cubit": [{"name": "coax_annular_sleeve", "volume": 38.451673891926546}],
    }

    row = reference[0]
    strict_summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-5)
    roundtrip_summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-4)

    assert row["volume"] == pytest.approx(38.453094079939064)
    assert _vol(sleeve) == pytest.approx(row["volume"])
    assert row["area"] == pytest.approx(86.70795723907828)
    assert row["terms"]["inner_void"] < 0.0
    assert row["parameters"] == {"inner_radius": 0.9, "outer_radius": 2.1, "height": 3.4}
    assert row["policy"] == "analytic_coax_annular_sleeve_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([4.2, 4.2, 3.4])
    assert strict_summary["status"] == "needs_attention"
    assert roundtrip_summary["status"] == "ok"
    assert roundtrip_summary["sources"] == ["build123d", "cubit"]
    assert roundtrip_summary["max_volume_rel_error"] == pytest.approx(3.693299710979559e-5)


def test_ribbed_busbar_heat_sink_reference_matches_volume_roundtrip_gate():
    reference = [ribbed_busbar_heat_sink_reference_row(
        8.0, 3.2, 0.35,
        4, 0.2, 0.75, 0.35,
        0.18, 3.2, 1.25,
    )]
    measured = {
        "build123d": [{"name": "ribbed_busbar_heat_sink_four_holes", "volume": reference[0]["volume"]}],
        "cubit": [{"name": "ribbed_busbar_heat_sink_four_holes", "volume": reference[0]["volume"]}],
    }

    row = reference[0]
    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert row["volume"] == pytest.approx(13.617497357233166)
    assert row["terms"]["base"] == pytest.approx(8.96)
    assert row["terms"]["straight_ribs"] == pytest.approx(4.8)
    assert row["terms"]["four_base_holes"] < 0.0
    assert row["fin_count"] == 4
    assert row["clearances"]["hole_to_fin_band"] > 0.0
    assert row["policy"] == "analytic_ribbed_busbar_heat_sink_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([8.0, 3.2, 1.1])
    assert summary["status"] == "ok"
    assert summary["sources"] == ["build123d", "cubit"]


def test_three_phase_busbar_snubber_plate_slot59_matches_volume_roundtrip_gate():
    reference = [three_phase_busbar_snubber_plate_reference_row(
        9.0, 3.6, 0.35,
        3, 1.0, 0.6, 0.45, 2.4,
        2, 1.2, 0.45, 0.25, 2.0,
        0.16, 3.9, 1.35,
        phase_tab_y0=0.15,
        snubber_pad_y0=-1.1,
    )]
    measured = {
        "build123d": [{"name": "three_phase_busbar_snubber_plate", "volume": 12.307405319295338}],
        "cubit": [{"name": "three_phase_busbar_snubber_plate", "volume": 12.307405319295343}],
    }

    row = reference[0]
    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert row["volume"] == pytest.approx(12.30740531929534)
    assert row["terms"]["base"] == pytest.approx(11.34)
    assert row["terms"]["three_phase_tabs"] == pytest.approx(0.81)
    assert row["terms"]["two_snubber_pads"] == pytest.approx(0.27)
    assert row["terms"]["four_mount_holes"] < 0.0
    assert row["counts"] == {"phase_tabs": 3, "snubber_pads": 2, "mount_holes": 4}
    assert row["clearances"]["mount_hole_to_y_edge"] == pytest.approx(0.29)
    assert row["clearances"]["phase_tab_gap"] == pytest.approx(1.4)
    assert row["policy"] == "analytic_three_phase_busbar_snubber_plate_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([9.0, 3.6, 0.8])
    assert summary["status"] == "ok"
    assert summary["sources"] == ["build123d", "cubit"]
    assert summary["max_volume_rel_error"] < 5.0e-16


def test_rcd_snubber_heat_spreader_slot75_matches_volume_roundtrip_gate():
    reference = [rcd_snubber_heat_spreader_reference_row(
        10.0, 4.0, 0.32,
        5, 8.0, 0.12, 0.45, 0.45,
        2, 1.25, 0.70, 0.38, 2.4,
        0.16, 4.2, 1.55,
        snubber_pad_y0=-1.55,
    )]
    measured = {
        "build123d": [{"name": "rcd_snubber_heat_spreader", "volume": 15.522056291927165}],
        "cubit": [{"name": "rcd_snubber_heat_spreader", "volume": 15.522056291927173}],
    }

    row = reference[0]
    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert row["volume"] == pytest.approx(15.52205629192717)
    assert row["terms"]["base"] == pytest.approx(12.8)
    assert row["terms"]["straight_ribs"] == pytest.approx(2.16)
    assert row["terms"]["snubber_pads"] == pytest.approx(0.665)
    assert row["terms"]["four_mount_holes"] < 0.0
    assert row["counts"] == {"ribs": 5, "snubber_pads": 2, "mount_holes": 4}
    assert row["clearances"]["snubber_pad_to_rib_band"] == pytest.approx(0.24)
    assert row["clearances"]["mount_hole_to_rib_band"] == pytest.approx(0.43)
    assert row["parameters"]["rib_pitch"] == pytest.approx(0.45)
    assert row["parameters"]["snubber_pad_y0"] == pytest.approx(-1.55)
    assert row["policy"] == "analytic_rcd_snubber_heat_spreader_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([10.0, 4.0, 0.77])
    assert summary["status"] == "ok"
    assert summary["sources"] == ["build123d", "cubit"]
    assert summary["max_volume_rel_error"] < 6.0e-16


def test_rcd_snubber_capacitance_sweep_keeps_cad_variant_provenance():
    rows = rcd_snubber_capacitance_sweep_rows(
        [0.047, 0.10, 0.22],
        [1.05, 1.25, 1.55],
    )

    assert [row["capacitance_uF"] for row in rows] == pytest.approx([0.047, 0.10, 0.22])
    assert [row["snubber_pad_x"] for row in rows] == pytest.approx([1.05, 1.25, 1.55])
    assert rows[0]["name"] == "rcd_snubber_heat_spreader_0.047uF"
    assert rows[-1]["terms"]["snubber_pads"] > rows[0]["terms"]["snubber_pads"]
    assert rows[-1]["volume"] > rows[0]["volume"]
    assert all(row["design_table_role"].startswith("RCD snubber capacitance") for row in rows)

    summary = shape_parameter_sweep_summary(
        rows,
        "capacitance_uF",
        metric_keys=("volume", "snubber_pad_volume"),
    )
    metrics = {row["metric"]: row for row in summary["metric_rows"]}

    assert summary["status"] == "ok"
    assert summary["parameter_strictly_increasing"] is True
    assert metrics["volume"]["monotonic_non_decreasing"] is True
    assert metrics["snubber_pad_volume"]["monotonic_non_decreasing"] is True
    assert metrics["snubber_pad_volume"]["delta_first_to_last"] == pytest.approx(
        rows[-1]["snubber_pad_volume"] - rows[0]["snubber_pad_volume"]
    )


def test_thermal_robin_cooling_plate_slot83_matches_volume_roundtrip_gate():
    reference = [thermal_robin_cooling_plate_reference_row(
        9.0, 4.5, 0.30,
        4, 7.2, 0.14, 0.55, 0.42,
        2, 1.6, 0.80, 0.32, 2.6,
        0.14, 3.8, 1.8,
        fin_y0=0.45,
        device_pad_y0=-1.45,
    )]
    measured = {
        "build123d": [{"name": "thermal_robin_cooling_plate", "volume": 15.112909740787572}],
        "cubit": [{"name": "thermal_robin_cooling_plate", "volume": 15.112909740787561}],
    }

    row = reference[0]
    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert row["volume"] == pytest.approx(15.11290974078757)
    assert row["terms"]["base"] == pytest.approx(12.15)
    assert row["terms"]["straight_cooling_fins"] == pytest.approx(2.2176)
    assert row["terms"]["device_pads"] == pytest.approx(0.8192)
    assert row["terms"]["four_mount_holes"] < 0.0
    assert row["counts"] == {"fins": 4, "device_pads": 2, "mount_holes": 4}
    assert row["clearances"]["device_pad_to_fin_band"] == pytest.approx(0.8)
    assert row["clearances"]["mount_hole_to_y_edge"] == pytest.approx(0.31)
    assert row["parameters"]["fin_y0"] == pytest.approx(0.45)
    assert row["parameters"]["device_pad_y0"] == pytest.approx(-1.45)
    assert row["policy"] == "analytic_thermal_robin_cooling_plate_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([9.0, 4.5, 0.85])
    assert summary["status"] == "ok"
    assert summary["sources"] == ["build123d", "cubit"]
    assert summary["max_volume_rel_error"] < 8.0e-16


def test_v_type_ipm_rotor_coupon_reference_matches_volume_roundtrip_gate():
    reference = [v_type_ipm_rotor_coupon_reference_row(
        8.0, 5.0, 0.35,
        2.2, 0.35, 28.0,
        1.65, 0.55,
        0.45,
    )]
    measured = {
        "build123d": [{"name": "v_type_ipm_rotor_coupon", "volume": reference[0]["volume"]}],
        "cubit": [{"name": "v_type_ipm_rotor_coupon", "volume": reference[0]["volume"]}],
    }

    row = reference[0]
    summary = shape_volume_crosscheck_summary(reference, measured, rtol=1.0e-12)

    assert row["volume"] == pytest.approx(13.238339620676824)
    assert row["terms"]["coupon"] == pytest.approx(14.0)
    assert row["terms"]["two_v_magnet_pockets"] == pytest.approx(-0.539)
    assert row["terms"]["central_bore"] == pytest.approx(-0.2226603793231766)
    assert row["counts"] == {"magnet_pockets": 2, "bore": 1}
    assert row["clearances"]["pocket_to_bore_x"] == pytest.approx(0.1466001243676493)
    assert row["clearances"]["mirrored_pocket_gap"] == pytest.approx(1.1932002487352986)
    assert row["parameters"]["magnet_slot_angle_deg"] == pytest.approx(28.0)
    assert row["policy"] == "analytic_v_type_ipm_rotor_coupon_volume_reference"
    assert row["bounding_box"]["size"] == pytest.approx([8.0, 5.0, 0.35])
    assert summary["status"] == "ok"
    assert summary["sources"] == ["build123d", "cubit"]


def test_build123d_volume_crosscheck_mcp_tool_dispatches_json():
    from radia_mcp.build123d.server import build123d_volume_crosscheck

    reference = [{"name": "box", "volume": 24.0}]
    measured = {"cubit": [{"name": "box", "volume": 24.0}]}
    payload = json.loads(build123d_volume_crosscheck(json.dumps(reference), json.dumps(measured)))

    assert payload["status"] == "ok"
    assert payload["n_sources"] == 1
    assert payload["comparison_sets"][0]["source"] == "cubit"
    assert payload["comparison_sets"][0]["rows"][0]["passed"] is True


def test_perforated_prism_roundtrip_requires_volume_and_boundary_topology():
    gate = shape_perforated_prism_roundtrip_gate(
        reference_volume=168472.3918002237,
        imported_volume=168472.39180031797,
        hole_count=625,
        hole_side_count=6,
        imported_surface_count=3756,
        imported_body_count=1,
        volume_rtol=1.0e-12,
    )
    assert gate["status"] == "ok"
    assert gate["expected_surface_count"] == 3756
    assert gate["volume_relative_error"] < 1.0e-12

    same_volume_missing_holes = shape_perforated_prism_roundtrip_gate(
        reference_volume=168472.3918002237,
        imported_volume=168472.3918002237,
        hole_count=625,
        hole_side_count=6,
        imported_surface_count=6,
        imported_body_count=1,
    )
    assert same_volume_missing_holes["status"] == "needs_attention"
    assert same_volume_missing_holes["checks"]["volume_agrees"] is True
    assert same_volume_missing_holes["checks"]["surface_topology_preserved"] is False


def test_perforated_prism_roundtrip_mcp_tool_dispatches_json():
    from radia_mcp.build123d.server import build123d_perforated_prism_roundtrip_gate

    payload = json.loads(build123d_perforated_prism_roundtrip_gate(
        168472.3918002237,
        168472.39180031797,
        625,
        6,
        3756,
        volume_rtol=1.0e-12,
    ))
    assert payload["status"] == "ok"
    assert payload["checks"]["surface_topology_preserved"] is True


def test_build123d_volume_crosscheck_with_units_normalizes_and_rejects_implicit_units():
    from radia_mcp.build123d.server import build123d_volume_crosscheck_with_units

    reference = [{"name": "box", "volume": 6000.0, "volume_unit": "mm^3"}]
    measured = {
        "external_cm": [{"name": "box", "volume": 6.0, "volume_unit": "cm^3"}],
        "external_m": [{"name": "box", "volume": 6.0e-6, "volume_unit": "m^3"}],
    }
    payload = json.loads(build123d_volume_crosscheck_with_units(
        json.dumps(reference),
        json.dumps(measured),
        rtol=1.0e-12,
    ))
    assert payload["status"] == "ok"
    assert payload["target_volume_unit"] == "mm^3"
    assert payload["normalized_measured_sets"]["external_cm"][0]["volume"] == pytest.approx(6000.0)
    assert payload["normalized_measured_sets"]["external_m"][0]["volume"] == pytest.approx(6000.0)

    missing = json.loads(build123d_volume_crosscheck_with_units(
        json.dumps([{"name": "box", "volume": 6000.0}]),
        json.dumps(measured),
    ))
    assert missing["status"] == "needs_attention"
    assert "volume_unit" in " ".join(missing["issues"])


def test_build123d_volume_crosscheck_source_gates_mcp_tool_dispatch_json():
    from radia_mcp.build123d.server import (
        build123d_volume_crosscheck,
        build123d_volume_crosscheck_source_coverage_gate,
        build123d_volume_crosscheck_source_identity_gate,
    )

    reference = [{"name": "box", "volume": 24.0}]
    measured = [
        {
            "source": "cubit",
            "rows": [{"name": "box", "volume": 24.0}],
            "measurement_method": "coreform_cubit_volume_command",
            "body_identity_key": "name",
            "source_artifact_id": "slot331_cubit_volume_rows_v1",
            "parameter_set_artifact_id": "slot391_cad_parameter_set_v1",
            "parameter_set_digest": "sha256:slot391-cad-parameter-set",
            "parameter_set_path": r"artifacts/build123d/slot391_cad_parameter_set.json",
            "objective_observable_id": "slot391_volume_quality_objective_v1",
            "objective_observable_family": "cad_volume_crosscheck",
        },
        {
            "source": "cst_import",
            "rows": [{"name": "box", "volume": 24.0}],
            "measurement_method": "cst_modeler_solid_volume_export",
            "body_identity_key": "name",
            "source_artifact_id": "slot331_cst_volume_rows_v1",
            "parameter_set_artifact_id": "slot391_cad_parameter_set_v1",
            "parameter_set_digest": "sha256:slot391-cad-parameter-set",
            "parameter_set_path": r"artifacts/build123d/slot391_cad_parameter_set.json",
            "objective_observable_id": "slot391_volume_quality_objective_v1",
            "objective_observable_family": "cad_volume_crosscheck",
        },
    ]
    summary_json = build123d_volume_crosscheck(json.dumps(reference), json.dumps(measured))
    coverage = json.loads(build123d_volume_crosscheck_source_coverage_gate(
        summary_json,
        required_sources_json=json.dumps(["cubit", "cst_import"]),
        max_allowed_volume_rel_error=1.0e-12,
    ))
    identity = json.loads(build123d_volume_crosscheck_source_identity_gate(
        summary_json,
        expected_measurement_methods_json=json.dumps({
            "cubit": "coreform_cubit_volume_command",
            "cst_import": "cst_modeler_solid_volume_export",
        }),
        expected_body_identity_keys_json=json.dumps({
            "cubit": "name",
            "cst_import": "name",
        }),
        expected_source_artifact_ids_json=json.dumps({
            "cubit": "slot331_cubit_volume_rows_v1",
            "cst_import": "slot331_cst_volume_rows_v1",
        }),
        expected_parameter_set_artifact_ids_json=json.dumps({
            "cubit": "slot391_cad_parameter_set_v1",
            "cst_import": "slot391_cad_parameter_set_v1",
        }),
        expected_parameter_set_digests_json=json.dumps({
            "cubit": "sha256:slot391-cad-parameter-set",
            "cst_import": "sha256:slot391-cad-parameter-set",
        }),
        expected_parameter_set_paths_json=json.dumps({
            "cubit": r"artifacts/build123d/slot391_cad_parameter_set.json",
            "cst_import": r"artifacts/build123d/slot391_cad_parameter_set.json",
        }),
        expected_objective_observable_ids_json=json.dumps({
            "cubit": "slot391_volume_quality_objective_v1",
            "cst_import": "slot391_volume_quality_objective_v1",
        }),
        expected_objective_observable_families_json=json.dumps({
            "cubit": "cad_volume_crosscheck",
            "cst_import": "cad_volume_crosscheck",
        }),
    ))

    assert coverage["status"] == "ok"
    assert coverage["checks"]["required_sources_present"] is True
    assert identity["status"] == "ok"
    assert identity["checks"]["expected_measurement_methods_match"] is True
    assert identity["checks"]["expected_source_artifact_ids_match"] is True
    assert identity["checks"]["expected_parameter_set_digests_match"] is True
    assert identity["checks"]["expected_objective_observable_families_match"] is True


def test_build123d_external_cad_volume_evidence_package_mcp_tool_dispatch_json():
    from radia_mcp.build123d.server import (
        build123d_external_cad_volume_evidence_package,
        build123d_volume_crosscheck,
    )

    shape_row = {
        "name": "box",
        "volume": 24.0,
        "area": 52.0,
        "is_valid": True,
        "geometry_id": "box_revB",
        "recipe_id": "slot339_box_recipe",
        "cad_kernel": "occt",
        "cad_kernel_version": "7.8-test",
        "script_path": r"artifacts/build123d/slot339_build123d.py",
        "export_id": "slot339_box_step",
        "authoring_source": "build123d_occt",
        "mesh_route": "cubit_hex_or_mixed_path",
        "length_unit": "m",
        "area_unit": "m^2",
        "volume_unit": "m^3",
        "bounding_box": {
            "min": [-1.0, -1.5, -2.0],
            "max": [1.0, 1.5, 2.0],
            "center": [0.0, 0.0, 0.0],
            "size": [2.0, 3.0, 4.0],
        },
    }
    measured = [
        {
            "source": "cubit",
            "rows": [{"name": "box", "volume": 24.0}],
            "measurement_method": "coreform_cubit_volume_command",
            "body_identity_key": "name",
            "source_artifact_id": "slot339_cubit_volume_rows_v1",
            "parameter_set_artifact_id": "slot391_cad_parameter_set_v1",
            "parameter_set_digest": "sha256:slot391-cad-parameter-set",
            "parameter_set_path": r"artifacts/build123d/slot391_cad_parameter_set.json",
            "objective_observable_id": "slot391_volume_quality_objective_v1",
            "objective_observable_family": "cad_volume_crosscheck",
        },
        {
            "source": "cst_import",
            "rows": [{"name": "box", "volume": 24.0}],
            "measurement_method": "cst_modeler_solid_volume_export",
            "body_identity_key": "name",
            "source_artifact_id": "slot339_cst_volume_rows_v1",
            "parameter_set_artifact_id": "slot391_cad_parameter_set_v1",
            "parameter_set_digest": "sha256:slot391-cad-parameter-set",
            "parameter_set_path": r"artifacts/build123d/slot391_cad_parameter_set.json",
            "objective_observable_id": "slot391_volume_quality_objective_v1",
            "objective_observable_family": "cad_volume_crosscheck",
        },
    ]
    volume_summary_json = build123d_volume_crosscheck(json.dumps([shape_row]), json.dumps(measured), rtol=1.0e-12)
    payload = json.loads(build123d_external_cad_volume_evidence_package(
        json.dumps([shape_row]),
        volume_summary_json,
        required_sources_json=json.dumps(["cubit", "cst_import"]),
        expected_measurement_methods_json=json.dumps({
            "cubit": "coreform_cubit_volume_command",
            "cst_import": "cst_modeler_solid_volume_export",
        }),
        expected_body_identity_keys_json=json.dumps({"cubit": "name", "cst_import": "name"}),
        expected_source_artifact_ids_json=json.dumps({
            "cubit": "slot339_cubit_volume_rows_v1",
            "cst_import": "slot339_cst_volume_rows_v1",
        }),
        expected_parameter_set_artifact_ids_json=json.dumps({
            "cubit": "slot391_cad_parameter_set_v1",
            "cst_import": "slot391_cad_parameter_set_v1",
        }),
        expected_parameter_set_digests_json=json.dumps({
            "cubit": "sha256:slot391-cad-parameter-set",
            "cst_import": "sha256:slot391-cad-parameter-set",
        }),
        expected_parameter_set_paths_json=json.dumps({
            "cubit": r"artifacts/build123d/slot391_cad_parameter_set.json",
            "cst_import": r"artifacts/build123d/slot391_cad_parameter_set.json",
        }),
        expected_objective_observable_ids_json=json.dumps({
            "cubit": "slot391_volume_quality_objective_v1",
            "cst_import": "slot391_volume_quality_objective_v1",
        }),
        expected_objective_observable_families_json=json.dumps({
            "cubit": "cad_volume_crosscheck",
            "cst_import": "cad_volume_crosscheck",
        }),
        expected_length_unit="m",
        expected_area_unit="m^2",
        expected_volume_unit="m^3",
        required_metadata_fields_json=json.dumps([
            "recipe_id",
            "cad_kernel",
            "cad_kernel_version",
            "script_path",
            "export_id",
        ]),
        max_allowed_volume_rel_error=1.0e-12,
    ))

    assert payload["status"] == "ok"
    assert payload["checks"]["coverage_gate_ok"] is True
    assert payload["checks"]["source_identity_gate_ok"] is True
    assert payload["checks"]["route_contract_gate_ok"] is True
    assert payload["source_identity_gate"]["checks"]["expected_parameter_set_digests_match"] is True


def test_build123d_mass_property_crosscheck_mcp_tool_dispatches_json():
    from radia_mcp.build123d.server import build123d_mass_property_crosscheck

    reference = [{
        "name": "box",
        "volume": 24.0,
        "area": 52.0,
        "is_valid": True,
        "bounding_box": {
            "min": [-1.0, -1.5, -2.0],
            "max": [1.0, 1.5, 2.0],
            "center": [0.0, 0.0, 0.0],
            "size": [2.0, 3.0, 4.0],
        },
    }]
    measured = {"cubit": [dict(reference[0])], "cst_import": [dict(reference[0])]}
    payload = json.loads(build123d_mass_property_crosscheck(json.dumps(reference), json.dumps(measured)))

    assert payload["status"] == "ok"
    assert payload["n_sources"] == 2
    assert payload["sources"] == ["cubit", "cst_import"]
    assert payload["comparison_sets"][0]["rows"][0]["passed"] is True


def test_build123d_cad_route_source_contract_mcp_tool_dispatches_json():
    from radia_mcp.build123d.server import build123d_cad_route_source_contract

    rows = [{
        "name": "box",
        "volume": 24.0,
        "area": 52.0,
        "is_valid": True,
        "geometry_id": "box_revA",
        "recipe_id": "slot243_box_recipe",
        "cad_kernel": "occt",
        "cad_kernel_version": "7.8-test",
        "script_path": r"artifacts/build123d/slot243_build123d.py",
        "export_id": "slot243_box_step",
        "source_kind": "build123d_occt",
        "mesh_route": "cubit_hex_or_mixed_path",
        "length_unit": "m",
        "area_unit": "m^2",
        "volume_unit": "m^3",
        "bounding_box": {
            "min": [-1.0, -1.5, -2.0],
            "max": [1.0, 1.5, 2.0],
            "center": [0.0, 0.0, 0.0],
            "size": [2.0, 3.0, 4.0],
        },
    }]
    crosscheck = shape_volume_crosscheck_summary(
        rows,
        {
            "coreform_cubit": [{"name": "box", "volume": 24.0}],
            "cst_import": [{"name": "box", "volume": 24.0}],
        },
        rtol=1.0e-12,
    )
    payload = json.loads(build123d_cad_route_source_contract(
        json.dumps(rows),
        json.dumps(crosscheck),
        expected_length_unit="m",
        expected_area_unit="m^2",
        expected_volume_unit="m^3",
        required_metadata_fields_json=json.dumps([
            "recipe_id",
            "cad_kernel",
            "cad_kernel_version",
            "script_path",
            "export_id",
        ]),
    ))

    assert payload["status"] == "ok"
    assert payload["external_crosscheck_sources"] == ["coreform_cubit", "cst_import"]
    assert payload["checks"]["required_source_groups_present"] is True
    assert payload["checks"]["shape_volume_unit_expected_ok"] is True
    assert payload["checks"]["required_shape_metadata_present"] is True


def test_build123d_cubit_solver_route_handoff_mcp_tool_dispatches_json():
    from radia_mcp.build123d.server import (
        build123d_cad_handoff_manifest,
        build123d_cubit_quality_ledger_handoff,
        build123d_cubit_solver_route_handoff,
    )

    rows = [{
        "name": "box",
        "volume": 24.0,
        "area": 52.0,
        "geometry_id": "box_revA",
        "mesh_route": "cubit_hex_or_mixed_path",
        "length_unit": "m",
        "area_unit": "m^2",
        "volume_unit": "m^3",
        "cad_measurement_convention": "occt_closed_solid_mass_properties",
        "cad_measurement_postprocess_row_convention_schema_id": (
            "build123d_occt_mass_property_row_convention_v1"
        ),
        "cad_measurement_component_basis_schema_id": (
            "build123d_occt_volume_area_bbox_component_basis_v1"
        ),
        "bounding_box": {
            "min": [-1.0, -1.5, -2.0],
            "max": [1.0, 1.5, 2.0],
            "center": [0.0, 0.0, 0.0],
            "size": [2.0, 3.0, 4.0],
        },
    }]
    external_volume = shape_volume_crosscheck_summary(
        rows,
        {"cubit": [{"name": "box", "volume": 24.0}]},
        rtol=1.0e-12,
    )
    solver_route_gate = cubit_mixed_solver_route_manifest_gate(
        {
            "volume_kind_counts": {"hex": 1, "pyramid": 1, "tet": 1},
            "surface_kind_counts": {"quad": 1, "triangle": 1},
            "routing_hint": "cubit_hex_or_mixed_path",
        },
        {
            "solver_route_package_id": "slot347_box_solver_route_v1",
            "routing_hint": "cubit_hex_or_mixed_path",
            "route_policy": "hex_primary_pyramid_transition_tet_compatibility",
            "downstream_solver": "NGSolve/radia-ngsolve",
            "solver_route_convention_schema_id": "coreform_mixed_hex_pyramid_tet_route_convention_v1",
            "tet_only_owner": "netgen_tri_tet_path",
            "no_implicit_tetization": True,
            "volume_routes": [
                {"volume_kind": "hex", "solver_role": "primary_volume_fem"},
                {"volume_kind": "pyramid", "solver_role": "transition_bridge", "not_primary_region": True},
                {"volume_kind": "tet", "solver_role": "compatibility_subregion_volume_fem"},
            ],
            "surface_routes": [
                {"surface_kind": "quad", "solver_role": "hex_boundary_trace"},
                {"surface_kind": "triangle", "solver_role": "tet_boundary_trace"},
            ],
        },
        expected_package_id="slot347_box_solver_route_v1",
        expected_solver_route_convention_schema_id="coreform_mixed_hex_pyramid_tet_route_convention_v1",
        require_solver_route_convention_schema=True,
    )
    solver_route_handoff = json.loads(build123d_cubit_solver_route_handoff(
        json.dumps(rows),
        json.dumps(solver_route_gate),
        expected_solver_route_package_id="slot347_box_solver_route_v1",
        expected_solver_route_convention_schema_id="coreform_mixed_hex_pyramid_tet_route_convention_v1",
        require_solver_route_convention_schema=True,
    ))
    quality_ledger_gate = {
        "policy": "cubit_mesh_quality_ledger_identity_gate",
        "status": "ok",
        "quality_artifact_id": "slot370_box_quality_ledger_v1",
        "quality_digest": "sha256:slot370-box-quality-ledger",
        "metric_set_id": "cubit_scaled_jacobian_hex_v1",
        "export_id": "slot370_box_hex_quality",
        "geometry_id": "box_revA",
        "mesh_artifact_id": "slot370_box_hex_vol_v1",
        "mesh_digest": "sha256:slot370-box-hex-vol",
        "routing_hint": "cubit_hex_or_mixed_path",
        "min_scaled_jacobian": 1.0,
        "negative_jacobian_count": 0,
        "element_type_counts": {"hex": 1},
        "inventory_is_tri_tet_only": False,
        "checks": {
            "quality_digest_recorded": True,
            "mesh_digest_recorded": True,
            "min_scaled_jacobian_above_threshold": True,
            "negative_jacobian_count_zero": True,
            "hex_or_mixed_volume_family_present": True,
            "not_tri_tet_only_for_cubit_quality_ledger": True,
        },
    }
    quality_ledger_handoff = json.loads(build123d_cubit_quality_ledger_handoff(
        json.dumps(rows),
        json.dumps(quality_ledger_gate),
        expected_quality_artifact_id="slot370_box_quality_ledger_v1",
        expected_quality_digest="sha256:slot370-box-quality-ledger",
        expected_metric_set_id="cubit_scaled_jacobian_hex_v1",
        expected_export_id="slot370_box_hex_quality",
        expected_mesh_artifact_id="slot370_box_hex_vol_v1",
        expected_mesh_digest="sha256:slot370-box-hex-vol",
    ))
    file_manifest = [
        {
            "kind": "step",
            "path": r"artifacts/build123d/slot347_box.step",
            "length_unit": "m",
            "area_unit": "m^2",
            "volume_unit": "m^3",
            "cad_measurement_convention": "occt_closed_solid_mass_properties",
            "cad_measurement_postprocess_row_convention_schema_id": (
                "build123d_occt_mass_property_row_convention_v1"
            ),
            "cad_measurement_component_basis_schema_id": (
                "build123d_occt_volume_area_bbox_component_basis_v1"
            ),
        },
        {
            "kind": "build123d_measurement_json",
            "path": r"artifacts/build123d/slot347_box_measure.json",
            "length_unit": "m",
            "area_unit": "m^2",
            "volume_unit": "m^3",
            "cad_measurement_convention": "occt_closed_solid_mass_properties",
            "cad_measurement_postprocess_row_convention_schema_id": (
                "build123d_occt_mass_property_row_convention_v1"
            ),
            "cad_measurement_component_basis_schema_id": (
                "build123d_occt_volume_area_bbox_component_basis_v1"
            ),
        },
        {
            "kind": "cubit_quality_ledger_json",
            "path": r"artifacts/cubit/slot370_box_quality_ledger.json",
            "length_unit": "m",
            "area_unit": "m^2",
            "volume_unit": "m^3",
            "cad_measurement_convention": "occt_closed_solid_mass_properties",
            "cad_measurement_postprocess_row_convention_schema_id": (
                "build123d_occt_mass_property_row_convention_v1"
            ),
            "cad_measurement_component_basis_schema_id": (
                "build123d_occt_volume_area_bbox_component_basis_v1"
            ),
        },
    ]
    handoff = json.loads(build123d_cad_handoff_manifest(
        json.dumps(rows),
        json.dumps(file_manifest),
        external_volume_summary_json=json.dumps(external_volume),
        cubit_quality_ledger_handoff_json=json.dumps(quality_ledger_handoff),
        cubit_solver_route_handoff_json=json.dumps(solver_route_handoff),
        required_file_kinds_json=json.dumps(["step", "build123d_measurement_json", "cubit_quality_ledger_json"]),
        expected_geometry_ids_json=json.dumps(["box_revA"]),
        expected_length_unit="m",
        expected_area_unit="m^2",
        expected_volume_unit="m^3",
        expected_measurement_convention="occt_closed_solid_mass_properties",
        expected_measurement_postprocess_row_convention_schema_id=(
            "build123d_occt_mass_property_row_convention_v1"
        ),
        expected_measurement_component_basis_schema_id=(
            "build123d_occt_volume_area_bbox_component_basis_v1"
        ),
        require_measurement_postprocess_row_convention_schema=True,
        require_measurement_component_basis_schema=True,
    ))

    assert solver_route_handoff["status"] == "ok"
    assert solver_route_handoff["checks"]["solver_route_pyramid_transition_role_recorded"] is True
    assert solver_route_handoff["checks"]["solver_route_no_implicit_tetization"] is True
    assert solver_route_handoff["checks"]["solver_route_convention_schema_id_recorded_when_required"] is True
    assert solver_route_handoff["checks"]["expected_solver_route_convention_schema_id_matches"] is True
    assert solver_route_handoff["solver_route_convention_schema_id"] == "coreform_mixed_hex_pyramid_tet_route_convention_v1"
    assert quality_ledger_handoff["status"] == "ok"
    assert quality_ledger_handoff["checks"]["expected_quality_digest_matches"] is True
    assert quality_ledger_handoff["checks"]["expected_mesh_digest_matches"] is True
    assert handoff["status"] == "ok"
    assert handoff["checks"]["expected_cad_measurement_postprocess_row_convention_schema_id_matches"] is True
    assert handoff["require_measurement_postprocess_row_convention_schema"] is True
    assert handoff["checks"]["expected_cad_measurement_component_basis_schema_id_matches"] is True
    assert handoff["require_measurement_component_basis_schema"] is True
    assert handoff["checks"]["cubit_quality_ledger_handoff_ok"] is True
    assert handoff["cubit_quality_ledger_handoff_policy"] == "build123d_cubit_quality_ledger_handoff_gate"
    assert handoff["checks"]["cubit_solver_route_handoff_ok"] is True
    assert handoff["cubit_solver_route_handoff_policy"] == "build123d_cubit_solver_route_handoff_gate"


def test_compare_shape_measurement_rows_compares_bbox_when_present():
    reference = [{
        "name": "box",
        "volume": 1.0,
        "area": 6.0,
        "bounding_box": {
            "min": [0.0, 0.0, 0.0],
            "max": [1.0, 1.0, 1.0],
            "center": [0.5, 0.5, 0.5],
            "size": [1.0, 1.0, 1.0],
        },
    }]
    measured = [{
        "name": "box",
        "volume": 1.0,
        "area": 6.0,
        "bounding_box": {
            "min": [0.0, 0.0, 0.0],
            "max": [1.0, 1.0, 1.0000002],
            "center": [0.5, 0.5, 0.5000001],
            "size": [1.0, 1.0, 1.0000002],
        },
    }]

    row = compare_shape_measurement_rows(reference, measured, bbox_atol=1.0e-5)[0]
    assert row["passed"]
    assert row["bbox_compared"]
    assert row["bbox_abs_error"] == pytest.approx(2.0e-7)

    measured[0]["bounding_box"]["max"] = [1.0, 1.0, 1.01]
    measured[0]["bounding_box"]["size"] = [1.0, 1.0, 1.01]
    bad = compare_shape_measurement_rows(reference, measured, bbox_atol=1.0e-5)[0]
    assert not bad["passed"]
    assert bad["reason"] == "bbox outside tolerance"


def test_shape_measurement_comparison_summary_compacts_errors():
    reference = [{"name": "box", "volume": 24.0, "area": 52.0}]
    measured = [{"name": "box", "volume": 24.0, "area": 52.0}]
    summary = shape_measurement_comparison_summary(reference, measured, measured_label="cubit")

    assert summary["measured_label"] == "cubit"
    assert summary["n_cases"] == 1
    assert summary["n_passed"] == 1
    assert summary["n_bbox_compared"] == 0
    assert summary["max_volume_rel_error"] == pytest.approx(0.0)
    assert summary["max_area_rel_error"] == pytest.approx(0.0)
    assert summary["max_bbox_abs_error"] == pytest.approx(0.0)


def test_shape_measurement_health_summary_reports_worst_mismatches():
    reference = [
        {"name": "ok", "volume": 10.0, "area": 20.0, "is_valid": True},
        {"name": "bad", "volume": 10.0, "area": 20.0, "is_valid": True},
        {"name": "missing", "volume": 1.0, "area": 2.0, "is_valid": True},
    ]
    measured = [
        {"name": "ok", "volume": 10.000001, "area": 20.000001},
        {"name": "bad", "volume": 10.2, "area": 20.0},
    ]

    health = shape_measurement_health_summary(
        reference,
        measured,
        rtol=1.0e-4,
        measured_label="external",
        worst_limit=2,
    )

    assert health["status"] == "needs_attention"
    assert health["ok_for_geometry_roundtrip"] is False
    assert health["checks"]["all_reference_shapes_valid"] is True
    assert health["checks"]["all_measurements_present_and_within_tolerance"] is False
    assert health["comparison_summary"]["n_passed"] == 1
    assert [row["name"] for row in health["worst_comparisons"]] == ["missing", "bad"]
    assert worst_shape_measurement_comparison_rows([], limit=0) == []
    with pytest.raises(ValueError, match="limit must be non-negative"):
        worst_shape_measurement_comparison_rows([], limit=-1)


def test_shape_mass_property_crosscheck_summary_handles_multiple_cad_sources():
    box = Box(2, 3, 5).solid()
    reference = [shape_measurement_row(box, name="cad_block")]
    measured = {
        "cubit": [dict(reference[0])],
        "cst_import": [dict(reference[0])],
    }

    summary = shape_mass_property_crosscheck_summary(
        reference,
        measured,
        rtol=1.0e-12,
        bbox_atol=1.0e-12,
    )

    assert summary["policy"] == "build123d_external_cad_volume_area_bbox_crosscheck"
    assert summary["status"] == "ok"
    assert summary["sources"] == ["cubit", "cst_import"]
    assert summary["n_sources"] == 2
    assert summary["n_failed_rows"] == 0
    assert summary["max_area_rel_error"] == pytest.approx(0.0)
    assert all(row["n_bbox_compared"] == 1 for row in summary["comparison_sets"])

    bad_cst = dict(reference[0])
    bad_cst["area"] *= 1.01
    bad_summary = shape_mass_property_crosscheck_summary(
        reference,
        {"cubit": [dict(reference[0])], "cst_import": [bad_cst]},
        rtol=1.0e-4,
        bbox_atol=1.0e-12,
    )
    by_source = {row["source"]: row for row in bad_summary["comparison_sets"]}

    assert bad_summary["status"] == "needs_attention"
    assert bad_summary["ok_for_cad_roundtrip_mass_properties"] is False
    assert by_source["cubit"]["status"] == "ok"
    assert by_source["cst_import"]["status"] == "needs_attention"
    assert by_source["cst_import"]["worst_comparisons"][0]["reason"] == "outside tolerance"


def test_segment_meshes_in_netgen():
    """CAE gate: the wedge tet-meshes cleanly through the build123d -> STEP -> Netgen pipeline
    (build123d's OCP kernel and Netgen's OCC binding are different, so STEP is the bridge)."""
    import tempfile
    from build123d import export_step
    from netgen.occ import OCCGeometry
    from ngsolve import Mesh
    seg = annular_segment(40, 55, 20, 0, 45)
    with tempfile.TemporaryDirectory() as d:
        f = os.path.join(d, "seg.step")
        export_step(seg, f)
        mesh = Mesh(OCCGeometry(f).GenerateMesh(maxh=5.0))
    assert mesh.ne > 50, f"wedge should tet-mesh (got {mesh.ne} elements)"


def main():
    test_annular_segment_volume_and_validity()
    test_tube_volume()
    test_racetrack_coil_builds()
    test_polar_array_full_ring_sums_to_annulus()
    test_polar_array_partial_fan()
    test_linear_array_count_and_spacing()
    test_mirrored_doubles_volume()
    test_assembly_keeps_regions_separate()
    test_segment_meshes_in_netgen()
    print("[OK] build123d modeling ops: annular wedge, tube, racetrack coil, polar/linear arrays, "
          "mirror, labelled assembly -- analytic volumes, OCCT-valid, region labels, Netgen-meshable.")


if __name__ == "__main__":
    from ngsolve import TaskManager

    with TaskManager():
        main()
