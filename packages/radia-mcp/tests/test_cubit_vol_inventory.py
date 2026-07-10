import json

from radia_mcp.cubit.server import (
    cubit_docs,
    cubit_journal_reproducibility_gate,
    cubit_mixed_order_series_gate,
    cubit_vol_inventory,
)
from radia_mcp.cubit.knowledge.netgen_workflow import get_netgen_documentation
from radia_mcp.cubit.vol_inventory import (
    cubit_bnd_area_interface_gate,
    cubit_curvilinear_handoff_manifest_gate,
    cubit_element_quality_gate,
    cubit_export_package_identity_gate,
    cubit_headless_batch_quality_package_gate,
    cubit_headless_installation_route_gate,
    cubit_hex_quality_gate,
    cubit_meshing_scheme_trace_gate,
    cubit_mesh_quality_ledger_identity_gate,
    cubit_mixed_interface_adjacency_gate,
    cubit_mixed_order_series_inventory_gate,
    cubit_mass_property_sidecar_gate,
    cubit_mixed_solver_route_manifest_gate,
    cubit_mixed_solver_ready_package_gate,
    cubit_mixed_transition_metadata_gate,
    cubit_quality_distribution_gate,
    cubit_release_feature_routing_gate,
    cubit_submodel_boundary_handoff_mesh_package_gate,
    cubit_vol_label_metadata_gate,
    summarize_netgen_vol_inventory,
)


MIXED_VOL = """
mesh3d
dimension
3

surfaceelements
2
1 1 1 0 3 1 2 3
2 2 1 0 4 1 2 3 4

volumeelements
3
1 8 1 2 3 4 5 6 7 8
2 5 1 2 3 4 5
3 4 1 2 3 5

points
8
0 0 0
1 0 0
1 1 0
0 1 0
0 0 1
1 0 1
1 1 1
0 1 1

materials
3
1 hex_core
2 pyramid_transition
3 tet_region

endmesh
"""


TRI_TET_VOL = """
mesh3d
dimension
3

surfaceelements
1
1 1 1 0 3 1 2 3

volumeelements
1
1 4 1 2 3 4

points
4
0 0 0
1 0 0
0 1 0
0 0 1

materials
1
1 domain

endmesh
"""


LABELLED_TRI_TET_VOL = """
mesh3d
dimension
3

surfaceelements
1
1 1 1 0 3 1 2 3

volumeelements
1
1 4 1 2 3 4

points
4
0 0 0
1 0 0
0 1 0
0 0 1

materials
1
1 core

bcnames
1
1 outer

endmesh
"""


CURVED_MIXED_VOL = MIXED_VOL + """
curvedelements
2
1 2 3 4
5 6 7 8
"""


LABELLED_MIXED_SUBMODEL_VOL = MIXED_VOL.replace(
    "endmesh",
    """bcnames
2
1 zoom_boundary_outer
2 far_boundary_guard

endmesh""",
)


SURFACE_UV_MIXED_VOL = """
mesh3d
dimension
3

surfaceelementsuv
2
1 1 1 0 3 1 2 3 0 0 1 0 0 1
2 2 1 0 4 1 2 3 4 0 0 1 0 1 1 0 1

volumeelements
3
1 8 1 2 3 4 5 6 7 8
2 5 1 2 3 4 5
3 4 1 2 3 5

points
8
0 0 0
1 0 0
1 1 0
0 1 0
0 0 1
1 0 1
1 1 1
0 1 1

materials
3
1 hex_core
2 pyramid_transition
3 tet_region

endmesh
"""


def test_cubit_vol_inventory_classifies_hex_pyramid_tet_mixed_mesh():
    inv = summarize_netgen_vol_inventory(MIXED_VOL, source="unit")

    assert inv["volume_kind_counts"] == {"hex": 1, "pyramid": 1, "tet": 1}
    assert inv["surface_kind_counts"] == {"quad": 1, "triangle": 1}
    assert inv["curvedelements_present"] is False
    assert inv["has_mixed_hex_transition"] is True
    assert inv["is_tri_tet_only"] is False
    assert inv["routing_hint"] == "cubit_hex_or_mixed_path"
    assert inv["materials"][2] == "pyramid_transition"
    assert "Cubit/Coreform owns hex-led" in inv["policy"]


def test_cubit_vol_inventory_keeps_tet_only_on_netgen_route():
    inv = summarize_netgen_vol_inventory(TRI_TET_VOL)

    assert inv["volume_kind_counts"] == {"tet": 1}
    assert inv["surface_kind_counts"] == {"triangle": 1}
    assert inv["curvedelements_present"] is False
    assert inv["has_mixed_hex_transition"] is False
    assert inv["is_tri_tet_only"] is True
    assert inv["routing_hint"] == "netgen_tri_tet_path"


def test_cubit_vol_label_metadata_gate_requires_materials_and_boundaries():
    inv = summarize_netgen_vol_inventory(LABELLED_TRI_TET_VOL)
    gate = cubit_vol_label_metadata_gate(
        inv,
        required_materials=("core",),
        required_boundaries=("outer",),
    )

    assert inv["materials"] == {1: "core"}
    assert inv["boundary_names"] == {1: "outer"}
    assert gate["policy"] == "cubit_vol_label_metadata_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["required_materials_present"] is True
    assert gate["checks"]["required_boundaries_present"] is True

    missing = cubit_vol_label_metadata_gate(
        summarize_netgen_vol_inventory(TRI_TET_VOL),
        required_materials=("core",),
        required_boundaries=("outer",),
    )
    assert missing["status"] == "needs_attention"
    assert missing["checks"]["required_materials_present"] is False
    assert missing["checks"]["boundary_names_recorded"] is False


def test_cubit_vol_inventory_keeps_curved_mixed_mesh_on_cubit_route():
    inv = summarize_netgen_vol_inventory(CURVED_MIXED_VOL)

    assert inv["curvedelements_present"] is True
    assert inv["volume_kind_counts"] == {"hex": 1, "pyramid": 1, "tet": 1}
    assert inv["surface_kind_counts"] == {"quad": 1, "triangle": 1}
    assert inv["routing_hint"] == "cubit_hex_or_mixed_path"
    assert "curvedelements section can grow with order" in inv["order_series_policy"]


def test_cubit_mixed_transition_metadata_gate_requires_pyramid_bridge_labels_and_routing():
    inv = summarize_netgen_vol_inventory(MIXED_VOL, source="unit")
    gate = cubit_mixed_transition_metadata_gate(inv)

    assert gate["policy"] == "cubit_mixed_transition_metadata_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["required_volume_kinds_present"] is True
    assert gate["checks"]["required_surface_kinds_present"] is True
    assert gate["checks"]["transition_kinds_present"] is True
    assert gate["checks"]["routing_hint_is_cubit_mixed"] is True
    assert gate["checks"]["transition_material_label_present"] is True
    assert gate["required_surface_kinds"] == ["quad", "triangle"]
    assert "zero sidecar material volume" in " ".join(gate["notes"])

    tet_only = cubit_mixed_transition_metadata_gate(summarize_netgen_vol_inventory(TRI_TET_VOL))
    assert tet_only["status"] == "needs_attention"
    assert tet_only["checks"]["contains_hex_region"] is False
    assert tet_only["checks"]["transition_kinds_present"] is False
    assert tet_only["checks"]["routing_hint_is_cubit_mixed"] is False

    missing_label_inv = dict(inv)
    missing_label_inv["materials"] = {1: "hex_core", 3: "tet_region"}
    missing_label = cubit_mixed_transition_metadata_gate(missing_label_inv)
    assert missing_label["status"] == "needs_attention"
    assert missing_label["checks"]["transition_material_label_present"] is False

    missing_quad_surface_inv = dict(inv)
    missing_quad_surface_inv["surface_kind_counts"] = {"triangle": 2}
    missing_quad_surface = cubit_mixed_transition_metadata_gate(missing_quad_surface_inv)
    assert missing_quad_surface["status"] == "needs_attention"
    assert missing_quad_surface["checks"]["required_surface_kinds_present"] is False


def test_cubit_mixed_order_series_inventory_gate_keeps_routing_topology_invariant():
    order1 = summarize_netgen_vol_inventory(MIXED_VOL, source="mixed_o1.vol")
    order3 = summarize_netgen_vol_inventory(CURVED_MIXED_VOL, source="mixed_o3.vol")
    gate = cubit_mixed_order_series_inventory_gate([
        {"order": 1, "inventory": order1},
        {"order": 3, "inventory": order3},
    ])

    assert gate["policy"] == "cubit_mixed_order_series_inventory_gate"
    assert gate["status"] == "ok"
    assert gate["orders"] == [1, 3]
    assert gate["checks"]["first_order_inventory_present"] is True
    assert gate["checks"]["first_order_inventory_not_curved"] is True
    assert gate["checks"]["volume_kind_counts_invariant"] is True
    assert gate["checks"]["surface_kind_counts_invariant"] is True
    assert gate["checks"]["routing_hint_is_cubit_mixed"] is True
    assert gate["checks"]["required_volume_kinds_present"] is True
    assert gate["checks"]["required_surface_kinds_present"] is True
    assert gate["series"][1]["curvedelements_present"] is True
    assert "curvedelements section may grow" in " ".join(gate["notes"])

    tet_only = summarize_netgen_vol_inventory(TRI_TET_VOL, source="tet_only_o3.vol")
    wrong_route = cubit_mixed_order_series_inventory_gate([
        {"order": 1, "inventory": order1},
        {"order": 3, "inventory": tet_only},
    ])
    assert wrong_route["status"] == "needs_attention"
    assert wrong_route["checks"]["volume_kind_counts_invariant"] is False
    assert wrong_route["checks"]["surface_kind_counts_invariant"] is False
    assert wrong_route["checks"]["routing_hint_is_cubit_mixed"] is False
    assert wrong_route["checks"]["has_mixed_hex_transition"] is False
    assert wrong_route["checks"]["not_tri_tet_only"] is False

    stale_topology = dict(order3)
    stale_topology["volume_kind_counts"] = {"hex": 1, "tet": 2}
    stale_gate = cubit_mixed_order_series_inventory_gate([
        {"order": 1, "inventory": order1},
        {"order": 3, "inventory": stale_topology},
    ])
    assert stale_gate["status"] == "needs_attention"
    assert stale_gate["checks"]["volume_kind_counts_invariant"] is False

    missing_first = cubit_mixed_order_series_inventory_gate([
        {"order": 2, "inventory": order3},
        {"order": 3, "inventory": order3},
    ])
    assert missing_first["status"] == "needs_attention"
    assert missing_first["checks"]["first_order_inventory_present"] is False


def test_cubit_mixed_order_series_mcp_tool_accepts_rows_and_rejects_drift():
    base = {
        "volume_kind_counts": {"hex": 4, "pyramid": 8, "tet": 12},
        "surface_kind_counts": {"quad": 6, "triangle": 10},
        "routing_hint": "cubit_hex_or_mixed_path",
        "is_tri_tet_only": False,
        "has_mixed_hex_transition": True,
    }
    rows = [
        {"order": 1, **base, "curvedelements_present": False},
        {"order": 2, **base, "curvedelements_present": True},
    ]
    ok = json.loads(cubit_mixed_order_series_gate(rows))
    assert ok["status"] == "ok"
    assert ok["orders"] == [1, 2]

    drift = [dict(row) for row in rows]
    drift[1] = {
        **drift[1],
        "volume_kind_counts": {"hex": 4, "pyramid": 0, "tet": 20},
    }
    bad = json.loads(cubit_mixed_order_series_gate(drift))
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["volume_kind_counts_invariant"] is False


def test_cubit_journal_reproducibility_gate_refuses_comment_only_root_cause():
    failed = """# NG label
reset
cylinder height 1 radius 1
mesh volume 1
"""
    passed = """# OK label

reset
cylinder   height 1   radius 1
mesh volume 1
"""
    same = json.loads(cubit_journal_reproducibility_gate(
        failed, passed, outcome_a="failed", outcome_b="passed"
    ))
    assert same["status"] == "needs_run_provenance"
    assert same["analysis_complete"] is True
    assert same["commands_equal"] is True
    assert same["command_digest_a"] == same["command_digest_b"]
    assert same["script_difference_explains_outcome"] is False
    assert "complete_solver_log" in same["required_run_provenance"]

    changed = json.loads(cubit_journal_reproducibility_gate(
        failed,
        passed.replace("mesh volume 1", "volume 1 scheme sweep\nmesh volume 1"),
        outcome_a="failed",
        outcome_b="passed",
    ))
    assert changed["status"] == "commands_differ"
    assert changed["commands_equal"] is False
    assert changed["differences"]


def test_cubit_vol_inventory_counts_surfaceelementsuv_from_stock_export_netgen():
    inv = summarize_netgen_vol_inventory(SURFACE_UV_MIXED_VOL)

    assert inv["surface_section"] == "surfaceelementsuv"
    assert inv["surface_kind_counts"] == {"quad": 1, "triangle": 1}
    assert inv["volume_kind_counts"] == {"hex": 1, "pyramid": 1, "tet": 1}
    assert inv["routing_hint"] == "cubit_hex_or_mixed_path"
    assert "surfaceelementsuv" in inv["policy"]


def test_cubit_vol_inventory_mcp_tool_dispatches_json():
    payload = json.loads(cubit_vol_inventory(text=MIXED_VOL))

    assert payload["status"] == "ok"
    assert payload["volume_kind_counts"]["hex"] == 1
    assert payload["volume_kind_counts"]["pyramid"] == 1
    assert payload["routing_hint"] == "cubit_hex_or_mixed_path"


def test_cubit_docs_route_tet_only_to_netgen_and_mixed_to_cubit():
    lab = cubit_docs("scripting_lab_policy")
    routing = cubit_docs("format_routing")
    pyramid = cubit_docs("scripting_pyramid_handling")

    assert "tet-only mesh request" in lab
    assert "prefer build123d/Netgen/OCC" in lab
    assert "hex + tet transition mesh" in lab
    assert "cubit_vol_inventory" in routing
    assert "cubit_hex_or_mixed_path" in routing
    assert "Inventory them explicitly" in routing
    assert "high-order `.vol` files" in routing
    assert "companion `.vol.json` material volume alone" in routing
    assert "quad and triangle surface families" in routing
    assert "mixed transition adjacency roles" in routing
    assert "hex_to_transition" in routing
    assert "transition_to_tet" in routing
    assert "Cubit meshing scheme traces" in routing
    assert "cubit_meshing_scheme_trace_gate" in routing
    assert "export netgen" in routing
    assert "Headless Cubit smoke tests can export a valid `.vol` before teardown warns" in routing
    assert "source`/`sink`/`sibc" in routing
    assert "hex-led and mixed hex+pyramid+tet lane" in pyramid


def test_netgen_workflow_records_o_grid_hex_sphere_gate():
    doc = get_netgen_documentation("overview")

    assert "O-grid hex sphere gate" in doc
    assert "volume 1 scheme sphere" in doc
    assert "56 hexes" in doc
    assert "order-3 rel err 0.00131" in doc
    assert "Do not call `mesh.Curve()`" in doc
    assert "coreform_cubit.com -nographics -batch" in doc
    assert "O-grid hex sphere eigenvalue gate" in doc
    assert "lambda_1=(pi/R)^2" in doc
    assert "4.493409457909064" in doc
    assert "lambda1=9.86268248871533" in doc
    assert "relative errors `7.01e-4` and `5.96e-4`" in doc


def test_netgen_workflow_records_mapped_hex_brick_area_gate():
    doc = get_netgen_documentation("overview")

    assert "mapped hex brick volume/area gate" in doc
    assert "volume 1 scheme map" in doc
    assert "192 hexes" in doc
    assert "surface-area rel err `2.05e-15`" in doc
    assert "Avoid multi-line `dict(...)`" in doc
    assert "Cubit mass-property sidecar gate" in doc
    assert "cubit_mass_property_sidecar_gate" in doc
    assert "get_volume_volume" in doc
    assert "Slot210 adds the unit contract" in doc
    assert "length_unit" in doc
    assert "mm^3" in doc
    assert "bbox size `[1.5, 2.0, 0.75]`" in doc
    assert "Cubit export package identity gate" in doc
    assert "cubit_export_package_identity_gate" in doc
    assert "`export_id`, `geometry_id`, order, and routing hint" in doc
    assert "old sidecar or a raw JSON from a different geometry" in doc
    assert "Slot411 adds schema identity" in doc
    assert "vol_sidecar_schema_id" in doc
    assert "require_vol_sidecar_schema=True" in doc
    assert "Slot418 adds solver-route convention schema identity" in doc
    assert "solver_route_convention_schema_id" in doc
    assert "require_solver_route_convention_schema=True" in doc


def test_cubit_hex_quality_gate_replays_mapped_hex_quality_slot():
    gate = cubit_hex_quality_gate(
        [1.0, 0.9999999999999999, 1.0],
        expected_hex_count=3,
        min_scaled_jacobian=0.99,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "cubit_hex_scaled_jacobian_quality_gate"
    assert gate["count"] == 3
    assert gate["min"] == 0.9999999999999999
    assert gate["bad_count"] == 0

    bad = cubit_hex_quality_gate([1.0, 0.1], expected_hex_count=3, min_scaled_jacobian=0.2)
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["min_scaled_jacobian_ok"] is False
    assert bad["checks"]["expected_hex_count_ok"] is False


def test_cubit_element_quality_gate_replays_tetra10_tri6_metrics():
    tet = cubit_element_quality_gate(
        [0.82, 0.76, 0.61],
        element_type="Tetra10",
        metric="scaled Jacobian",
        expected_count=3,
        min_value=0.5,
    )

    assert tet["status"] == "ok"
    assert tet["policy"] == "cubit_element_quality_metric_gate"
    assert tet["element_type"] == "tetra10"
    assert tet["metric"] == "scaled_jacobian"
    assert tet["bad_count"] == 0
    assert "Tetra10/Tri6 Jacobian" in tet["version_note"]

    tri = cubit_element_quality_gate(
        [0.91, 0.64, 0.18],
        element_type="Tri6",
        metric="Jacobian",
        expected_count=4,
        min_value=0.25,
    )
    assert tri["status"] == "needs_attention"
    assert tri["checks"]["min_value_ok"] is False
    assert tri["checks"]["expected_count_ok"] is False
    assert tri["bad_count"] == 1


def test_cubit_quality_distribution_gate_records_quantiles_and_histogram():
    gate = cubit_quality_distribution_gate(
        [1.0, 0.98, 0.97, 0.90, 0.82],
        element_type="hex",
        metric="scaled jacobian",
        expected_count=5,
        min_value=0.8,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "cubit_quality_distribution_gate"
    assert gate["finite_count"] == 5
    assert gate["nonfinite_count"] == 0
    assert gate["quantiles"]["p00"] == 0.82
    assert gate["quantiles"]["p50"] == 0.97
    assert gate["quantiles"]["p100"] == 1.0
    assert gate["histogram"]["underflow"] == 0
    assert gate["histogram"]["overflow"] == 0
    assert sum(bin_record["count"] for bin_record in gate["histogram"]["bins"]) == 5
    assert "distribution" in gate["version_note"]

    bad = cubit_quality_distribution_gate(
        [1.0, 0.75, 0.1, float("nan")],
        expected_count=5,
        min_value=0.2,
    )
    assert bad["status"] == "needs_attention"
    assert bad["nonfinite_count"] == 1
    assert bad["low_count"] == 1
    assert bad["bad_count"] == 2
    assert bad["checks"]["all_finite"] is False
    assert bad["checks"]["expected_count_ok"] is False


def test_cubit_mesh_quality_ledger_identity_gate_binds_quality_to_mesh_artifact():
    quality = cubit_quality_distribution_gate(
        [1.0, 0.99, 0.98, 0.97],
        element_type="hex",
        metric="scaled jacobian",
        expected_count=4,
        min_value=0.95,
    )
    quality.update({
        "export_id": "slot369_hex_quality_ledger",
        "geometry_id": "unit_brick_quality_ledger_v1",
        "quality_component_basis_schema_id": (
            "coreform_hex_scaled_jacobian_element_component_basis_v1"
        ),
    })
    export_package = {
        "policy": "cubit_export_package_identity_gate",
        "status": "ok",
        "export_ids": ["slot369_hex_quality_ledger"],
        "geometry_ids": ["unit_brick_quality_ledger_v1"],
        "export_output_artifact_id": "slot369_hex_quality_ledger_vol_v1",
        "export_output_digest": "sha256:vol369",
    }
    headless = {
        "policy": "cubit_headless_batch_quality_package_gate",
        "status": "ok",
        "export_id": "slot369_hex_quality_ledger",
        "geometry_id": "unit_brick_quality_ledger_v1",
        "quality_count": 4,
        "qualityComponentBasisSchemaId": (
            "coreform_hex_scaled_jacobian_element_component_basis_v1"
        ),
    }
    inventory = {
        "source": r"artifacts/cubit/coreform_slot369_hex_quality_ledger.vol",
        "volume_kind_counts": {"hex": 4},
        "surface_kind_counts": {"quad": 24},
        "routing_hint": "cubit_hex_or_mixed_path",
        "is_tri_tet_only": False,
    }
    ledger = {
        "mesh_quality_artifact_id": "slot369_hex_quality_ledger_json_v1",
        "mesh_quality_digest": "sha256:quality369",
        "mesh_quality_path": r"artifacts/cubit/coreform_slot369_hex_quality_ledger.json",
        "quality_metric_set_id": "cubit_scaled_jacobian_hex_v1",
        "mesh_quality_postprocess_row_convention_schema_id": (
            "coreform_scaled_jacobian_hex_quality_row_convention_v1"
        ),
        "mesh_quality_component_basis_schema_id": (
            "coreform_hex_scaled_jacobian_element_component_basis_v1"
        ),
        "export_id": "slot369_hex_quality_ledger",
        "geometry_id": "unit_brick_quality_ledger_v1",
        "mesh_artifact_id": "slot369_hex_quality_ledger_vol_v1",
        "mesh_digest": "sha256:vol369",
        "metric": "scaled_jacobian",
        "min_scaled_jacobian": 0.97,
        "negative_jacobian_count": 0,
        "element_type_counts": {"hex": 4},
        "routing_hint": "cubit_hex_or_mixed_path",
        "version": "2025.12",
        "created_at_utc": "2026-07-01T10:07:11Z",
        "elapsed_s": 0.010,
        "parameter_set_artifact_id": "coreform_slot390_hex_quality_parameter_set_v1",
        "parameter_set_digest": "sha256:coreform-slot390-hex-quality-parameter-set-v1",
        "parameter_set_path": r"artifacts/cubit/slot390_hex_quality_parameter_set.json",
        "objective_observable_id": "coreform_slot390_min_scaled_jacobian_objective_v1",
        "objective_observable_family": "mesh_quality_min_scaled_jacobian_objective",
        "parameter_set": {
            "parameterSetArtifactId": "coreform_slot390_hex_quality_parameter_set_v1",
            "parameterSetDigest": "sha256:coreform-slot390-hex-quality-parameter-set-v1",
            "parameterSetPath": r"artifacts/cubit/slot390_hex_quality_parameter_set.json",
        },
        "objective": {
            "objectiveObservableId": "coreform_slot390_min_scaled_jacobian_objective_v1",
            "objectiveObservableFamily": "mesh_quality_min_scaled_jacobian_objective",
        },
        "timing_breakdown_s": {
            "headless_startup": 0.001,
            "mesh_quality_query": 0.003,
            "vol_inventory": 0.002,
            "write_json": 0.001,
        },
    }

    gate = cubit_mesh_quality_ledger_identity_gate(
        ledger,
        quality_gate=quality,
        export_package_gate=export_package,
        headless_batch_quality_gate=headless,
        inventory=inventory,
        expected_quality_artifact_id="slot369_hex_quality_ledger_json_v1",
        expected_quality_digest="sha256:quality369",
        expected_metric_set_id="cubit_scaled_jacobian_hex_v1",
        expected_export_id="slot369_hex_quality_ledger",
        expected_geometry_id="unit_brick_quality_ledger_v1",
        expected_mesh_artifact_id="slot369_hex_quality_ledger_vol_v1",
        expected_mesh_digest="sha256:vol369",
        expected_version="2025.12",
        expected_parameter_set_artifact_id="coreform_slot390_hex_quality_parameter_set_v1",
        expected_parameter_set_digest="sha256:coreform-slot390-hex-quality-parameter-set-v1",
        expected_parameter_set_path=r"artifacts/cubit/slot390_hex_quality_parameter_set.json",
        expected_objective_observable_id="coreform_slot390_min_scaled_jacobian_objective_v1",
        expected_objective_observable_family="mesh_quality_min_scaled_jacobian_objective",
        expected_quality_postprocess_row_convention_schema_id=(
            "coreform_scaled_jacobian_hex_quality_row_convention_v1"
        ),
        expected_quality_component_basis_schema_id=(
            "coreform_hex_scaled_jacobian_element_component_basis_v1"
        ),
        expected_element_type_counts={"hex": 4},
        min_scaled_jacobian=0.95,
        require_execution_metadata=True,
        require_parameter_set_artifact=True,
        require_quality_postprocess_row_convention_schema=True,
        require_quality_component_basis_schema=True,
        require_timing_breakdown=True,
    )

    assert gate["policy"] == "cubit_mesh_quality_ledger_identity_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["expected_quality_digest_matches"] is True
    assert gate["checks"]["quality_distribution_gate_ok"] is True
    assert gate["checks"]["mesh_digest_matches_export_package_when_present"] is True
    assert gate["checks"]["inventory_counts_match_ledger_counts"] is True
    assert gate["checks"]["not_tri_tet_only_for_cubit_quality_ledger"] is True
    assert gate["checks"]["created_at_utc_parseable_when_present"] is True
    assert gate["checks"]["expected_version_matches"] is True
    assert gate["checks"]["parameter_set_artifact_id_consistent_when_present"] is True
    assert gate["checks"]["parameter_set_digest_consistent_when_present"] is True
    assert gate["checks"]["parameter_set_path_consistent_when_present"] is True
    assert gate["checks"]["objective_observable_id_consistent_when_present"] is True
    assert gate["checks"]["objective_observable_family_consistent_when_present"] is True
    assert gate["checks"]["parameter_set_artifact_id_recorded_when_required"] is True
    assert gate["checks"]["parameter_set_digest_recorded_when_required"] is True
    assert gate["checks"]["parameter_set_path_recorded_when_required"] is True
    assert gate["checks"]["parameter_set_artifact_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_parameter_set_artifact_id_matches"] is True
    assert gate["checks"]["parameter_set_digest_recorded_when_expected"] is True
    assert gate["checks"]["expected_parameter_set_digest_matches"] is True
    assert gate["checks"]["parameter_set_path_recorded_when_expected"] is True
    assert gate["checks"]["expected_parameter_set_path_matches"] is True
    assert gate["checks"]["objective_observable_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_objective_observable_id_matches"] is True
    assert gate["checks"]["objective_observable_family_recorded_when_expected"] is True
    assert gate["checks"]["expected_objective_observable_family_matches"] is True
    assert gate["checks"]["quality_postprocess_row_convention_schema_id_consistent_when_present"] is True
    assert gate["checks"]["quality_postprocess_row_convention_schema_id_recorded_when_required"] is True
    assert gate["checks"]["quality_postprocess_row_convention_schema_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_quality_postprocess_row_convention_schema_id_matches"] is True
    assert gate["checks"]["quality_component_basis_schema_id_consistent_when_present"] is True
    assert gate["checks"]["quality_component_basis_schema_id_recorded_when_required"] is True
    assert gate["checks"]["quality_component_basis_schema_id_recorded_when_expected"] is True
    assert gate["checks"]["expected_quality_component_basis_schema_id_matches"] is True
    assert gate["checks"]["timing_breakdown_has_required_stage_count"] is True
    assert gate["checks"]["timing_breakdown_total_within_elapsed_when_present"] is True
    assert gate["parameter_set_artifact_id"] == "coreform_slot390_hex_quality_parameter_set_v1"
    assert gate["parameter_set_digest"] == "sha256:coreform-slot390-hex-quality-parameter-set-v1"
    assert gate["parameter_set_path"] == r"artifacts/cubit/slot390_hex_quality_parameter_set.json"
    assert gate["objective_observable_id"] == "coreform_slot390_min_scaled_jacobian_objective_v1"
    assert gate["objective_observable_family"] == "mesh_quality_min_scaled_jacobian_objective"
    assert gate["quality_postprocess_row_convention_schema_id"] == (
        "coreform_scaled_jacobian_hex_quality_row_convention_v1"
    )
    assert gate["quality_postprocess_row_convention_schema_ids"] == [
        "coreform_scaled_jacobian_hex_quality_row_convention_v1"
    ]
    assert gate["quality_postprocess_row_convention_schema_required"] is True
    assert gate["quality_component_basis_schema_id"] == (
        "coreform_hex_scaled_jacobian_element_component_basis_v1"
    )
    assert gate["quality_component_basis_schema_ids"] == [
        "coreform_hex_scaled_jacobian_element_component_basis_v1"
    ]
    assert gate["quality_component_basis_schema_required"] is True

    stale_digest = cubit_mesh_quality_ledger_identity_gate(
        {**ledger, "mesh_quality_digest": "sha256:old"},
        expected_quality_digest="sha256:quality369",
    )
    assert stale_digest["status"] == "needs_attention"
    assert stale_digest["checks"]["expected_quality_digest_matches"] is False

    inverted = cubit_mesh_quality_ledger_identity_gate(
        {**ledger, "negative_jacobian_count": 2, "min_scaled_jacobian": -0.1},
        min_scaled_jacobian=0.95,
    )
    assert inverted["status"] == "needs_attention"
    assert inverted["checks"]["negative_jacobian_count_zero"] is False
    assert inverted["checks"]["min_scaled_jacobian_above_threshold"] is False

    tet_inventory = {
        **inventory,
        "volume_kind_counts": {"tet": 4},
        "surface_kind_counts": {"triangle": 24},
        "routing_hint": "netgen_tri_tet_path",
        "is_tri_tet_only": True,
    }
    wrong_inventory = cubit_mesh_quality_ledger_identity_gate(
        ledger,
        inventory=tet_inventory,
    )
    assert wrong_inventory["status"] == "needs_attention"
    assert wrong_inventory["inventory_is_tri_tet_only"] is True
    assert wrong_inventory["checks"]["inventory_routing_hint_matches_ledger"] is False
    assert wrong_inventory["checks"]["inventory_counts_match_ledger_counts"] is False
    assert wrong_inventory["checks"]["not_tri_tet_only_for_cubit_quality_ledger"] is False

    stale_version = cubit_mesh_quality_ledger_identity_gate(
        {**ledger, "version": "2024.8"},
        expected_version="2025.12",
        require_execution_metadata=True,
        require_timing_breakdown=True,
    )
    assert stale_version["status"] == "needs_attention"
    assert stale_version["checks"]["expected_version_matches"] is False

    bad_timestamp = cubit_mesh_quality_ledger_identity_gate(
        {**ledger, "created_at_utc": "not-a-date"},
        require_execution_metadata=True,
        require_timing_breakdown=True,
    )
    assert bad_timestamp["status"] == "needs_attention"
    assert bad_timestamp["checks"]["created_at_utc_parseable_when_present"] is False

    sparse_timing = cubit_mesh_quality_ledger_identity_gate(
        {**ledger, "timing_breakdown_s": {"mesh_quality_query": 0.003}},
        require_execution_metadata=True,
        require_timing_breakdown=True,
    )
    assert sparse_timing["status"] == "needs_attention"
    assert sparse_timing["checks"]["timing_breakdown_has_required_stage_count"] is False

    impossible_timing = cubit_mesh_quality_ledger_identity_gate(
        {**ledger, "elapsed_s": 0.001},
        require_execution_metadata=True,
        require_timing_breakdown=True,
    )
    assert impossible_timing["status"] == "needs_attention"
    assert impossible_timing["checks"]["timing_breakdown_total_within_elapsed_when_present"] is False

    stale_parameter_set_digest = {
        **ledger,
        "parameter_set_digest": "sha256:old-parameter-set",
        "parameter_set": {
            **ledger["parameter_set"],
            "parameterSetDigest": "sha256:old-parameter-set",
        },
    }
    stale_parameter_set_digest_gate = cubit_mesh_quality_ledger_identity_gate(
        stale_parameter_set_digest,
        expected_parameter_set_digest="sha256:coreform-slot390-hex-quality-parameter-set-v1",
        require_parameter_set_artifact=True,
    )
    assert stale_parameter_set_digest_gate["status"] == "needs_attention"
    assert stale_parameter_set_digest_gate["checks"]["parameter_set_digest_consistent_when_present"] is True
    assert stale_parameter_set_digest_gate["checks"]["expected_parameter_set_digest_matches"] is False

    missing_parameter_set_path = {
        key: value for key, value in ledger.items() if key != "parameter_set_path"
    }
    missing_parameter_set_path["parameter_set"] = dict(ledger["parameter_set"])
    missing_parameter_set_path["parameter_set"].pop("parameterSetPath")
    missing_parameter_set_path_gate = cubit_mesh_quality_ledger_identity_gate(
        missing_parameter_set_path,
        require_parameter_set_artifact=True,
    )
    assert missing_parameter_set_path_gate["status"] == "needs_attention"
    assert missing_parameter_set_path_gate["checks"]["parameter_set_artifact_id_recorded_when_required"] is True
    assert missing_parameter_set_path_gate["checks"]["parameter_set_digest_recorded_when_required"] is True
    assert missing_parameter_set_path_gate["checks"]["parameter_set_path_recorded_when_required"] is False

    wrong_objective_family = {
        **ledger,
        "objective_observable_family": "surface_area_objective",
        "objective": {
            **ledger["objective"],
            "objectiveObservableFamily": "surface_area_objective",
        },
    }
    wrong_objective_family_gate = cubit_mesh_quality_ledger_identity_gate(
        wrong_objective_family,
        expected_objective_observable_family="mesh_quality_min_scaled_jacobian_objective",
    )
    assert wrong_objective_family_gate["status"] == "needs_attention"
    assert wrong_objective_family_gate["checks"]["objective_observable_family_consistent_when_present"] is True
    assert wrong_objective_family_gate["checks"]["expected_objective_observable_family_matches"] is False

    stale_quality_row_convention_schema = {
        **ledger,
        "mesh_quality_postprocess_row_convention_schema_id": (
            "coreform_quality_min_scalar_row_v0"
        ),
    }
    stale_quality_row_convention_schema_gate = cubit_mesh_quality_ledger_identity_gate(
        stale_quality_row_convention_schema,
        expected_metric_set_id="cubit_scaled_jacobian_hex_v1",
        expected_quality_postprocess_row_convention_schema_id=(
            "coreform_scaled_jacobian_hex_quality_row_convention_v1"
        ),
        require_quality_postprocess_row_convention_schema=True,
    )
    assert stale_quality_row_convention_schema_gate["status"] == "needs_attention"
    assert stale_quality_row_convention_schema_gate["checks"]["expected_metric_set_id_matches"] is True
    assert (
        stale_quality_row_convention_schema_gate["checks"][
            "quality_postprocess_row_convention_schema_id_recorded_when_required"
        ]
        is True
    )
    assert (
        stale_quality_row_convention_schema_gate["checks"][
            "expected_quality_postprocess_row_convention_schema_id_matches"
        ]
        is False
    )

    missing_quality_row_convention_schema = {
        key: value
        for key, value in ledger.items()
        if key != "mesh_quality_postprocess_row_convention_schema_id"
    }
    missing_quality_row_convention_schema_gate = cubit_mesh_quality_ledger_identity_gate(
        missing_quality_row_convention_schema,
        expected_quality_postprocess_row_convention_schema_id=(
            "coreform_scaled_jacobian_hex_quality_row_convention_v1"
        ),
        require_quality_postprocess_row_convention_schema=True,
    )
    assert missing_quality_row_convention_schema_gate["status"] == "needs_attention"
    assert (
        missing_quality_row_convention_schema_gate["checks"][
            "quality_postprocess_row_convention_schema_id_recorded_when_required"
        ]
        is False
    )
    assert (
        missing_quality_row_convention_schema_gate["checks"][
            "quality_postprocess_row_convention_schema_id_recorded_when_expected"
        ]
        is False
    )

    stale_quality_component_basis_schema = {
        **ledger,
        "mesh_quality_component_basis_schema_id": (
            "coreform_quality_scalar_value_component_basis_v0"
        ),
    }
    stale_quality_component_basis_schema_gate = cubit_mesh_quality_ledger_identity_gate(
        stale_quality_component_basis_schema,
        expected_metric_set_id="cubit_scaled_jacobian_hex_v1",
        expected_quality_component_basis_schema_id=(
            "coreform_hex_scaled_jacobian_element_component_basis_v1"
        ),
        require_quality_component_basis_schema=True,
    )
    assert stale_quality_component_basis_schema_gate["status"] == "needs_attention"
    assert stale_quality_component_basis_schema_gate["checks"]["expected_metric_set_id_matches"] is True
    assert (
        stale_quality_component_basis_schema_gate["checks"][
            "quality_component_basis_schema_id_recorded_when_required"
        ]
        is True
    )
    assert (
        stale_quality_component_basis_schema_gate["checks"][
            "expected_quality_component_basis_schema_id_matches"
        ]
        is False
    )

    missing_quality_component_basis_schema = {
        key: value
        for key, value in ledger.items()
        if key != "mesh_quality_component_basis_schema_id"
    }
    missing_quality_component_basis_schema_gate = cubit_mesh_quality_ledger_identity_gate(
        missing_quality_component_basis_schema,
        expected_quality_component_basis_schema_id=(
            "coreform_hex_scaled_jacobian_element_component_basis_v1"
        ),
        require_quality_component_basis_schema=True,
    )
    assert missing_quality_component_basis_schema_gate["status"] == "needs_attention"
    assert (
        missing_quality_component_basis_schema_gate["checks"][
            "quality_component_basis_schema_id_recorded_when_required"
        ]
        is False
    )
    assert (
        missing_quality_component_basis_schema_gate["checks"][
            "quality_component_basis_schema_id_recorded_when_expected"
        ]
        is False
    )


def test_cubit_release_feature_routing_gate_maps_2026_6_to_lab_lanes():
    features = [
        {
            "feature_key": "anisotropic_tetrahedral_meshing",
            "category": "Meshing",
            "lab_route": "advanced tet reference, not the default tet-only education path",
            "validation_note": "record anisotropy metadata before comparing against Netgen tri/tet examples",
        },
        {
            "feature_key": "cohesive_element_generation",
            "category": "Meshing",
            "lab_route": "interface/crack teaching examples with explicit block identity",
            "validation_note": "require cohesive block and unmerged surface provenance",
        },
        {
            "feature_key": "higher_order_quality_metrics",
            "category": "Meshing",
            "lab_route": "Tetra10/Tri6 Jacobian quality replay",
            "validation_note": "archive raw metric lists and replay cubit_element_quality_gate",
        },
        {
            "feature_key": "tri_tet_meshing_robustness",
            "category": "Meshing",
            "lab_route": "composite-surface robustness watchlist",
            "validation_note": "keep tet-only solver examples on Netgen/OCC unless Cubit adds value",
        },
        {
            "feature_key": "sculpt_refinement_memory",
            "category": "Meshing",
            "lab_route": "large hex refinement validation lane",
            "validation_note": "record memory and element-count envelopes",
        },
        {
            "feature_key": "solver_io_compatibility",
            "category": "Input/Output",
            "lab_route": "Exodus/Abaqus/I-DEAS import-export compatibility notes",
            "validation_note": "record 64-bit ID and degenerate-element handling separately from .vol",
        },
        {
            "feature_key": "python_312_runtime",
            "category": "Miscellaneous",
            "lab_route": "headless scripting compatibility",
            "validation_note": "keep Python script syntax simple for batch execution",
        },
    ]

    gate = cubit_release_feature_routing_gate(
        features,
        release_version="2026.6",
        source_url="https://coreform.com/coreform-cubit/release-notes/v2026-6/",
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "cubit_release_feature_routing_gate"
    assert gate["feature_count"] == 7
    assert gate["checks"]["source_url_is_coreform"] is True
    assert gate["checks"]["required_features_present"] is True
    assert gate["checks"]["cubit_role_is_hex_or_mixed"] is True
    assert gate["checks"]["tet_only_owner_is_netgen"] is True

    missing = cubit_release_feature_routing_gate(
        features[:-1],
        release_version="2026.6",
        source_url="https://coreform.com/coreform-cubit/release-notes/v2026-6/",
    )
    assert missing["status"] == "needs_attention"
    assert missing["missing_features"] == ["python_312_runtime"]


def test_cubit_headless_installation_route_gate_separates_installed_version_from_release_watchlist():
    gate = cubit_headless_installation_route_gate(
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

    assert gate["policy"] == "cubit_headless_installation_route_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["required_headless_flags_present"] is True
    assert gate["checks"]["binary_path_is_console_com"] is True
    assert gate["checks"]["gui_daemon_disabled_by_default"] is True
    assert gate["checks"]["live_claim_matches_installed_version"] is True
    assert gate["checks"]["release_note_watchlist_not_live_claim"] is True
    assert gate["checks"]["license_status_recorded"] is True
    assert gate["checks"]["license_status_allows_headless_probe"] is True
    assert gate["checks"]["version_probe_is_synchronous_console"] is True
    assert gate["checks"]["version_probe_uses_recorded_binary"] is True
    assert gate["checks"]["version_probe_summary_records_installed_version"] is True
    assert gate["checks"]["version_probe_summary_records_license_status"] is True

    overclaim = cubit_headless_installation_route_gate(
        {
            "installed_version": "2025.12",
            "binary_path": "C:/Program Files/Coreform Cubit 2025.12/bin/coreform_cubit.com",
            "binary_exists": True,
            "headless_flags": ["-nographics", "-batch"],
            "gui_policy": "headless_no_gui_daemon_by_default",
            "release_note_version": "2026.6",
            "release_note_status": "installed",
            "live_claimed_release_version": "2026.6",
        }
    )
    assert overclaim["status"] == "needs_attention"
    assert overclaim["checks"]["live_claim_matches_installed_version"] is False
    assert overclaim["checks"]["release_note_watchlist_not_live_claim"] is False

    gui_binary = cubit_headless_installation_route_gate(
        {
            "installed_version": "2025.12",
            "binary_path": "C:/Program Files/Coreform Cubit 2025.12/bin/coreform_cubit.exe",
            "binary_exists": True,
            "headless_flags": ["-nographics", "-batch"],
            "gui_policy": "headless_no_gui_daemon_by_default",
            "live_claimed_release_version": "2025.12",
            "version_probe_command": "C:/Program Files/Coreform Cubit 2025.12/bin/coreform_cubit.com -version",
        }
    )
    assert gui_binary["status"] == "needs_attention"
    assert gui_binary["checks"]["binary_path_is_console_com"] is False
    assert gui_binary["checks"]["version_probe_is_synchronous_console"] is True
    assert gui_binary["checks"]["version_probe_uses_recorded_binary"] is False

    stale_probe_summary = cubit_headless_installation_route_gate(
        {
            "installed_version": "2025.12",
            "binary_path": "C:/Program Files/Coreform Cubit 2025.12/bin/coreform_cubit.com",
            "binary_exists": True,
            "headless_flags": ["-nographics", "-batch"],
            "gui_policy": "headless_no_gui_daemon_by_default",
            "live_claimed_release_version": "2025.12",
            "license_status": "ValidStudent",
            "version_probe_summary": {
                "license_status": "ValidStudent",
                "version_line": "Coreform Cubit Version 2026.6 Build future",
            },
        }
    )
    assert stale_probe_summary["status"] == "needs_attention"
    assert stale_probe_summary["checks"]["version_probe_summary_records_installed_version"] is False
    assert stale_probe_summary["checks"]["version_probe_summary_records_license_status"] is True


def test_cubit_curvilinear_handoff_manifest_gate_keeps_import_geometry_order_and_quality_together():
    manifest = {
        "mesh_id": "slot170_imported_hex_cylinder",
        "export_id": "slot170_hex_curved_o3",
        "source_mesh": {
            "kind": "third_party_mesh",
            "path": r"artifacts/cubit/slot170_imported_hex_cylinder.inp",
            "volume_kinds": ["hex"],
            "surface_kinds": ["quad"],
        },
        "geometry_association": {
            "cad_source": "cylinder.step",
            "projection_policy": "project_boundary_nodes_to_cad_curves_and_surfaces",
            "boundary_ids_preserved": True,
            "projection_quality": {
                "max_distance": 2.0e-8,
                "tolerance": 1.0e-6,
                "normal_deviation_max": 0.002,
            },
        },
        "curved_export": {
            "format": "netgen_vol",
            "order": 3,
            "path": r"artifacts/cubit/slot170_hex_curved_o3.vol",
            "routing_hint": "cubit_hex_or_mixed_path",
            "implicit_element_conversion": False,
        },
        "quality": {
            "metric": "scaled_jacobian",
            "min": 0.995,
            "count": 224,
            "negative_jacobian_count": 0,
        },
        "provenance": {
            "literature_note": "third-party curvilinear handoff requires mesh-to-CAD association before high-order export",
        },
    }

    gate = cubit_curvilinear_handoff_manifest_gate(
        manifest,
        expected_mesh_id="slot170_imported_hex_cylinder",
        expected_export_id="slot170_hex_curved_o3",
        min_order=2,
        min_scaled_jacobian=0.2,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "cubit_curvilinear_handoff_manifest_gate"
    assert gate["order"] == 3
    assert gate["quality_min"] == 0.995
    assert gate["checks"]["hex_or_mixed_cubit_route"] is True
    assert gate["checks"]["geometry_association_recorded"] is True
    assert gate["checks"]["projection_error_recorded"] is True
    assert gate["checks"]["projection_error_within_tolerance"] is True
    assert gate["checks"]["no_implicit_element_conversion"] is True
    assert gate["checks"]["negative_jacobian_count_recorded"] is True
    assert gate["checks"]["negative_jacobian_count_zero"] is True

    missing_geometry = dict(manifest)
    missing_geometry["geometry_association"] = {"cad_source": "", "boundary_ids_preserved": False}
    missing_gate = cubit_curvilinear_handoff_manifest_gate(missing_geometry)
    assert missing_gate["status"] == "needs_attention"
    assert missing_gate["checks"]["geometry_association_recorded"] is False
    assert missing_gate["checks"]["boundary_ids_preserved"] is False

    first_order = dict(manifest)
    first_order["curved_export"] = {**manifest["curved_export"], "order": 1}
    first_order_gate = cubit_curvilinear_handoff_manifest_gate(first_order)
    assert first_order_gate["status"] == "needs_attention"
    assert first_order_gate["checks"]["curved_export_order_ok"] is False

    poor_projection = dict(manifest)
    poor_projection["geometry_association"] = {
        **manifest["geometry_association"],
        "projection_quality": {"max_distance": 3.0e-4, "tolerance": 1.0e-6},
    }
    poor_projection_gate = cubit_curvilinear_handoff_manifest_gate(poor_projection)
    assert poor_projection_gate["status"] == "needs_attention"
    assert poor_projection_gate["checks"]["projection_error_within_tolerance"] is False

    inverted = dict(manifest)
    inverted["quality"] = {**manifest["quality"], "negative_jacobian_count": 2}
    inverted_gate = cubit_curvilinear_handoff_manifest_gate(inverted)
    assert inverted_gate["status"] == "needs_attention"
    assert inverted_gate["checks"]["negative_jacobian_count_zero"] is False

    tet_only = dict(manifest)
    tet_only["source_mesh"] = {**manifest["source_mesh"], "volume_kinds": ["tet"], "surface_kinds": ["triangle"]}
    tet_only_gate = cubit_curvilinear_handoff_manifest_gate(tet_only)
    assert tet_only_gate["status"] == "needs_attention"
    assert tet_only_gate["checks"]["hex_or_mixed_cubit_route"] is False


def test_cubit_bnd_area_interface_gate_matches_external_plus_interface_once():
    gate = cubit_bnd_area_interface_gate(
        external_area=57.0,
        material_interface_area=12.0,
        ngsolve_bnd_area=69.0,
        rel_tol=1e-12,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "cubit_ngsolve_bnd_area_includes_material_interfaces_once"
    assert gate["expected_bnd_area"] == 69.0
    assert gate["matches_external_only"] is False
    assert gate["checks"]["matches_external_plus_interface"] is True

    bad = cubit_bnd_area_interface_gate(
        external_area=57.0,
        material_interface_area=12.0,
        ngsolve_bnd_area=57.0,
        rel_tol=1e-12,
    )
    assert bad["status"] == "needs_attention"
    assert bad["matches_external_only"] is True
    assert bad["checks"]["matches_external_plus_interface"] is False


def test_cubit_mass_property_sidecar_gate_records_volume_area_and_bbox():
    rows = [
        {
            "name": "hex_brick",
            "volume": 2.25,
            "area": 11.25,
            "length_unit": "m",
            "volume_unit": "m^3",
            "area_unit": "m^2",
            "bounding_box": {"size": [1.5, 2.0, 0.75]},
        }
    ]
    gate = cubit_mass_property_sidecar_gate(
        rows,
        expected_total_volume=2.25,
        expected_total_area=11.25,
        expected_bbox_size=[1.5, 2.0, 0.75],
        expected_length_unit="m",
        expected_volume_unit="m^3",
        expected_area_unit="m^2",
        rel_tol=1e-12,
        abs_tol=1e-12,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "cubit_mass_property_sidecar_gate"
    assert gate["checks"]["total_volume_expected_ok"] is True
    assert gate["checks"]["total_area_expected_ok"] is True
    assert gate["checks"]["bbox_size_expected_ok"] is True
    assert gate["checks"]["length_unit_expected_ok"] is True
    assert gate["checks"]["volume_unit_expected_ok"] is True
    assert gate["checks"]["area_unit_expected_ok"] is True
    assert gate["units"] == {"length": ["m"], "volume": ["m^3"], "area": ["m^2"]}
    assert gate["volume_rel_error"] == 0.0
    assert gate["area_rel_error"] == 0.0
    assert "Volume alone is not enough" in " ".join(gate["notes"])
    assert "mm^3 and m^3" in " ".join(gate["notes"])

    bad = cubit_mass_property_sidecar_gate(
        [
            {
                "name": "hex_brick",
                "volume": 2.25,
                "area": 10.0,
                "units": {"length": "m", "volume": "mm^3", "area": "m^2"},
                "bounding_box": {"size": [1.5, 2.0, 0.7]},
            }
        ],
        expected_total_volume=2.25,
        expected_total_area=11.25,
        expected_bbox_size=[1.5, 2.0, 0.75],
        expected_volume_unit="m^3",
        rel_tol=1e-12,
        abs_tol=1e-12,
    )
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["total_volume_expected_ok"] is True
    assert bad["checks"]["total_area_expected_ok"] is False
    assert bad["checks"]["bbox_size_expected_ok"] is False
    assert bad["checks"]["volume_unit_expected_ok"] is False


def test_cubit_mass_property_sidecar_gate_binds_material_labels_and_allows_zero_transition():
    rows = [
        {
            "name": "hex_core",
            "volume": 2.0,
            "area": 10.0,
            "length_unit": "m",
            "volume_unit": "m^3",
            "area_unit": "m^2",
            "bounding_box": {"size": [2.0, 1.0, 1.0]},
        },
        {
            "name": "pyramid_transition",
            "volume": 0.0,
            "area": 0.0,
            "length_unit": "m",
            "volume_unit": "m^3",
            "area_unit": "m^2",
        },
        {
            "name": "tet_region",
            "volume": 0.5,
            "area": 3.0,
            "length_unit": "m",
            "volume_unit": "m^3",
            "area_unit": "m^2",
        },
    ]

    gate = cubit_mass_property_sidecar_gate(
        rows,
        expected_total_volume=2.5,
        expected_total_area=13.0,
        expected_bbox_size=[2.0, 1.0, 1.0],
        expected_length_unit="m",
        expected_volume_unit="m^3",
        expected_area_unit="m^2",
        expected_material_names=("hex_core", "pyramid_transition", "tet_region"),
        allow_zero_measurement_names=("pyramid_transition",),
        rel_tol=1e-12,
        abs_tol=1e-12,
    )

    assert gate["status"] == "ok"
    assert gate["checks"]["expected_material_names_present"] is True
    assert gate["checks"]["all_volumes_positive"] is True
    assert gate["checks"]["all_areas_positive"] is True
    assert gate["zero_volume_row_names"] == ["pyramid_transition"]
    assert gate["zero_area_row_names"] == ["pyramid_transition"]
    assert gate["disallowed_zero_volume_row_names"] == []

    missing_transition = cubit_mass_property_sidecar_gate(
        [row for row in rows if row["name"] != "pyramid_transition"],
        expected_total_volume=2.5,
        expected_total_area=13.0,
        expected_material_names=("hex_core", "pyramid_transition", "tet_region"),
        allow_zero_measurement_names=("pyramid_transition",),
    )
    assert missing_transition["status"] == "needs_attention"
    assert missing_transition["missing_material_names"] == ["pyramid_transition"]
    assert missing_transition["checks"]["expected_material_names_present"] is False

    zero_not_named = cubit_mass_property_sidecar_gate(
        rows,
        expected_material_names=("hex_core", "pyramid_transition", "tet_region"),
    )
    assert zero_not_named["status"] == "needs_attention"
    assert zero_not_named["disallowed_zero_volume_row_names"] == ["pyramid_transition"]
    assert zero_not_named["checks"]["all_volumes_positive"] is False


def test_cubit_export_package_identity_gate_pairs_vol_sidecar_raw_and_ids():
    vol_path = "artifacts/slot146_hex_brick_o3.vol"
    artifacts = [
        {
            "kind": "vol",
            "path": vol_path,
            "export_id": "slot146_hex_brick_o3",
            "geometry_id": "hex_brick_v1",
            "order": 3,
        },
        {
            "kind": "vol_sidecar",
            "path": vol_path + ".json",
            "export_id": "slot146_hex_brick_o3",
            "geometry_id": "hex_brick_v1",
            "order": 3,
            "n_elements": 12,
            "n_points": 13,
            "vol_sidecar_schema_id": "coreform_netgen_vol_sidecar_inventory_v1",
        },
        {
            "kind": "raw_result",
            "path": "artifacts/slot146_hex_brick_raw.json",
            "export_id": "slot146_hex_brick_o3",
            "geometry_id": "hex_brick_v1",
        },
        {
            "kind": "mass_property_sidecar",
            "path": "artifacts/slot146_hex_brick_mass.json",
            "export_id": "slot146_hex_brick_o3",
            "geometry_id": "hex_brick_v1",
        },
    ]
    for row in artifacts:
        row["export_observable_id"] = "slot306_netgen_vol_inventory_v1"
        row["export_observable_family"] = "netgen_vol_inventory"
    inventory = {
        "source": vol_path,
        "routing_hint": "cubit_hex_or_mixed_path",
        "volume_elements": 12,
        "points": 13,
        "inventory_observable_id": "slot306_netgen_vol_inventory_v1",
        "inventory_observable_family": "netgen_vol_inventory",
    }

    gate = cubit_export_package_identity_gate(
        artifacts,
        expected_export_id="slot146_hex_brick_o3",
        expected_geometry_id="hex_brick_v1",
        expected_order=3,
        expected_routing_hint="cubit_hex_or_mixed_path",
        expected_export_observable_id="slot306_netgen_vol_inventory_v1",
        expected_export_observable_family="netgen_vol_inventory",
        require_export_observable=True,
        expected_vol_sidecar_schema_id="coreform_netgen_vol_sidecar_inventory_v1",
        require_vol_sidecar_schema=True,
        require_vol_sidecar_inventory_counts=True,
        inventory=inventory,
    )

    assert gate["policy"] == "cubit_export_package_identity_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["vol_sidecar_pairs_vol"] is True
    assert gate["checks"]["geometry_id_unique"] is True
    assert gate["checks"]["inventory_source_matches_vol"] is True
    assert gate["checks"]["inventory_routing_hint_matches_expected"] is True
    assert gate["checks"]["vol_sidecar_inventory_counts_recorded_when_required"] is True
    assert gate["checks"]["vol_sidecar_schema_id_recorded_when_required"] is True
    assert gate["checks"]["expected_vol_sidecar_schema_id_matches"] is True
    assert gate["checks"]["vol_sidecar_element_count_matches_inventory"] is True
    assert gate["checks"]["vol_sidecar_point_count_matches_inventory"] is True
    assert gate["checks"]["vol_sidecar_order_matches_expected"] is True
    assert gate["vol_sidecar_element_counts"] == [12]
    assert gate["vol_sidecar_point_counts"] == [13]
    assert gate["vol_sidecar_schema_id"] == "coreform_netgen_vol_sidecar_inventory_v1"
    assert gate["checks"]["export_observable_id_recorded_when_required"] is True
    assert gate["checks"]["expected_export_observable_id_matches"] is True
    assert gate["checks"]["expected_export_observable_family_matches"] is True
    assert gate["export_observable_id"] == "slot306_netgen_vol_inventory_v1"
    assert gate["export_observable_family"] == "netgen_vol_inventory"

    wrong_sidecar = [dict(row) for row in artifacts]
    wrong_sidecar[1]["path"] = "artifacts/slot146_other.vol.json"
    bad_sidecar = cubit_export_package_identity_gate(
        wrong_sidecar,
        expected_export_id="slot146_hex_brick_o3",
        expected_geometry_id="hex_brick_v1",
        expected_order=3,
        expected_routing_hint="cubit_hex_or_mixed_path",
        inventory=inventory,
    )
    assert bad_sidecar["status"] == "needs_attention"
    assert bad_sidecar["checks"]["vol_sidecar_pairs_vol"] is False

    stale_sidecar_count = [dict(row) for row in artifacts]
    stale_sidecar_count[1]["n_elements"] = 11
    bad_sidecar_count = cubit_export_package_identity_gate(
        stale_sidecar_count,
        expected_export_id="slot146_hex_brick_o3",
        expected_geometry_id="hex_brick_v1",
        expected_order=3,
        expected_routing_hint="cubit_hex_or_mixed_path",
        require_vol_sidecar_inventory_counts=True,
        inventory=inventory,
    )
    assert bad_sidecar_count["status"] == "needs_attention"
    assert bad_sidecar_count["checks"]["vol_sidecar_element_count_matches_inventory"] is False

    missing_sidecar_count = [dict(row) for row in artifacts]
    missing_sidecar_count[1].pop("n_points")
    bad_missing_sidecar_count = cubit_export_package_identity_gate(
        missing_sidecar_count,
        expected_export_id="slot146_hex_brick_o3",
        expected_geometry_id="hex_brick_v1",
        expected_order=3,
        expected_routing_hint="cubit_hex_or_mixed_path",
        require_vol_sidecar_inventory_counts=True,
        inventory=inventory,
    )
    assert bad_missing_sidecar_count["status"] == "needs_attention"
    assert bad_missing_sidecar_count["checks"]["vol_sidecar_inventory_counts_recorded_when_required"] is False

    stale_sidecar_schema = [dict(row) for row in artifacts]
    stale_sidecar_schema[1]["vol_sidecar_schema_id"] = "coreform_legacy_material_volume_sidecar_v0"
    bad_sidecar_schema = cubit_export_package_identity_gate(
        stale_sidecar_schema,
        expected_export_id="slot146_hex_brick_o3",
        expected_geometry_id="hex_brick_v1",
        expected_order=3,
        expected_vol_sidecar_schema_id="coreform_netgen_vol_sidecar_inventory_v1",
        require_vol_sidecar_schema=True,
        inventory=inventory,
    )
    assert bad_sidecar_schema["status"] == "needs_attention"
    assert bad_sidecar_schema["checks"]["expected_vol_sidecar_schema_id_matches"] is False

    missing_sidecar_schema = [dict(row) for row in artifacts]
    missing_sidecar_schema[1].pop("vol_sidecar_schema_id")
    bad_missing_sidecar_schema = cubit_export_package_identity_gate(
        missing_sidecar_schema,
        expected_export_id="slot146_hex_brick_o3",
        expected_geometry_id="hex_brick_v1",
        expected_order=3,
        require_vol_sidecar_schema=True,
        inventory=inventory,
    )
    assert bad_missing_sidecar_schema["status"] == "needs_attention"
    assert bad_missing_sidecar_schema["checks"]["vol_sidecar_schema_id_recorded_when_required"] is False

    wrong_geometry = [dict(row) for row in artifacts]
    wrong_geometry[-1]["geometry_id"] = "hex_brick_v2"
    bad_geometry = cubit_export_package_identity_gate(wrong_geometry)
    assert bad_geometry["status"] == "needs_attention"
    assert bad_geometry["checks"]["geometry_id_unique"] is False

    wrong_observable = [dict(row) for row in artifacts]
    wrong_observable[1]["export_observable_id"] = "slot305_quality_distribution_v1"
    bad_observable = cubit_export_package_identity_gate(
        wrong_observable,
        expected_export_observable_id="slot306_netgen_vol_inventory_v1",
        require_export_observable=True,
        inventory=inventory,
    )
    assert bad_observable["status"] == "needs_attention"
    assert bad_observable["checks"]["export_observable_id_consistent_when_present"] is False
    assert bad_observable["checks"]["expected_export_observable_id_matches"] is False

    wrong_observable_family = [dict(row) for row in artifacts]
    wrong_observable_family[0]["export_observable_family"] = "artifact_quality_distribution"
    bad_observable_family = cubit_export_package_identity_gate(
        wrong_observable_family,
        expected_export_observable_family="netgen_vol_inventory",
        require_export_observable=True,
        inventory=inventory,
    )
    assert bad_observable_family["status"] == "needs_attention"
    assert bad_observable_family["checks"]["export_observable_family_consistent_when_present"] is False
    assert bad_observable_family["checks"]["expected_export_observable_family_matches"] is False


def test_cubit_mixed_solver_ready_package_gate_binds_transition_export_bnd_and_quality():
    vol_path = r"artifacts/cubit/slot178_mixed_hex_pyramid_tet.vol"
    inv = summarize_netgen_vol_inventory(MIXED_VOL, source=vol_path)
    transition = cubit_mixed_transition_metadata_gate(inv)
    export_package = cubit_export_package_identity_gate(
        [
            {
                "kind": "vol",
                "path": vol_path,
                "export_id": "slot178_mixed_hex_pyramid_tet",
                "geometry_id": "mixed_transition_fixture",
                "order": 2,
                "export_output_artifact_id": "slot298_mixed_export_package_v1",
                "export_output_digest": "sha256:slot298_mixed_export_package_v1",
                "export_output_path": r"artifacts/cubit/slot298_mixed_export_package.json",
            },
            {
                "kind": "vol_sidecar",
                "path": vol_path + ".json",
                "export_id": "slot178_mixed_hex_pyramid_tet",
                "geometry_id": "mixed_transition_fixture",
                "order": 2,
                "export_output_artifact_id": "slot298_mixed_export_package_v1",
                "export_output_digest": "sha256:slot298_mixed_export_package_v1",
                "export_output_path": r"artifacts/cubit/slot298_mixed_export_package.json",
            },
            {
                "kind": "raw_result",
                "path": r"artifacts/cubit/slot178_raw.json",
                "export_id": "slot178_mixed_hex_pyramid_tet",
                "geometry_id": "mixed_transition_fixture",
                "export_output_artifact_id": "slot298_mixed_export_package_v1",
                "export_output_digest": "sha256:slot298_mixed_export_package_v1",
                "export_output_path": r"artifacts/cubit/slot298_mixed_export_package.json",
            },
        ],
        expected_export_id="slot178_mixed_hex_pyramid_tet",
        expected_geometry_id="mixed_transition_fixture",
        expected_order=2,
        expected_routing_hint="cubit_hex_or_mixed_path",
        expected_export_output_artifact_id="slot298_mixed_export_package_v1",
        expected_export_output_digest="sha256:slot298_mixed_export_package_v1",
        require_export_output_artifact=True,
        inventory=inv,
    )
    bnd = cubit_bnd_area_interface_gate(
        external_area=6.0,
        material_interface_area=2.0,
        ngsolve_bnd_area=8.0,
        rel_tol=1.0e-12,
    )
    quality = cubit_quality_distribution_gate(
        [0.82, 0.91, 0.76, 0.88],
        element_type="mixed_hex_pyramid_tet",
        expected_count=4,
        min_value=0.2,
    )
    batch = {
        "pass": True,
        "export_id": "slot178_mixed_hex_pyramid_tet",
        "geometry_id": "mixed_transition_fixture",
        "command_line": r"coreform_cubit.com -nographics -batch artifacts/cubit/slot330_mixed_solver_ready.py",
        "process_mode": "headless_batch",
        "batch_script": r"artifacts/cubit/slot330_mixed_solver_ready.py",
        "journal_policy": "batch_script_archived_no_gui_daemon",
        "gui_daemon": False,
        "exit_code": 0,
        "solver_ready_claimed": True,
        "version": "Coreform Cubit 2025.12",
        "output_paths": [vol_path],
    }
    headless_batch = cubit_headless_batch_quality_package_gate(
        batch,
        {**quality, "export_id": "slot178_mixed_hex_pyramid_tet", "geometry_id": "mixed_transition_fixture"},
        expected_export_id="slot178_mixed_hex_pyramid_tet",
        expected_geometry_id="mixed_transition_fixture",
        expected_element_type=None,
        export_inventory=inv,
    )
    interface = cubit_mixed_interface_adjacency_gate(
        [
            {
                "surface_id": 101,
                "role": "hex_to_transition",
                "surface_kind": "quad",
                "adjacent_material_names": ["hex_core", "pyramid_transition"],
                "adjacent_volume_kinds": ["hex", "pyramid"],
                "boundary_name": "hex_pyramid_interface",
            },
            {
                "surface_id": 102,
                "role": "transition_to_tet",
                "surface_kind": "triangle",
                "adjacent_material_names": ["pyramid_transition", "tet_region"],
                "adjacent_volume_kinds": ["pyramid", "tet"],
                "boundary_name": "pyramid_tet_interface",
            },
        ],
    )
    scheme_trace = cubit_meshing_scheme_trace_gate(
        {
            "trace_id": "slot290_mixed_scheme_trace",
            "command_digest": "sha256:slot290-map-tetmesh-export",
            "commands": [
                "imprint all",
                "merge all",
                "volume 1 scheme map",
                "volume 2 scheme tetmesh",
                "export netgen \"slot290_mixed.vol\" order 2 overwrite",
            ],
            "volume_schemes": {
                "1": "map",
                "2": "tetmesh",
            },
            "export_order": 2,
            "export_output_artifact_id": "slot298_mixed_export_package_v1",
            "export_output_digest": "sha256:slot298_mixed_export_package_v1",
            "export_output_path": r"artifacts/cubit/slot298_mixed_export_package.json",
        },
        expected_trace_id="slot290_mixed_scheme_trace",
        expected_command_digest="sha256:slot290-map-tetmesh-export",
        expected_volume_schemes={"1": "map", "2": "tetmesh"},
        expected_export_order=2,
        expected_export_output_artifact_id="slot298_mixed_export_package_v1",
        expected_export_output_digest="sha256:slot298_mixed_export_package_v1",
        expected_export_output_path=r"artifacts/cubit/slot298_mixed_export_package.json",
    )
    routing_policy = cubit_release_feature_routing_gate(
        [
            {
                "feature_key": "tri_tet_meshing_robustness",
                "category": "Meshing",
                "lab_route": "Cubit is reserved for hex-led or mixed transition meshes",
                "validation_note": "tet-only education stays on Netgen/OCC unless the slot explicitly overrides it",
            }
        ],
        release_version="2026.6",
        source_url="https://coreform.com/coreform-cubit/release-notes/v2026-6/",
        required_features=("tri_tet_meshing_robustness",),
    )
    curvilinear = cubit_curvilinear_handoff_manifest_gate(
        {
            "mesh_id": "slot338_mixed_imported_mesh",
            "export_id": "slot178_mixed_hex_pyramid_tet",
            "source_mesh": {
                "kind": "third_party_mesh",
                "volume_kinds": ["hex", "pyramid", "tet"],
                "surface_kinds": ["quad", "triangle"],
            },
            "geometry_association": {
                "cad_source": "mixed_transition_fixture.step",
                "projection_policy": "project_boundary_nodes_to_cad_curves_and_surfaces",
                "boundary_ids_preserved": True,
                "projection_quality": {"max_distance": 2.0e-8, "tolerance": 1.0e-6},
            },
            "curved_export": {
                "format": "netgen_vol",
                "order": 2,
                "routing_hint": "cubit_hex_or_mixed_path",
                "implicit_element_conversion": False,
            },
            "quality": {
                "metric": "scaled_jacobian",
                "min": 0.76,
                "count": 4,
                "negative_jacobian_count": 0,
            },
            "provenance": {
                "literature_note": "third-party curvilinear handoff needs CAD association and projection validation before reuse",
            },
        },
        expected_export_id="slot178_mixed_hex_pyramid_tet",
    )
    solver_route = cubit_mixed_solver_route_manifest_gate(
        inv,
        {
            "solver_route_package_id": "slot346_mixed_solver_route_v1",
            "routing_hint": "cubit_hex_or_mixed_path",
            "route_policy": "hex_primary_pyramid_transition_tet_compatibility",
            "solver_route_convention_schema_id": "coreform_mixed_hex_pyramid_tet_route_convention_v1",
            "downstream_solver": "NGSolve/radia-ngsolve",
            "downstream_solver_contract_artifact_id": "slot383_ngsolve_mixed_element_reader_contract_v1",
            "downstream_solver_contract_digest": "sha256:slot383-ngsolve-mixed-element-reader-contract-v1",
            "downstream_solver_contract_path": r"artifacts/cubit/slot383_ngsolve_mixed_element_reader_contract.json",
            "tet_only_owner": "netgen_tri_tet_path",
            "no_implicit_tetization": True,
            "volume_routes": [
                {
                    "volume_kind": "hex",
                    "material": "hex_core",
                    "solver_role": "primary_volume_fem",
                    "space_family": ["H1", "HCurl"],
                },
                {
                    "volume_kind": "pyramid",
                    "material": "pyramid_transition",
                    "solver_role": "transition_bridge",
                    "space_family": ["H1"],
                    "not_primary_region": True,
                },
                {
                    "volume_kind": "tet",
                    "material": "tet_region",
                    "solver_role": "compatibility_subregion_volume_fem",
                    "space_family": ["H1", "HCurl"],
                },
            ],
            "surface_routes": [
                {"surface_kind": "quad", "solver_role": "hex_boundary_trace"},
                {"surface_kind": "triangle", "solver_role": "tet_boundary_trace"},
            ],
        },
        expected_package_id="slot346_mixed_solver_route_v1",
        expected_solver_contract_artifact_id="slot383_ngsolve_mixed_element_reader_contract_v1",
        expected_solver_contract_digest="sha256:slot383-ngsolve-mixed-element-reader-contract-v1",
        expected_solver_contract_path=r"artifacts/cubit/slot383_ngsolve_mixed_element_reader_contract.json",
        expected_solver_route_convention_schema_id="coreform_mixed_hex_pyramid_tet_route_convention_v1",
        require_solver_contract_artifact=True,
        require_solver_route_convention_schema=True,
    )

    gate = cubit_mixed_solver_ready_package_gate(
        inv,
        transition,
        export_package,
        bnd,
        quality,
        routing_policy_gate=routing_policy,
        interface_adjacency_gate=interface,
        scheme_trace_gate=scheme_trace,
        headless_batch_quality_gate=headless_batch,
        curvilinear_handoff_gate=curvilinear,
        solver_route_gate=solver_route,
    )

    assert gate["policy"] == "cubit_mixed_solver_ready_package_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["inventory_is_mixed_hex_pyramid_tet"] is True
    assert gate["checks"]["transition_gate_ok"] is True
    assert gate["checks"]["export_package_gate_ok"] is True
    assert gate["checks"]["bnd_area_gate_ok"] is True
    assert gate["checks"]["quality_distribution_gate_ok"] is True
    assert gate["checks"]["routing_policy_gate_ok"] is True
    assert gate["checks"]["interface_adjacency_gate_ok"] is True
    assert gate["checks"]["scheme_trace_gate_ok"] is True
    assert gate["checks"]["headless_batch_quality_gate_ok"] is True
    assert gate["checks"]["curvilinear_handoff_gate_ok"] is True
    assert gate["checks"]["solver_route_manifest_gate_ok"] is True
    assert gate["headless_batch_quality_policy"] == "cubit_headless_batch_quality_package_gate"
    assert gate["curvilinear_handoff_policy"] == "cubit_curvilinear_handoff_manifest_gate"
    assert gate["solver_route_policy"] == "cubit_mixed_solver_route_manifest_gate"
    assert solver_route["checks"]["hex_primary_volume_role_recorded"] is True
    assert solver_route["checks"]["pyramid_transition_role_recorded"] is True
    assert solver_route["checks"]["tet_compatibility_role_recorded"] is True
    assert solver_route["checks"]["no_implicit_tetization_recorded"] is True
    assert solver_route["checks"]["solver_route_convention_schema_id_recorded_when_required"] is True
    assert solver_route["checks"]["expected_solver_route_convention_schema_id_matches"] is True
    assert solver_route["checks"]["solver_contract_artifact_id_recorded_when_required"] is True
    assert solver_route["checks"]["solver_contract_digest_recorded_when_required"] is True
    assert solver_route["checks"]["solver_contract_path_recorded_when_required"] is True
    assert solver_route["checks"]["expected_solver_contract_artifact_id_matches"] is True
    assert solver_route["checks"]["expected_solver_contract_digest_matches"] is True
    assert solver_route["checks"]["expected_solver_contract_path_matches"] is True
    assert solver_route["solver_route_convention_schema_id"] == "coreform_mixed_hex_pyramid_tet_route_convention_v1"
    assert solver_route["solver_contract_artifact_id"] == "slot383_ngsolve_mixed_element_reader_contract_v1"
    assert solver_route["solver_contract_digest"] == "sha256:slot383-ngsolve-mixed-element-reader-contract-v1"
    assert solver_route["solver_contract_path"] == r"artifacts/cubit/slot383_ngsolve_mixed_element_reader_contract.json"
    assert export_package["checks"]["export_output_artifact_id_recorded_when_required"] is True
    assert export_package["checks"]["export_output_digest_recorded_when_required"] is True
    assert export_package["checks"]["export_output_path_recorded_when_required"] is True
    assert export_package["checks"]["expected_export_output_artifact_id_matches"] is True
    assert export_package["checks"]["expected_export_output_digest_matches"] is True
    assert export_package["export_output_artifact_id"] == "slot298_mixed_export_package_v1"
    assert export_package["export_output_digest"] == "sha256:slot298_mixed_export_package_v1"
    assert export_package["export_output_path"] == r"artifacts/cubit/slot298_mixed_export_package.json"
    assert scheme_trace["checks"]["expected_export_output_artifact_id_matches"] is True
    assert scheme_trace["checks"]["expected_export_output_digest_matches"] is True
    assert scheme_trace["checks"]["expected_export_output_path_matches"] is True
    assert scheme_trace["export_output_artifact_id"] == "slot298_mixed_export_package_v1"
    assert scheme_trace["export_output_digest"] == "sha256:slot298_mixed_export_package_v1"
    assert scheme_trace["export_output_path"] == r"artifacts/cubit/slot298_mixed_export_package.json"

    stale_solver_contract = cubit_mixed_solver_route_manifest_gate(
        inv,
        {
            **solver_route,
            "solver_contract_digest": "sha256:stale-reader-contract",
        },
        expected_package_id="slot346_mixed_solver_route_v1",
        expected_solver_contract_digest="sha256:slot383-ngsolve-mixed-element-reader-contract-v1",
        require_solver_contract_artifact=True,
    )
    assert stale_solver_contract["status"] == "needs_attention"
    assert stale_solver_contract["checks"]["expected_solver_contract_digest_matches"] is False

    stale_route_convention = cubit_mixed_solver_route_manifest_gate(
        inv,
        {
            **solver_route,
            "solver_route_convention_schema_id": "coreform_value_only_mixed_route_v0",
        },
        expected_package_id="slot346_mixed_solver_route_v1",
        expected_solver_route_convention_schema_id="coreform_mixed_hex_pyramid_tet_route_convention_v1",
        require_solver_route_convention_schema=True,
    )
    assert stale_route_convention["status"] == "needs_attention"
    assert (
        stale_route_convention["checks"]["expected_solver_route_convention_schema_id_matches"]
        is False
    )
    assert stale_route_convention["checks"]["hex_primary_volume_role_recorded"] is True
    assert stale_route_convention["checks"]["pyramid_transition_role_recorded"] is True

    missing_route_convention = dict(solver_route)
    missing_route_convention.pop("solver_route_convention_schema_id")
    missing_route_convention_gate = cubit_mixed_solver_route_manifest_gate(
        inv,
        missing_route_convention,
        expected_package_id="slot346_mixed_solver_route_v1",
        require_solver_route_convention_schema=True,
    )
    assert missing_route_convention_gate["status"] == "needs_attention"
    assert (
        missing_route_convention_gate["checks"]["solver_route_convention_schema_id_recorded_when_required"]
        is False
    )

    missing_solver_contract_path = cubit_mixed_solver_route_manifest_gate(
        inv,
        {
            **solver_route,
            "solver_contract_path": "",
        },
        expected_package_id="slot346_mixed_solver_route_v1",
        require_solver_contract_artifact=True,
    )
    assert missing_solver_contract_path["status"] == "needs_attention"
    assert missing_solver_contract_path["checks"]["solver_contract_path_recorded_when_required"] is False

    bad_bnd = cubit_bnd_area_interface_gate(
        external_area=6.0,
        material_interface_area=2.0,
        ngsolve_bnd_area=6.0,
        rel_tol=1.0e-12,
    )
    bad = cubit_mixed_solver_ready_package_gate(inv, transition, export_package, bad_bnd, quality)
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["bnd_area_gate_ok"] is False

    bad_route = cubit_release_feature_routing_gate(
        [
            {
                "feature_key": "tri_tet_meshing_robustness",
                "category": "Meshing",
                "lab_route": "tet-only Cubit default",
                "validation_note": "bad route for this mixed package",
            }
        ],
        release_version="2026.6",
        source_url="https://coreform.com/coreform-cubit/release-notes/v2026-6/",
        lab_cubit_role="tet_only",
        tet_only_owner="cubit",
        required_features=("tri_tet_meshing_robustness",),
    )
    bad_package = cubit_mixed_solver_ready_package_gate(
        inv,
        transition,
        export_package,
        bnd,
        quality,
        routing_policy_gate=bad_route,
    )
    assert bad_package["status"] == "needs_attention"
    assert bad_package["checks"]["routing_policy_gate_ok"] is False

    bad_interface = cubit_mixed_interface_adjacency_gate(
        [
            {
                "surface_id": 101,
                "role": "hex_to_transition",
                "surface_kind": "triangle",
                "adjacent_material_names": ["hex_core", "pyramid_transition"],
                "adjacent_volume_kinds": ["hex", "pyramid"],
            },
            {
                "surface_id": 102,
                "role": "transition_to_tet",
                "surface_kind": "triangle",
                "adjacent_material_names": ["pyramid_transition", "tet_region"],
                "adjacent_volume_kinds": ["pyramid", "tet"],
            },
        ],
    )
    bad_interface_package = cubit_mixed_solver_ready_package_gate(
        inv,
        transition,
        export_package,
        bnd,
        quality,
        interface_adjacency_gate=bad_interface,
    )
    assert bad_interface_package["status"] == "needs_attention"
    assert bad_interface_package["checks"]["interface_adjacency_gate_ok"] is False

    stale_scheme = cubit_meshing_scheme_trace_gate(
        {
            "trace_id": "slot290_mixed_scheme_trace",
            "command_digest": "sha256:old-journal",
            "commands": [
                "imprint all",
                "merge all",
                "volume 1 scheme tetmesh",
                "volume 2 scheme tetmesh",
                "export netgen \"slot290_mixed.vol\" order 2 overwrite",
            ],
            "volume_schemes": {
                "1": "tetmesh",
                "2": "tetmesh",
            },
            "export_order": 2,
        },
        expected_trace_id="slot290_mixed_scheme_trace",
        expected_command_digest="sha256:slot290-map-tetmesh-export",
        expected_volume_schemes={"1": "map", "2": "tetmesh"},
        expected_export_order=2,
    )
    assert stale_scheme["status"] == "needs_attention"
    assert stale_scheme["checks"]["expected_command_digest_matches"] is False
    assert stale_scheme["checks"]["expected_volume_schemes_match"] is False
    bad_scheme_package = cubit_mixed_solver_ready_package_gate(
        inv,
        transition,
        export_package,
        bnd,
        quality,
        scheme_trace_gate=stale_scheme,
    )
    assert bad_scheme_package["status"] == "needs_attention"
    assert bad_scheme_package["checks"]["scheme_trace_gate_ok"] is False

    gui_headless_batch = cubit_headless_batch_quality_package_gate(
        {**batch, "command_line": r"coreform_cubit.exe artifacts/cubit/slot330_mixed_solver_ready.py", "process_mode": "gui", "gui_daemon": True},
        {**quality, "export_id": "slot178_mixed_hex_pyramid_tet", "geometry_id": "mixed_transition_fixture"},
        expected_export_id="slot178_mixed_hex_pyramid_tet",
        expected_geometry_id="mixed_transition_fixture",
        expected_element_type=None,
        export_inventory=inv,
    )
    assert gui_headless_batch["status"] == "needs_attention"
    assert gui_headless_batch["checks"]["process_mode_is_headless_batch"] is False
    assert gui_headless_batch["checks"]["nographics_flag_present"] is False
    bad_headless_package = cubit_mixed_solver_ready_package_gate(
        inv,
        transition,
        export_package,
        bnd,
        quality,
        headless_batch_quality_gate=gui_headless_batch,
    )
    assert bad_headless_package["status"] == "needs_attention"
    assert bad_headless_package["checks"]["headless_batch_quality_gate_ok"] is False

    stale_curvilinear = cubit_curvilinear_handoff_manifest_gate(
        {
            "mesh_id": "slot338_mixed_imported_mesh",
            "export_id": "slot178_mixed_hex_pyramid_tet",
            "source_mesh": {
                "kind": "third_party_mesh",
                "volume_kinds": ["hex", "pyramid", "tet"],
                "surface_kinds": ["quad", "triangle"],
            },
            "geometry_association": {
                "cad_source": "mixed_transition_fixture.step",
                "projection_policy": "project_boundary_nodes_to_cad_curves_and_surfaces",
                "boundary_ids_preserved": True,
                "projection_quality": {"max_distance": 3.0e-4, "tolerance": 1.0e-6},
            },
            "curved_export": {
                "format": "netgen_vol",
                "order": 2,
                "routing_hint": "cubit_hex_or_mixed_path",
                "implicit_element_conversion": False,
            },
            "quality": {
                "metric": "scaled_jacobian",
                "min": 0.76,
                "count": 4,
                "negative_jacobian_count": 0,
            },
            "provenance": {
                "literature_note": "stale projection row for negative-control package test",
            },
        },
        expected_export_id="slot178_mixed_hex_pyramid_tet",
    )
    assert stale_curvilinear["status"] == "needs_attention"
    assert stale_curvilinear["checks"]["projection_error_within_tolerance"] is False
    bad_curvilinear_package = cubit_mixed_solver_ready_package_gate(
        inv,
        transition,
        export_package,
        bnd,
        quality,
        curvilinear_handoff_gate=stale_curvilinear,
    )
    assert bad_curvilinear_package["status"] == "needs_attention"
    assert bad_curvilinear_package["checks"]["curvilinear_handoff_gate_ok"] is False

    stale_solver_route = cubit_mixed_solver_route_manifest_gate(
        inv,
        {
            "solver_route_package_id": "slot346_mixed_solver_route_v1",
            "routing_hint": "cubit_hex_or_mixed_path",
            "route_policy": "silently_tetize_mixed_mesh",
            "downstream_solver": "NGSolve/radia-ngsolve",
            "tet_only_owner": "cubit",
            "no_implicit_tetization": False,
            "volume_routes": [
                {"volume_kind": "hex", "solver_role": "primary_volume_fem"},
                {"volume_kind": "pyramid", "solver_role": "primary_volume_fem"},
                {"volume_kind": "tet", "solver_role": "primary_volume_fem"},
            ],
            "surface_routes": [
                {"surface_kind": "quad", "solver_role": "hex_boundary_trace"},
                {"surface_kind": "triangle", "solver_role": "tet_boundary_trace"},
            ],
        },
        expected_package_id="slot346_mixed_solver_route_v1",
    )
    assert stale_solver_route["status"] == "needs_attention"
    assert stale_solver_route["checks"]["pyramid_transition_role_recorded"] is False
    assert stale_solver_route["checks"]["no_implicit_tetization_recorded"] is False
    bad_solver_route_package = cubit_mixed_solver_ready_package_gate(
        inv,
        transition,
        export_package,
        bnd,
        quality,
        solver_route_gate=stale_solver_route,
    )
    assert bad_solver_route_package["status"] == "needs_attention"
    assert bad_solver_route_package["checks"]["solver_route_manifest_gate_ok"] is False

    missing_export = cubit_meshing_scheme_trace_gate(
        {
            "trace_id": "slot290_mixed_scheme_trace",
            "command_digest": "sha256:slot290-map-tetmesh-export",
            "commands": [
                "imprint all",
                "merge all",
                "volume 1 scheme map",
                "volume 2 scheme tetmesh",
            ],
            "volume_schemes": {
                "1": "map",
                "2": "tetmesh",
            },
            "export_order": 2,
        },
        expected_trace_id="slot290_mixed_scheme_trace",
        expected_command_digest="sha256:slot290-map-tetmesh-export",
        expected_volume_schemes={"1": "map", "2": "tetmesh"},
        expected_export_order=2,
    )
    assert missing_export["status"] == "needs_attention"
    assert missing_export["checks"]["required_command_fragments_present"] is False

    stale_scheme_output = cubit_meshing_scheme_trace_gate(
        {
            "trace_id": "slot290_mixed_scheme_trace",
            "command_digest": "sha256:slot290-map-tetmesh-export",
            "commands": [
                "imprint all",
                "merge all",
                "volume 1 scheme map",
                "volume 2 scheme tetmesh",
                "export netgen \"slot290_mixed.vol\" order 2 overwrite",
            ],
            "volume_schemes": {
                "1": "map",
                "2": "tetmesh",
            },
            "export_order": 2,
            "export_output_artifact_id": "slot298_mixed_export_package_v1",
            "export_output_digest": "sha256:old_export_package",
            "export_output_path": r"artifacts/cubit/slot298_mixed_export_package.json",
        },
        expected_trace_id="slot290_mixed_scheme_trace",
        expected_command_digest="sha256:slot290-map-tetmesh-export",
        expected_volume_schemes={"1": "map", "2": "tetmesh"},
        expected_export_order=2,
        expected_export_output_artifact_id="slot298_mixed_export_package_v1",
        expected_export_output_digest="sha256:slot298_mixed_export_package_v1",
        expected_export_output_path=r"artifacts/cubit/slot298_mixed_export_package.json",
    )
    assert stale_scheme_output["status"] == "needs_attention"
    assert stale_scheme_output["checks"]["expected_command_digest_matches"] is True
    assert stale_scheme_output["checks"]["expected_volume_schemes_match"] is True
    assert stale_scheme_output["checks"]["expected_export_output_digest_matches"] is False
    bad_scheme_output_package = cubit_mixed_solver_ready_package_gate(
        inv,
        transition,
        export_package,
        bnd,
        quality,
        scheme_trace_gate=stale_scheme_output,
    )
    assert bad_scheme_output_package["status"] == "needs_attention"
    assert bad_scheme_output_package["checks"]["scheme_trace_gate_ok"] is False

    stale_output_artifact_rows = [
        dict(row)
        for row in [
            {
                "kind": "vol",
                "path": vol_path,
                "export_id": "slot178_mixed_hex_pyramid_tet",
                "geometry_id": "mixed_transition_fixture",
                "order": 2,
                "export_output_artifact_id": "slot298_mixed_export_package_v1",
                "export_output_digest": "sha256:slot298_mixed_export_package_v1",
                "export_output_path": r"artifacts/cubit/slot298_mixed_export_package.json",
            },
            {
                "kind": "vol_sidecar",
                "path": vol_path + ".json",
                "export_id": "slot178_mixed_hex_pyramid_tet",
                "geometry_id": "mixed_transition_fixture",
                "order": 2,
                "export_output_artifact_id": "slot298_old_export_package",
                "export_output_digest": "sha256:slot298_mixed_export_package_v1",
                "export_output_path": r"artifacts/cubit/slot298_mixed_export_package.json",
            },
            {
                "kind": "raw_result",
                "path": r"artifacts/cubit/slot178_raw.json",
                "export_id": "slot178_mixed_hex_pyramid_tet",
                "geometry_id": "mixed_transition_fixture",
                "export_output_artifact_id": "slot298_mixed_export_package_v1",
                "export_output_digest": "sha256:slot298_mixed_export_package_v1",
                "export_output_path": r"artifacts/cubit/slot298_mixed_export_package.json",
            },
        ]
    ]
    stale_output_artifact = cubit_export_package_identity_gate(
        stale_output_artifact_rows,
        expected_export_id="slot178_mixed_hex_pyramid_tet",
        expected_geometry_id="mixed_transition_fixture",
        expected_export_output_artifact_id="slot298_mixed_export_package_v1",
        require_export_output_artifact=True,
        inventory=inv,
    )
    assert stale_output_artifact["status"] == "needs_attention"
    assert stale_output_artifact["checks"]["export_output_artifact_id_consistent_when_present"] is False
    assert stale_output_artifact["checks"]["expected_export_output_artifact_id_matches"] is True

    stale_output_digest_rows = [dict(row) for row in stale_output_artifact_rows]
    stale_output_digest_rows[1]["export_output_artifact_id"] = "slot298_mixed_export_package_v1"
    stale_output_digest_rows[1]["export_output_digest"] = "sha256:slot298_old_digest"
    stale_output_digest = cubit_export_package_identity_gate(
        stale_output_digest_rows,
        expected_export_output_digest="sha256:slot298_mixed_export_package_v1",
        require_export_output_artifact=True,
        inventory=inv,
    )
    assert stale_output_digest["status"] == "needs_attention"
    assert stale_output_digest["checks"]["export_output_digest_consistent_when_present"] is False
    assert stale_output_digest["checks"]["expected_export_output_digest_matches"] is True

    missing_output_path_rows = [dict(row) for row in stale_output_artifact_rows]
    for row in missing_output_path_rows:
        row.pop("export_output_path")
    missing_output_path = cubit_export_package_identity_gate(
        missing_output_path_rows,
        require_export_output_artifact=True,
        inventory=inv,
    )
    assert missing_output_path["status"] == "needs_attention"
    assert missing_output_path["checks"]["export_output_path_recorded_when_required"] is False


def test_cubit_mixed_interface_adjacency_gate_keeps_transition_faces_named():
    rows = [
        {
            "surface_id": 101,
            "role": "hex_to_transition",
            "surface_kind": "quad",
            "adjacent_material_names": ["hex_core", "pyramid_transition"],
            "adjacent_volume_kinds": ["hex", "pyramid"],
            "boundary_name": "hex_pyramid_interface",
        },
        {
            "surface_id": 102,
            "role": "transition_to_tet",
            "surface_kind": "triangle",
            "adjacent_material_names": ["pyramid_transition", "tet_region"],
            "adjacent_volume_kinds": ["pyramid", "tet"],
            "boundary_name": "pyramid_tet_interface",
        },
    ]

    gate = cubit_mixed_interface_adjacency_gate(rows)

    assert gate["policy"] == "cubit_mixed_interface_adjacency_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["required_roles_present"] is True
    assert gate["checks"]["transition_material_touches_every_interface"] is True
    assert gate["checks"]["role_surface_kinds_match"] is True
    assert gate["checks"]["role_volume_kind_pairs_match"] is True
    assert gate["roles_present"] == ["hex_to_transition", "transition_to_tet"]

    missing_transition = cubit_mixed_interface_adjacency_gate(
        [
            {**rows[0], "adjacent_material_names": ["hex_core", "air"]},
            rows[1],
        ]
    )
    assert missing_transition["status"] == "needs_attention"
    assert missing_transition["checks"]["transition_material_touches_every_interface"] is False
    assert missing_transition["transition_missing_rows"] == [1]

    wrong_surface_kind = cubit_mixed_interface_adjacency_gate(
        [
            {**rows[0], "surface_kind": "triangle"},
            rows[1],
        ]
    )
    assert wrong_surface_kind["status"] == "needs_attention"
    assert wrong_surface_kind["checks"]["role_surface_kinds_match"] is False
    assert wrong_surface_kind["role_surface_mismatch_rows"] == [1]

    wrong_adjacency = cubit_mixed_interface_adjacency_gate(
        [
            rows[0],
            {**rows[1], "adjacent_volume_kinds": ["hex", "pyramid"]},
        ]
    )
    assert wrong_adjacency["status"] == "needs_attention"
    assert wrong_adjacency["checks"]["role_volume_kind_pairs_match"] is False
    assert wrong_adjacency["role_volume_mismatch_rows"] == [2]

    missing_role = cubit_mixed_interface_adjacency_gate([rows[0]])
    assert missing_role["status"] == "needs_attention"
    assert missing_role["checks"]["required_roles_present"] is False


def test_cubit_submodel_boundary_handoff_mesh_package_gate_binds_vol_labels_to_handoff():
    inv = summarize_netgen_vol_inventory(
        LABELLED_MIXED_SUBMODEL_VOL,
        source=r"artifacts/cubit/slot250_submodel_hex_mixed.vol",
    )
    handoff = {
        "parent_model_id": "slot249_global_plate_bending_coarse_v1",
        "parent_mesh_id": "mesh_global_h2",
        "submodel_region_id": "slot250_zoom_region_tip_01",
        "local_mesh_id": "slot250_submodel_hex_mixed",
        "zoom_boundary_id": "zoom_boundary_outer",
        "boundary_trace_id": "trace_parent_to_cubit_zoom_01",
        "boundary_transfer_quantity": "displacement+slope",
        "boundary_transfer_error_estimate": 0.018,
        "boundary_transfer_error_unit": "relative",
        "local_refinement_rule": "Cubit hex-led local refinement with recorded boundary trace",
        "transition_policy": "keep pyramid bridge as an explicit conformal hex-to-tet transition",
        "target_observable_id": "tip_bending_moment",
    }

    gate = cubit_submodel_boundary_handoff_mesh_package_gate(
        inv,
        handoff,
        expected_boundary_name="zoom_boundary_outer",
        max_boundary_transfer_error=0.02,
    )

    assert gate["policy"] == "cubit_submodel_boundary_handoff_mesh_package_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["inventory_source_is_vol"] is True
    assert gate["checks"]["volume_kind_counts_recorded"] is True
    assert gate["checks"]["expected_volume_kinds_present"] is True
    assert gate["checks"]["hex_family_present_for_cubit_submodel"] is True
    assert gate["checks"]["transition_policy_recorded_when_present"] is True
    assert gate["present_volume_kinds"] == ["hex", "pyramid", "tet"]
    assert gate["checks"]["routing_hint_matches_expected"] is True
    assert gate["checks"]["not_tri_tet_only_for_cubit_submodel"] is True
    assert gate["checks"]["zoom_boundary_present_in_vol_inventory"] is True
    assert gate["checks"]["boundary_transfer_error_within_limit"] is True
    assert gate["checks"]["parent_local_mesh_identity_separated"] is True
    assert gate["checks"]["boundary_handoff_not_value_only"] is True

    wrong_boundary = cubit_submodel_boundary_handoff_mesh_package_gate(
        inv,
        {**handoff, "zoom_boundary_id": "missing_zoom_boundary"},
        max_boundary_transfer_error=0.02,
    )
    assert wrong_boundary["status"] == "needs_attention"
    assert wrong_boundary["checks"]["zoom_boundary_present_in_vol_inventory"] is False

    over_budget = cubit_submodel_boundary_handoff_mesh_package_gate(
        inv,
        {**handoff, "boundary_transfer_error_estimate": 0.031},
        expected_boundary_name="zoom_boundary_outer",
        max_boundary_transfer_error=0.02,
    )
    assert over_budget["status"] == "needs_attention"
    assert over_budget["checks"]["boundary_transfer_error_within_limit"] is False

    missing_transition_policy = cubit_submodel_boundary_handoff_mesh_package_gate(
        inv,
        {key: value for key, value in handoff.items() if key != "transition_policy"},
        expected_boundary_name="zoom_boundary_outer",
        max_boundary_transfer_error=0.02,
    )
    assert missing_transition_policy["status"] == "needs_attention"
    assert missing_transition_policy["checks"]["transition_policy_recorded_when_present"] is False

    missing_hex_family = cubit_submodel_boundary_handoff_mesh_package_gate(
        {**inv, "volume_kind_counts": {"pyramid": 1, "tet": 2}},
        handoff,
        expected_boundary_name="zoom_boundary_outer",
        max_boundary_transfer_error=0.02,
    )
    assert missing_hex_family["status"] == "needs_attention"
    assert missing_hex_family["checks"]["expected_volume_kinds_present"] is False
    assert missing_hex_family["checks"]["hex_family_present_for_cubit_submodel"] is False

    tet_only = cubit_submodel_boundary_handoff_mesh_package_gate(
        summarize_netgen_vol_inventory(TRI_TET_VOL, source="slot250_tet_only.vol"),
        handoff,
        max_boundary_transfer_error=0.02,
    )
    assert tet_only["status"] == "needs_attention"
    assert tet_only["checks"]["not_tri_tet_only_for_cubit_submodel"] is False
    assert tet_only["checks"]["boundary_names_recorded"] is False


def test_cubit_headless_batch_quality_package_gate_keeps_raw_and_quality_together():
    batch = {
        "pass": True,
        "export_id": "slot154_hex_batch_A",
        "geometry_id": "three_block_hex_quality_v1",
        "command_line": r"coreform_cubit.com -nographics -batch artifacts/cubit/slot154.py",
        "process_mode": "headless_batch",
        "batch_script": r"artifacts/cubit/slot154.py",
        "journal_policy": "batch_script_archived_no_gui_daemon",
        "gui_daemon": False,
        "exit_code": 0,
        "version": "Coreform Cubit 2025.12",
        "output_paths": [
            r"artifacts/cubit/slot154_raw.json",
            r"artifacts/cubit/slot154_quality.json",
            r"artifacts/cubit/slot154_hex.vol",
        ],
    }
    quality = {
        "policy": "cubit_quality_distribution_gate",
        "status": "ok",
        "export_id": "slot154_hex_batch_A",
        "geometry_id": "three_block_hex_quality_v1",
        "element_type": "hex",
        "count": 216,
    }
    inventory = {
        "source": r"artifacts/cubit/slot154_hex.vol",
        "volume_elements": 216,
        "volume_kind_counts": {"hex": 216},
        "surface_kind_counts": {"quad": 384},
        "routing_hint": "cubit_hex_or_mixed_path",
    }

    gate = cubit_headless_batch_quality_package_gate(
        batch,
        quality,
        expected_export_id="slot154_hex_batch_A",
        expected_geometry_id="three_block_hex_quality_v1",
        export_inventory=inventory,
    )
    assert gate["policy"] == "cubit_headless_batch_quality_package_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["headless_command_recorded"] is True
    assert gate["checks"]["process_mode_is_headless_batch"] is True
    assert gate["checks"]["nographics_flag_present"] is True
    assert gate["checks"]["batch_flag_present"] is True
    assert gate["checks"]["gui_daemon_disabled"] is True
    assert gate["checks"]["batch_script_appears_in_command"] is True
    assert gate["checks"]["process_exit_code_success_or_documented"] is True
    assert gate["checks"]["journal_policy_records_batch_not_gui_daemon"] is True
    assert gate["checks"]["export_id_matches_quality"] is True
    assert gate["checks"]["quality_count_positive"] is True
    assert gate["checks"]["export_inventory_source_in_output_paths"] is True
    assert gate["checks"]["export_inventory_volume_elements_positive"] is True
    assert gate["checks"]["export_inventory_count_matches_quality"] is True
    assert gate["checks"]["export_inventory_contains_quality_element"] is True
    assert gate["checks"]["export_inventory_not_tri_tet_only_for_cubit_hex_route"] is True
    assert gate["export_inventory_is_tri_tet_only"] is False

    stale_quality = {**quality, "geometry_id": "old_geometry"}
    stale_gate = cubit_headless_batch_quality_package_gate(batch, stale_quality)
    assert stale_gate["status"] == "needs_attention"
    assert stale_gate["checks"]["geometry_id_matches_quality"] is False

    gui_batch = {**batch, "command_line": r"coreform_cubit.exe artifacts/cubit/slot154.py", "headless": False}
    gui_gate = cubit_headless_batch_quality_package_gate(gui_batch, quality)
    assert gui_gate["status"] == "needs_attention"
    assert gui_gate["checks"]["headless_command_recorded"] is False
    assert gui_gate["checks"]["nographics_flag_present"] is False
    assert gui_gate["checks"]["batch_flag_present"] is False

    documented_nonzero_batch = {
        **batch,
        "exit_code": 3,
        "process_exit_note": "Known headless startup warning after valid batch artifacts were written.",
    }
    documented_nonzero_gate = cubit_headless_batch_quality_package_gate(documented_nonzero_batch, quality)
    assert documented_nonzero_gate["status"] == "ok"
    assert documented_nonzero_gate["checks"]["process_exit_code_success_or_documented"] is True

    artifact_evidence_batch = {
        **batch,
        "exit_code": 1,
        "process_exit_note": "Batch process ended nonzero after valid artifact evidence was written; do not claim solver-ready from exit status alone.",
        "process_exit_policy": "artifact_evidence_over_process_exit",
        "solver_ready_claimed": False,
    }
    artifact_evidence_gate = cubit_headless_batch_quality_package_gate(
        artifact_evidence_batch,
        quality,
        export_inventory=inventory,
        expected_process_exit_policy="artifact_evidence_over_process_exit",
    )
    assert artifact_evidence_gate["status"] == "ok"
    assert artifact_evidence_gate["process_exit_policy"] == "artifact_evidence_over_process_exit"
    assert artifact_evidence_gate["solver_ready_claimed"] is False
    assert artifact_evidence_gate["checks"]["process_exit_policy_recorded_when_nonzero"] is True
    assert artifact_evidence_gate["checks"]["expected_process_exit_policy_matches"] is True
    assert artifact_evidence_gate["checks"]["nonzero_exit_does_not_claim_solver_ready"] is True
    assert artifact_evidence_gate["checks"]["process_exit_code_success_or_documented"] is True

    solver_ready_from_nonzero_exit = {
        **artifact_evidence_batch,
        "solver_ready_claimed": True,
    }
    solver_ready_from_nonzero_gate = cubit_headless_batch_quality_package_gate(
        solver_ready_from_nonzero_exit,
        quality,
        export_inventory=inventory,
        expected_process_exit_policy="artifact_evidence_over_process_exit",
    )
    assert solver_ready_from_nonzero_gate["status"] == "needs_attention"
    assert solver_ready_from_nonzero_gate["checks"]["nonzero_exit_does_not_claim_solver_ready"] is False
    assert solver_ready_from_nonzero_gate["checks"]["process_exit_code_success_or_documented"] is False

    wrong_exit_policy_gate = cubit_headless_batch_quality_package_gate(
        {**artifact_evidence_batch, "process_exit_policy": "ignore_exit_code"},
        quality,
        export_inventory=inventory,
        expected_process_exit_policy="artifact_evidence_over_process_exit",
    )
    assert wrong_exit_policy_gate["status"] == "needs_attention"
    assert wrong_exit_policy_gate["checks"]["expected_process_exit_policy_matches"] is False
    assert wrong_exit_policy_gate["checks"]["process_exit_code_success_or_documented"] is False

    daemon_batch = {**batch, "gui_daemon": True, "exit_code": 3}
    daemon_gate = cubit_headless_batch_quality_package_gate(daemon_batch, quality)
    assert daemon_gate["status"] == "needs_attention"
    assert daemon_gate["checks"]["gui_daemon_disabled"] is False
    assert daemon_gate["checks"]["process_exit_code_success_or_documented"] is False

    empty_quality = {**quality, "count": 0}
    empty_gate = cubit_headless_batch_quality_package_gate(batch, empty_quality)
    assert empty_gate["status"] == "needs_attention"
    assert empty_gate["checks"]["quality_count_positive"] is False

    stale_inventory = {**inventory, "volume_kind_counts": {"hex": 128}}
    stale_inventory_gate = cubit_headless_batch_quality_package_gate(
        batch,
        quality,
        export_inventory=stale_inventory,
    )
    assert stale_inventory_gate["status"] == "needs_attention"
    assert stale_inventory_gate["checks"]["export_inventory_count_matches_quality"] is False

    tri_tet_inventory = {
        **inventory,
        "volume_kind_counts": {"tet": 216},
        "surface_kind_counts": {"triangle": 384},
        "routing_hint": "netgen_tri_tet_path",
        "is_tri_tet_only": True,
    }
    tri_tet_gate = cubit_headless_batch_quality_package_gate(
        batch,
        quality,
        export_inventory=tri_tet_inventory,
    )
    assert tri_tet_gate["status"] == "needs_attention"
    assert tri_tet_gate["export_inventory_is_tri_tet_only"] is True
    assert tri_tet_gate["checks"]["export_inventory_contains_quality_element"] is False
    assert tri_tet_gate["checks"]["export_inventory_not_tri_tet_only_for_cubit_hex_route"] is False


def test_netgen_workflow_records_mapped_hex_quality_replay_gate():
    doc = get_netgen_documentation("overview")

    assert "mapped hex quality replay gate" in doc
    assert "384 hexes" in doc
    assert "scaled-Jacobian min `0.9999999999999999`" in doc
    assert "cubit_hex_quality_gate" in doc
    assert "quality distribution replay gate" in doc
    assert "cubit_quality_distribution_gate" in doc
    assert "p05/p50/p95" in doc
    assert "surfaceelementsuv" in doc
    assert ".vol label metadata gate" in doc
    assert "cubit_vol_label_metadata_gate" in doc
    assert "headless batch quality package gate" in doc
    assert "cubit_headless_batch_quality_package_gate" in doc
    assert "zero quality count" in doc
    assert "plugin freshness warning" in doc
    assert "Keep the batch process exit status separate" in doc
    assert "plugin-specific export evidence" in doc
    assert "export_inventory" in doc
    assert "slot194 lesson" in doc
    assert "volume-kind" in doc
    assert "Slot314 tightens the noisy-exit path" in doc
    assert "process_exit_policy=artifact_evidence_over_process_exit" in doc
    assert "Slot369 adds a mesh-quality ledger identity gate" in doc
    assert "cubit_mesh_quality_ledger_identity_gate" in doc
    assert "mesh_quality_artifact_id" in doc
    assert "quality_metric_set_id" in doc
    assert "negative-Jacobian count" in doc
    assert "stale quality digest" in doc
    assert "tri/tet-only inventory paired with a Cubit\nhex-led quality row" in doc
    assert "Slot330 binds that headless/process evidence to the mixed solver-ready package" in doc
    assert "headless_batch_quality_gate" in doc
    assert "disabled GUI daemon" in doc
    assert "Slot362 binds Cubit scheme traces to the actual exported mesh artifact" in doc
    assert "A fresh journal\ntrace with an old `.vol` digest" in doc
    assert "Slot346 adds the downstream solver-route manifest" in doc
    assert "cubit_mixed_solver_route_manifest_gate" in doc
    assert "no_implicit_tetization=true" in doc
    assert "Slot383 adds downstream solver-reader contract identity" in doc
    assert "solver_contract_artifact_id" in doc
    assert "downstream_solver_contract_*" in doc
    assert "stale solver-reader digest" in doc
    assert "Slot418 adds solver-route convention schema identity" in doc
    assert "hex\nprimary, pyramid transition, tet compatibility/subregion" in doc
    assert "value-only route convention" in doc
    assert "Slot404 binds the emitted Netgen `.vol` file" in doc
    assert "require_vol_sidecar_inventory_counts=True" in doc
    assert "sidecar `n_elements`" in doc
    assert "12 volume elements, 13\npoints, and order 1" in doc
    assert "do not infer\norder from the filename" in doc
    assert "omitting the\n`order` argument" in doc
    assert "sidecar\nwith `order=2`" in doc
    assert "Slot390 extends the same quality ledger" in doc
    assert "parameter_set_artifact_id" in doc
    assert "parameter_set_digest" in doc
    assert "parameter_set_path" in doc
    assert "objective_observable_id" in doc
    assert "objective_observable_family" in doc
    assert "stale parameter-set digests" in doc
    assert "wrong\nobjective families" in doc
    assert "Slot425 adds postprocess-row convention schema identity" in doc
    assert "mesh_quality_postprocess_row_convention_schema_id" in doc
    assert "require_quality_postprocess_row_convention_schema=True" in doc
    assert "stale scalar-row conventions" in doc
    assert "Slot432 adds component-basis schema identity" in doc
    assert "mesh_quality_component_basis_schema_id" in doc
    assert "require_quality_component_basis_schema=True" in doc
    assert "stale scalar-value component bases" in doc
    assert "solver_ready_claimed=false" in doc
    assert "nonzero exit with `solver_ready_claimed=true`" in doc
    assert "Slot234 adds the route-separation check" in doc
    assert "not be tri/tet-only" in doc
    assert "submodel boundary handoff mesh package gate" in doc
    assert "cubit_submodel_boundary_handoff_mesh_package_gate" in doc
    assert "zoom_boundary_id" in doc
    assert "quality element kind" in doc


def test_netgen_workflow_records_coreform_2026_6_release_routing_gate():
    doc = get_netgen_documentation("overview")

    assert "Coreform Cubit 2026.6 release routing gate" in doc
    assert "released on 1 June 2026" in doc
    assert "anisotropic tetrahedral meshing" in doc
    assert "cohesive element generation" in doc
    assert "Tetra10/Tri6 Jacobian" in doc
    assert "64-bit Exodus IDs" in doc
    assert "Python 3.12 runtime" in doc
    assert "cubit_release_feature_routing_gate" in doc
    assert "hex-led and mixed hex+pyramid+tet lane" in doc
    assert "tet-only\neducation stays on Netgen/OCC" in doc
    assert "rather than `.vol` parser\nrelaxation" in doc
    assert "close a Coreform/Cubit slot only when the headless batch\ncommand" in doc
    assert "result_artifact_id" in doc
    assert "mesh-quality evidence" in doc
    assert "third-party curvilinear handoff manifest gate" in doc
    assert "cubit_curvilinear_handoff_manifest_gate" in doc
    assert "no implicit tetization" in doc
    assert "projection_quality.max_distance <= projection_quality.tolerance" in doc
    assert "negative_jacobian_count = 0" in doc
    assert "hex-led or mixed curvilinear\nhandoff lane" in doc
    assert "Slot226 rechecked the installed-version lane" in doc
    assert "coreform_cubit.com -version" in doc
    assert "status: ValidStudent" in doc
    assert "Coreform Cubit Version 2025.12 Build 3d8d3af7" in doc
    assert "Slot354 tightens installed-version evidence" in doc
    assert "coreform_cubit.exe` stub" in doc
    assert "`.com -version` probe command to use that same recorded binary" in doc


def test_netgen_workflow_records_mixed_order_series_gate():
    doc = get_netgen_documentation("overview")

    assert "mixed hex+pyramid+tet order-series gate" in doc
    assert "cubit_mixed_order_series_inventory_gate" in doc
    assert "1 hex, 1 pyramid, 10 tets" in doc
    assert "curvedelements" in doc
    assert "cubit_hex_or_mixed_path" in doc
    assert "zero material volume" in doc
    assert "mixed transition metadata gate" in doc
    assert "cubit_mixed_transition_metadata_gate" in doc
    assert "pyramid transition block" in doc
    assert "Slot282 adds the interface-adjacency ledger" in doc
    assert "cubit_mixed_interface_adjacency_gate" in doc
    assert "hex-pyramid interface" in doc
    assert "mixed solver-ready package gate" in doc
    assert "cubit_mixed_solver_ready_package_gate" in doc
    assert "cubit_ngsolve_bnd_area_includes_material_interfaces_once" in doc
    assert "cubit_curvilinear_handoff_manifest_gate" in doc
    assert "projection tolerance" in doc
    assert "Slot338 adds the literature-driven curvilinear handoff row" in doc
    assert "Pyramid cells are not display-only\nmesh noise" in doc
    assert "Slot202 adds the routing-policy row" in doc
    assert "routing_policy_gate" in doc
    assert "tet-only default route" in doc


def test_netgen_workflow_records_live_mixed_bnd_surfaceelementsuv_gate():
    doc = get_netgen_documentation("overview")

    assert "live mixed hex+pyramid+tet NGSolve BND gate" in doc
    assert "surfaceelementsuv" in doc
    assert "volume matched `2e-9`" in doc
    assert "matched `11e-6`, not the\nexternal brick area `10e-6`" in doc
    assert "split material interface of area\n`1e-6` is included once" in doc


def test_netgen_workflow_records_two_block_interface_bnd_gate():
    doc = get_netgen_documentation("overview")

    assert "two-block hex interface gate" in doc
    assert "144 hexes" in doc
    assert "external area 42" in doc
    assert "material interface of area 6" in doc
    assert "Integrate(1, mesh, BND)" in doc
    assert "returns 48, not 42" in doc
    assert "get_relatives" in doc
    assert "slot90 replay" in doc
    assert "Coreform Cubit `2025.12`" in doc
    assert "Coreform 2026.6 release notes" in doc
    assert "anisotropic tetrahedral meshing" in doc
    assert "cohesive element generation" in doc
    assert "higher-order Tetra10/Tri6 Jacobian metrics" in doc
    assert "216 hexes, 270 quad surface records" in doc
    assert "BND equals external area plus the shared\ninterface once" in doc
    assert "avoid top-level\nmulti-line dict literals" in doc


def test_netgen_workflow_records_three_block_hex_quality_gate():
    doc = get_netgen_documentation("overview")

    assert "three-block hex quality gate" in doc
    assert "216 hexes" in doc
    assert "external area 57" in doc
    assert "total area is 12" in doc
    assert "scaled Jacobian `0.9999999999999999`" in doc
    assert "returns 69" in doc
    assert "cubit_bnd_area_interface_gate" in doc


def test_netgen_workflow_records_curved_hex_cylinder_order_series_gate():
    doc = get_netgen_documentation("overview")

    assert "curved hex cylinder order-series gate" in doc
    assert "volume all scheme auto" in doc
    assert "Trouble finding logical box" in doc
    assert "224 hexes" in doc
    assert "1.57e-7" in doc
    assert "6.99e-8" in doc
    assert "material-interface area" in doc


def test_netgen_workflow_records_annular_hex_tube_field_gate():
    doc = get_netgen_documentation("overview")

    assert "annular hex tube capacitance field gate" in doc
    assert "volume all scheme sweep" in doc
    assert "864 hexes" in doc
    assert "min scaled Jacobian `0.9951847266721953`" in doc
    assert "orders 1, 3,\nand 5" in doc
    assert "capacitance per length `2*pi/log(b/a)`" in doc
    assert "capacitance-per-length rel err `1.00e-6`" in doc
    assert "volume rel err `4.11e-9`" in doc
    assert "BND-area rel err `2.27e-9`" in doc
    assert "capacitance-per-length rel err stayed\nat `1.02e-6`" in doc
    assert "abs(n dot rhat) > 0.9" in doc
    assert "C++ plugin is out of\ndate" in doc
    assert "curved hex export preserves a field quantity" in doc
