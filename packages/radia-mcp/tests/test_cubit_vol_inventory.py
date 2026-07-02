import json

from radia_mcp.cubit.server import cubit_docs, cubit_vol_inventory
from radia_mcp.cubit.knowledge.netgen_workflow import get_netgen_documentation
from radia_mcp.cubit.vol_inventory import (
    cubit_bnd_area_interface_gate,
    cubit_element_quality_gate,
    cubit_export_package_identity_gate,
    cubit_hex_quality_gate,
    cubit_mass_property_sidecar_gate,
    cubit_mixed_transition_metadata_gate,
    cubit_quality_distribution_gate,
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
    assert gate["checks"]["transition_kinds_present"] is True
    assert gate["checks"]["routing_hint_is_cubit_mixed"] is True
    assert gate["checks"]["transition_material_label_present"] is True
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
    assert "bbox size `[1.5, 2.0, 0.75]`" in doc
    assert "Cubit export package identity gate" in doc
    assert "cubit_export_package_identity_gate" in doc
    assert "`export_id`, `geometry_id`, order, and routing hint" in doc
    assert "old sidecar or a raw JSON from a different geometry" in doc


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
            "bounding_box": {"size": [1.5, 2.0, 0.75]},
        }
    ]
    gate = cubit_mass_property_sidecar_gate(
        rows,
        expected_total_volume=2.25,
        expected_total_area=11.25,
        expected_bbox_size=[1.5, 2.0, 0.75],
        rel_tol=1e-12,
        abs_tol=1e-12,
    )

    assert gate["status"] == "ok"
    assert gate["policy"] == "cubit_mass_property_sidecar_gate"
    assert gate["checks"]["total_volume_expected_ok"] is True
    assert gate["checks"]["total_area_expected_ok"] is True
    assert gate["checks"]["bbox_size_expected_ok"] is True
    assert gate["volume_rel_error"] == 0.0
    assert gate["area_rel_error"] == 0.0
    assert "Volume alone is not enough" in " ".join(gate["notes"])

    bad = cubit_mass_property_sidecar_gate(
        [{"name": "hex_brick", "volume": 2.25, "area": 10.0, "bounding_box": {"size": [1.5, 2.0, 0.7]}}],
        expected_total_volume=2.25,
        expected_total_area=11.25,
        expected_bbox_size=[1.5, 2.0, 0.75],
        rel_tol=1e-12,
        abs_tol=1e-12,
    )
    assert bad["status"] == "needs_attention"
    assert bad["checks"]["total_volume_expected_ok"] is True
    assert bad["checks"]["total_area_expected_ok"] is False
    assert bad["checks"]["bbox_size_expected_ok"] is False


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
    inventory = {"source": vol_path, "routing_hint": "cubit_hex_or_mixed_path"}

    gate = cubit_export_package_identity_gate(
        artifacts,
        expected_export_id="slot146_hex_brick_o3",
        expected_geometry_id="hex_brick_v1",
        expected_order=3,
        expected_routing_hint="cubit_hex_or_mixed_path",
        inventory=inventory,
    )

    assert gate["policy"] == "cubit_export_package_identity_gate"
    assert gate["status"] == "ok"
    assert gate["checks"]["vol_sidecar_pairs_vol"] is True
    assert gate["checks"]["geometry_id_unique"] is True
    assert gate["checks"]["inventory_source_matches_vol"] is True
    assert gate["checks"]["inventory_routing_hint_matches_expected"] is True

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

    wrong_geometry = [dict(row) for row in artifacts]
    wrong_geometry[-1]["geometry_id"] = "hex_brick_v2"
    bad_geometry = cubit_export_package_identity_gate(wrong_geometry)
    assert bad_geometry["status"] == "needs_attention"
    assert bad_geometry["checks"]["geometry_id_unique"] is False


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


def test_netgen_workflow_records_mixed_order_series_gate():
    doc = get_netgen_documentation("overview")

    assert "mixed hex+pyramid+tet order-series gate" in doc
    assert "1 hex, 1 pyramid, 10 tets" in doc
    assert "curvedelements" in doc
    assert "cubit_hex_or_mixed_path" in doc
    assert "zero material volume" in doc
    assert "mixed transition metadata gate" in doc
    assert "cubit_mixed_transition_metadata_gate" in doc
    assert "pyramid transition block" in doc


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
