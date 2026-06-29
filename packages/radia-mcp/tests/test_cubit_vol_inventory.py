import json

from radia_mcp.cubit.server import cubit_docs, cubit_vol_inventory
from radia_mcp.cubit.knowledge.netgen_workflow import get_netgen_documentation
from radia_mcp.cubit.vol_inventory import summarize_netgen_vol_inventory


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


CURVED_MIXED_VOL = MIXED_VOL + """
curvedelements
2
1 2 3 4
5 6 7 8
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


def test_cubit_vol_inventory_keeps_curved_mixed_mesh_on_cubit_route():
    inv = summarize_netgen_vol_inventory(CURVED_MIXED_VOL)

    assert inv["curvedelements_present"] is True
    assert inv["volume_kind_counts"] == {"hex": 1, "pyramid": 1, "tet": 1}
    assert inv["surface_kind_counts"] == {"quad": 1, "triangle": 1}
    assert inv["routing_hint"] == "cubit_hex_or_mixed_path"
    assert "curvedelements section can grow with order" in inv["order_series_policy"]


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


def test_netgen_workflow_records_mapped_hex_brick_area_gate():
    doc = get_netgen_documentation("overview")

    assert "mapped hex brick volume/area gate" in doc
    assert "volume 1 scheme map" in doc
    assert "192 hexes" in doc
    assert "surface-area rel err `2.05e-15`" in doc
    assert "Avoid multi-line `dict(...)`" in doc


def test_netgen_workflow_records_mixed_order_series_gate():
    doc = get_netgen_documentation("overview")

    assert "mixed hex+pyramid+tet order-series gate" in doc
    assert "1 hex, 1 pyramid, 10 tets" in doc
    assert "curvedelements" in doc
    assert "cubit_hex_or_mixed_path" in doc
    assert "zero material volume" in doc


def test_netgen_workflow_records_two_block_interface_bnd_gate():
    doc = get_netgen_documentation("overview")

    assert "two-block hex interface gate" in doc
    assert "144 hexes" in doc
    assert "external area 42" in doc
    assert "material interface of area 6" in doc
    assert "Integrate(1, mesh, BND)" in doc
    assert "returns 48, not 42" in doc
    assert "get_relatives" in doc


def test_netgen_workflow_records_three_block_hex_quality_gate():
    doc = get_netgen_documentation("overview")

    assert "three-block hex quality gate" in doc
    assert "216 hexes" in doc
    assert "external area 57" in doc
    assert "total area is 12" in doc
    assert "scaled Jacobian `0.9999999999999999`" in doc
    assert "returns 69" in doc


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
    assert "capacitance per length `2*pi/log(b/a)`" in doc
    assert "capacitance-per-length rel err `1.00e-6`" in doc
    assert "curved hex export preserves a field quantity" in doc
