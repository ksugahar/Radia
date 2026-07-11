import copy
import json

from radia_mcp.build123d.brep_roundtrip_gate import brep_mass_topology_roundtrip_gate
from radia_mcp.build123d.server import build123d_brep_mass_topology_roundtrip_gate as mcp_gate
from radia_mcp.cubit.api_reference import get_api_reference


REFERENCE = {
    "center_semantics": "mass_centroid",
    "volume": 94.77361530411152,
    "surface_area": 211.51619773911386,
    "bbox_min": [-11.000000096893165, -1.5090618034959675, -1.0040745362639427e-7],
    "bbox_max": [11.000000097053237, 1.5090618034959618, 5.643719262371525],
    "center_of_mass": [-2.0326094477530603e-8, 2.5561918787657665e-12, 3.950937975761567],
    "vertex_count": 18,
    "edge_count": 27,
    "face_count": 11,
    "solid_count": 1,
    "boundary_euler_characteristic": 2,
}


def _measured(mode):
    return {
        "import_mode": mode,
        "center_semantics": "mass_centroid",
        "volume": 94.77656568071431,
        "surface_area": 211.86923296460708,
        "bbox_min": [-11.0971548062, -1.511335072149, -0.001000014901],
        "bbox_max": [11.097154807135, 1.511335072149, 5.651410736144],
        "center_of_mass": [-0.0001055111683, -2.95796e-9, 3.951004521621],
        "vertex_count": 18,
        "edge_count": 27,
        "face_count": 11,
        "solid_count": 1,
        "boundary_euler_characteristic": 2,
    }


def test_brep_roundtrip_accepts_heal_noheal_with_mass_centroid_and_euler_invariant():
    result = brep_mass_topology_roundtrip_gate(
        REFERENCE,
        [_measured("heal"), _measured("noheal")],
        expected_volume=94.77361455046953,
    )
    assert result["status"] == "ok"
    assert result["reference_boundary_euler_characteristic"] == 2
    assert max(row["volume_relative_error"] for row in result["comparisons"]) < 4.0e-5


def test_brep_roundtrip_rejects_bbox_center_mislabeled_as_mass_centroid_and_genus_drift():
    bad = _measured("noheal")
    bad["center_semantics"] = "bbox_center"
    bad["center_of_mass"] = [0.0, 0.0, 2.825]
    bad["boundary_euler_characteristic"] = 0
    result = brep_mass_topology_roundtrip_gate(REFERENCE, [_measured("heal"), bad])
    assert result["status"] == "needs_attention"
    assert result["checks"]["measured_centers_are_mass_centroids"] is False
    assert result["checks"]["all_euler_characteristics_match"] is False


def test_brep_roundtrip_mcp_dispatches_json():
    result = json.loads(mcp_gate(
        json.dumps(REFERENCE),
        json.dumps([_measured("heal"), _measured("noheal")]),
        expected_volume=94.77361455046953,
    ))
    assert result["status"] == "ok"


def test_cubit_geometry_reference_distinguishes_center_from_mass_centroid():
    docs = get_api_reference("geometry_queries")
    assert "representative geometric center" in docs
    assert "volume_id).centroid()" in docs
    assert "CenterOf.MASS" in docs
