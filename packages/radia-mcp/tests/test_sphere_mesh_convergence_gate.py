import copy
import json
import math

from radia_mcp.radia_ngsolve.server import (
    linear_sphere_geometry_convergence_gate as mcp_gate,
)
from radia_mcp.radia_ngsolve.sphere_mesh_convergence_gate import (
    linear_sphere_geometry_convergence_gate,
)


def _rows():
    volume_exact = 4.0 * math.pi / 3.0
    area_exact = 4.0 * math.pi
    rows = []
    for level, volume_error, area_error in (
        (0, 0.40, 0.24),
        (1, 0.12, 0.07),
        (2, 0.032, 0.018),
        (3, 0.008, 0.0045),
        (4, 0.002, 0.0011),
    ):
        rows.append(
            {
                "level": level,
                "points": 10 * 4**level + 3,
                "triangles": 20 * 4**level,
                "tets": 20 * 4**level,
                "volume": volume_exact * (1.0 - volume_error),
                "surface_area": area_exact * (1.0 - area_error),
                "boundary_orientation": "outward",
                "maximum_surface_radius_error": 2.0e-16,
                "volume_reader_relative_error": 2.0e-14,
                "surface_reader_relative_error": 7.0e-15,
            }
        )
    return rows


def _replay(rows):
    row = rows[3]
    return {
        "level": row["level"],
        "points": row["points"],
        "triangles": row["triangles"],
        "tets": row["tets"],
        "volume": row["volume"],
        "surface_area": row["surface_area"],
    }


def test_accepts_outward_second_order_sphere_family():
    rows = _rows()
    result = linear_sphere_geometry_convergence_gate(
        rows,
        analytic_volume=4.0 * math.pi / 3.0,
        analytic_surface_area=4.0 * math.pi,
        replay=_replay(rows),
    )
    assert result["status"] == "ok"
    assert result["metrics"]["volume_observed_orders"][-1] > 1.8


def test_rejects_inward_stale_and_nonconvergent_family():
    rows = _rows()
    rows[2]["boundary_orientation"] = "inward"
    rows[3]["volume_reader_relative_error"] = 1.0e-4
    rows[4]["surface_area"] = rows[3]["surface_area"] * 0.99
    replay = _replay(_rows())
    replay["volume"] *= 1.01
    result = linear_sphere_geometry_convergence_gate(
        rows,
        analytic_volume=4.0 * math.pi / 3.0,
        analytic_surface_area=4.0 * math.pi,
        replay=replay,
    )
    assert result["status"] == "needs_attention"
    assert result["checks"]["all_boundaries_are_outward"] is False
    assert result["checks"]["independent_readers_agree"] is False
    assert result["checks"]["independent_replay_is_exact"] is False


def test_mcp_dispatch_and_bad_shape():
    rows = _rows()
    result = json.loads(
        mcp_gate(rows, 4.0 * math.pi / 3.0, 4.0 * math.pi, _replay(rows))
    )
    assert result["status"] == "ok"
    bad = json.loads(mcp_gate([], 1.0, 1.0, {}))
    assert bad["status"] == "invalid_input"
