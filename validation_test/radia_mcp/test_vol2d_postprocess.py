"""Independent numerical and Gmsh checks for the 2-D postprocess contract."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

gmsh = pytest.importorskip("gmsh")
pytest.importorskip("ngsolve")

from radia_mcp.radia_ngsolve.vol2d_circuit import MU0, write_structured_rect_vol
from radia_mcp.radia_ngsolve.vol2d_dynamics import assemble_vol2d_dynamics
from radia_mcp.radia_ngsolve.vol2d_postprocess import analyze_vol2d_postprocess


def _request() -> dict:
    path = Path(r"C:\temp\radia_vol2d_postprocess_validation.vol")
    write_structured_rect_vol(
        path,
        x0=0.0,
        x1=1.0,
        y0=0.0,
        y1=1.0,
        nx=10,
        ny=10,
        material="coil",
    )
    return {
        "operation": "solve",
        "vol_text": path.read_text(encoding="utf-8"),
        "source_name": path.name,
        "element_family": "P1",
        "formulation": "planar",
        "dirichlet_boundaries": ["bottom", "right", "top", "left"],
        "materials": {
            "coil": {
                "permeability_h_per_m": MU0,
                "conductivity_s_per_m": 0.0,
            }
        },
        "branches": [{"name": "winding", "material": "coil", "turns": 1.0}],
        "current_rows_a": [[1.75]],
        "point_probes": [
            {"name": "center", "coordinates_m": [0.5, 0.5]},
            {"name": "off_axis", "coordinates_m": [0.6, 0.55]},
        ],
        "path_probes": [
            {
                "name": "forward",
                "points_m": [[0.2, 0.35], [0.8, 0.35]],
                "samples_per_segment": 20,
            },
            {
                "name": "reverse",
                "points_m": [[0.8, 0.35], [0.2, 0.35]],
                "samples_per_segment": 20,
            },
        ],
        "region_materials": ["coil"],
        "model_depth_m": 0.12,
        "export_basename": "vol2d_post_validation",
    }


def _axis_request() -> dict:
    path = Path(r"C:\temp\radia_vol2d_postprocess_validation_axi.vol")
    write_structured_rect_vol(
        path,
        x0=0.1,
        x1=0.3,
        y0=-0.1,
        y1=0.1,
        nx=8,
        ny=8,
        material="coil",
    )
    return {
        "operation": "solve",
        "vol_text": path.read_text(encoding="utf-8"),
        "source_name": path.name,
        "element_family": "P1",
        "formulation": "axisymmetric_henrotte",
        "dirichlet_boundaries": ["bottom", "right", "top", "left"],
        "materials": {
            "coil": {
                "permeability_h_per_m": MU0,
                "conductivity_s_per_m": 0.0,
            }
        },
        "branches": [{"name": "winding", "material": "coil", "turns": 1.0}],
        "current_rows_a": [[1.0], [2.0]],
        "point_probes": [
            {"name": "center", "coordinates_m": [0.2, 0.0]},
            {"name": "upper", "coordinates_m": [0.2, 0.04]},
        ],
        "path_probes": [
            {
                "name": "open_path",
                "points_m": [[0.14, -0.04], [0.26, -0.04], [0.26, 0.04]],
                "samples_per_segment": 4,
            }
        ],
        "region_materials": ["coil"],
        "export_basename": "vol2d_post_validation_axi",
    }


def test_energy_identity_and_oriented_path_reversal() -> None:
    request = _request()
    result = analyze_vol2d_postprocess(request)
    contract = result["result_contract"]
    row = contract["sweep_rows"][0]
    operators = assemble_vol2d_dynamics(request)
    state = np.asarray(contract["field_state_rows"][0], dtype=float)
    stiffness = np.asarray(operators["assembly"]["field_matrix"], dtype=float)
    matrix_energy = 0.5 * float(state @ stiffness @ state) * request["model_depth_m"]
    assert row["total_magnetic_energy_j"] == pytest.approx(matrix_energy, rel=2.0e-11)
    forward, reverse = row["path_probes"]
    assert reverse["b_tangent_line_integral_t_m"] == pytest.approx(
        -forward["b_tangent_line_integral_t_m"], rel=2.0e-11, abs=1.0e-18
    )
    assert reverse["h_tangent_line_integral_a"] == pytest.approx(
        -forward["h_tangent_line_integral_a"], rel=2.0e-11, abs=1.0e-12
    )


def test_gmsh_v41_reopens_with_three_field_views() -> None:
    result = analyze_vol2d_postprocess(_request())
    target = Path(r"C:\temp\vol2d_post_validation_reopen.msh")
    target.write_text(result["exports"]["gmsh_msh"]["content"], encoding="utf-8")
    gmsh.initialize()
    try:
        gmsh.open(str(target))
        assert len(gmsh.model.mesh.getNodes()[0]) == 121
        assert sum(len(tags) for tags in gmsh.model.mesh.getElements(2)[1]) == 200
        assert len(gmsh.view.getTags()) == 3
    finally:
        gmsh.finalize()


def test_negative_current_preserves_energy_and_reverses_b() -> None:
    request = _request()
    request["current_rows_a"] = [[1.0], [-1.0]]
    result = analyze_vol2d_postprocess(request)["result_contract"]["sweep_rows"]
    assert result[1]["total_magnetic_energy_j"] == pytest.approx(
        result[0]["total_magnetic_energy_j"], rel=2.0e-11
    )
    positive = result[0]["point_probes"][1]["b_t"]
    negative = result[1]["point_probes"][1]["b_t"]
    assert negative == pytest.approx([-positive[0], -positive[1]], rel=2.0e-11)
    assert math.isfinite(result[0]["maximum_point_probe_b_t"])


def test_axisymmetric_full_revolution_volume_and_energy_scaling() -> None:
    result = analyze_vol2d_postprocess(_axis_request())["result_contract"]
    first, second = result["sweep_rows"]
    volume = first["region_integrals"][0]["full_revolution_volume_m3"]
    expected = math.pi * (0.3**2 - 0.1**2) * 0.2
    assert volume == pytest.approx(expected, rel=2.0e-12)
    assert second["total_magnetic_energy_j"] == pytest.approx(
        4.0 * first["total_magnetic_energy_j"], rel=2.0e-11
    )
    assert result["request_contract"]["coordinate_contract"]["field_components"] == [
        "B_r",
        "B_z",
    ]


def test_axisymmetric_probe_at_r_zero_fails_closed() -> None:
    request = _axis_request()
    request["point_probes"][0]["coordinates_m"] = [0.0, 0.0]
    with pytest.raises(ValueError, match="r >"):
        analyze_vol2d_postprocess(request)
