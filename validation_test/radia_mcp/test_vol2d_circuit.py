from __future__ import annotations

import asyncio
import json
from pathlib import Path

import numpy as np
import pytest

from radia_mcp.motor.server import motor_vol2d_circuit_analysis
from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol
from radia_mcp.radia_ngsolve.vol2d_circuit import (
    MESH_SCHEMA,
    analyze_vol2d_circuit,
    assemble_vol2d_field,
    inspect_netgen_2d_vol,
    write_structured_rect_vol,
)


def _rect_text(tmp_path: Path, *, quads: bool) -> str:
    path = tmp_path / ("rect_quad.vol" if quads else "rect_tri.vol")
    write_structured_rect_vol(
        path,
        x0=0.1,
        x1=1.1,
        y0=-0.5,
        y1=0.5,
        nx=2,
        ny=2,
        quads=quads,
        material="coil",
    )
    return path.read_text(encoding="utf-8")


def _curved_triangle_text(tmp_path: Path) -> str:
    from netgen.geom2d import SplineGeometry
    from ngsolve import Mesh

    geometry = SplineGeometry()
    geometry.AddCircle((1.1, 0.0), 1.0, leftdomain=1, bc="outer")
    geometry.SetMaterial(1, "coil")
    mesh = Mesh(geometry.GenerateMesh(maxh=0.55))
    mesh.Curve(2)
    path = tmp_path / "curved_tri.vol"
    mesh.ngmesh.Save(str(path))
    return path.read_text(encoding="utf-8")


def _assembly_request(text: str, family: str, *, formulation: str) -> dict:
    boundaries = ["outer"] if family == "P2_curved" else ["bottom", "right", "top", "left"]
    return {
        "vol_text": text,
        "source_name": f"generated_{family}.vol",
        "element_family": family,
        "formulation": formulation,
        "dirichlet_boundaries": boundaries,
        "permeability_h_per_m": 1.0,
        "branches": [{"name": "winding", "material": "coil", "turns": 12.0}],
    }


def test_dimension_two_inventory_is_neutral_and_does_not_weaken_3d_policy(tmp_path):
    text = _rect_text(tmp_path, quads=True)
    contract = inspect_netgen_2d_vol(text, source_name="generated.vol")

    assert contract["schema"] == MESH_SCHEMA
    assert contract["dimension"] == 2
    assert contract["points"] == 9
    assert contract["cells"] == 4
    assert contract["triangles"] == 0
    assert contract["quadrilaterals"] == 4
    assert contract["material_names"] == ["coil"]
    assert contract["material_areas_m2"]["coil"] == pytest.approx(1.0)
    assert len(contract["sha256"]) == 64
    assert len(contract["contract_sha256"]) == 64
    assert contract["source_name"] == "generated.vol"
    assert not (tmp_path / "tracked-fixture.vol").exists()

    tri_text = _rect_text(tmp_path, quads=False)
    with pytest.raises(ValueError, match="dimension 3"):
        parse_netgen_tri_tet_vol(tri_text)


@pytest.mark.parametrize(
    ("family", "mesh_kind"),
    [
        ("P1", "tri"),
        ("Q1", "quad"),
        ("P2", "tri"),
        ("Q2", "quad"),
        ("P2_curved", "curved_tri"),
        ("Q2_curved", "quad"),
    ],
)
def test_all_six_axifem_families_assemble_from_real_vol(tmp_path, family, mesh_kind):
    if mesh_kind == "curved_tri":
        text = _curved_triangle_text(tmp_path)
    else:
        text = _rect_text(tmp_path, quads=mesh_kind == "quad")
    result = assemble_vol2d_field(
        _assembly_request(text, family, formulation="axisymmetric_henrotte")
    )
    matrix = np.asarray(result["field_matrix"])
    source = np.asarray(result["source_matrix"])

    assert result["status"] == "assembled"
    assert result["backend"] == "radia-axifem-h1henrotte"
    assert result["element_family"] == family
    assert result["generated_mesh_git_required"] is False
    assert matrix.shape == (result["free_field_dofs"],) * 2
    assert source.shape == (result["free_field_dofs"], 1)
    assert np.allclose(matrix, matrix.T)
    assert np.linalg.eigvalsh(matrix).min() > 0.0
    assert result["source_column_norms"][0] > 0.0


@pytest.mark.parametrize(("family", "quads"), [("P1", False), ("Q1", True)])
def test_planar_p1_q1_assembly_and_signed_turns(tmp_path, family, quads):
    text = _rect_text(tmp_path, quads=quads)
    request = _assembly_request(text, family, formulation="planar")
    request["branches"] = [
        {"name": "forward", "material": "coil", "turns": 8.0},
        {"name": "return", "material": "coil", "turns": -8.0},
    ]
    result = assemble_vol2d_field(request)
    source = np.asarray(result["source_matrix"])

    assert result["backend"] == "ngsolve-h1"
    assert result["branch_turns"] == [8.0, -8.0]
    assert source[:, 0] == pytest.approx(-source[:, 1])


def test_real_vol_assembly_feeds_augmented_parallel_circuit(tmp_path):
    text = _rect_text(tmp_path, quads=False)
    request = _assembly_request(text, "P1", formulation="planar")
    request["branches"] = [
        {"name": "left", "material": "coil", "turns": 6.0},
        {"name": "right", "material": "coil", "turns": -4.0},
    ]
    request["circuit"] = {
        "operation": "parallel",
        "branch_impedance_ohm": [1.0, 2.0],
        "total_current_a": 3.0,
        "frequency_hz": 50.0,
    }
    result = analyze_vol2d_circuit(request)

    assert result["status"] == "solved"
    assert result["execution_version"]["radia_mcp"]
    assert result["solution"]["status"] == "solved"
    assert result["solution"]["residual"]["maximum"] < 1.0e-9


def test_dimension_region_family_and_constraint_errors_fail_closed(tmp_path):
    tri_text = _rect_text(tmp_path, quads=False)
    request = _assembly_request(tri_text, "Q1", formulation="planar")
    with pytest.raises(ValueError, match="all-quad"):
        assemble_vol2d_field(request)

    request = _assembly_request(tri_text, "P1", formulation="planar")
    request["branches"][0]["material"] = "missing"
    with pytest.raises(ValueError, match="absent"):
        assemble_vol2d_field(request)

    request = _assembly_request(tri_text, "P1", formulation="planar")
    request["dirichlet_boundaries"] = ["missing"]
    with pytest.raises(ValueError, match="absent"):
        assemble_vol2d_field(request)

    request = _assembly_request(tri_text, "P1", formulation="planar")
    request["maximum_dense_dofs"] = 0
    with pytest.raises(ValueError, match="limited"):
        assemble_vol2d_field(request)


def test_mcp_tool_calls_real_vol_assembly_and_rejects_bad_json(tmp_path):
    text = _rect_text(tmp_path, quads=False)
    request = _assembly_request(text, "P1", formulation="planar")
    request["circuit"] = {
        "operation": "series",
        "branch_impedance_ohm": [0.5],
        "circuit_current_a": 2.0,
        "frequency_hz": 0.0,
    }

    response = json.loads(asyncio.run(motor_vol2d_circuit_analysis(json.dumps(request))))
    invalid = json.loads(asyncio.run(motor_vol2d_circuit_analysis("{")))

    assert response["status"] == "solved"
    assert response["assembly"]["mesh_contract"]["cells"] == 8
    assert invalid["status"] == "invalid_input"
