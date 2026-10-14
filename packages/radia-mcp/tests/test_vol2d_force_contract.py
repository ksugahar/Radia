from __future__ import annotations

import asyncio
import json

import pytest

from radia_mcp.motor.server import motor_vol2d_force_analysis
from radia_mcp.radia_ngsolve.vol2d_force import (
    normalize_vol2d_force_target,
    vol2d_force_refinement_gate,
    vol2d_force_virtual_work_gate,
)
from radia_mcp.radia_ngsolve.vol2d_circuit import (
    Netgen2DBoundaryEdge,
    Netgen2DCell,
    Netgen2DVolMesh,
)


def _mesh_view() -> Netgen2DVolMesh:
    return Netgen2DVolMesh(
        points=(
            (-0.1, -0.1, 0.0),
            (0.1, -0.1, 0.0),
            (0.1, 0.1, 0.0),
            (-0.1, 0.1, 0.0),
            (-1.0, -1.0, 0.0),
            (1.0, -1.0, 0.0),
            (1.0, 1.0, 0.0),
            (-1.0, 1.0, 0.0),
        ),
        cells=(
            Netgen2DCell(2, (1, 2, 3)),
            Netgen2DCell(2, (1, 3, 4)),
            Netgen2DCell(1, (5, 6, 2)),
            Netgen2DCell(1, (5, 2, 1)),
            Netgen2DCell(1, (6, 7, 3)),
            Netgen2DCell(1, (6, 3, 2)),
            Netgen2DCell(1, (7, 8, 4)),
            Netgen2DCell(1, (7, 4, 3)),
            Netgen2DCell(1, (8, 5, 1)),
            Netgen2DCell(1, (8, 1, 4)),
        ),
        boundary_edges=(
            Netgen2DBoundaryEdge(1, (5, 6)),
            Netgen2DBoundaryEdge(1, (6, 7)),
            Netgen2DBoundaryEdge(1, (7, 8)),
            Netgen2DBoundaryEdge(1, (8, 5)),
        ),
        materials={1: "air", 2: "target"},
        boundary_names={1: "outer"},
        content_sha256="a" * 64,
        has_curved_geometry=False,
        source_name="generated.vol",
    )


def _target() -> dict:
    mesh = _mesh_view()
    return {
        "target_material": "target",
        "air_material": "air",
        "center_m": [0.0, 0.0],
        "inner_radius_m": 0.2,
        "outer_radius_m": 0.5,
        "outer_boundary_names": ["outer"],
        "method": "weighted_stress",
        "model_depth_m": 0.01,
        "force_component_frame": "global_cartesian_xy",
        "mesh_contract_sha256": mesh.contract()["contract_sha256"],
        "material_contract_sha256": "b" * 64,
        "operator_sha256": "c" * 64,
    }


def test_force_target_binds_air_band_boundaries_units_and_digests() -> None:
    result = normalize_vol2d_force_target(
        _mesh_view(),
        {"contract_sha256": "b" * 64},
        "c" * 64,
        _target(),
        formulation="planar",
        dirichlet_boundaries=["outer"],
    )
    target = result["target"]
    assert target["force_unit_basis"] == "N_per_m_out_of_plane_and_depth_integrated_N"
    assert target["band_air_cell_count"] > 0
    assert len(target["boundary_contract_sha256"]) == 64
    assert len(target["force_request_sha256"]) == 64


@pytest.mark.parametrize(
    "mutation,match",
    [
        (lambda row: row.__setitem__("operator_sha256", "d" * 64), "operator_sha256"),
        (lambda row: row.pop("model_depth_m"), "model_depth_m"),
        (lambda row: row.__setitem__("inner_radius_m", 0.1), "strictly enclose"),
        (
            lambda row: row.__setitem__("force_component_frame", "local_polar"),
            "force_component_frame",
        ),
    ],
)
def test_force_target_lineage_and_geometry_fail_closed(mutation, match: str) -> None:
    target = _target()
    mutation(target)
    with pytest.raises(ValueError, match=match):
        normalize_vol2d_force_target(
            _mesh_view(),
            {"contract_sha256": "b" * 64},
            "c" * 64,
            target,
            formulation="planar",
            dirichlet_boundaries=["outer"],
        )


def test_virtual_work_gate_accepts_centered_fixed_current_derivative() -> None:
    result = vol2d_force_virtual_work_gate(
        {
            "rows": [
                {
                    "displacement_m": -0.01,
                    "coenergy": 0.9,
                    "coenergy_unit": "J",
                    "physics_identity_sha256": "a" * 64,
                    "fixed_current_sha256": "b" * 64,
                },
                {
                    "displacement_m": 0.0,
                    "coenergy": 1.0,
                    "coenergy_unit": "J",
                    "physics_identity_sha256": "a" * 64,
                    "fixed_current_sha256": "b" * 64,
                },
                {
                    "displacement_m": 0.01,
                    "coenergy": 1.1,
                    "coenergy_unit": "J",
                    "physics_identity_sha256": "a" * 64,
                    "fixed_current_sha256": "b" * 64,
                },
            ],
            "weighted_stress_force": 10.0,
            "relative_tolerance": 1.0e-12,
        }
    )
    assert result["status"] == "ok"
    assert result["force_from_coenergy_derivative"] == pytest.approx(10.0)


def test_virtual_work_gate_rejects_stale_current_identity() -> None:
    rows = [
        {
            "displacement_m": value,
            "coenergy": 1.0 + value,
            "coenergy_unit": "J_per_m",
            "physics_identity_sha256": "a" * 64,
            "fixed_current_sha256": ("b" if value <= 0 else "c") * 64,
        }
        for value in (-0.01, 0.0, 0.01)
    ]
    with pytest.raises(ValueError, match="fixed_current_sha256"):
        vol2d_force_virtual_work_gate(
            {"rows": rows, "weighted_stress_force": 1.0}
        )


def test_refinement_gate_accepts_nonzero_convergent_force() -> None:
    result = vol2d_force_refinement_gate(
        {
            "rows": [
                {
                    "mesh_cells": 100,
                    "force_vector": [1.0, 0.0],
                    "force_unit": "N_per_m",
                    "physics_identity_sha256": "d" * 64,
                },
                {
                    "mesh_cells": 400,
                    "force_vector": [1.05, 0.0],
                    "force_unit": "N_per_m",
                    "physics_identity_sha256": "d" * 64,
                },
                {
                    "mesh_cells": 1600,
                    "force_vector": [1.06, 0.0],
                    "force_unit": "N_per_m",
                    "physics_identity_sha256": "d" * 64,
                },
            ],
            "terminal_relative_tolerance": 0.02,
        }
    )
    assert result["status"] == "ok"
    assert result["terminal_relative_change"] < 0.02


@pytest.mark.parametrize(
    "payload,match",
    [
        (
            {
                "rows": [
                    {
                        "mesh_cells": 100,
                        "force_vector": [0.0, 0.0],
                        "force_unit": "N",
                        "physics_identity_sha256": "e" * 64,
                    },
                    {
                        "mesh_cells": 200,
                        "force_vector": [0.0, 0.0],
                        "force_unit": "N",
                        "physics_identity_sha256": "e" * 64,
                    },
                ]
            },
            "near-zero",
        ),
        (
            {
                "rows": [
                    {
                        "mesh_cells": 100,
                        "force_vector": [1.0, 0.0],
                        "force_unit": "N",
                        "physics_identity_sha256": "e" * 64,
                    },
                    {
                        "mesh_cells": 200,
                        "force_vector": [-1.0, 0.0],
                        "force_unit": "N",
                        "physics_identity_sha256": "e" * 64,
                    },
                ],
                "terminal_relative_tolerance": 0.1,
            },
            "not converged",
        ),
    ],
)
def test_refinement_gate_fails_closed(payload: dict, match: str) -> None:
    with pytest.raises(ValueError, match=match):
        vol2d_force_refinement_gate(payload)


def test_force_mcp_rejects_invalid_json_before_worker_launch() -> None:
    result = json.loads(asyncio.run(motor_vol2d_force_analysis("{")))
    assert result["schema"] == "radia.vol2d-force-analysis.v1"
    assert result["status"] == "invalid_input"
