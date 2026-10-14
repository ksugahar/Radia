"""Heavy `.vol` force validation: conductor, passive target, and axisymmetry."""

from __future__ import annotations

import hashlib
import json
import math

from radia_mcp.radia_ngsolve.vol2d_circuit import write_structured_material_rect_vol
from radia_mcp.radia_ngsolve.vol2d_dynamics import assemble_vol2d_dynamics
from radia_mcp.radia_ngsolve.vol2d_force import (
    solve_vol2d_force,
    vol2d_force_virtual_work_gate,
)


MU0 = 4.0e-7 * math.pi
BOUNDARIES = ["bottom", "right", "top", "left"]


def _sha(value) -> str:
    return hashlib.sha256(
        json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")
    ).hexdigest()


def _bind(request: dict, target: dict) -> dict:
    operators = assemble_vol2d_dynamics(request)
    result = dict(request)
    result["force_target"] = {
        **target,
        "outer_boundary_names": BOUNDARIES,
        "mesh_contract_sha256": operators["assembly"]["mesh_contract"]["contract_sha256"],
        "material_contract_sha256": operators["material_contract"]["contract_sha256"],
        "operator_sha256": operators["operator_sha256"],
    }
    return result


def _planar_request(text: str, *, passive_target: bool) -> dict:
    materials = {
        "air": {"permeability_h_per_m": MU0},
        "source": {"permeability_h_per_m": MU0},
        "target": {
            "permeability_h_per_m": 600.0 * MU0 if passive_target else MU0
        },
    }
    branches = [{"name": "source_branch", "material": "source", "turns": 1.0}]
    currents = [1000.0]
    if not passive_target:
        branches.append({"name": "target_branch", "material": "target", "turns": 1.0})
        currents.append(1000.0)
    return {
        "vol_text": text,
        "source_name": "generated-force.vol",
        "element_family": "P1",
        "formulation": "planar",
        "dirichlet_boundaries": BOUNDARIES,
        "branches": branches,
        "materials": materials,
        "branch_current_a": currents,
        "maximum_dense_dofs": 900,
    }


def _planar_mesh(path, target_x: float) -> str:
    write_structured_material_rect_vol(
        path,
        x0=-1.5,
        x1=1.5,
        y0=-1.0,
        y1=1.0,
        nx=24,
        ny=16,
        rectangles=[
            {"name": "source", "x0": -0.75, "x1": -0.5, "y0": -0.125, "y1": 0.125},
            {
                "name": "target",
                "x0": target_x - 0.125,
                "x1": target_x + 0.125,
                "y0": -0.125,
                "y1": 0.125,
            },
        ],
    )
    return path.read_text(encoding="utf-8")


def test_planar_dual_force_and_passive_virtual_work(tmp_path) -> None:
    text = _planar_mesh(tmp_path / "dual.vol", 0.5)
    dual = solve_vol2d_force(
        _bind(
            _planar_request(text, passive_target=False),
            {
                "target_material": "target",
                "air_material": "air",
                "center_m": [0.5, 0.0],
                "inner_radius_m": 0.19,
                "outer_radius_m": 0.32,
                "method": "dual_lorentz_weighted_stress",
                "target_branch": "target_branch",
                "agreement_relative_tolerance": 0.05,
                "model_depth_m": 0.0254,
                "force_component_frame": "global_cartesian_xy",
            },
        )
    )
    assert dual["force"]["dual_method_relative_disagreement"] < 0.05
    assert dual["force"]["weighted_stress_n_per_m"][0] < 0.0

    runs = []
    positions = [0.375, 0.5, 0.625]
    for position in positions:
        text = _planar_mesh(tmp_path / f"passive-{position}.vol", position)
        runs.append(
            solve_vol2d_force(
                _bind(
                    _planar_request(text, passive_target=True),
                    {
                        "target_material": "target",
                        "air_material": "air",
                        "center_m": [position, 0.0],
                        "inner_radius_m": 0.19,
                        "outer_radius_m": 0.32,
                        "method": "weighted_stress",
                        "model_depth_m": 0.0254,
                        "force_component_frame": "global_cartesian_xy",
                    },
                )
            )
        )
    identity = _sha({"geometry": "translated passive target", "current_a": 1000.0})
    fixed_current = _sha({"source_branch_a": 1000.0})
    gate = vol2d_force_virtual_work_gate(
        {
            "rows": [
                {
                    "displacement_m": position,
                    "coenergy": run["field"]["coenergy"],
                    "coenergy_unit": "J_per_m",
                    "physics_identity_sha256": identity,
                    "fixed_current_sha256": fixed_current,
                }
                for position, run in zip(positions, runs)
            ],
            "weighted_stress_force": runs[1]["force"]["weighted_stress_n_per_m"][0],
            "relative_tolerance": 0.07,
        }
    )
    assert gate["relative_disagreement"] < 0.07
    assert runs[1]["force"]["lorentz_scope"] == "not_applicable_to_passive_magnetic_target"


def test_axisymmetric_force_is_full_revolution_axial_only(tmp_path) -> None:
    path = tmp_path / "axisymmetric.vol"
    write_structured_material_rect_vol(
        path,
        x0=0.1,
        x1=1.6,
        y0=-1.2,
        y1=1.2,
        nx=15,
        ny=24,
        quads=True,
        rectangles=[
            {"name": "source", "x0": 0.7, "x1": 0.9, "y0": -0.5, "y1": -0.3},
            {"name": "target", "x0": 0.7, "x1": 0.9, "y0": 0.3, "y1": 0.5},
        ],
    )
    text = path.read_text(encoding="utf-8")
    request = {
        "vol_text": text,
        "source_name": "generated-axisymmetric.vol",
        "element_family": "Q1",
        "formulation": "axisymmetric_henrotte",
        "dirichlet_boundaries": BOUNDARIES,
        "branches": [
            {"name": "source_branch", "material": "source", "turns": 20.0},
        ],
        "materials": {
            "air": {"permeability_h_per_m": MU0},
            "source": {"permeability_h_per_m": MU0},
            "target": {"permeability_h_per_m": 600.0 * MU0},
        },
        "branch_current_a": [10.0],
        "maximum_dense_dofs": 900,
    }
    result = solve_vol2d_force(
        _bind(
            request,
            {
                "target_material": "target",
                "air_material": "air",
                "center_m": [0.8, 0.4],
                "inner_radius_m": 0.15,
                "outer_radius_m": 0.24,
                "method": "weighted_stress",
                "force_component_frame": "meridional_rz_axial_resultant_only",
            },
        )
    )
    assert result["force"]["weighted_stress_full_revolution_axial_n"] < 0.0
    assert result["force"]["net_radial_force_n"] == 0.0
    assert result["force"]["toroidal_weight"] == "2*pi*r"
