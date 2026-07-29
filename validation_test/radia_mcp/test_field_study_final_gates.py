"""Production-class analytical gates for the complete 2-D Field Study surface."""

from __future__ import annotations

import math
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("gmsh")

from ngsolve import Mesh
from netgen.geom2d import SplineGeometry

from radia_mcp.radia_ngsolve.vol2d_circuit import write_structured_rect_vol
from radia_mcp.radia_ngsolve.vol2d_dynamics import solve_vol2d_harmonic
from radia_mcp.radia_ngsolve.vol2d_electrostatic import (
    solve_vol2d_electrostatic_system,
)
from radia_mcp.radia_ngsolve.vol2d_thermal import solve_vol2d_transient_heat
from radia_mcp.radia_ngsolve.vol2d_scalar import solve_vol2d_scalar


MU0 = 4.0e-7 * math.pi
EPS0 = 8.8541878128e-12


def _rectangle(tmp_path: Path, *, quads: bool = False) -> str:
    path = tmp_path / ("field_quad.vol" if quads else "field_tri.vol")
    write_structured_rect_vol(
        path,
        x0=0.0,
        x1=1.0,
        y0=0.0,
        y1=1.0,
        nx=6,
        ny=5,
        quads=quads,
        material="domain",
    )
    return path.read_text(encoding="utf-8")


def _three_conductor_mesh(tmp_path: Path) -> str:
    geometry = SplineGeometry()
    geometry.AddRectangle((0.0, 0.0), (4.0, 3.0), leftdomain=1, bc="outer")
    geometry.AddCircle((1.25, 1.5), 0.32, leftdomain=0, rightdomain=1, bc="left_conductor")
    geometry.AddCircle((2.75, 1.5), 0.32, leftdomain=0, rightdomain=1, bc="right_conductor")
    geometry.SetMaterial(1, "dielectric")
    mesh = Mesh(geometry.GenerateMesh(maxh=0.18))
    path = tmp_path / "three_conductor.vol"
    mesh.ngmesh.Save(str(path))
    return path.read_text(encoding="utf-8")


def test_uniform_transient_heat_matches_exact_energy_ramp(tmp_path: Path) -> None:
    result = solve_vol2d_transient_heat(
        {
            "physics": "transient_heat",
            "vol_text": _rectangle(tmp_path, quads=True),
            "source_name": "generated_uniform_heat.vol",
            "element_family": "Q1",
            "formulation": "planar",
            "model_depth_m": 0.25,
            "dirichlet_values": {},
            "robin_boundaries": {},
            "materials": {
                "domain": {
                    "coefficient_si": 2.0,
                    "volumetric_source_si": 8.0,
                    "volumetric_heat_capacity_j_per_m3_k": 4.0,
                }
            },
            "initial_temperature_k": 300.0,
            "time_s": [0.0, 0.1, 0.2],
            "theta": 1.0,
            "export_basename": "uniform_transient_heat",
        }
    )
    contract = result["result_contract"]

    assert contract["minimum_temperature_history_k"][-1] == pytest.approx(300.4, abs=2.0e-11)
    assert contract["maximum_temperature_history_k"][-1] == pytest.approx(300.4, abs=2.0e-11)
    assert contract["maximum_step_residual_inf"] < 1.0e-10
    assert result["exports"]["gmsh_msh"]["content"].startswith("$MeshFormat")


def test_steady_heat_conduction_convection_closes(tmp_path: Path) -> None:
    depth = 0.2
    result = solve_vol2d_scalar(
        {
            "physics": "steady_heat",
            "vol_text": _rectangle(tmp_path),
            "source_name": "generated_steady_heat.vol",
            "element_family": "P1",
            "formulation": "planar",
            "model_depth_m": depth,
            "dirichlet_values": {"left": 400.0},
            "robin_boundaries": {
                "right": {"transfer_w_per_m2_k": 10.0, "ambient_k": 300.0}
            },
            "materials": {
                "domain": {"coefficient_si": 5.0, "volumetric_source_si": 0.0}
            },
            "export_basename": "steady_heat_production_gate",
        }
    )
    observables = result["result_contract"]["observables"]
    expected = 100.0 / (1.0 / (5.0 * depth) + 1.0 / (10.0 * depth))
    assert observables["convection_outflow_w"] == pytest.approx(expected, rel=1.0e-10)
    assert observables["heat_balance_residual_w"] == pytest.approx(0.0, abs=1.0e-9)


@pytest.mark.parametrize("frequency_hz", [0.0, 1000.0])
def test_current_flow_terminal_admittance_and_power_close(
    tmp_path: Path, frequency_hz: float
) -> None:
    depth = 0.2
    material = {"coefficient_si": 5.0, "volumetric_source_si": 0.0}
    if frequency_hz:
        material["relative_permittivity"] = 3.0
    result = solve_vol2d_scalar(
        {
            "physics": "current_flow",
            "frequency_hz": frequency_hz,
            "vol_text": _rectangle(tmp_path),
            "source_name": "generated_current_flow.vol",
            "element_family": "P1",
            "formulation": "planar",
            "model_depth_m": depth,
            "dirichlet_values": {"left": 0.0, "right": 10.0},
            "terminal_pair": {
                "positive_boundary": "right",
                "negative_boundary": "left",
            },
            "materials": {"domain": material},
            "export_basename": "current_flow_production_gate",
        }
    )
    contract = result["result_contract"]
    observables = contract["observables"]
    actual = complex(*observables["admittance_s"])
    expected = 5.0 * depth
    if frequency_hz:
        expected += 1j * 2.0 * math.pi * frequency_hz * EPS0 * 3.0 * depth
    assert actual == pytest.approx(expected, rel=1.0e-10)
    assert complex(*contract["terminal"]["reaction_closure"]) == pytest.approx(
        0.0, abs=1.0e-9
    )
    if frequency_hz == 0.0:
        assert observables["conduction_power_w"] == pytest.approx(100.0, rel=1.0e-10)


def test_three_conductor_capacitance_force_and_energy_close(tmp_path: Path) -> None:
    result = solve_vol2d_electrostatic_system(
        {
            "physics": "electrostatic_system",
            "vol_text": _three_conductor_mesh(tmp_path),
            "source_name": "generated_three_conductor.vol",
            "element_family": "P2",
            "formulation": "planar",
            "model_depth_m": 0.1,
            "materials": {
                "dielectric": {
                    "coefficient_si": 2.5 * EPS0,
                    "volumetric_source_si": 0.0,
                }
            },
            "conductors": ["left_conductor", "right_conductor", "outer"],
            "applied_voltages_v": [1.0, -1.0, 0.0],
            "force_boundaries": ["left_conductor", "right_conductor", "outer"],
            "export_basename": "three_conductor_electrostatic",
        }
    )
    contract = result["result_contract"]
    matrix = np.asarray(contract["capacitance_matrix_f"])
    forces = contract["electrostatic_force_on_conductor_n"]

    assert all(contract["checks"].values())
    assert np.allclose(matrix, matrix.T, rtol=1.0e-10, atol=1.0e-20)
    assert contract["energy_relative_error"] < 1.0e-9
    assert forces["left_conductor"][0] == pytest.approx(
        -forces["right_conductor"][0], rel=2.0e-2
    )
    total_force = np.sum(np.asarray(list(forces.values())), axis=0)
    assert np.linalg.norm(total_force) < 5.0e-2 * max(
        np.linalg.norm(forces["left_conductor"]), 1.0e-30
    )


def test_parallel_plate_maxwell_force_matches_closed_form(tmp_path: Path) -> None:
    gap = 0.4
    depth = 0.2
    voltage = 20.0
    path = tmp_path / "parallel_plate.vol"
    write_structured_rect_vol(
        path,
        x0=0.0,
        x1=1.0,
        y0=0.0,
        y1=gap,
        nx=20,
        ny=8,
        material="dielectric",
    )
    result = solve_vol2d_scalar(
        {
            "physics": "electrostatic",
            "vol_text": path.read_text(encoding="utf-8"),
            "source_name": path.name,
            "element_family": "P2",
            "formulation": "planar",
            "model_depth_m": depth,
            "dirichlet_values": {"bottom": 0.0, "top": voltage},
            "materials": {
                "dielectric": {
                    "coefficient_si": EPS0,
                    "volumetric_source_si": 0.0,
                }
            },
            "force_boundaries": ["top"],
            "export_basename": "parallel_plate_force",
        }
    )
    dielectric_force = result["result_contract"]["observables"][
        "electrostatic_force_by_boundary_n"
    ]["top"]
    expected = 0.5 * EPS0 * (voltage / gap) ** 2 * depth

    assert dielectric_force[0] == pytest.approx(0.0, abs=1.0e-18)
    assert dielectric_force[1] == pytest.approx(expected, rel=2.0e-10)


def test_nonlinear_harmonic_peak_secant_converges_and_closes_power(tmp_path: Path) -> None:
    request = {
        "vol_text": _rectangle(tmp_path),
        "source_name": "generated_nonlinear_harmonic.vol",
        "element_family": "P1",
        "formulation": "planar",
        "dirichlet_boundaries": ["bottom", "right", "top", "left"],
        "branches": [{"name": "coil", "material": "domain", "turns": 12.0}],
        "materials": {
            "domain": {
                "bh_curve": [
                    {"b_t": 0.0, "h_a_per_m": 0.0},
                    {"b_t": 0.2, "h_a_per_m": 80.0},
                    {"b_t": 0.8, "h_a_per_m": 600.0},
                    {"b_t": 1.2, "h_a_per_m": 4000.0},
                    {"b_t": 1.6, "h_a_per_m": 30000.0},
                ],
                "conductivity_s_per_m": 3.0,
            }
        },
        "frequency_hz": 200.0,
        "branch_current_a": [[100.0, -25.0]],
        "relaxation": 0.2,
        "relative_tolerance": 1.0e-10,
        "maximum_iterations": 400,
        "export_basename": "nonlinear_harmonic",
    }
    result = solve_vol2d_harmonic(request)

    assert result["converged"] is True
    assert result["iterations"] >= 3
    assert result["material_model"] == "single_frequency_peak_secant_BH_fixed_point"
    assert result["harmonics_resolved"] == [1]
    assert result["hysteresis_resolved"] is False
    assert result["nonlinear_operator_relative_change_from_initial"] > 1.0e-4
    assert result["residual_inf"] < 1.0e-4
    assert abs(result["power_closure_error_w"]) < 1.0e-7 * max(1.0, result["eddy_loss_w"])
    assert result["cycle_secant_energy_proxy_j"] > 0.0
