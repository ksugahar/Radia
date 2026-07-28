"""Analytical validation for the portable dimension-2 scalar-PDE artifact."""

from __future__ import annotations

import math
from pathlib import Path

import pytest

gmsh = pytest.importorskip("gmsh")
pytest.importorskip("ngsolve")

from radia_mcp.radia_ngsolve.vol2d_circuit import write_structured_rect_vol
from radia_mcp.radia_ngsolve.vol2d_scalar import analyze_vol2d_scalar


def _mesh(name: str, *, x0: float = 0.0, x1: float = 1.0, y1: float = 1.0, nx: int = 16, ny: int = 8) -> str:
    path = Path(r"C:\temp") / name
    write_structured_rect_vol(
        path,
        x0=x0,
        x1=x1,
        y0=0.0,
        y1=y1,
        nx=nx,
        ny=ny,
        material="domain",
    )
    return path.read_text(encoding="utf-8")


def _base(physics: str, coefficient: float) -> dict:
    return {
        "operation": "solve",
        "physics": physics,
        "vol_text": _mesh("radia_vol2d_scalar_validation.vol"),
        "source_name": "radia_vol2d_scalar_validation.vol",
        "element_family": "P1",
        "formulation": "planar",
        "model_depth_m": 0.2,
        "dirichlet_values": {"left": 0.0, "right": 10.0},
        "terminal_pair": {"positive_boundary": "right", "negative_boundary": "left"},
        "materials": {"domain": {"coefficient_si": coefficient, "volumetric_source_si": 0.0}},
        "export_basename": f"vol2d_scalar_{physics}",
    }


def test_planar_electrostatic_capacitance_and_replay() -> None:
    result = analyze_vol2d_scalar(_base("electrostatic", 2.0))
    observables = result["result_contract"]["observables"]
    assert observables["capacitance_f"] == pytest.approx(2.0 * 0.2, rel=2.0e-12)
    assert observables["electric_energy_j"] == pytest.approx(20.0, rel=2.0e-12)
    replay = analyze_vol2d_scalar({"operation": "replay_gate", "replay_artifact": result})
    assert replay["status"] == "accepted"
    tampered = result.copy()
    tampered["result_contract"] = dict(result["result_contract"])
    tampered["result_contract"]["factorization_sha256"] = "0" * 64
    assert analyze_vol2d_scalar({"operation": "replay_gate", "replay_artifact": tampered})["status"] == "rejected"


def test_planar_dc_and_ac_current_terminal_observables() -> None:
    dc = analyze_vol2d_scalar(_base("current_flow", 5.0))["result_contract"]["observables"]
    assert dc["resistance_ohm"] == pytest.approx(1.0, rel=2.0e-12)
    assert dc["conduction_power_w"] == pytest.approx(100.0, rel=2.0e-12)

    request = _base("current_flow", 5.0)
    request["frequency_hz"] = 1000.0
    request["materials"]["domain"]["relative_permittivity"] = 3.0
    ac = analyze_vol2d_scalar(request)["result_contract"]["observables"]
    admittance = complex(*ac["admittance_s"])
    expected = 1.0 + 1j * 2.0 * math.pi * 1000.0 * 8.8541878128e-12 * 3.0 * 0.2
    assert admittance == pytest.approx(expected, rel=3.0e-11)


def test_heat_dirichlet_robin_and_uniform_source_balance() -> None:
    request = _base("steady_heat", 5.0)
    request.update({
        "dirichlet_values": {"left": 400.0},
        "terminal_pair": None,
        "robin_boundaries": {"right": {"transfer_w_per_m2_k": 10.0, "ambient_k": 300.0}},
    })
    result = analyze_vol2d_scalar(request)["result_contract"]["observables"]
    expected_heat = 100.0 / (1.0 / (5.0 * 0.2) + 1.0 / (10.0 * 0.2))
    assert result["convection_outflow_w"] == pytest.approx(expected_heat, rel=3.0e-11)
    assert result["heat_balance_residual_w"] == pytest.approx(0.0, abs=2.0e-10)

    source_request = _base("steady_heat", 5.0)
    source_request.update({
        "dirichlet_values": {},
        "terminal_pair": None,
        "robin_boundaries": {
            name: {"transfer_w_per_m2_k": 10.0, "ambient_k": 300.0}
            for name in ("left", "right", "top", "bottom")
        },
    })
    source_request["materials"]["domain"]["volumetric_source_si"] = 100.0
    source = analyze_vol2d_scalar(source_request)["result_contract"]["observables"]
    assert source["generated_heat_w"] == pytest.approx(20.0, rel=2.0e-12)
    assert source["convection_outflow_w"] == pytest.approx(20.0, rel=2.0e-11)
    assert source["heat_balance_residual_w"] == pytest.approx(0.0, abs=2.0e-10)


@pytest.mark.parametrize(
    ("physics", "coefficient", "observable"),
    [
        ("electrostatic", 2.0, "capacitance_f"),
        ("current_flow", 5.0, "admittance_s"),
        ("steady_heat", 7.0, "field_gradient_quadratic_half"),
    ],
)
def test_axisymmetric_full_revolution_logarithmic_solution(physics: str, coefficient: float, observable: str) -> None:
    a, b, axial_length = 0.2, 1.0, 0.5
    request = _base(physics, coefficient)
    request.update({
        "vol_text": _mesh("radia_vol2d_scalar_axi_validation.vol", x0=a, x1=b, y1=axial_length, nx=48, ny=6),
        "source_name": "radia_vol2d_scalar_axi_validation.vol",
        "formulation": "axisymmetric",
        "model_depth_m": None,
        "export_basename": f"vol2d_scalar_axi_{physics}",
    })
    result = analyze_vol2d_scalar(request)["result_contract"]["observables"]
    lumped = 2.0 * math.pi * coefficient * axial_length / math.log(b / a)
    if physics == "electrostatic":
        assert result[observable] == pytest.approx(lumped, rel=8.0e-4)
    elif physics == "current_flow":
        assert complex(*result[observable]).real == pytest.approx(lumped, rel=8.0e-4)
    else:
        expected_half = 0.5 * lumped * 10.0**2
        assert result[observable] == pytest.approx(expected_half, rel=8.0e-4)


def test_gmsh_v41_export_reopens() -> None:
    result = analyze_vol2d_scalar(_base("electrostatic", 2.0))
    target = Path(r"C:\temp\vol2d_scalar_validation_reopen.msh")
    target.write_text(result["exports"]["gmsh_msh"]["content"], encoding="utf-8")
    gmsh.initialize()
    try:
        gmsh.open(str(target))
        assert len(gmsh.model.mesh.getNodes()[0]) > 0
        assert len(gmsh.view.getTags()) == 3
    finally:
        gmsh.finalize()
