"""Pure contract tests for the portable dimension-2 scalar solver."""

from __future__ import annotations

from pathlib import Path

import pytest

from radia_mcp.radia_ngsolve.vol2d_circuit import write_structured_rect_vol
from radia_mcp.radia_ngsolve.vol2d_scalar import _prepare, scalar_replay_gate


def _vol(*, x0: float = 0.0) -> str:
    path = Path(r"C:\temp\radia_vol2d_scalar_contract.vol")
    write_structured_rect_vol(
        path,
        x0=x0,
        x1=1.0,
        y0=0.0,
        y1=1.0,
        nx=2,
        ny=2,
        material="domain",
    )
    return path.read_text(encoding="utf-8")


def _request(physics: str = "electrostatic") -> dict:
    return {
        "physics": physics,
        "vol_text": _vol(),
        "source_name": "contract.vol",
        "element_family": "P1",
        "formulation": "planar",
        "model_depth_m": 0.25,
        "dirichlet_values": {"left": 0.0, "right": 1.0},
        "materials": {
            "domain": {
                "coefficient_si": [2.0, 3.0],
                "volumetric_source_si": 0.0,
            }
        },
    }


@pytest.mark.parametrize("physics", ["electrostatic", "current_flow", "steady_heat"])
def test_closed_world_physics_and_diagonal_tensor(physics: str) -> None:
    prepared, mesh, materials, boundaries = _prepare(_request(physics))
    assert prepared["physics"] == physics
    assert prepared["model_depth_m"] == pytest.approx(0.25)
    assert mesh.contract()["dimension"] == 2
    assert materials["materials"]["domain"]["coefficient_si"] == [2.0, 3.0]
    assert boundaries["dirichlet_values"]["right"] == [1.0, 0.0]


def test_axisymmetric_measure_rejects_depth_and_negative_radius() -> None:
    request = _request()
    request.update({"formulation": "axisymmetric", "model_depth_m": None})
    prepared, *_ = _prepare(request)
    assert prepared["model_depth_m"] is None
    request["model_depth_m"] = 1.0
    with pytest.raises(ValueError, match="must not specify"):
        _prepare(request)
    request["model_depth_m"] = None
    request["vol_text"] = _vol(x0=-0.1)
    with pytest.raises(ValueError, match="nonnegative radius"):
        _prepare(request)


def test_mixed_units_nonpositive_tensor_and_nullspace_fail_closed() -> None:
    request = _request()
    request["materials"]["domain"]["conductivity_ms_per_m"] = 2.0
    with pytest.raises(ValueError, match="mixed-unit"):
        _prepare(request)
    request = _request()
    request["materials"]["domain"]["coefficient_si"] = [2.0, 0.0]
    with pytest.raises(ValueError, match="positive"):
        _prepare(request)
    request = _request()
    request["dirichlet_values"] = {}
    with pytest.raises(ValueError, match="nullspace"):
        _prepare(request)


def test_robin_is_heat_only_and_disjoint_from_dirichlet() -> None:
    request = _request("current_flow")
    request["robin_boundaries"] = {
        "top": {"transfer_w_per_m2_k": 5.0, "ambient_k": 300.0}
    }
    with pytest.raises(ValueError, match="only for heat studies"):
        _prepare(request)
    request = _request("steady_heat")
    request["robin_boundaries"] = {
        "left": {"transfer_w_per_m2_k": 5.0, "ambient_k": 300.0}
    }
    with pytest.raises(ValueError, match="both Dirichlet and Robin"):
        _prepare(request)


def test_replay_rejects_stale_or_incomplete_digest_surface() -> None:
    result = scalar_replay_gate(
        {
            "result_contract": {
                "schema": "radia.vol2d-scalar-analysis.v1",
                "status": "solved",
                "request_contract": {},
                "request_contract_sha256": "stale",
            },
            "exports": {},
        }
    )
    assert result["status"] == "rejected"
    assert result["checks"]["request_contract_sha256"] is False

