from __future__ import annotations

import asyncio
import json

import pytest

from radia_mcp.motor.server import motor_vol2d_dynamic_analysis
from radia_mcp.radia_ngsolve.vol2d_dynamics import normalize_vol2d_materials


def _mesh_contract() -> dict:
    return {
        "material_names": ["air", "iron"],
        "contract_sha256": "a" * 64,
    }


def _materials() -> dict:
    return {
        "air": {"permeability_h_per_m": 1.2566370614359173e-6},
        "iron": {
            "bh_curve": [
                {"b_t": 0.0, "h_a_per_m": 0.0},
                {"b_t": 1.0, "h_a_per_m": 100.0},
                {"b_t": 1.5, "h_a_per_m": 1000.0},
            ],
            "conductivity_s_per_m": 2.0e6,
        },
    }


def test_material_contract_is_complete_si_and_digest_bound() -> None:
    result = normalize_vol2d_materials(_mesh_contract(), _materials())

    assert result["schema"] == "radia.vol2d-material-contract.v1"
    assert result["units"]["conductivity"] == "S/m"
    assert result["materials"]["iron"]["kind"] == "nonlinear_bh"
    assert result["materials"]["iron"]["initial_permeability_h_per_m"] == pytest.approx(0.01)
    assert len(result["contract_sha256"]) == 64
    assert len(result["materials"]["iron"]["material_sha256"]) == 64


@pytest.mark.parametrize(
    "mutation,match",
    [
        (lambda row: row.pop("air"), "cover mesh materials exactly"),
        (
            lambda row: row["iron"]["bh_curve"].__setitem__(
                2, {"b_t": 0.8, "h_a_per_m": 1000.0}
            ),
            "B values must be strictly increasing",
        ),
        (
            lambda row: row["iron"].__setitem__("conductivity_s_per_m", -1.0),
            "must be nonnegative",
        ),
    ],
)
def test_material_contract_fails_closed(mutation, match: str) -> None:
    materials = _materials()
    mutation(materials)
    with pytest.raises(ValueError, match=match):
        normalize_vol2d_materials(_mesh_contract(), materials)


def test_dynamic_mcp_rejects_invalid_json_before_worker_launch() -> None:
    result = json.loads(asyncio.run(motor_vol2d_dynamic_analysis("{")))
    assert result["schema"] == "radia.vol2d-dynamic-analysis.v1"
    assert result["status"] == "invalid_input"
