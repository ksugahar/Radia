from __future__ import annotations

import copy
import json

import pytest

from radia_mcp.radia_ngsolve.material_contrast_force_gate import (
    material_contrast_force_gate,
)
from radia_mcp.radia_ngsolve.server import material_contrast_force_gate as mcp_gate


def _cases() -> list[dict]:
    common = {"force_unit": "N", "coordinate_frame": "cartesian"}
    return [
        {**common, "role": "background", "force_n": [-4.0e-7, 0.0, 0.0]},
        {**common, "role": "attractive", "force_n": [-1.5e-5, 1.0e-12, 0.0]},
        {**common, "role": "repulsive_low", "force_n": [4.3e-5, 0.0, 1.0e-12]},
        {**common, "role": "repulsive_high", "force_n": [1.25e-4, 1.0e-12, 0.0]},
    ]


def test_material_contrast_force_gate_accepts_sign_and_strength_order() -> None:
    result = material_contrast_force_gate(_cases())
    assert result["status"] == "ok"
    assert result["metrics"]["stronger_repulsion_ratio"] > 2.0


def test_material_contrast_force_gate_rejects_lost_attraction() -> None:
    cases = copy.deepcopy(_cases())
    cases[1]["force_n"][0] *= -1.0
    result = material_contrast_force_gate(cases)
    assert result["status"] == "needs_attention"
    assert result["checks"]["material_contrast_reverses_axial_force"] is False


def test_material_contrast_force_gate_rejects_large_background() -> None:
    cases = copy.deepcopy(_cases())
    cases[0]["force_n"][0] = 2.0e-5
    result = material_contrast_force_gate(cases)
    assert result["status"] == "needs_attention"
    assert result["checks"]["background_force_is_small"] is False


def test_material_contrast_force_mcp_tool_dispatches() -> None:
    result = json.loads(mcp_gate(_cases()))
    assert result["status"] == "ok"


def test_material_contrast_force_gate_rejects_nonordering_ratio() -> None:
    with pytest.raises(ValueError, match="greater than one"):
        material_contrast_force_gate(_cases(), min_stronger_repulsion_ratio=1.0)
