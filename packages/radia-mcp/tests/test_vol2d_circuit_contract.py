from __future__ import annotations

import asyncio
import json

import pytest

from radia_mcp.motor.server import motor_vol2d_circuit_analysis
from radia_mcp.radia_ngsolve.netgen_vol import parse_netgen_tri_tet_vol
from radia_mcp.radia_ngsolve.vol2d_circuit import (
    assemble_vol2d_field,
    inspect_netgen_2d_vol,
)


TRIANGLE_VOL = """mesh3d
dimension
2
surfaceelements
1
1 1 0 0 3 1 2 3
volumeelements
0
edgesegmentsgi2
3
1 0 1 2
2 0 2 3
3 0 3 1
points
3
0 0 0
1 0 0
0 1 0
materials
1
1 coil
bcnames
3
1 bottom
2 diagonal
3 left
endmesh
"""


def test_dimension_two_contract_is_stable_without_solver_imports():
    first = inspect_netgen_2d_vol(TRIANGLE_VOL, source_name="generated.vol")
    second = inspect_netgen_2d_vol(TRIANGLE_VOL, source_name="replayed.vol")

    assert first["dimension"] == 2
    assert first["triangles"] == 1
    assert first["quadrilaterals"] == 0
    assert first["material_areas_m2"]["coil"] == pytest.approx(0.5)
    assert first["sha256"] == second["sha256"]
    assert first["contract_sha256"] == second["contract_sha256"]
    assert first["source_name"] != second["source_name"]


def test_two_and_three_dimensional_vol_contracts_fail_closed_across_lanes():
    with pytest.raises(ValueError, match="requires dimension 3"):
        parse_netgen_tri_tet_vol(TRIANGLE_VOL)

    dimension_three = TRIANGLE_VOL.replace("dimension\n2", "dimension\n3", 1)
    with pytest.raises(ValueError, match="requires dimension 2"):
        inspect_netgen_2d_vol(dimension_three)


def test_vol2d_mcp_rejects_invalid_json_before_worker_launch():
    result = json.loads(asyncio.run(motor_vol2d_circuit_analysis("{")))

    assert result["schema"] == "radia.vol2d-circuit-analysis.v1"
    assert result["status"] == "invalid_input"


def test_inventory_shorthand_is_rejected_before_native_assembly():
    with pytest.raises(ValueError, match="facedescriptors"):
        assemble_vol2d_field(
            {
                "vol_text": TRIANGLE_VOL,
                "element_family": "P1",
                "formulation": "planar",
                "dirichlet_boundaries": ["bottom"],
                "permeability_h_per_m": 1.2566370614359173e-6,
                "branches": [{"name": "coil", "material": "coil", "turns": 1.0}],
            }
        )
