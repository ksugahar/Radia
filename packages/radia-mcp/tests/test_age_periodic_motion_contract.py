from __future__ import annotations

import math
from types import SimpleNamespace

import pytest

from radia_mcp.radia_ngsolve.age_periodic_motion import (
    _angle_grid,
    _component_contract,
    age_sector_torque_gate,
    normalize_periodic_sector,
)
from radia_mcp.radia_ngsolve.vol2d_circuit import (
    Netgen2DBoundaryEdge,
    Netgen2DCell,
)


def _sector(**updates):
    value = {
        "slots": 12,
        "poles": 4,
        "sector_count": 4,
        "sector_angle_deg": 90.0,
        "boundary": "anti-periodic",
        "boundary_phase": -1,
    }
    value.update(updates)
    return value


def _rows():
    identity = {
        "mesh_contract_sha256": "a" * 64,
        "material_contract_sha256": "b" * 64,
        "operator_sha256": "c" * 64,
        "age_factorization_sha256": "d" * 64,
        "excitation_sha256": "e" * 64,
    }
    return [
        {
            "rotor_angle_rad": index * 0.1,
            "sector_torque_nm": torque,
            "full_machine_torque_nm": 4.0 * torque,
            **identity,
        }
        for index, torque in enumerate((0.5, 0.1, -0.4, -0.2))
    ]


def test_periodic_sector_is_derived_from_slot_pole_identity():
    result = normalize_periodic_sector(_sector())
    assert result["sectors"] == 4
    assert result["poles_per_sector"] == 1
    assert result["boundary"] == "anti-periodic"
    assert result["boundary_phase"] == -1.0
    assert len(result["periodicity_contract_sha256"]) == 64
    assert normalize_periodic_sector(result) == result

    periodic = normalize_periodic_sector(
        _sector(poles=8, boundary="periodic", boundary_phase=1)
    )
    assert periodic["poles_per_sector"] == 2
    assert periodic["boundary"] == "periodic"


@pytest.mark.parametrize(
    "updates",
    [
        {"sector_count": 3},
        {"sector_angle_deg": 45.0},
        {"boundary": "periodic"},
        {"boundary_phase": 1},
    ],
)
def test_periodic_sector_rejects_inconsistent_sign_or_geometry(updates):
    with pytest.raises(ValueError):
        normalize_periodic_sector(_sector(**updates))


def test_sector_torque_gate_requires_scaling_and_fixed_lineage():
    result = age_sector_torque_gate(
        {"periodic_sector": _sector(), "rows": _rows()}
    )
    assert result["status"] == "ok"
    assert result["identity_reused_without_remesh"] is True

    stale = _rows()
    stale[2]["operator_sha256"] = "f" * 64
    with pytest.raises(ValueError, match="reuse"):
        age_sector_torque_gate({"periodic_sector": _sector(), "rows": stale})

    wrong_scale = _rows()
    wrong_scale[1]["full_machine_torque_nm"] = 99.0
    with pytest.raises(ValueError, match="multiplier"):
        age_sector_torque_gate(
            {"periodic_sector": _sector(), "rows": wrong_scale}
        )


def test_angle_grid_rejects_repeated_endpoint_and_nonuniform_sampling():
    valid = _angle_grid([index * math.pi / 4.0 for index in range(4)], math.pi)
    assert valid["endpoint_policy"] == "exclude_repeated_period_endpoint"
    with pytest.raises(ValueError, match="exclude"):
        _angle_grid([index * math.pi / 4.0 for index in range(5)], math.pi)
    with pytest.raises(ValueError, match="uniform"):
        _angle_grid([0.0, 0.2, 0.5, 0.8], 1.0)


def test_age_topology_rejects_connected_rotor_and_stator_regions():
    cells = (
        Netgen2DCell(1, (1, 2, 3)),
        Netgen2DCell(2, (3, 4, 5)),
    )
    mesh = SimpleNamespace(
        cells=cells,
        boundary_edges=(
            Netgen2DBoundaryEdge(1, (1, 2)),
            Netgen2DBoundaryEdge(2, (4, 5)),
        ),
        boundary_names={1: "rotor_ring", 2: "stator_ring"},
        material_name=lambda number: {1: "rotor", 2: "stator"}[number],
    )
    names = {
        "rotor_ring": "rotor_ring",
        "stator_ring": "stator_ring",
        "rotor_material": "rotor",
        "stator_material": "stator",
    }
    with pytest.raises(ValueError, match="exactly two disconnected"):
        _component_contract(mesh, names)
