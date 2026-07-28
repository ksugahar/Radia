from __future__ import annotations

import math

import pytest
from netgen.geom2d import SplineGeometry
from ngsolve import Mesh

from radia_mcp.radia_ngsolve.age_periodic_motion import (
    age_sector_torque_gate,
    solve_age_periodic_motion,
)


def _vol_text(tmp_path) -> str:
    path = tmp_path / "age_two_region.vol"
    geometry = SplineGeometry()
    geometry.AddCircle((0, 0), 0.20, leftdomain=2, rightdomain=0, bc="outer")
    geometry.AddCircle((0, 0), 0.11, leftdomain=0, rightdomain=2, bc="stator_ring")
    geometry.AddCircle((0, 0), 0.10, leftdomain=1, rightdomain=0, bc="rotor_ring")
    geometry.AddCircle((0, 0), 0.05, leftdomain=0, rightdomain=1, bc="rotor_inner")
    geometry.SetMaterial(1, "rotor")
    geometry.SetMaterial(2, "stator")
    mesh = Mesh(geometry.GenerateMesh(maxh=0.02))
    mesh.ngmesh.Save(str(path))
    return path.read_text(encoding="utf-8")


def _request(vol_text: str, frequency_hz: float = 0.0) -> dict:
    return {
        "operation": "solve",
        "vol_text": vol_text,
        "source_name": "generated_age_two_region.vol",
        "airgap": {
            "inner_radius_m": 0.10,
            "outer_radius_m": 0.11,
            "rotor_ring": "rotor_ring",
            "stator_ring": "stator_ring",
            "rotor_inner": "rotor_inner",
            "outer": "outer",
            "rotor_material": "rotor",
            "stator_material": "stator",
            "harmonics": [2],
        },
        "materials": {
            "rotor": {"relative_permeability": 1.0, "conductivity_s_per_m": 0.0},
            "stator": {
                "relative_permeability": 1.0,
                "conductivity_s_per_m": 2.0e6 if frequency_hz else 0.0,
            },
        },
        "periodic_sector": {
            "slots": 12,
            "poles": 4,
            "sector_count": 4,
            "sector_angle_deg": 90.0,
            "boundary": "anti-periodic",
            "boundary_phase": -1,
        },
        "excitation": {"2": {"rotor_amplitude": 1.0, "stator_amplitude": 0.3}},
        "rotor_angles_rad": [index * math.pi / 8.0 for index in range(8)],
        "axial_length_m": 0.05,
        "frequency_hz": frequency_hz,
        "element_order": 2,
    }


def test_real_and_eddy_age_sweeps_reuse_one_mesh_operator_and_factor(tmp_path):
    vol_text = _vol_text(tmp_path)
    static = solve_age_periodic_motion(_request(vol_text))
    eddy = solve_age_periodic_motion(_request(vol_text, frequency_hz=200.0))

    for result in (static, eddy):
        summary = result["torque_summary"]
        assert result["status"] == "solved"
        assert summary["peak_to_peak_nm"] > 0.0
        assert summary["closure_relative_error"] < 1.0e-8
        assert result["mesh_reused_all_angles"] is True
        assert result["operator_reused_all_angles"] is True
        assert result["factorization_reused_all_angles"] is True
        assert len({row["operator_sha256"] for row in result["torque_rows"]}) == 1
        assert len({row["age_factorization_sha256"] for row in result["torque_rows"]}) == 1

    assert static["torque_summary"]["phase_sign_reversal_observed"] is True
    assert static["torque_summary"]["phase_sign_reversal_required"] is True
    assert eddy["torque_summary"]["phase_sign_reversal_required"] is False

    assert static["operator_sha256"] != eddy["operator_sha256"]

    multiplier = static["periodicity_contract"]["whole_machine_multiplier"]
    rows = [
        {
            **row,
            "sector_torque_nm": row["torque_nm"] / multiplier,
            "full_machine_torque_nm": row["torque_nm"],
        }
        for row in static["torque_rows"]
    ]
    gate = age_sector_torque_gate(
        {"periodic_sector": _request(vol_text)["periodic_sector"], "rows": rows}
    )
    assert gate["status"] == "ok"

    incomplete = _request(vol_text)
    del incomplete["airgap"]["rotor_material"]
    with pytest.raises(ValueError, match="four boundaries and two materials"):
        solve_age_periodic_motion(incomplete)

    wrong_radius = _request(vol_text)
    wrong_radius["airgap"]["inner_radius_m"] = 0.101
    with pytest.raises(ValueError, match="declared radius"):
        solve_age_periodic_motion(wrong_radius)

    duplicate_harmonic = _request(vol_text)
    duplicate_harmonic["airgap"]["harmonics"] = [2, 2]
    with pytest.raises(ValueError, match="unique positive orders"):
        solve_age_periodic_motion(duplicate_harmonic)
