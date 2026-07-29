from __future__ import annotations

import json
from pathlib import Path

import numpy as np
import pytest

from radia_mcp.radia_ngsolve import field_study_cli, vol2d_dynamics
from radia_mcp.radia_ngsolve.vol2d_circuit import _runtime_vol_path
from radia_mcp.radia_ngsolve.vol2d_dynamics import _solve_harmonic_matrices


def test_harmonic_matrix_kernel_closes_branch_power_and_joule_loss() -> None:
    frequency_hz = 50.0
    omega = 2.0 * np.pi * frequency_hz
    solved = _solve_harmonic_matrices(
        np.array([[4.0]]),
        np.array([[3.0]]),
        np.array([[2.0]]),
        frequency_hz=frequency_hz,
        branch_current_a=[[1.0, -0.25]],
    )

    state = complex(solved["state"][0])
    expected_loss = omega**2 * 3.0 * abs(state) ** 2
    expected_energy = 0.5 * 4.0 * abs(state) ** 2
    assert solved["eddy_loss"] == pytest.approx(expected_loss)
    assert solved["magnetic_energy"] == pytest.approx(expected_energy)
    assert solved["apparent_power"].real == pytest.approx(expected_loss)
    assert solved["power_error"] == pytest.approx(0.0, abs=1.0e-12)
    assert np.linalg.norm(solved["residual"], ord=np.inf) < 1.0e-12


@pytest.mark.parametrize(
    ("frequency", "current", "match"),
    [
        (0.0, [[1.0, 0.0]], "positive"),
        (50.0, [[0.0, 0.0]], "nonzero excitation"),
        (50.0, [[1.0, 0.0], [2.0, 0.0]], "one phasor per branch"),
    ],
)
def test_harmonic_matrix_kernel_rejects_invalid_studies(
    frequency: float, current: list[list[float]], match: str
) -> None:
    with pytest.raises(ValueError, match=match):
        _solve_harmonic_matrices(
            np.array([[4.0]]),
            np.array([[3.0]]),
            np.array([[2.0]]),
            frequency_hz=frequency,
            branch_current_a=current,
        )


def test_field_study_cli_materializes_owned_gmsh_bundle(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
) -> None:
    msh = "$MeshFormat\n4.1 0 8\n$EndMeshFormat\n"
    exports = {
        "gmsh_msh": {"content": msh},
        "gmsh_geo": {"content": 'Merge "field.msh";\n'},
        "gmsh_geo_opt": {"content": "General.Axes = 1;\n"},
        "gmsh_msh_opt": {"content": "Mesh.SurfaceEdges = 1;\n"},
    }
    monkeypatch.setattr(
        field_study_cli,
        "analyze_vol2d_scalar",
        lambda request: {
            "schema": "radia.vol2d-scalar-analysis.v1",
            "status": "solved",
            "exports": exports,
            "request": request,
        },
    )
    target = tmp_path / "field.msh"
    result = field_study_cli.run(
        {"physics": "electrostatic", "export_basename": "ignored"},
        msh_output=target,
    )

    assert result["status"] == "solved"
    assert target.read_text(encoding="utf-8").startswith("$MeshFormat\n4.1")
    assert target.with_suffix(".geo").is_file()
    assert Path(str(target.with_suffix(".geo")) + ".opt").is_file()
    assert Path(str(target) + ".opt").is_file()


def test_field_study_cli_rejects_unknown_physics(tmp_path: Path) -> None:
    with pytest.raises(ValueError, match="physics must be"):
        field_study_cli.run({"physics": "magical"}, msh_output=tmp_path / "x.msh")


def test_harmonic_dispatch_records_execution_version(
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    monkeypatch.setattr(
        vol2d_dynamics,
        "solve_vol2d_harmonic",
        lambda request: {"status": "solved", "request": dict(request)},
    )
    result = vol2d_dynamics.analyze_vol2d_dynamics({"operation": "harmonic"})

    assert result["status"] == "solved"
    assert result["executed_at_utc"]
    assert result["execution_version"]["radia_mcp"]
    assert result["execution_version"]["ngsolve"]


def test_runtime_vol_normalizes_matlab_crlf(monkeypatch: pytest.MonkeyPatch, tmp_path: Path) -> None:
    monkeypatch.setenv("RADIA_MCP_TEMP", str(tmp_path))
    path = _runtime_vol_path("materials\r\n1\r\n1 domain\r\n", "a" * 64)

    assert b"\r" not in path.read_bytes()
    assert path.read_text(encoding="utf-8") == "materials\n1\n1 domain\n"
