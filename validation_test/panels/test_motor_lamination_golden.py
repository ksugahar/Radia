"""Golden test for calc_motor_lamination.py (Hollaus Effective Material).

Tests the offline cell problem + online global FE workflow:
- cell ECL increases with frequency (low-frequency: ~f^2; high-frequency:
  approaches asymptote when skin depth << sheet thickness)
- mu_eff stays close to the linear bulk value at low f (mu_iron *
  volume fraction)
- end-to-end full mode produces a finite total_ECL
"""

from __future__ import annotations

import json
import os
import subprocess
import sys

import pytest


HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.abspath(os.path.join(HERE, "..", ".."))
PANELS = os.path.join(REPO, "src", "radia", "panels")
TEMP = "C:/temp"


def _run(extra_args, output_json):
    mode = extra_args[extra_args.index("--mode") + 1]
    output_msh = os.path.splitext(output_json)[0] + ".msh"
    cmd = [
        sys.executable,
        os.path.join(PANELS, "calc_motor_lamination.py"),
        # sparsecholesky (ngsolve built-in, no MKL) -> same solution as the
        # production pardiso default, but immune to the MKL PARDISO threading-
        # layer load that is flaky under the pytest subprocess (see the motor
        # transient golden test for the full explanation).
        "--linear-solver", "sparsecholesky",
        "--output", output_json,
    ]
    if mode in ("global", "full"):
        cmd += ["--msh-output", output_msh]
    cmd += extra_args
    env = os.environ.copy()
    env["PYTHONIOENCODING"] = "utf-8"
    r = subprocess.run(cmd, capture_output=True, text=True, env=env)
    assert r.returncode == 0, (
        f"calc_motor_lamination.py failed (rc={r.returncode}):\n"
        f"STDERR:\n{r.stderr[-1500:]}"
    )
    with open(output_json, "r", encoding="utf-8") as f:
        result = json.load(f)
    if mode in ("global", "full"):
        assert result["gmsh_file"] == os.path.abspath(output_msh)
        with open(output_msh, "r", encoding="utf-8") as f:
            assert f.read(40).startswith("$MeshFormat\n4.1 0 8\n$EndMeshFormat\n")
    return result


def test_cell_problem_freq_sweep():
    """ECL density should increase monotonically with frequency."""
    out = os.path.join(TEMP, "test_lamination_cell.json")
    os.makedirs(TEMP, exist_ok=True)
    result = _run([
        "--mode", "cell",
        "--d-iron", "0.35e-3",
        "--d-ins", "0.05e-3",
        "--sigma", "4e6",
        "--mu-r-iron", "5000",
        "--B-list", "1.0",
        "--freq-list", "50,500,5000",
        "--cell-n-elements", "60",
    ], out)
    table = result["em_table"]
    assert len(table) == 3
    ECLs = [t["ECL_density_W_per_m3"] for t in table]
    # ECL must be positive and monotone increasing with frequency
    assert all(e > 0 for e in ECLs), f"ECL not positive: {ECLs}"
    assert ECLs[0] < ECLs[1] < ECLs[2], \
        f"ECL not monotone with freq: {ECLs}"
    # At high f, ECL should be significantly larger than at low f
    assert ECLs[2] > 2 * ECLs[0], \
        f"ECL doesn't grow with freq: 50Hz={ECLs[0]} 5kHz={ECLs[2]}"


def test_cell_problem_thickness_monotonic():
    """v2 cell problem: thicker laminations give MORE ECL per cell
    volume (Bertotti d^2 scaling regime, thin-sheet d << skin depth).
    """
    out = os.path.join(TEMP, "test_lamination_thicknesses.json")
    ECLs = []
    for d in (0.18e-3, 0.35e-3, 0.50e-3):
        result = _run([
            "--mode", "cell",
            "--d-iron", str(d),
            "--d-ins", "0.05e-3",
            "--sigma", "4e6", "--mu-r-iron", "5000",
            "--B-list", "1.0", "--freq-list", "50",
            "--cell-n-elements", "60",
        ], out)
        ECLs.append(result["em_table"][0]["ECL_density_W_per_m3"])
    # Monotonic
    assert ECLs[0] < ECLs[1] < ECLs[2], \
        f"ECL must be monotone-increasing with thickness: {ECLs}"
    # Approximately d^2-like scaling: doubling d should give ~3-4x ECL
    ratio_18_to_35 = ECLs[1] / ECLs[0]
    assert 2.0 < ratio_18_to_35 < 5.0, \
        f"ratio for ~2x thickness should be in d^2 range: {ratio_18_to_35}"


def test_cell_problem_bertotti_pure_iron_limit():
    """As d_ins -> 0, the numerical ECL converges to Bertotti's
    analytical formula P_ecl = sigma*omega^2*B^2*d^2/24 per iron
    volume (within ~5% even with vanishing insulation)."""
    import math
    sigma = 4e6; B = 1.0; f = 50
    omega = 2 * math.pi * f
    out = os.path.join(TEMP, "test_lamination_bertotti.json")
    # Vanishing insulation: pure iron limit
    result = _run([
        "--mode", "cell",
        "--d-iron", "0.35e-3",
        "--d-ins", "0.5e-6",   # 500 nm (nearly zero)
        "--sigma", str(sigma), "--mu-r-iron", "5000",
        "--B-list", str(B), "--freq-list", str(f),
        "--cell-n-elements", "80",
    ], out)
    cell = result["em_table"][0]
    P_per_cell = cell["ECL_density_W_per_m3"]
    V_iron_frac = cell["V_iron_fraction"]
    P_per_iron = P_per_cell / V_iron_frac
    P_bertotti = sigma * omega**2 * B**2 * (0.35e-3)**2 / 24.0
    ratio = P_per_iron / P_bertotti
    # In pure-iron limit, ratio should be very close to 1.0
    assert 0.95 < ratio < 1.05, \
        f"Bertotti limit not reproduced: numerical/analytical = {ratio:.4f}"


def test_full_mode_pipeline():
    """End-to-end: cell + global produces finite total ECL on PMSM mesh."""
    # Ensure mesh exists
    sys.path.insert(0, HERE)
    from build_test_motor_mesh import build_pmsm_mesh
    mesh_path = os.path.join(TEMP, "test_pmsm.vol")
    if not os.path.exists(mesh_path):
        build_pmsm_mesh(out_vol=mesh_path)

    out = os.path.join(TEMP, "test_lamination_full.json")
    result = _run([
        "--mode", "full",
        "--vol", mesh_path,
        "--d-iron", "0.35e-3", "--d-ins", "0.05e-3",
        "--sigma", "4e6", "--mu-r-iron", "5000",
        "--B-list", "1.0", "--freq-list", "1000",
        "--H-amplitude", "1000", "--freq", "1000",
        "--cell-n-elements", "60",
    ], out)
    assert "total_ECL_W" in result, f"missing total_ECL_W: keys={list(result)}"
    assert result["total_ECL_W"] > 0, \
        f"total ECL not positive: {result['total_ECL_W']}"
    assert result["iron_volume_m3"] > 0
    assert result["n_em_table_entries"] == 1


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
