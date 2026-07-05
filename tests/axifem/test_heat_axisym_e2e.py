"""End-to-end test for `calc_heat_axisym.py` (FEMM-canonical standard
H1 + 2 pi r weighting; matches FEMM 4.2 hsolv/prob1big.cpp).  Runs the
production CLI as a subprocess on the structured-grid cylinder fixture
and asserts the temperature rise is in the analytically-expected band.

Geometry: Cu-disk-style steel cylinder R = 25 mm, H = 25 mm, structured
9x9 quad grid (NR x NZ = 9 x 9 = 81 quads, 100 nodes).
Source: uniform q_surf = 1e5 W/m^2 on the outer (r = R) edge.
Time: backward-euler, dt = 1 s, t_end = 10 s.
Material: steel (rho * cp = 3.64e6 J/(m^3 K), k = 46.6 W/(m K)).

Analytical expectation:
  Heat input rate    = q * 2 pi R H = 1e5 * 2 pi * 0.025 * 0.025 ~ 392.7 W
  Heat capacity      = rho cp * pi R^2 H = 3.64e6 * pi * 0.025^2 * 0.025 ~ 178.8 J/K
  Average dT/dt      = 392.7 / 178.8 ~ 2.20 K/s
  After 10 s, T_avg  ~ 25 + 22 ~ 47 C
  T_max (outer surf) > T_avg; T_min (axis) < T_avg.
  Expected band: T_max in [55, 65] C, T_min in [30, 40] C.
"""
import json
import os
import subprocess
import sys

import pytest


REPO = os.path.normpath(
    os.path.join(os.path.dirname(__file__), "..", ".."))
FIXTURE_DIR = os.path.join(REPO, "validation_test", "panels", "fixtures")
FIXTURE = os.path.join(FIXTURE_DIR,
                       "heat_workpiece_cylinder_R25_H25_axisym.vol")
SCRIPT = os.path.join(REPO, "src", "radia", "panels",
                       "calc_heat_axisym.py")


@pytest.fixture(scope="module")
def regenerate_fixture():
    """Ensure the structured-grid .vol exists (regen if missing)."""
    if not os.path.isfile(FIXTURE):
        gen = os.path.join(FIXTURE_DIR, "generate_heat_cylinder_axisym.py")
        subprocess.run([sys.executable, gen], check=True)
    return FIXTURE


def test_calc_heat_axisym_temperature_band(regenerate_fixture, tmp_path):
    """Run calc_heat_axisym.py end-to-end and check T_max / T_min are
    in the analytically-expected band."""
    cmd = [sys.executable, SCRIPT,
           "--wp-vol", regenerate_fixture,
           "--surface-label", "outer",
           "--q-uniform", "1e5",
           "--material", "steel",
           "--t-end", "10",
           "--dt", "1",
           "--t-initial", "25",
           "--t-ext", "25",
           "--h-conv", "0",
           # sparsecholesky (ngsolve built-in, no MKL): identical direct-solver
           # solution to pardiso (temperature band unchanged), but immune to the
           # MKL PARDISO threading-layer load that FATALs under the pytest
           # subprocess on this machine.  Verified standalone: pardiso gives the
           # same T_max=59.8 C / T_min=34.6 C.  Production keeps pardiso default.
           "--linear-solver", "sparsecholesky",
           "--fes-order", "1",
           "--msh-output", str(tmp_path / "heat_axisym_T.msh")]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
    assert proc.returncode == 0, (
        f"calc_heat_axisym exited {proc.returncode}\n"
        f"STDERR:\n{proc.stderr}\nSTDOUT tail:\n{proc.stdout[-500:]}")
    # Last line of stdout is the JSON summary
    json_line = proc.stdout.strip().splitlines()[-1]
    result = json.loads(json_line)
    assert "T_max_C" in result, f"missing T_max_C in result: {result}"

    T_max = float(result["T_max_C"])
    T_min = float(result["T_min_C"])
    assert 55.0 <= T_max <= 65.0, (
        f"T_max_C = {T_max} outside expected band [55, 65] for "
        f"q=1e5 W/m^2 * 10 s on steel cylinder R=25mm")
    assert 30.0 <= T_min <= 40.0, (
        f"T_min_C = {T_min} outside expected band [30, 40] for "
        f"q=1e5 W/m^2 * 10 s on steel cylinder R=25mm")
    # Heat input + surface area check (geometry parity)
    Q_input = float(result["Q_input_J"])
    expected_Q = 1e5 * 2 * 3.14159265 * 0.025 * 0.025 * 10  # q*A*t
    assert abs(Q_input - expected_Q) / expected_Q < 0.01, (
        f"Q_input = {Q_input} J vs expected {expected_Q} J")

    # Confirm we're using the standard H1 FESpace (FEMM-canonical;
    # matches FEMM 4.2 hsolv/prob1big.cpp).
    assert "FES:H1" in proc.stderr, (
        f"calc_heat_axisym did not log H1 FESpace usage; current stderr "
        f"shows: {proc.stderr[:500]}")
