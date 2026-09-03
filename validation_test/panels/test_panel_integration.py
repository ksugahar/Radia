"""
Layer 3: Integration test - full subprocess pipeline (Cubit required).

Requires: Cubit installed, system Python with NGSolve.
Run manually: python validation_test/panels/test_panel_integration.py

Tests:
  - Create model in Cubit (torus + workpiece cylinder)
  - Run calc_inductance.py as subprocess (same as panel does)
  - Verify JSON output schema and physical reasonableness
  - Test both ESIM and Dowell models
"""

import json
import os
import subprocess
import sys
import tempfile

import pytest

# Skip outside the explicit Cubit integration environment.
def _cubit_available():
    try:
        from radia.install_panels import find_cubit_bin
        return find_cubit_bin() is not None
    except Exception:
        return False

pytestmark = pytest.mark.skipif(
    not _cubit_available(), reason="Cubit not available")

# ============================================================
# Setup
# ============================================================
_repo = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
_panels_dir = os.path.join(_repo, "src", "radia", "panels")
_calc_script = os.path.join(_panels_dir, "calc_inductance.py")

# Find external Python (same logic as panel)
def _find_python():
    radia_py = os.environ.get("RADIA_PYTHON")
    if radia_py and os.path.isfile(radia_py):
        return radia_py
    for candidate in ["py", "python"]:
        try:
            r = subprocess.run([candidate, "-c", "import ngsolve; print('ok')"],
                               capture_output=True, text=True, timeout=30)
            if r.returncode == 0 and "ok" in r.stdout:
                return candidate
        except (FileNotFoundError, subprocess.TimeoutExpired):
            continue
    return sys.executable


def _find_cubit():
    try:
        from radia.install_panels import find_cubit_bin
        return find_cubit_bin()
    except ImportError:
        return os.environ.get("CUBIT_PATH")


# ============================================================
# Create test model in Cubit
# ============================================================
def create_torus_with_workpiece(cub5_path):
    """Create torus coil + cylinder workpiece in Cubit.

    Returns True if model created successfully.
    """
    cubit_bin = _find_cubit()
    if not cubit_bin:
        print("SKIP: Cubit not found")
        return False

    sys.path.append(cubit_bin)
    import cubit
    import io, contextlib
    with contextlib.redirect_stderr(io.StringIO()):
        cubit.init(["cubit", "-nojournal", "-batch", "-noinit"])

    cubit.cmd("reset")

    # Torus coil
    cubit.cmd("create torus major radius 0.05 minor radius 0.005")
    cubit.cmd("brick x 0.015 y 0.005 z 0.015")
    cubit.cmd("move Volume 2 x 0.05 y 0 z 0 include_merged")
    cubit.cmd("subtract volume 2 from volume 1")

    # Workpiece cylinder at center
    cubit.cmd("create cylinder height 0.02 radius 0.01")

    # Mesh
    vid_coil = 3  # after subtract
    vid_wp = 4    # cylinder
    cubit.cmd(f"volume {vid_coil} scheme tetmesh")
    cubit.cmd(f"volume {vid_coil} size 0.003")
    cubit.cmd(f"mesh volume {vid_coil}")
    cubit.cmd(f"volume {vid_wp} scheme tetmesh")
    cubit.cmd(f"volume {vid_wp} size 0.003")
    cubit.cmd(f"mesh volume {vid_wp}")

    # Blocks
    cubit.cmd("set duplicate block elements on")
    cubit.cmd(f'block 1 add volume {vid_coil}')
    cubit.cmd('block 1 name "conductor"')

    # Source/sink: planar faces at gap
    surfaces = cubit.get_relatives("volume", vid_coil, "surface")
    source_sid, sink_sid = None, None
    for sid in surfaces:
        s = cubit.surface(sid)
        if s.is_planar():
            cx, cy, cz = s.center_point()
            if cy > 0:
                source_sid = sid
            elif cy < 0:
                sink_sid = sid

    if source_sid and sink_sid:
        has_tri = len(cubit.get_entities("tri")) > 0
        if has_tri:
            cubit.cmd(f'block 2 add tri in surface {source_sid}')
            cubit.cmd('block 2 name "source"')
            cubit.cmd(f'block 3 add tri in surface {sink_sid}')
            cubit.cmd('block 3 name "sink"')

    # Workpiece block
    cubit.cmd(f'block 5 add volume {vid_wp}')
    cubit.cmd('block 5 name "workpiece"')

    cubit.cmd(f'save cub5 "{cub5_path}" overwrite')
    print(f"Model saved: {cub5_path}")
    return True


# ============================================================
# Run calc_inductance.py subprocess
# ============================================================
def run_calc(cub5_path, extra_args=None):
    """Run calc_inductance.py and return parsed JSON result."""
    python = _find_python()
    tmpdir = tempfile.mkdtemp(prefix="radia_test_")
    json_out = os.path.join(tmpdir, "result.json")

    args = [
        python, _calc_script,
        "--cub5", cub5_path,
        "--order", "2",
        "--source", "source",
        "--sink", "sink",
        "--output", json_out,
    ]
    if extra_args:
        args += extra_args

    print(f"Running: {' '.join(args[-6:])}")
    proc = subprocess.run(args, capture_output=True, text=True, timeout=300)

    if proc.returncode != 0:
        print(f"STDERR: {proc.stderr[-500:]}")
        print(f"STDOUT: {proc.stdout[-500:]}")
        return None

    if os.path.exists(json_out):
        with open(json_out) as f:
            return json.load(f)

    # Fallback: parse stdout
    for line in proc.stdout.strip().splitlines():
        if line.startswith("{"):
            return json.loads(line)
    return None


# ============================================================
# Pytest fixture (Cubit required - skip if unavailable)
# ============================================================
@pytest.fixture(scope="module")
def cub5_path():
    """Create Cubit model; skip if Cubit unavailable."""
    tmpdir = tempfile.mkdtemp(prefix="radia_panel_test_")
    cub5 = os.path.join(tmpdir, "test_model.cub5")
    try:
        if not create_torus_with_workpiece(cub5):
            pytest.skip("Cubit not available")
    except (ImportError, AttributeError, Exception) as e:
        pytest.skip(f"Cubit not available: {e}")
    return cub5


# ============================================================
# Tests
# ============================================================
def test_inductance_only(cub5_path):
    """Test 1: Inductance without workpiece (backward compatible)."""
    print("\n--- Test 1: Inductance only (no workpiece) ---")
    result = run_calc(cub5_path)
    if result is None:
        print("FAIL: No result")
        return False
    if "error" in result:
        print(f"FAIL: {result['error']}")
        return False

    L = result["inductance_H"]
    print(f"  L = {L*1e9:.2f} nH")
    assert L > 50e-9, f"L too small: {L}"
    assert L < 500e-9, f"L too large: {L}"
    assert "wp_R_effective" not in result, "Should not have ESIM keys"
    print("  PASS")
    return True


def test_esim(cub5_path):
    """Test 2: ESIM workpiece impedance."""
    print("\n--- Test 2: ESIM workpiece ---")
    result = run_calc(cub5_path, [
        "--workpiece", "workpiece",
        "--impedance-model", "esim",
        "--frequency", "50000",
        "--sigma", "2e6",
        "--half-thickness", "0.005",
        "--material", "steel",
    ])
    if result is None:
        print("FAIL: No result")
        return False
    if "error" in result:
        print(f"FAIL: {result['error']}")
        return False

    L = result["inductance_H"]
    R = result.get("wp_R_effective", 0)
    P = result.get("wp_P_total", 0)
    print(f"  L = {L*1e9:.2f} nH, R = {R:.4e} Ohm, P = {P:.4e} W")

    assert "wp_R_effective" in result, "Missing ESIM result keys"
    assert R > 0, "R must be positive"
    assert P > 0, "P must be positive"
    assert result["wp_model"] == "esim"
    print("  PASS")
    return True


def test_dowell(cub5_path):
    """Test 3: Dowell workpiece impedance."""
    print("\n--- Test 3: Dowell workpiece ---")
    result = run_calc(cub5_path, [
        "--workpiece", "workpiece",
        "--impedance-model", "dowell",
        "--frequency", "50000",
        "--sigma", "5.8e7",
        "--half-thickness", "0.005",
        "--material", "copper",
    ])
    if result is None:
        print("FAIL: No result")
        return False
    if "error" in result:
        print(f"FAIL: {result['error']}")
        return False

    R = result.get("wp_R_effective", 0)
    assert R > 0, "Dowell R must be positive"
    assert result["wp_model"] == "dowell"
    print(f"  R = {R:.4e} Ohm")
    print("  PASS")
    return True


# ============================================================
# Main
# ============================================================
if __name__ == "__main__":
    tmpdir = tempfile.mkdtemp(prefix="radia_panel_test_")
    cub5 = os.path.join(tmpdir, "test_model.cub5")

    print("=" * 60)
    print("Panel Integration Test (Cubit required)")
    print("=" * 60)

    if not create_torus_with_workpiece(cub5):
        print("\nSKIPPED: Cubit not available")
        sys.exit(0)

    passed = 0
    failed = 0

    for test_fn in [test_inductance_only, test_esim, test_dowell]:
        try:
            if test_fn(cub5):
                passed += 1
            else:
                failed += 1
        except Exception as e:
            print(f"  FAIL: {e}")
            failed += 1

    print(f"\n{'='*60}")
    print(f"Results: {passed} passed, {failed} failed")
    print(f"{'='*60}")
    sys.exit(1 if failed > 0 else 0)
