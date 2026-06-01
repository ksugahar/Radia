"""End-to-end test: Cubit peec block → STEP + .vol → FEM + PEEC filament.

Tests the full Cubit pipeline:
  1. Cubit: create torus coil + air sphere
  2. Subtract coil from air (hole), mesh air
  3. add_kelvin_cubit() → Kelvin sphere with copy mesh + periodic BC
  4. export netgen → .vol, export step → coil.step
  5. calc_fem_kelvin.solve_fem(peec_step=...) → L
  6. Verify L vs Grover analytical (±15%)

Requires: Cubit license, NGSolve, cubit_mesh_export.ccm plugin.
"""

import math
import os
import sys
import tempfile

import pytest

HERE = os.path.dirname(os.path.abspath(__file__))
REPO_ROOT = os.path.dirname(os.path.dirname(HERE))
SRC = os.path.join(REPO_ROOT, "src", "radia")
PANELS = os.path.join(SRC, "panels")
for p in (SRC, PANELS):
    if p in sys.path:
        sys.path.remove(p)
    sys.path.insert(0, p)

MU_0 = 4e-7 * math.pi

# Geometry parameters
R_COIL = 0.030
A_COIL = 0.003
R_AIR = 0.060
MAXH_AIR = 0.012
MAXH_KELVIN = 0.020


def _grover_L(R, a):
    """Grover self-inductance (DC, circular cross-section)."""
    return MU_0 * R * (math.log(8 * R / a) - 2 + 0.25)


def _setup_cubit():
    """Initialize Cubit with plugin."""
    from install_panels import find_cubit_bin
    cubit_path = find_cubit_bin()
    if cubit_path and cubit_path not in sys.path:
        sys.path.append(cubit_path)
    plugin_dir = os.path.join(os.path.dirname(cubit_path), "bin", "plugins") \
        if cubit_path else None
    if plugin_dir and os.path.isdir(plugin_dir):
        os.environ["CUBIT_PLUGIN_DIR"] = plugin_dir

    import cubit
    cubit.init(["cubit", "-nojournal", "-batch",
                "-commandplugindir", plugin_dir or ""])
    return cubit


@pytest.fixture(scope="module")
def peec_cubit_model():
    """Build coil + air + Kelvin in Cubit, export .vol + coil.step."""
    cubit = _setup_cubit()
    tmpdir = tempfile.mkdtemp(prefix="radia_test_peec_cubit_")

    cubit.cmd("reset")

    # --- 1. Air sphere (volume 1) ---
    cubit.cmd(f"create sphere radius {R_AIR}")

    # --- 2. Coil torus (volume 2) ---
    cubit.cmd(f"create torus major radius {R_COIL} minor radius {A_COIL}")

    # --- 3. Subtract coil from air (hole approach) ---
    # After subtract: volume 1 = air-with-hole, volume 2 = coil (kept)
    cubit.cmd("subtract volume 2 from volume 1 keep")
    # volume 1 is now the air with a toroidal hole

    # --- 4. Webcut air for Kelvin symmetry (equator cut) ---
    cubit.cmd("webcut volume 1 with plane zplane")
    # Now: volume 1 = air top half, volume 3 = air bottom half

    # --- 5. Imprint + merge air halves ---
    cubit.cmd("imprint volume 1 3")
    cubit.cmd("merge volume 1 3")

    # --- 6. Mesh air ---
    cubit.cmd(f"volume 1 3 size {MAXH_AIR}")
    cubit.cmd("volume 1 3 scheme tetmesh")
    cubit.cmd("mesh volume 1 3")

    # --- 7. Assign air block ---
    cubit.cmd("block 1 add volume 1 3")
    cubit.cmd('block 1 name "air"')

    # --- 8. Add Kelvin sphere via add_kelvin_cubit() ---
    from add_kelvin import add_kelvin_cubit
    kelvin_info = add_kelvin_cubit(
        R=R_AIR, air_block="air", symmetry=["z"])

    # --- 9. Export ---
    vol_path = os.path.join(tmpdir, "model.vol")
    step_path = os.path.join(tmpdir, "coil.step")

    cubit.cmd(f'export netgen "{vol_path}" order 1 overwrite')
    cubit.cmd(f'export step "{step_path}" volume 2 overwrite')

    return {
        "vol_path": vol_path,
        "step_path": step_path,
        "tmpdir": tmpdir,
        "R_coil": R_COIL,
        "a_coil": A_COIL,
        "L_analytical": _grover_L(R_COIL, A_COIL),
        "kelvin_info": kelvin_info,
    }


class TestFEMPEECCubitIntegration:
    """End-to-end Cubit → .vol + coil.step → FEM + PEEC."""

    def test_vol_exists(self, peec_cubit_model):
        assert os.path.isfile(peec_cubit_model["vol_path"])

    def test_step_exists(self, peec_cubit_model):
        assert os.path.isfile(peec_cubit_model["step_path"])

    def test_vol_has_kelvin(self, peec_cubit_model):
        """Mesh should have air + kelvin with periodic identification."""
        from ngsolve import Mesh
        mesh = Mesh(peec_cubit_model["vol_path"])
        mats = set(mesh.GetMaterials())
        assert "air" in mats, f"Missing 'air', got {mats}"
        assert "kelvin" in mats, f"Missing 'kelvin', got {mats}"
        n_ident = mesh.ngmesh.GetNrIdentifications()
        assert n_ident > 0, "No periodic identification in .vol"

    def test_solve_fem_peec_inductance(self, peec_cubit_model):
        """Solve FEM with PEEC filament source, verify L vs Grover."""
        import importlib
        # Force source-tree modules (not site-packages)
        sys.path[:] = [p for p in sys.path
                       if not (os.sep + "site-packages" + os.sep in p
                               and "panels" in p.lower())]
        sys.path.insert(0, PANELS)
        sys.path.insert(0, SRC)
        for mod_name in list(sys.modules):
            if "calc_fem_kelvin" in mod_name or "calc_common" in mod_name:
                del sys.modules[mod_name]
        import calc_fem_kelvin
        solve_fem = calc_fem_kelvin.solve_fem

        result = solve_fem(
            vol_file=peec_cubit_model["vol_path"],
            fes_order=1,
            frequency=0,
            peec_step=peec_cubit_model["step_path"],
            peec_sigma=5.8e7,
            I_total=1.0,
            solver="pardiso",
            reg=1e-6,
            nthreads=0,
        )

        assert "error" not in result, f"Solver error: {result.get('error')}"
        L = result["L"]
        assert L > 0, f"Inductance must be positive, got {L}"

        L_ana = peec_cubit_model["L_analytical"]
        err = abs(L - L_ana) / L_ana
        print(f"\nL_fem  = {L*1e9:.2f} nH")
        print(f"L_ana  = {L_ana*1e9:.2f} nH")
        print(f"Error  = {err*100:.1f}%")
        assert err < 0.15, \
            f"L={L*1e9:.1f} nH vs analytical {L_ana*1e9:.1f} nH ({err*100:.1f}%)"
