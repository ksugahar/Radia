"""End-to-end test: Cubit peec block → STEP + .vol → FEM + PEEC filament.

Tests the full Cubit pipeline:
  1. Cubit: create air sphere (solid) + torus coil (STEP only)
  2. Mesh the air; the coil filament lies inside it
  3. add_kelvin_cubit() → Kelvin sphere with copy mesh
  4. export vol (with the Kelvin pairs) → .vol, export step → coil.step
  5. calc_fem_kelvin.solve_fem(peec_step=...) → L
  6. Verify L vs Grover analytical (±15%)

Requires: Cubit license, NGSolve, the cubit_to_vol_export.dll plugin
(RADIA_CUBIT_PLUGIN_DIR; skipped visibly when unset).
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
ORACLE_DIR = os.path.join(REPO_ROOT, "validation_test", "cubit")
for p in (SRC, PANELS, ORACLE_DIR):
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
    """Start headless Cubit with the cubit_to_vol_export plugin.

    Same load order as the production export worker
    (radia_mcp.cubit.vol_export): the Netgen runtime first, then Cubit with
    CUBIT_PLUGIN_DIR.  Skips visibly without RADIA_CUBIT_PLUGIN_DIR.
    """
    plugin_dir = os.environ.get("RADIA_CUBIT_PLUGIN_DIR")
    if not plugin_dir:
        pytest.skip("RADIA_CUBIT_PLUGIN_DIR (directory holding "
                    "cubit_to_vol_export.dll) is not set; needs licensed Cubit")
    if not os.path.isfile(os.path.join(plugin_dir, "cubit_to_vol_export.dll")):
        raise FileNotFoundError(f"cubit_to_vol_export.dll not found in {plugin_dir}")
    import ngsolve  # noqa: F401
    import netgen
    from radia.cubit_locator import find_cubit_bin

    cubit_bin = find_cubit_bin()
    if cubit_bin is None:
        raise RuntimeError("Cubit not found: set CUBIT_BIN_DIR")
    os.add_dll_directory(os.path.dirname(netgen.__file__))
    # validation_test/ is on sys.path for this package, and its cubit/
    # folder would shadow Cubit's module: put Cubit first and drop a shadow.
    if cubit_bin in sys.path:
        sys.path.remove(cubit_bin)
    sys.path.insert(0, cubit_bin)
    if not hasattr(sys.modules.get("cubit"), "init"):
        sys.modules.pop("cubit", None)
    os.add_dll_directory(cubit_bin)
    os.environ["CUBIT_PLUGIN_DIR"] = os.path.abspath(plugin_dir)
    import cubit

    cubit.init(["cubit", "-nographics", "-nojournal", "-batch", "-noinitfile"])
    return cubit


@pytest.fixture(scope="module")
def peec_cubit_model():
    """Build coil + air + Kelvin in Cubit, export .vol + coil.step."""
    cubit = _setup_cubit()
    tmpdir = tempfile.mkdtemp(prefix="radia_test_peec_cubit_")

    cubit.cmd("reset")

    # --- 1. Solid air sphere; the coil torus is only the STEP source ---
    # calc_fem_kelvin's PEEC route integrates the filament line source over
    # the air elements, so the coil must lie inside the air mesh: no hole.
    # The torus is not meshed and not in a block, so export vol skips it.
    cubit.cmd(f"create sphere radius {R_AIR}")
    air = cubit.get_last_id("volume")
    cubit.cmd(f"create torus major radius {R_COIL} minor radius {A_COIL}")
    coil = cubit.get_last_id("volume")

    # --- 2. Webcut the air at the equator (Kelvin z symmetry) ---
    before = set(cubit.parse_cubit_list("volume", "all"))
    cubit.cmd(f"webcut volume {air} with plane zplane")
    air_vols = sorted({air} | (set(cubit.parse_cubit_list("volume", "all")) - before))
    assert len(air_vols) == 2, air_vols
    air_str = " ".join(map(str, air_vols))

    # --- 3. Imprint + merge air halves, mesh the air ---
    cubit.cmd(f"imprint volume {air_str}")
    cubit.cmd(f"merge volume {air_str}")
    cubit.cmd(f"volume {air_str} size {MAXH_AIR}")
    cubit.cmd(f"volume {air_str} scheme tetmesh")
    cubit.cmd(f"mesh volume {air_str}")

    # --- 4. Assign air block ---
    cubit.cmd(f"block 1 add volume {air_str}")
    cubit.cmd('block 1 name "air"')

    # --- 5. Add Kelvin sphere via add_kelvin_cubit() ---
    from radia_mcp.cubit.kelvin import add_kelvin_cubit, export_arguments
    kelvin_info = add_kelvin_cubit(
        R=R_AIR, air_block="air", symmetry=["z"])

    # --- 6. Export ---
    vol_path = os.path.join(tmpdir, "model.vol")
    step_path = os.path.join(tmpdir, "coil.step")

    cubit.cmd(f'export vol "{vol_path}" order 1 {export_arguments(kelvin_info)} overwrite')
    cubit.cmd(f'export step "{step_path}" volume {coil} overwrite')

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
        """The exported .vol carries the complete Kelvin pairs.

        Oracle: the same .vol with its identification sections removed, paired
        after load by add_kelvin_identification (offset detected from the
        mesh, assignment by distance). The exporter's pairs must equal that
        set, and the detected offset must equal the offset passed to export.
        """
        from ngsolve import Mesh
        import kelvin_pair_oracle
        vol = peec_cubit_model["vol_path"]
        mesh = Mesh(vol)
        mats = set(mesh.GetMaterials())
        assert "air" in mats, f"Missing 'air', got {mats}"
        assert "kelvin" in mats, f"Missing 'kelvin', got {mats}"
        bnds = set(mesh.GetBoundaries())
        assert {"kelvin_int", "kelvin_ext"} <= bnds, f"Kelvin labels missing: {bnds}"
        info = kelvin_pair_oracle.compare(vol, peec_cubit_model["tmpdir"])
        assert info["groups"] == 1 and info["covers_kelvin"], info
        assert info["bare_groups_before_pairing"] == 0, info
        assert info["posthoc_unmatched"] == 0 and info["posthoc_max_dist"] < 1e-10, info
        assert info["exported_pairs"] == info["posthoc_pairs"] > 0, info
        assert info["differing_pairs"] == 0, info
        center = peec_cubit_model["kelvin_info"]["center"]
        assert max(abs(a - b) for a, b in zip(info["posthoc_offset"], center)) < 1e-9, (
            info["posthoc_offset"], center)

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
            solver="auto",
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
