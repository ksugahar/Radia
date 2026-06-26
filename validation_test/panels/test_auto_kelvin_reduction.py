"""Cubit-headless end-to-end test for reduction-mode Auto-Kelvin.

Locks the full pipeline:

    reduced air block (1/4 quadrant: x>=0, z>=0)
      --> kelvin_reduction={"x": "ht=0", "z": "bn=0"} in launcher config
         --> add_kelvin.add_kelvin_cubit(reduction=...)
            --> sym_ht=0_x / sym_bn=0_z sidesets on air+Kelvin cut faces
               --> export netgen
                  --> .vol + .vol.json with the sym_*_* bcnames

Verifies (2026-04-25 regression guard):

  1. Launcher config `kelvin_reduction` dict is forwarded to
     add_kelvin.add_kelvin_cubit(reduction=...).
  2. Cubit accepts `sym_bn=0_<axis>` / `sym_ht=0_<axis>` sideset names
     verbatim (including the `=` character).
  3. `ExportNetgenCommand.cpp`'s Kelvin auto-detect does NOT overwrite
     user-set sideset names (reduction-mode regression added 2026-04-25).
  4. Kelvin periodic pair count is 1:1 (copy-mesh matched perfectly).
  5. NGSolve loads the resulting .vol without a NgException.
  6. The .vol.json companion lists sym_bn=0_z / sym_ht=0_x with CAD
     areas.

Runtime ~15 s on LAB.  Skipped when Cubit is absent (CI / mdx).
"""
from __future__ import annotations

import json
import os
import subprocess
from pathlib import Path

import pytest


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
# Auto-Kelvin moved to cubit-mesh-export 2026-05-05.
CUBIT_HELPERS = ROOT / "packages" / "cubit-mesh-export" / "src" / \
                "cubit_mesh_export" / "cubit_helpers"
ENTRY_SCRIPT = CUBIT_HELPERS / "auto_kelvin_entry.py"


def _find_cubit() -> Path | None:
    p = Path(r"C:\Program Files\Coreform Cubit 2025.12\bin\coreform_cubit.exe")
    return p if p.is_file() else None


CUBIT = _find_cubit()


# Create a 1/4 sphere (x>=0, z>=0) air block via boolean intersection
# with a brick (avoids stale-topology issues from webcut+delete).
MAKE_QUARTER_AIR_PY = r"""
import cubit

cubit.cmd("reset")

R = 0.05
cubit.cmd("create sphere radius " + str(R))
sphere = cubit.get_last_id("volume")

BIG = 10 * R
cubit.cmd("brick x " + str(BIG) + " y " + str(BIG) + " z " + str(BIG))
brick = cubit.get_last_id("volume")
cubit.cmd("move volume " + str(brick) + " x " + str(BIG/2) + " y 0 z " + str(BIG/2))

cubit.cmd("intersect volume " + str(sphere) + " " + str(brick) + " keep")
cubit.cmd("delete volume " + str(sphere) + " " + str(brick))
vols = list(cubit.parse_cubit_list("volume", "all"))
assert len(vols) == 1
quarter = vols[0]

cubit.cmd("volume " + str(quarter) + " scheme tetmesh")
cubit.cmd("volume " + str(quarter) + " size 0.020")
cubit.cmd("mesh volume " + str(quarter))
cubit.cmd("block 1 add volume " + str(quarter))
cubit.cmd('block 1 name "air"')
"""


def _write_probe(probe_py: Path, probe_json: Path):
    """Single-line-per-statement probe (Cubit play has quirks with
    multi-line expressions)."""
    probe_py.write_text(
        "import json, cubit\n"
        "blocks = {bid: cubit.get_exodus_entity_name('block', bid) for bid in cubit.get_block_id_list()}\n"
        "sidesets = {sid: cubit.get_exodus_entity_name('sideset', sid) for sid in cubit.get_sideset_id_list()}\n"
        "nodesets = {nid: cubit.get_exodus_entity_name('nodeset', nid) for nid in cubit.get_nodeset_id_list()}\n"
        "open(r'" + str(probe_json) + "', 'w').write("
        "json.dumps({'blocks': blocks, 'sidesets': sidesets, 'nodesets': nodesets}))\n",
        encoding="ascii",
    )


@pytest.mark.skipif(CUBIT is None,
                    reason="Coreform Cubit not installed")
@pytest.mark.slow
def test_reduction_quarter_xz_emits_sym_labels(tmp_path):
    """A launcher config with kelvin_reduction={x: ht=0, z: bn=0}
    must produce sym_ht=0_x + sym_bn=0_z sidesets, and the resulting
    .vol must load in NGSolve with those exact boundary labels.
    """
    make_py = tmp_path / "make_geom.py"
    make_py.write_text(MAKE_QUARTER_AIR_PY, encoding="ascii")

    probe_py = tmp_path / "probe.py"
    probe_json = tmp_path / "probe.json"
    _write_probe(probe_py, probe_json)

    cfg = tmp_path / "launcher.json"
    cfg.write_text(json.dumps({
        "add_kelvin":        True,
        "kelvin_air_block":  "air",
        "kelvin_block_name": "kelvin",
        "kelvin_mesh_size":  None,
        "kelvin_reduction":  {"x": "ht=0", "z": "bn=0"},
    }), encoding="utf-8")

    vol_path = tmp_path / "quarter.vol"
    volj_path = tmp_path / "quarter.vol.json"

    jou = tmp_path / "master.jou"
    jou.write_text(
        f'play "{make_py.as_posix()}"\n'
        f'play "{ENTRY_SCRIPT.as_posix()}"\n'
        f'play "{probe_py.as_posix()}"\n'
        f'export netgen "{vol_path.as_posix()}" order 1 overwrite\n',
        encoding="ascii",
    )

    env = os.environ.copy()
    env["RADIA_LAUNCHER_CONFIG"] = str(cfg)
    env["CUBIT_HELPERS_DIR"] = str(CUBIT_HELPERS)

    proc = subprocess.run(
        [str(CUBIT), "-batch", "-nographics", "-nojournal", str(jou)],
        env=env, capture_output=True, timeout=600)
    stdout = (proc.stdout or b"").decode("utf-8", errors="replace")
    stderr = (proc.stderr or b"").decode("utf-8", errors="replace")

    # Probe the Cubit-side state (confirm sidesets were created with
    # the exact sym_<bc>_<axis> names, including the '=' character).
    assert probe_json.is_file(), (
        f"Probe JSON not produced.  Stdout tail:\n"
        + "\n".join(stdout.splitlines()[-30:])
        + "\nStderr tail:\n"
        + "\n".join(stderr.splitlines()[-20:]))
    with open(probe_json, "r", encoding="utf-8") as f:
        probe = json.load(f)
    ss_names = set(probe["sidesets"].values())
    assert "sym_ht=0_x" in ss_names, (
        f"Missing sym_ht=0_x sideset.  Got: {ss_names}")
    assert "sym_bn=0_z" in ss_names, (
        f"Missing sym_bn=0_z sideset.  Got: {ss_names}")
    assert "kelvin_int" in ss_names and "kelvin_ext" in ss_names

    # Verify the .vol + companion JSON were written.
    assert vol_path.is_file(), f".vol missing.  Stdout tail:\n{stdout[-2000:]}"
    assert volj_path.is_file(), f".vol.json missing."

    with open(volj_path, "r", encoding="utf-8") as f:
        volj = json.load(f)
    vol_bnd = set(volj["boundaries"].keys())
    assert "sym_ht=0_x" in vol_bnd, (
        f".vol.json missing sym_ht=0_x.  Got: {vol_bnd}")
    assert "sym_bn=0_z" in vol_bnd, (
        f".vol.json missing sym_bn=0_z.  Got: {vol_bnd}")

    # Load with NGSolve -- this catches the "Ask for unused hash-value"
    # failure mode from the original 2026-04-25 C++ auto-detect bug
    # (kelvin_ext overwrote sym_*_* names, causing periodic-pair mismatch).
    from ngsolve import Mesh
    with TaskManager():
        m = Mesh(str(vol_path))
        mats = list(m.GetMaterials())
        bnds = set(m.GetBoundaries())
        bbbs = list(m.GetBBBoundaries())

        assert "air" in mats and "kelvin" in mats, mats
        assert "sym_ht=0_x" in bnds, bnds
        assert "sym_bn=0_z" in bnds, bnds
        assert "kelvin_int" in bnds and "kelvin_ext" in bnds
        assert "GND" in bbbs, bbbs

        # Confirm the Kelvin periodic pairs are complete (0 unmatched).
        # This is the main signal that auto-detect respected the user's
        # sym_* sidesets and did NOT mislabel flat cut faces as kelvin_ext.
        assert "Kelvin periodic" in stdout, \
            "Kelvin periodic diagnostic missing"
        # Expected format: "Kelvin periodic: N/N pairs (max_dist=X, 0 unmatched)"
        import re
        m_pair = re.search(
            r"Kelvin periodic: (\d+)/(\d+) pairs \(max_dist=([\d.e+-]+), "
            r"(\d+) unmatched\)", stdout)
        assert m_pair, f"Kelvin periodic line not found in stdout"
        n_paired, n_total, max_dist, n_unmatched = (
            int(m_pair.group(1)), int(m_pair.group(2)),
            float(m_pair.group(3)), int(m_pair.group(4)))
        assert n_paired == n_total, (
            f"Kelvin periodic mismatch: {n_paired}/{n_total} "
            f"({n_unmatched} unmatched, max_dist={max_dist:.2e}).  "
            f"Expected 1:1 pairing (copy-mesh produces matching vertex sets).")
        assert n_unmatched == 0, (
            f"{n_unmatched} Kelvin vertices unpaired")


# Create a 1/8 sphere (x>=0, y>=0, z>=0 octant) air block.
MAKE_OCTANT_AIR_PY = r"""
import cubit
cubit.cmd("reset")
R = 0.05
cubit.cmd("create sphere radius " + str(R))
sphere = cubit.get_last_id("volume")
BIG = 10 * R
cubit.cmd("brick x " + str(BIG) + " y " + str(BIG) + " z " + str(BIG))
brick = cubit.get_last_id("volume")
cubit.cmd("move volume " + str(brick) + " x " + str(BIG/2) + " y " + str(BIG/2) + " z " + str(BIG/2))
cubit.cmd("intersect volume " + str(sphere) + " " + str(brick) + " keep")
cubit.cmd("delete volume " + str(sphere) + " " + str(brick))
vols = list(cubit.parse_cubit_list("volume", "all"))
assert len(vols) == 1
oct_vol = vols[0]
cubit.cmd("volume " + str(oct_vol) + " scheme tetmesh")
cubit.cmd("volume " + str(oct_vol) + " size 0.020")
cubit.cmd("mesh volume " + str(oct_vol))
cubit.cmd("block 1 add volume " + str(oct_vol))
cubit.cmd('block 1 name "air"')
"""


@pytest.mark.skipif(CUBIT is None,
                    reason="Coreform Cubit not installed")
@pytest.mark.slow
def test_reduction_eighth_emits_sym_and_kelvin_far(tmp_path):
    """A launcher config with kelvin_reduction={x: ht=0, y: ht=0,
    z: bn=0} (ELF '-x-y+z' 1/8 convention) must produce all three
    sym_*_* sidesets PLUS a `kelvin_far` sideset for the offset-axis
    cut plane.  The .vol must load in NGSolve.  Locks the 1/8 path
    introduced 2026-04-25 (single-axis Kelvin offset + 3-plane cut +
    kelvin_far Dirichlet anchor).
    """
    make_py = tmp_path / "make_geom.py"
    make_py.write_text(MAKE_OCTANT_AIR_PY, encoding="ascii")

    probe_py = tmp_path / "probe.py"
    probe_json = tmp_path / "probe.json"
    _write_probe(probe_py, probe_json)

    cfg = tmp_path / "launcher.json"
    cfg.write_text(json.dumps({
        "add_kelvin":        True,
        "kelvin_air_block":  "air",
        "kelvin_block_name": "kelvin",
        "kelvin_mesh_size":  None,
        "kelvin_reduction":  {"x": "ht=0", "y": "ht=0", "z": "bn=0"},
    }), encoding="utf-8")

    vol_path = tmp_path / "octant.vol"

    jou = tmp_path / "master.jou"
    jou.write_text(
        f'play "{make_py.as_posix()}"\n'
        f'play "{ENTRY_SCRIPT.as_posix()}"\n'
        f'play "{probe_py.as_posix()}"\n'
        f'export netgen "{vol_path.as_posix()}" order 1 overwrite\n',
        encoding="ascii",
    )

    env = os.environ.copy()
    env["RADIA_LAUNCHER_CONFIG"] = str(cfg)
    env["CUBIT_HELPERS_DIR"] = str(CUBIT_HELPERS)

    proc = subprocess.run(
        [str(CUBIT), "-batch", "-nographics", "-nojournal", str(jou)],
        env=env, capture_output=True, timeout=600)
    stdout = (proc.stdout or b"").decode("utf-8", errors="replace")
    stderr = (proc.stderr or b"").decode("utf-8", errors="replace")

    assert probe_json.is_file(), \
        f"Probe JSON not produced.  Stdout tail:\n{stdout[-2000:]}"
    with open(probe_json, "r", encoding="utf-8") as f:
        probe = json.load(f)
    ss_names = set(probe["sidesets"].values())
    # All three sym sidesets, one per reduction axis.
    assert "sym_ht=0_x" in ss_names, ss_names
    assert "sym_ht=0_y" in ss_names, ss_names
    assert "sym_bn=0_z" in ss_names, ss_names
    # kelvin_far is the new label for the 1/8 case.
    assert "kelvin_far" in ss_names, ss_names

    assert vol_path.is_file()

    from ngsolve import Mesh
    from ngsolve import TaskManager
    m = Mesh(str(vol_path))
    bnds = set(m.GetBoundaries())
    bbbs = list(m.GetBBBoundaries())
    assert "sym_ht=0_x" in bnds and "sym_ht=0_y" in bnds, bnds
    assert "sym_bn=0_z" in bnds, bnds
    assert "kelvin_far" in bnds, bnds
    assert "kelvin_int" in bnds and "kelvin_ext" in bnds
    assert "GND" in bbbs, bbbs


@pytest.mark.skipif(CUBIT is None,
                    reason="Coreform Cubit not installed")
@pytest.mark.slow
def test_reduction_eighth_all_bn_rejected_at_runtime(tmp_path):
    """The all-bn=0 1/8 configuration is physically impossible (B
    parallel to 3 mutually perpendicular planes => B = 0).  When
    passed via the launcher config, add_kelvin_cubit should raise
    ValueError; auto_kelvin_entry catches it and prints a WARN line
    so the launcher continues without Kelvin.
    """
    make_py = tmp_path / "make_geom.py"
    make_py.write_text(MAKE_OCTANT_AIR_PY, encoding="ascii")

    cfg = tmp_path / "launcher.json"
    cfg.write_text(json.dumps({
        "add_kelvin":        True,
        "kelvin_air_block":  "air",
        "kelvin_block_name": "kelvin",
        "kelvin_mesh_size":  None,
        "kelvin_reduction":  {"x": "bn=0", "y": "bn=0", "z": "bn=0"},
    }), encoding="utf-8")

    jou = tmp_path / "master.jou"
    jou.write_text(
        f'play "{make_py.as_posix()}"\n'
        f'play "{ENTRY_SCRIPT.as_posix()}"\n',
        encoding="ascii",
    )

    env = os.environ.copy()
    env["RADIA_LAUNCHER_CONFIG"] = str(cfg)
    env["CUBIT_HELPERS_DIR"] = str(CUBIT_HELPERS)

    proc = subprocess.run(
        [str(CUBIT), "-batch", "-nographics", "-nojournal", str(jou)],
        env=env, capture_output=True, timeout=300)
    stdout = (proc.stdout or b"").decode("utf-8", errors="replace")
    # The error message text must surface in stdout (auto_kelvin_entry
    # catches the ValueError and prints "[auto_kelvin_entry] ERROR: ...").
    assert "physically impossible" in stdout, (
        f"Expected the all-bn=0 ValueError to surface as 'physically "
        f"impossible' in auto_kelvin_entry output.  Stdout tail:\n"
        f"{stdout[-2000:]}")
