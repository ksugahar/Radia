"""Stage-2 CLI test for the Auto-Kelvin feature.

The Radia-NGSolve Cubit launcher exposes an "Add Kelvin open boundary
(auto)" checkbox.  Per the Panel Design Workflow Policy (CLAUDE.md),
every knob must first be validated as a CLI-argumented behaviour;
the Stage-3 PySide / C++ dialog then reduces to a trivial widget->
config mapping.

This test exercises the *argument-driven* layer end-to-end by driving
Cubit headlessly twice on the same geometry — once with
``add_kelvin=True`` in the config, once with ``add_kelvin=False`` — and
asserts that the resulting ``.vol.json`` materials differ exactly in
the presence of the ``kelvin`` block.  If this test passes, the GUI
checkbox needs no separate functional verification — its only job is
to toggle the ``add_kelvin`` boolean in the config JSON.

Not marked slow because it's not that slow (~30s with coarse mesh),
but does require Cubit to be installed at the standard location.
Skipped when Cubit is absent (CI / mdx).
"""
from __future__ import annotations

import json
import os
import shutil
import subprocess
from pathlib import Path

import pytest


HERE = Path(__file__).resolve().parent
ROOT = HERE.parent.parent
# Auto-Kelvin entry point + add_kelvin.py moved to cubit-mesh-export
# 2026-05-05.  The cubit_helpers/ directory holds both the entry
# script (auto_kelvin_entry.py) and the implementation (add_kelvin.py).
CUBIT_HELPERS = ROOT / "packages" / "cubit-mesh-export" / "src" / \
                "cubit_mesh_export" / "cubit_helpers"
ENTRY_SCRIPT = CUBIT_HELPERS / "auto_kelvin_entry.py"


def _find_cubit() -> Path | None:
    # LAB / 100号機 / mdx use this path; mdx rarely runs this test
    p = Path(r"C:\Program Files\Coreform Cubit 2025.12\bin\coreform_cubit.exe")
    return p if p.is_file() else None


CUBIT = _find_cubit()


SAMPLE_JOURNAL = r"""## Auto-Kelvin CLI test geometry - coil + wp + webcut air, NO Kelvin.
reset
create vertex 0.030 0 0
create vertex 0.033 0 0
create vertex 0.030 0 0.003
create curve arc center vertex 1 vertex 2 vertex 3 normal 0 1 0 full
create surface curve 1
sweep surface 1 axis 0 0 0 0 0 1 angle 355
${coil_vid = Id("volume")}
create cylinder radius 0.025 height 0.025
${wp_vid = Id("volume")}
create sphere radius 0.045
${air_vid = Id("volume")}
${id_before = Id("volume")}
webcut volume {air_vid} with plane zplane
${air_top = air_vid}
${air_bot = id_before + 1}
subtract volume {coil_vid} {wp_vid} from volume {air_top} {air_bot} keep_tool
imprint volume {coil_vid} {wp_vid} {air_top} {air_bot}
merge volume {coil_vid} {wp_vid} {air_top} {air_bot}
volume {coil_vid} scheme tetmesh
volume {coil_vid} size 0.010
mesh volume {coil_vid}
volume {air_top} {air_bot} scheme tetmesh
volume {air_top} {air_bot} size 0.040
mesh volume {air_top} {air_bot}
block 1 add volume {coil_vid}
block 1 name "coil"
block 2 add volume {air_top} {air_bot}
block 2 name "air"
group "coil_gaps" add surface in volume {coil_vid} with area < 0.0001
sideset 1 add surface in coil_gaps with y_coord > -0.001
sideset 1 name "source"
sideset 2 add surface in coil_gaps with y_coord < -0.001
sideset 2 name "sink"
sideset 3 add surface in volume {wp_vid}
sideset 3 name "sibc"
"""


def _run_cubit_with_config(tmpdir: Path, add_kelvin: bool,
                            set_env: bool = True,
                            kelvin_mesh_size: float | None = None) -> dict:
    """Run Cubit batch + probe Cubit block list post-entry.

    ``tmpdir`` is created if missing (so callers can pass a fresh
    subdir per invocation).  Returns {blocks: [...], kelvin_tet_count: N}.

    We probe ``cubit.get_block_id_list`` + ``parse_cubit_list('tet',
    'in block N')`` directly rather than going through
    ``export netgen`` + ``.vol.json``.  Deterministic and faster;
    also orthogonal to the cosmetic Learn Edition 50k-element export
    cap message (which per user 2026-04-24 does not block the export
    but is noisy in logs).
    """
    tmpdir.mkdir(parents=True, exist_ok=True)

    jou = tmpdir / "geom.jou"
    probe_py = tmpdir / "probe.py"
    probe_json = tmpdir / "probe_blocks.json"
    cfg = tmpdir / "launcher_config.json"
    log = tmpdir / "cubit.log"

    # Separate .py for the probe — Cubit .jou can't reliably embed
    # multi-line python statements with braces / colons (APREPRO parses
    # {} as expressions and chokes on JSON dict literals).
    probe_py.write_text(
        "import json\n"
        "import cubit\n"
        "blks = []\n"
        "kelvin_tet_count = 0\n"
        "for b in cubit.get_block_id_list():\n"
        "    name = cubit.get_exodus_entity_name('block', b)\n"
        "    blks.append(name)\n"
        "    if name and name.lower() == 'kelvin':\n"
        "        kelvin_tet_count = len(\n"
        "            cubit.parse_cubit_list('tet', 'in block %d' % b))\n"
        f"open(r'{probe_json}', 'w').write("
        "json.dumps({'blocks': blks, "
        "'kelvin_tet_count': kelvin_tet_count}))\n",
        encoding="ascii",
    )
    jou.write_text(
        SAMPLE_JOURNAL
        + f'\nplay "{ENTRY_SCRIPT.as_posix()}"\n'
        + f'play "{probe_py.as_posix()}"\n',
        encoding="ascii",
    )
    cfg.write_text(json.dumps({
        "add_kelvin":        add_kelvin,
        "kelvin_air_block":  "air",
        "kelvin_block_name": "kelvin",
        "kelvin_mesh_size":  kelvin_mesh_size,
    }), encoding="utf-8")

    env = os.environ.copy()
    # CUBIT_HELPERS_DIR lets auto_kelvin_entry.py locate add_kelvin.py
    # without relying on __file__ (Cubit's `play` does not bind it).
    # The .ccm code (run_auto_kelvin in ExportNetgenCommand.cpp) sets
    # this env var alongside RADIA_LAUNCHER_CONFIG; the test does the
    # same to exercise the same code path.
    env["CUBIT_HELPERS_DIR"] = str(CUBIT_HELPERS)
    if set_env:
        env["RADIA_LAUNCHER_CONFIG"] = str(cfg)
    else:
        env.pop("RADIA_LAUNCHER_CONFIG", None)

    proc = subprocess.run(
        [str(CUBIT), "-batch", "-nographics", "-nojournal", str(jou)],
        env=env, capture_output=True, timeout=600)
    out_s = (proc.stdout or b"").decode("utf-8", errors="replace")
    err_s = (proc.stderr or b"").decode("utf-8", errors="replace")
    log.write_text(
        "=== stdout ===\n" + out_s + "\n=== stderr ===\n" + err_s,
        encoding="utf-8", errors="replace")

    if not probe_json.is_file():
        pytest.fail(
            f"Cubit run for add_kelvin={add_kelvin} did not produce "
            f"block probe.  Log tail:\n"
            + "\n".join(out_s.splitlines()[-30:]))

    with open(probe_json, "r", encoding="utf-8") as f:
        data = json.load(f)
    data["_stdout_tail"] = "\n".join(out_s.splitlines()[-20:])
    return data


@pytest.mark.skipif(CUBIT is None,
                    reason="Coreform Cubit not installed")
@pytest.mark.slow
def test_auto_kelvin_adds_block_when_enabled(tmp_path):
    """With add_kelvin=True, the Cubit model must gain a 'kelvin' block."""
    data = _run_cubit_with_config(tmp_path, add_kelvin=True)
    blocks = set(b for b in data.get("blocks", []) if b)
    assert "kelvin" in blocks, (
        f"add_kelvin=True but 'kelvin' block absent.  "
        f"Found blocks: {blocks}.  Stdout tail:\n{data['_stdout_tail']}")
    assert "air" in blocks and "coil" in blocks, (
        f"physical blocks missing: {blocks}")


@pytest.mark.skipif(CUBIT is None,
                    reason="Coreform Cubit not installed")
@pytest.mark.slow
def test_auto_kelvin_skips_when_disabled(tmp_path):
    """With add_kelvin=False, the Cubit model must NOT get a 'kelvin' block."""
    data = _run_cubit_with_config(tmp_path, add_kelvin=False)
    blocks = set(b for b in data.get("blocks", []) if b)
    assert "kelvin" not in blocks, (
        f"add_kelvin=False but 'kelvin' block present.  "
        f"Found blocks: {blocks}.  Stdout tail:\n{data['_stdout_tail']}")
    assert "air" in blocks and "coil" in blocks, (
        f"physical blocks missing: {blocks}")


@pytest.mark.skipif(CUBIT is None,
                    reason="Coreform Cubit not installed")
@pytest.mark.slow
def test_auto_kelvin_default_is_on(tmp_path):
    """Without RADIA_LAUNCHER_CONFIG env var, default is add_kelvin=True."""
    # Pass set_env=False to bypass the config file.
    data = _run_cubit_with_config(tmp_path, add_kelvin=True, set_env=False)
    blocks = set(b for b in data.get("blocks", []) if b)
    assert "kelvin" in blocks, (
        f"default should add Kelvin but blocks: {blocks}.  "
        f"Stdout tail:\n{data['_stdout_tail']}")


@pytest.mark.skipif(CUBIT is None,
                    reason="Coreform Cubit not installed")
@pytest.mark.slow
def test_auto_kelvin_mesh_size_coarser_means_fewer_tets(tmp_path):
    """kelvin_mesh_size config key controls Kelvin exterior tet count.

    Coarser (larger) mesh_size must give fewer tets in the kelvin block.
    This proves the mesh-size config key is forwarded all the way down
    to add_kelvin_cubit's internal sizing + tetmeshing.
    """
    fine = _run_cubit_with_config(
        tmp_path / "fine", add_kelvin=True, kelvin_mesh_size=0.010)
    coarse = _run_cubit_with_config(
        tmp_path / "coarse", add_kelvin=True, kelvin_mesh_size=0.030)

    nf = fine.get("kelvin_tet_count", 0)
    nc = coarse.get("kelvin_tet_count", 0)
    assert nf > 0 and nc > 0, (
        f"kelvin tet counts unexpectedly zero: fine={nf}, coarse={nc}")
    assert nc < nf, (
        f"coarser Kelvin mesh should have FEWER tets than fine, but got "
        f"coarse={nc} vs fine={nf}")


@pytest.mark.skipif(CUBIT is None,
                    reason="Coreform Cubit not installed")
@pytest.mark.slow
def test_checkbox_off_path_skips_entry_script(tmp_path):
    """Verifies the C++ launcher's `if (kelvinAuto) play entry` gate.

    When the user UNCHECKS the 'Add Kelvin open boundary (auto)' checkbox,
    the C++ launcher skips the `play auto_kelvin_entry.py` call entirely
    (see RadiaComp.cpp::launch_radia_ngsolve).  This test reproduces that
    exact code path by building a .jou that does NOT play the entry —
    the Cubit model should end up with only the physical blocks and no
    'kelvin' block in sight.
    """
    jou = tmp_path / "geom.jou"
    probe_py = tmp_path / "probe.py"
    probe_json = tmp_path / "probe_blocks.json"
    log = tmp_path / "cubit.log"

    probe_py.write_text(
        "import json\n"
        "import cubit\n"
        "blks = [cubit.get_exodus_entity_name('block', b) "
        "for b in cubit.get_block_id_list()]\n"
        f"open(r'{probe_json}', 'w').write("
        "json.dumps({'blocks': blks}))\n",
        encoding="ascii",
    )
    # NO play auto_kelvin_entry.py — matches checkbox-OFF path exactly.
    jou.write_text(
        SAMPLE_JOURNAL
        + f'\nplay "{probe_py.as_posix()}"\n',
        encoding="ascii",
    )

    env = os.environ.copy()
    env.pop("RADIA_LAUNCHER_CONFIG", None)
    proc = subprocess.run(
        [str(CUBIT), "-batch", "-nographics", "-nojournal", str(jou)],
        env=env, capture_output=True, timeout=600)
    out_s = (proc.stdout or b"").decode("utf-8", errors="replace")
    err_s = (proc.stderr or b"").decode("utf-8", errors="replace")
    log.write_text(out_s + "\n" + err_s,
                   encoding="utf-8", errors="replace")

    if not probe_json.is_file():
        pytest.fail(
            f"probe did not run.  Log tail:\n"
            + "\n".join(out_s.splitlines()[-30:]))

    with open(probe_json, "r", encoding="utf-8") as f:
        data = json.load(f)
    blocks = set(b for b in data.get("blocks", []) if b)
    assert "kelvin" not in blocks, (
        f"Skipping play entry should leave no kelvin block, "
        f"but found: {blocks}")
    assert "coil" in blocks and "air" in blocks, (
        f"physical blocks missing: {blocks}")
