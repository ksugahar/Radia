"""Regression test for the export phantom-block bug.

BUG (2026-05-30, reported by keiko on 100号機 with a 6-turn loft coil
journal at W:/31_Go-Tech/.../2026_05_21_6turn_coil_loft/):
  Running `export netgen` against a journal with K user-defined
  blocks produces:
    * the error "ERROR: No block with ID K+1 was found", and
    * a phantom block (id K+1) gets created in the Cubit session and
      survives the export.

Cause: MeshExportInterface::get_block_list returns a "default block"
handle covering elements not assigned to any user block.  Calling
iface->id_from_handle() on that handle MATERIALISES the phantom block
in Cubit's database as a side effect, with id ~ K+1.  Downstream
parse_cubit_list("volume", "in block K+1") then errors because the
phantom block has no volume membership.

Fix: src/cubit_plugin/MeshData.cpp:extract_elements -- take a snapshot
of CubitInterface::parse_cubit_list("block", "all") BEFORE the
side-effecting iface->get_block_list call, and skip block IDs not in
the snapshot.

This test exercises the user-visible symptoms:
  (a) the export must not print "No block with ID N was found"
  (b) the Cubit block count must NOT change across the export

The first assert (a) requires capturing the export's stdout.  We assert
(b) directly via parse_cubit_list -- robust regardless of how Cubit
reports messages.

Skipped automatically when Cubit is not on the LAB / 100号機 / mdx path
(running locally without Cubit installed).
"""

import os
import sys
import pytest

# --- Setup Cubit (bypassing conftest.py's auto-detect because that
#     scanning has missed the cubit path on LAB before).
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))), 'src'))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(
    os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))

_cubit_path = None
try:
    from radia.install_panels import find_cubit_bin
    _cubit_path = find_cubit_bin()
except Exception:
    try:
        from install_panels import find_cubit_bin
        _cubit_path = find_cubit_bin()
    except Exception:
        _cubit_path = None

if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)

_plugin_dir = os.path.join(os.path.dirname(_cubit_path), 'bin', 'plugins') \
    if _cubit_path else None
if _plugin_dir and os.path.isdir(_plugin_dir):
    os.environ['CUBIT_PLUGIN_DIR'] = _plugin_dir

_cubit_import_err = None
try:
    import cubit
    _cubit_available = True
except ImportError as e:
    _cubit_available = False
    _cubit_import_err = repr(e)

pytestmark = pytest.mark.skipif(
    not _cubit_available,
    reason=(f"Cubit not importable: {_cubit_import_err} "
            f"(cubit_path={_cubit_path}, plugin_dir={_plugin_dir})")
)


# Initialise Cubit once for the module.  Each test resets state explicitly.
def _ensure_cubit_inited():
    if not getattr(_ensure_cubit_inited, "_done", False):
        cubit.init(['cubit', '-nojournal', '-batch', '-nographics',
                    '-commandplugindir', _plugin_dir or ''])
        _ensure_cubit_inited._done = True


@pytest.fixture
def fresh_cubit():
    """A clean Cubit session: reset model + verify no leftover blocks."""
    _ensure_cubit_inited()
    cubit.cmd("reset")
    assert cubit.parse_cubit_list("block", "all") == [], \
        "fresh_cubit fixture: model not clean after reset"
    yield


def _snapshot_cubit_state():
    """Capture every entity collection that export should NOT modify."""
    return {
        "blocks":   tuple(cubit.parse_cubit_list("block",   "all")),
        "sidesets": tuple(cubit.parse_cubit_list("sideset", "all")),
        "nodesets": tuple(cubit.parse_cubit_list("nodeset", "all")),
        "volumes":  tuple(cubit.parse_cubit_list("volume",  "all")),
        "surfaces": tuple(cubit.parse_cubit_list("surface", "all")),
        "curves":   tuple(cubit.parse_cubit_list("curve",   "all")),
    }


def test_radia_export_netgen_does_not_create_phantom_block(fresh_cubit, tmp_path):
    """Reproduces the keiko 2026-05-30 100号機 bug.

    A journal with exactly 1 user-defined block (the workpiece sphere)
    must not leave a phantom block in the Cubit session after the export.
    """
    # Build a sphere with EXACTLY one user-defined block.  The phantom
    # block in the bug took id 2 because the user had explicitly created
    # block 1 -- the reproducer needs the same setup.
    cubit.cmd("create sphere radius 0.01")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add volume 1")
    cubit.cmd("block 1 name 'workpiece'")

    snapshot_before = _snapshot_cubit_state()
    assert snapshot_before["blocks"] == (1,), \
        f"setup: expected exactly block 1, got {snapshot_before['blocks']}"

    vol_path = str(tmp_path / "out.vol")
    # The bug message: ERROR: No block with ID 2 was found.
    cubit.cmd(f'export netgen "{vol_path}" order 1')

    snapshot_after = _snapshot_cubit_state()

    # PRIMARY ASSERTION: no Cubit state drift.
    drifted = {k: (snapshot_before[k], snapshot_after[k])
               for k in snapshot_before
               if snapshot_before[k] != snapshot_after[k]}
    assert not drifted, (
        f"export netgen drifted Cubit state -- likely a phantom-block "
        f"side effect.  Drifted entries: {drifted}"
    )

    # SECONDARY: the export must have produced a non-empty .vol.
    assert os.path.isfile(vol_path), f"export did not produce {vol_path}"
    assert os.path.getsize(vol_path) > 0, f"empty .vol at {vol_path}"


def test_radia_export_netgen_is_state_idempotent(fresh_cubit, tmp_path):
    """Calling export N times back-to-back must keep Cubit state
    pinned to the first snapshot.

    This catches the broader "any side effect during export" class of
    bugs even if the phantom-block specific symptom were patched in a
    way that leaves another residue (e.g. a leftover transient
    nodeset).  Idempotence is the contract.
    """
    cubit.cmd("create sphere radius 0.01")
    cubit.cmd("mesh volume 1")
    cubit.cmd("block 1 add volume 1")
    cubit.cmd("block 1 name 'workpiece'")

    # First export -- treat its post-state as the reference.
    cubit.cmd(f'export netgen "{tmp_path / "out_0.vol"}" order 1')
    s0 = _snapshot_cubit_state()

    # Subsequent exports must not drift.
    for i in range(1, 3):
        cubit.cmd(f'export netgen "{tmp_path / f"out_{i}.vol"}" order 1')
        si = _snapshot_cubit_state()
        drifted = {k: (s0[k], si[k]) for k in s0 if s0[k] != si[k]}
        assert not drifted, (
            f"export iteration {i} drifted from iteration 0: {drifted}"
        )
