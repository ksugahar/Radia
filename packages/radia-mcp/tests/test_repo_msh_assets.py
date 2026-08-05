"""Standing audit: every tracked .msh in the durable lanes is sound.

Runs the structural audit (no Jacobians -- pure Python, CI-fast) over
docs/, validation_test/, and src/radia/panels/samples on every test
run, so a v2.2 regression or corrupt export is caught at commit time.

The allowlist below names the KNOWN legacy v2.2 assets.  It may only
SHRINK: a companion test fails the moment an allowlisted file becomes
v4.1 (or disappears), forcing the entry to be removed.
"""

from pathlib import Path

import pytest

from radia_mcp.gmsh.msh_inspect import audit_msh_directory, validate_msh

_REPO = Path(__file__).resolve().parents[3]
_LANES = ("docs", "validation_test", "src/radia/panels/samples")

# Known legacy v2.2 assets. Never ADD to this list -- new .msh must be
# v4.1 (repo format policy 2026-04-15).
_KNOWN_LEGACY = {
    # docs/cubit_mesh_export: intentional historical references (README).
    "docs/cubit_mesh_export/gmsh/cube_2nd_order.msh",
    "docs/cubit_mesh_export/gmsh/cube_v2.msh",
    "docs/cubit_mesh_export/gmsh/gmsh_2nd_order/brick_layer.msh",
    "docs/cubit_mesh_export/gmsh/gmsh_2nd_order/hex20.msh",
    "docs/cubit_mesh_export/gmsh/gmsh_2nd_order/wedge15.msh",
    # (peec_integration was migrated to v4.1 on 2026-08-05.)
}


def _existing_lanes():
    lanes = [_REPO / lane for lane in _LANES if (_REPO / lane).is_dir()]
    if not lanes:
        pytest.skip("repo durable lanes not present (installed-package run)")
    return lanes


def test_durable_lane_msh_assets_are_structurally_sound():
    unexpected = []
    for lane in _existing_lanes():
        audit = audit_msh_directory(lane)
        assert audit["ok"] is True
        for issue in audit["issues"]:
            rel = (lane.relative_to(_REPO).as_posix() + "/"
                   + issue["path"].replace("\\", "/"))
            if rel in _KNOWN_LEGACY:
                continue
            unexpected.append((rel, issue["failed_checks"],
                               issue["errors"][:1]))
    assert unexpected == [], (
        "unsound .msh assets outside the known-legacy allowlist "
        f"(new .msh must be valid v4.1): {unexpected}")


def test_known_legacy_allowlist_only_shrinks():
    _existing_lanes()
    stale_entries = []
    for rel in sorted(_KNOWN_LEGACY):
        path = _REPO / rel
        if not path.is_file():
            stale_entries.append((rel, "file no longer exists"))
            continue
        result = validate_msh(path)
        if result["checks"].get("format_is_v41", False):
            stale_entries.append((rel, "migrated to v4.1"))
    assert stale_entries == [], (
        "allowlist entries are stale -- remove them from _KNOWN_LEGACY: "
        f"{stale_entries}")
