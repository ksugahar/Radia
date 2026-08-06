"""Standing audit: every tracked .msh in the durable lanes is sound.

Runs the structural audit (no Jacobians -- pure Python, CI-fast) over
docs/, validation_test/, and src/radia/panels/samples on every test
run, so a v2.2 regression or corrupt export is caught at commit time.

The allowlist below names the KNOWN legacy v2.2 assets.  It may only
SHRINK: a companion test fails the moment an allowlisted file becomes
v4.1 (or disappears), forcing the entry to be removed.
"""

import subprocess
from pathlib import Path

import pytest

from radia_mcp.gmsh.msh_inspect import validate_msh

_REPO = Path(__file__).resolve().parents[3]
_LANES = ("docs", "validation_test", "src/radia/panels/samples")

# Known legacy v2.2 assets. Never ADD to this list -- new .msh must be
# v4.1 (repo format policy 2026-04-15).
_KNOWN_LEGACY = set()


def _existing_lanes():
    lanes = [_REPO / lane for lane in _LANES if (_REPO / lane).is_dir()]
    if not lanes:
        pytest.skip("repo durable lanes not present (installed-package run)")
    return lanes


def _tracked_lane_assets():
    """Return tracked .msh paths only, independent of ignored local assets."""
    try:
        result = subprocess.run(
            ["git", "ls-files", "--", "*.msh"],
            cwd=_REPO,
            check=True,
            capture_output=True,
            text=True,
        )
    except (OSError, subprocess.CalledProcessError):
        return None
    lane_prefixes = tuple(f"{lane}/" for lane in _LANES)
    return {
        rel.replace("\\", "/")
        for rel in result.stdout.splitlines()
        if rel.replace("\\", "/").startswith(lane_prefixes)
    }


def test_durable_lane_msh_assets_are_structurally_sound():
    _existing_lanes()
    tracked = _tracked_lane_assets()
    if tracked is None:
        pytest.skip("git tracked-file inventory unavailable")
    unexpected = []
    for rel in sorted(tracked):
        result = validate_msh(_REPO / rel)
        if result["ok"] or rel in _KNOWN_LEGACY:
            continue
        failed = [name for name, passed in result["checks"].items()
                  if not passed]
        unexpected.append((rel, failed, result["errors"][:1]))
    assert unexpected == [], (
        "unsound .msh assets outside the known-legacy allowlist "
        f"(new .msh must be valid v4.1): {unexpected}")


def test_known_legacy_allowlist_only_shrinks():
    _existing_lanes()
    tracked = _tracked_lane_assets()
    if tracked is None:
        pytest.skip("git tracked-file inventory unavailable")
    stale_entries = []
    for rel in sorted(_KNOWN_LEGACY):
        if rel not in tracked:
            stale_entries.append((rel, "file is not tracked"))
            continue
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
