"""Standing audit: every tracked post-processing MSH is valid ASCII 4.2.

Inspect durable docs/validation/sample artifacts; ignored local files never
satisfy this contract. Unsupported versions receive no exceptions.
"""

import subprocess
from pathlib import Path

import pytest

from radia_mcp.gmsh.msh_inspect import validate_msh

_REPO = Path(__file__).resolve().parents[2]
_LANES = ("docs", "validation_test", "src/radia/panels/samples")

_KNOWN_GOOD_MSH42_SAMPLE = 'validation_test/gmsh_post/fixtures/unit_cube_post.msh'


def _existing_lanes():
    lanes = [_REPO / lane for lane in _LANES if (_REPO / lane).is_dir()]
    if not lanes:
        pytest.fail("repository integration lane requires durable artifact directories")
    return lanes


def _tracked_lane_assets():
    """Return tracked .msh paths only, independent of ignored local assets."""
    try:
        result = subprocess.run(
            ["git", "-c", f"safe.directory={_REPO.as_posix()}", "ls-files", "--", "*.msh"],
            cwd=_REPO,
            check=True,
            capture_output=True,
            text=True, encoding="utf-8",
        )
    except (OSError, subprocess.CalledProcessError) as exc:
        pytest.fail(f"required tracked-file inventory unavailable: {exc}")
    lane_prefixes = tuple(f"{lane}/" for lane in _LANES)
    return {
        rel.replace("\\", "/")
        for rel in result.stdout.splitlines()
        if rel.replace("\\", "/").startswith(lane_prefixes)
    }


def test_durable_lane_msh_assets_are_structurally_sound():
    _existing_lanes()
    tracked = _tracked_lane_assets()
    assert tracked, "expected committed MSH assets in the repository durable lanes"
    assert _KNOWN_GOOD_MSH42_SAMPLE in tracked, 'required known-good post-only MSH 4.2 sample'
    unexpected = []
    for rel in sorted(tracked):
        result = validate_msh(_REPO / rel)
        if result["ok"]:
            continue
        failed = [name for name, passed in result["checks"].items()
                  if not passed]
        unexpected.append((rel, failed, result["errors"][:1]))
    assert unexpected == [], (
        "unsound post-processing .msh assets "
        f"(new .msh must be valid v4.2): {unexpected}")
