"""Version-consistency gate: pyproject.toml `version` MUST equal the
package `__init__.py` `__version__`, for every package in the monorepo.

WHY (bug pattern `init-py-version-mismatch-vs-pyproject`): the two lines
are read by different tools.  radia-mcp 0.36.6 once silently published a
wheel whose pyproject said 0.36.6 while __init__ said 0.34.2 -- the wheel
'Verify' step at PyPI publish time was the ONLY thing that caught it, and
only AFTER a version number was burned.  This test makes the mismatch a
CI failure on every push (it is a pure file read -- no radia / ngsolve /
mcp import needed, so it runs everywhere, including the self-hosted
"Run basic tests" step and any minimal-dep environment).

Mirrors the `version` gate in tools/ci_preflight.py (run that locally
BEFORE push).
"""
import os
import re

import pytest

_HERE = os.path.dirname(os.path.abspath(__file__))
_REPO = os.path.dirname(_HERE)

# (label, pyproject path, __init__ path) -- relative to repo root.
_PACKAGES = [
    ("radia", "pyproject.toml", "src/radia/__init__.py"),
    ("radia-mcp",
     "packages/radia-mcp/pyproject.toml",
     "packages/radia-mcp/src/radia_mcp/__init__.py"),
    ("cubit-mesh-export",
     "packages/cubit-mesh-export/pyproject.toml",
     "packages/cubit-mesh-export/src/cubit_mesh_export/__init__.py"),
]

_VER_RE = re.compile(r'(?:^version|^__version__)\s*=\s*"([^"]+)"', re.M)


def _read_version(rel_path):
    path = os.path.join(_REPO, rel_path)
    if not os.path.exists(path):
        return None
    with open(path, encoding="utf-8") as f:
        m = _VER_RE.search(f.read())
    return m.group(1) if m else None


@pytest.mark.parametrize("label,pyproject,init", _PACKAGES,
                         ids=[p[0] for p in _PACKAGES])
def test_pyproject_matches_init(label, pyproject, init):
    v_pyproject = _read_version(pyproject)
    v_init = _read_version(init)
    assert v_pyproject is not None, f"{label}: no version in {pyproject}"
    assert v_init is not None, f"{label}: no __version__ in {init}"
    assert v_pyproject == v_init, (
        f"{label}: VERSION MISMATCH -- {pyproject} says {v_pyproject!r} but "
        f"{init} says {v_init!r}.  Bump BOTH in lockstep (release-quad "
        f"Phase 2).")
