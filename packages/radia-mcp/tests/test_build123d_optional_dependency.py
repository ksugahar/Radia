"""Contracts for build123d metadata gates under the minimal MCP install."""

from __future__ import annotations

import os
from pathlib import Path
import subprocess
import sys


def test_mass_property_gate_does_not_require_build123d():
    package_root = Path(__file__).resolve().parents[1]
    env = os.environ.copy()
    env["PYTHONPATH"] = os.pathsep.join(
        [str(package_root / "src"), env.get("PYTHONPATH", "")]
    ).rstrip(os.pathsep)
    script = r"""
import json
import sys

sys.modules["build123d"] = None

from radia_mcp.build123d.modeling import tube
from radia_mcp.build123d.server import build123d_mass_property_crosscheck

reference = [{"name": "box", "volume": 1.0, "area": 6.0}]
measured = {"external_cad": [{"name": "box", "volume": 1.0, "area": 6.0}]}
result = json.loads(
    build123d_mass_property_crosscheck(json.dumps(reference), json.dumps(measured))
)
assert result["status"] == "ok", result

try:
    tube(1.0, 2.0, 3.0)
except ImportError as exc:
    assert "optional 'build123d' dependency" in str(exc)
else:
    raise AssertionError("geometry helper unexpectedly ran without build123d")
"""
    completed = subprocess.run(
        [sys.executable, "-c", script],
        cwd=package_root,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert completed.returncode == 0, completed.stdout + completed.stderr
