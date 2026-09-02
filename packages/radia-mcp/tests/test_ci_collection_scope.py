"""Exercise pytest's impact-selection hook without optional dependencies."""

from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys


PACKAGE_ROOT = Path(__file__).resolve().parents[1]


def test_node_selection_keeps_only_requested_contract():
    selected = (
        "tests/test_meta_health.py::test_meta_catalog_has_at_least_30_servers"
    )
    env = os.environ.copy()
    env["RADIA_MCP_CI_SELECTION_JSON"] = json.dumps([selected])
    result = subprocess.run(
        [
            sys.executable,
            "-m",
            "pytest",
            "tests/test_meta_health.py",
            "--collect-only",
            "-q",
            "-p",
            "no:cacheprovider",
        ],
        cwd=PACKAGE_ROOT,
        env=env,
        capture_output=True,
        text=True,
        encoding="utf-8",
        errors="replace",
        timeout=60,
    )
    output = result.stdout + result.stderr
    assert result.returncode == 0, output
    assert "test_meta_catalog_has_at_least_30_servers" in output
    assert "test_meta_health_all_subpackages_import" not in output
