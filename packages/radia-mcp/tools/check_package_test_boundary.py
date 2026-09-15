"""Exercise the refactored contracts with only the package tree available.

This is a source-package boundary probe, not wheel or full-suite acceptance.
The temporary copy has no monorepo siblings and cannot use the live editable.
"""

from __future__ import annotations

import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
from xml.etree import ElementTree


TESTS = (
    "test_axisymmetric_point_source_contract.py",
    "test_ci_selection.py",
    "test_matlab_optuna_policy_section.py",
    "test_matlab_radia_mex_contract.py",
    "test_mcp_sdk_dependency_contract.py",
    "test_md2html.py",
    "test_radiation_pressure.py",
    "test_surface_triangle_maxwell_force.py",
    "test_time_harmonic_maxwell_stress.py",
    "test_validation_lane_separation.py",
)


def main() -> int:
    source = Path(__file__).resolve().parents[1]
    scratch = Path("C:/temp") if os.name == "nt" else None
    with tempfile.TemporaryDirectory(prefix="mcp-test-boundary-", dir=scratch) as temp:
        root = Path(temp) / "standalone"
        for folder in ("src", "tools"):
            shutil.copytree(source / folder, root / folder,
                            ignore=shutil.ignore_patterns("__pycache__", "*.pyc", "*.egg-info"))
        # Copy all tests so selection/boundary checks see the real test inventory.
        shutil.copytree(source / "tests", root / "tests",
                        ignore=shutil.ignore_patterns("__pycache__", "*.pyc"))
        shutil.copy2(source / "pyproject.toml", root / "pyproject.toml")
        shutil.copy2(source / "README.md", root / "README.md")
        env = os.environ.copy()
        env.pop("RADIA_MCP_CI_SELECTION_JSON", None)
        env["PYTHONPATH"] = str(root / "src")
        env["RADIA_REPO_ROOT"] = str(root / "no-monorepo")
        env.pop("RADIA_MCP_FORCE_MINIMAL", None)
        env["PYTEST_DISABLE_PLUGIN_AUTOLOAD"] = "1"
        print(f"Testing package-only source at {root}", flush=True)
        report = Path(temp) / "results.xml"
        result = subprocess.run(
            [sys.executable, "-m", "pytest", "-v", "-ra", "--confcutdir=tests",
             f"--junitxml={report}",
             *(f"tests/{name}" for name in TESTS)], cwd=root, env=env,
        )
        if result.returncode:
            return result.returncode
        cases = ElementTree.parse(report).iter("testcase")
        executed = {
            case.attrib.get("classname", "").split(".")[-1]
            for case in cases if case.find("skipped") is None
        }
        missing = {Path(name).stem for name in TESTS} - executed
        if missing:
            print(f"Required package boundary modules did not execute: {sorted(missing)}")
            return 1
        return 0


if __name__ == "__main__":
    raise SystemExit(main())
