"""Native CI dependency setup must not mutate another task's host runtime."""

from pathlib import Path
import shutil
import subprocess

import pytest
import yaml


def dependency_step():
    path = Path(__file__).resolve().parents[1] / ".github/workflows/build-test.yml"
    workflow = yaml.safe_load(path.read_text(encoding="utf-8"))
    return next(step["run"] for step in workflow["jobs"]["build-test"]["steps"]
                if step.get("name") == "Install test dependencies")


def test_dependencies_use_explicit_venv_without_host_cleanup():
    script = dependency_step()
    assert "& $ciPython -m pip install" in script
    assert "sys.prefix != sys.base_prefix" in script
    assert "os.path.samefile(sys.prefix, os.environ['RADIA_CI_VENV'])" in script
    for forbidden in ("Stop-Process", "taskkill", "Remove-Item", "site.getsitepackages"):
        assert forbidden.lower() not in script.lower()


@pytest.mark.skipif(shutil.which("pwsh") is None, reason="PowerShell runtime required")
def test_missing_venv_fails_before_any_install():
    script = "$env:RADIA_CI_VENV = ''\n" + dependency_step()
    result = subprocess.run(["pwsh", "-NoProfile", "-NonInteractive", "-Command", script],
                            capture_output=True, text=True, timeout=20)
    assert result.returncode != 0
    assert "RADIA_CI_VENV is required" in result.stderr


@pytest.mark.skipif(shutil.which("pwsh") is None, reason="PowerShell runtime required")
def test_missing_venv_python_fails_before_any_install(tmp_path):
    path = str(tmp_path).replace("'", "''")
    script = f"$env:RADIA_CI_VENV = '{path}'\n" + dependency_step()
    result = subprocess.run(["pwsh", "-NoProfile", "-NonInteractive", "-Command", script],
                            capture_output=True, text=True, timeout=20)
    assert result.returncode != 0
    assert "Isolated CI Python is missing" in result.stderr
