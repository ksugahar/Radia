"""Keep the formerly omitted compiled module in the parity inventory."""
import json
import yaml
import subprocess
import shutil
import os
import runpy
from types import SimpleNamespace
import pytest
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]


def test_engine_timeout_retains_evidence_when_cleanup_fails(tmp_path, monkeypatch):
    runner = runpy.run_path(str(ROOT/"validation_test/ngsolve_matlab_parity/run_sparsesolv_parity.py"))
    output = tmp_path/"nested/failure.json"

    def fail_cleanup(command, **kwargs):
        assert command == ["taskkill", "/PID", "12345", "/T", "/F"]
        assert json.loads(output.read_text())["passed"] is False
        raise OSError("access denied")

    monkeypatch.setattr(subprocess, "run", fail_cleanup)
    runner["record_timeout"](SimpleNamespace(pid=12345), output)
    record = json.loads(output.read_text())
    assert record["passed"] is False
    assert record["cleanup_error"] == "access denied"


def test_sparsesolv_has_native_and_explicit_fallback_owners():
    manifest = json.loads((ROOT / "matlab/python_api_parity_manifest.json").read_text())
    entry = next(e for e in manifest["binary_extensions"] if e["python"] == "sparsesolv_ngsolve.pyd")
    assert entry["status"] == "focused-native-commands"
    assert (ROOT / entry["fallback"]).is_file()
    source = (ROOT / "src/matlab/radia_mex.cpp").read_text()
    for command in entry["native_commands"]:
        assert f'command == "{command}"' in source
    for owner in entry["matlab"].split("; "):
        assert (ROOT / owner).is_file()


def test_ams_setup_uses_ngsolve_gradient_without_taskmanager():
    source = (ROOT / "src/matlab/radia_mex.cpp").read_text()
    body = source.split("void SparseSolvAMS(")[1].split("void SparseSolvIC(")[0]
    assert "hc->CreateGradient()" in body
    assert "RegionTaskManager" not in body
    assert "a.fespace != space.fespace" in body


def test_matlab_lane_runs_engine_and_retains_json_on_mdx():
    workflow = yaml.safe_load((ROOT/".github/workflows/sparsesolv.yml").read_text())
    job = workflow["jobs"]["ams-regression"]
    assert "mdx" in job["runs-on"]
    step = next(s for s in job["steps"] if s.get("name") == "Build and verify native MATLAB parity")
    assert step["if"] == "steps.matlab-impact.outputs.required == 'true'"
    assert "-MatlabMexOnly" in step["run"]
    assert "run_sparsesolv_parity.py --output" in step["run"]
    assert "sparsesolv-matlab.json" in job["steps"][-1]["with"]["path"]


def test_missing_diff_base_selects_matlab_without_failing_step(tmp_path):
    pwsh = shutil.which("pwsh")
    if not pwsh:
        pytest.skip("PowerShell runner contract")
    workflow = yaml.safe_load((ROOT/".github/workflows/sparsesolv.yml").read_text())
    step = next(s for s in workflow["jobs"]["ams-regression"]["steps"] if s.get("id") == "matlab-impact")
    script = step["run"].replace("${{ github.event_name }}", "pull_request")
    output = tmp_path/"output"
    result = subprocess.run([pwsh, "-NoProfile", "-Command",
        "function git { $global:LASTEXITCODE=1 }; " + script],
        env={**os.environ, "GITHUB_OUTPUT": str(output)}, capture_output=True, text=True)
    assert result.returncode == 0, result.stderr
    assert output.read_text().strip() == "required=true"
