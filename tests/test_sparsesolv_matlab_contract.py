"""Retain historical source checks and enforce independent GPU CI ownership."""
import json
import os
import runpy
import shutil
import subprocess
from pathlib import Path
from types import SimpleNamespace
import pytest
import yaml

ROOT = Path(__file__).resolve().parents[1]

def test_owned_engine_does_not_inspect_or_attach_to_caller_sessions(monkeypatch):
    runner = runpy.run_path(str(ROOT/"validation_test/ngsolve_matlab_parity/run_sparsesolv_parity.py"))
    def forbidden(*args, **kwargs):
        raise AssertionError("Caller-owned MATLAB sessions must remain untouched")
    monkeypatch.setattr(subprocess, "run", forbidden)
    started = []
    api = SimpleNamespace(find_matlab=forbidden, connect_matlab=forbidden,
                          start_matlab=lambda options: started.append(options) or "owned")
    assert runner["start_owned_engine"](api) == "owned"
    assert started == ["-nodesktop -nosplash -singleCompThread"]

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


@pytest.mark.parametrize("argv,expected", [([], 900), (["--timeout", "1200"], 1200)])
def test_parity_runner_waits_for_the_worker_with_the_configured_timeout(
        tmp_path, monkeypatch, argv, expected):
    runner = runpy.run_path(str(ROOT/"validation_test/ngsolve_matlab_parity/run_sparsesolv_parity.py"))
    calls = {}

    class Child:
        pid = 4242

        def wait(self, timeout=None):
            calls["timeout"] = timeout
            return 0

    def popen(command, **kwargs):
        calls["command"] = command
        return Child()

    monkeypatch.setattr(subprocess, "Popen", popen)
    monkeypatch.setattr("sys.argv", ["run_sparsesolv_parity.py", "--output",
                                     str(tmp_path/"r.json"), *argv])
    with pytest.raises(SystemExit) as stop:
        runner["main"]()
    assert stop.value.code == 0
    assert calls["timeout"] == expected
    assert "--worker" in calls["command"] and "--timeout" not in calls["command"]


def test_sparsesolv_is_cpu_python_and_cuda_pcg_has_a_separate_owner():
    manifest = json.loads((ROOT / "packages/radia-gpu/matlab/python_api_parity_manifest.json").read_text())
    entry = next(e for e in manifest["binary_extensions"] if e["python"] == "sparsesolv_ngsolve.pyd")
    assert entry["classification"] == "cpu-python-only"
    assert entry["matlab"] == []
    assert manifest["cpu_mex_supported"] is False
    feature = next(e for e in manifest["gpu_features"] if e["id"] == "csr-pcg")
    assert feature["cpu_fallback"] is False
    assert (ROOT / feature["source"]).is_file()


def test_ams_setup_uses_ngsolve_gradient_without_taskmanager():
    source = (ROOT / "packages/radia-optuna/native/optuna_mex.cpp").read_text()
    body = source.split("void SparseSolvAMS(")[1].split("void SparseSolvIC(")[0]
    assert "hc->CreateGradient()" in body
    assert "RegionTaskManager" not in body
    assert "a.fespace != space.fespace" in body



def test_cpu_build_rejects_retired_mex_switch():
    pwsh = shutil.which("pwsh")
    if not pwsh:
        pytest.skip("PowerShell build entry contract")
    result = subprocess.run([pwsh, "-NoProfile", "-File", str(ROOT/"Build.ps1"),
                             "-MatlabMexOnly"], capture_output=True,
                            encoding="utf-8", errors="replace", timeout=20)
    assert result.returncode != 0
    assert "MatlabMexOnly" in result.stderr


def test_gpu_ci_has_its_own_build_and_true_residual_acceptance():
    cpu = yaml.safe_load((ROOT/".github/workflows/sparsesolv.yml").read_text())
    gpu = yaml.safe_load((ROOT/".github/workflows/radia-gpu.yml").read_text())
    steps = cpu["jobs"]["ams-regression"]["steps"]
    assert all(step.get("id") != "matlab-impact" for step in steps)
    assert "MatlabMexOnly" not in str(steps)
    job = gpu["jobs"]["cuda-mex"]
    assert "100" in job["runs-on"]
    commands = "\n".join(step.get("run", "") for step in job["steps"])
    assert "-CudaMexOnly -RequireNativeProvenance" in commands
    assert "run_poisson.py" in commands
    assert "--maxh 0.24" in commands
    assert "prepare_ci_evidence.py --poisson" in commands
    uploads = [step for step in job["steps"] if step.get("uses") == "actions/upload-artifact@v4"]
    evidence = next(step for step in uploads if step["with"]["name"] == "radia-gpu-cuda-evidence")
    assert "radia-gpu-public-" in evidence["with"]["path"]
