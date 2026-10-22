"""GPU selection preserves native coverage without scheduling documentation builds."""
import importlib.util
import json
from pathlib import Path

import pytest
import yaml

ROOT = Path(__file__).resolve().parents[2]
spec = importlib.util.spec_from_file_location("gpu_impact", ROOT / "tools/sparsesolv_ci_impact.py")
impact = importlib.util.module_from_spec(spec)
spec.loader.exec_module(impact)


def decide(paths, blobs=None, event="pull_request", raw=None):
    def git(*args):
        if args[0] == "merge-base":
            return "base"
        if args[0] == "diff":
            return raw if raw is not None else "".join("M\0" + p + "\0" for p in paths)
        if args[0] == "show":
            return (blobs or {})[args[1]]
        raise AssertionError(args)
    payload = {"pull_request": {"base": {"sha": "base"}}, "before": "base"}
    return impact.gpu_required_for_event(event, payload, git, impact.manifest())[0]


@pytest.mark.parametrize("path", [
    "packages/radia-simulink/README.md", "docs/gpu/tutorial.md",
    "packages/radia-simulink/src/radia_simulink/cli.py",
    "packages/radia-simulink/src/radia_simulink/ltspice/parser/cir_to_schemdraw.py",
    "packages/radia-simulink/tests/ltspice/test_cli.py",
    "packages/radia-mcp/src/radia_mcp/ltspice/server.py",
])
def test_non_gpu_changes_do_not_acquire_gpu_runner(path):
    assert not decide([path])


@pytest.mark.parametrize("path", [
    "packages/radia-simulink/native/radia_ams_cuda.cu",
    "packages/radia-simulink/native/shared.cuh",
    "packages/radia-simulink/tests/native/CMakeLists.txt",
    "packages/radia-simulink/matlab/+radia/cudaGateway.m",
    "packages/radia-simulink/tools/build_ams_cuda.py",
    "packages/radia-simulink/src/radia_simulink/_cuda_runtime.py",
    "packages/radia-simulink/src/radia_simulink/resident_amg.py",
    "packages/radia-simulink/src/radia_simulink/hdiv.py",
    "packages/radia-simulink/tests/test_ams_acceptance.py",
    "packages/radia-simulink/native_payload.py",
    "validation_test/ams_cuda/run.py", "validation_test/ic_cuda/input.vol",
    "validation_test/ngsolve_cuda_mex/run_poisson.py",
    "validation_test/hmatrix_cuda/run.py", "validation_test/fmm_cuda/run.py",
    "validation_test/radia_simulink/run.py",
    "src/core/rad_hdiv_vim.h", "src/ext/HACApK/kernel.c",
    "src/ext/sparsesolv/include/sparsesolv/preconditioners/compact_ams.hpp",
    "tests/cpp/test_radia_ih_thermal_acceptance.cpp",
    "Build.ps1", "CMakeLists.txt", "tools/ci_gpu_artifact.py",
    "tools/ci_tree_receipt.py", "tools/sparsesolv_ci_impact.py",
    ".github/workflows/radia-simulink.yml",
])
def test_required_gpu_contracts_still_acquire_gpu_runner(path):
    assert decide([path])


@pytest.mark.parametrize("path", ["pyproject.toml", "packages/radia-simulink/pyproject.toml"])
def test_version_only_metadata_skips_but_dependency_changes_select(path):
    old = '[project]\nversion="1.0"\ndependencies=["numpy"]\n'
    new = old.replace('1.0', '1.1')
    blobs = {"base:" + path: old, "HEAD:" + path: new}
    assert not decide([path], blobs)
    blobs["HEAD:" + path] = new.replace('numpy', 'numpy>=2')
    assert decide([path], blobs)
    blobs["HEAD:" + path] = '[project]\ndependencies=["numpy"]\n'
    assert decide([path], blobs)  # Removing a version is not a version bump.
    assert decide([path])  # Missing/deleted/malformed input cannot exempt a build.


def test_python_version_only_exemption_does_not_hide_initialization_changes():
    path = "packages/radia-simulink/src/radia_simulink/__init__.py"
    old = '__version__ = "1.0"\nfrom . import cli\n'
    blobs = {"base:" + path: old, "HEAD:" + path: old.replace('1.0', '1.1')}
    assert not decide([path], blobs)
    blobs["HEAD:" + path] += 'from . import _cuda_runtime\n'
    assert decide([path], blobs)
    blobs["HEAD:" + path] = 'from . import cli\n'
    assert decide([path], blobs)


def test_mixed_changes_and_both_rename_sides_keep_native_gate():
    native = "packages/radia-simulink/native/kernel.cu"
    assert decide(["packages/radia-simulink/README.md", native])
    assert decide([], raw="R100\0" + native + "\0archive/kernel.txt\0")
    assert decide([], raw="R100\0archive/kernel.txt\0" + native + "\0")
    assert decide([], raw="D\0" + native + "\0")


@pytest.mark.parametrize("event", ["workflow_dispatch", "schedule", "release", "unknown"])
def test_qualification_and_unknown_events_never_skip_gpu(event):
    assert decide([], event=event)


def test_bad_event_and_comparison_fail_closed():
    assert impact.gpu_required_for_event("pull_request", {}, None, {})[0]
    assert decide([], raw="R100\0truncated\0")


def test_workflow_consumes_gpu_decision_and_cpu_ams_is_unchanged():
    workflow = yaml.load((ROOT / ".github/workflows/radia-simulink.yml").read_text(), Loader=yaml.BaseLoader)
    assert workflow["jobs"]["ci-impact"]["outputs"]["gpu_required"] == "${{ steps.route.outputs.gpu_required }}"
    assert "outputs.gpu_required == 'true'" in workflow["jobs"]["cuda-mex"]["if"]
    cpu = yaml.load((ROOT / ".github/workflows/sparsesolv.yml").read_text(), Loader=yaml.BaseLoader)
    assert "gpu-rtx" not in cpu["jobs"]["ams-regression"]["runs-on"]


def test_cli_emits_separate_gpu_output(tmp_path, monkeypatch):
    event = tmp_path / "event.json"
    event.write_text(json.dumps({}))
    output = tmp_path / "output"
    monkeypatch.setenv("GITHUB_EVENT_PATH", str(event))
    monkeypatch.setenv("GITHUB_EVENT_NAME", "pull_request")
    monkeypatch.setenv("GITHUB_OUTPUT", str(output))
    monkeypatch.setattr(impact, "required_for_event", lambda *args: (True, "distribution selected"))
    monkeypatch.setattr(impact, "gpu_required_for_event", lambda *args: (False, "no GPU inputs"))
    impact.main(["--lane", "radia-simulink"])
    values = dict(line.split("=", 1) for line in output.read_text().splitlines())
    assert values["required"] == "true"
    assert values["gpu_required"] == "false"


@pytest.mark.parametrize("event", ["pull_request", "push"])
@pytest.mark.parametrize("extra", [None, "src/radia/ih_thermal.py", "src/core/rad_hdiv_vim.h", "CMakeLists.txt"])
def test_ams_version_exemption_cannot_hide_numerical_or_build_changes(extra, event):
    path = "packages/radia-simulink/pyproject.toml"
    before = '[project]\nversion="0.1.5rc1"\ndependencies=["numpy"]\n'
    after = before.replace('0.1.5rc1', '0.1.5')
    def git(*args):
        if args[0] == "merge-base": return "base"
        if args[0] == "diff": return "M\0" + path + "\0" + ("M\0" + extra + "\0" if extra else "")
        if args[0] == "show": return before if args[1].startswith("base:") else after
        raise AssertionError(args)
    selected = impact.required_for_event("sparsesolv", event,
        {"before": "base", "pull_request": {"base": {"sha": "base"}}}, git, impact.manifest())[0]
    assert selected is bool(extra)


@pytest.mark.parametrize("after", [
    '[project]\nversion="2"\ndependencies=["numpy>=2"]\n',
    '[project]\nversion="2"\ndependencies=["numpy"]\n[build-system]\nrequires=["setuptools>=90"]\n',
    '[project]\ndependencies=["numpy"]\n', 'not valid toml',
])
def test_ams_dependency_build_or_invalid_metadata_still_selects_full_gate(after):
    path = "packages/radia-simulink/pyproject.toml"
    def git(*args):
        if args[0] == "merge-base": return "base"
        if args[0] == "diff": return "M\0" + path + "\0"
        if args[0] == "show":
            return '[project]\nversion="1"\ndependencies=["numpy"]\n' if args[1].startswith("base:") else after
        raise AssertionError(args)
    assert impact.required_for_event("sparsesolv", "pull_request",
        {"pull_request": {"base": {"sha": "base"}}}, git, impact.manifest())[0]


def test_ams_qualification_still_runs_full_gate():
    for event in ("workflow_dispatch", "schedule", "release"):
        assert impact.required_for_event("sparsesolv", event, {}, None, impact.manifest())[0]
