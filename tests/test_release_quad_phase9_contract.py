"""Fail-closed solver-only cross-host probe contracts; no SSH required."""

import ast
import importlib.util
import subprocess
from pathlib import Path

import pytest


spec = importlib.util.spec_from_file_location(
    "quad_phase9", Path(__file__).resolve().parents[1] / "tools/release_quad.py")
quad = importlib.util.module_from_spec(spec)
spec.loader.exec_module(quad)


def probe():
    values = {"VER radia": "4.95.91"}
    for name in ("simulink/application.py", "panels/calc_inductance.py",
                 "panels/calc_fem_kelvin.py", "panels/calc_fem_coilmesh.py"):
        values["SHA radia/" + name] = "0123456789ab"
    return "\n".join(f"{key} = {value}" for key, value in values.items())


def test_required_keys_match_both_actual_probe_file_lists():
    for script in (quad.CROSS_MACHINE_PROBE, quad.CROSS_MACHINE_PROBE_LAB):
        tree = ast.parse(script)
        paths = next(ast.literal_eval(node.iter) for node in ast.walk(tree)
                     if isinstance(node, ast.For) and isinstance(node.iter, ast.List))
        assert {"SHA radia/" + path for path in paths} == {
            key for key in quad._PHASE9_FIELDS if key.startswith("SHA radia/")}
        assert "cubit_mesh_export" not in script
        assert "radia-mcp" not in script
    assert set(quad._parse_phase9_probe("LAB", probe())) == set(quad._PHASE9_FIELDS)


def test_phase9_accepts_complete_shuffled_output(monkeypatch):
    monkeypatch.setattr(
        quad, "_probe", lambda _label, *_args: "\n".join(reversed(probe().splitlines())))
    assert quad.cmd_phase9(None) == 0


@pytest.mark.parametrize("corrupt", [
    lambda text: "\n".join(text.splitlines()[:-1]),
    lambda text: text + "\n" + text.splitlines()[0],
    lambda text: text + "\nUNRECOGNIZED = 1",
    lambda text: text.replace("VER radia =", "VER radia "),
    lambda text: text.replace("4.95.91", ""),
    lambda text: text.replace("0123456789ab", "MISSING"),
    lambda text: text.replace("4.95.91", "MISSING"),
    lambda text: text.replace("4.95.91", "1garbage"),
    lambda text: text.replace("4.95.91", "N/A"),
    lambda text: text.replace("0123456789ab", "not-a-hash!!"),
])
def test_phase9_rejects_malformed_probe_even_when_all_hosts_agree(
        monkeypatch, corrupt):
    monkeypatch.setattr(quad, "_probe", lambda _label, *_args: corrupt(probe()))
    assert quad.cmd_phase9(None) == 4


def test_independent_package_drift_is_not_part_of_solver_probe(monkeypatch):
    assert "VER cubit-mesh-export" not in quad._PHASE9_FIELDS
    assert "VER radia-mcp" not in quad._PHASE9_FIELDS
    monkeypatch.setattr(quad, "_probe", lambda _label, *_args: probe())
    assert quad.cmd_phase9(None) == 0


@pytest.mark.parametrize("failure", [None, ""])
def test_failed_host_cannot_be_omitted(monkeypatch, failure):
    monkeypatch.setattr(
        quad, "_probe",
        lambda label, *_args: failure if label == "mdx2" else probe())
    assert quad.cmd_phase9(None) == 4


@pytest.mark.parametrize(
    "error", [OSError("no executable"), subprocess.TimeoutExpired("probe", 120)])
def test_probe_launch_failure_is_not_acceptance(monkeypatch, error):
    def run(*_args, **_kwargs):
        raise error
    monkeypatch.setattr(quad.subprocess, "run", run)
    assert quad._probe("LAB", ["python", "-"]) is None


def test_failed_exit_rejects_even_complete_stdout(monkeypatch):
    monkeypatch.setattr(
        quad.subprocess, "run",
        lambda *args, **_kwargs: subprocess.CompletedProcess(
            args, 1, stdout=probe(), stderr="failed after output"))
    assert quad._probe("LAB", ["python", "-"]) is None


def test_missing_version_parser_fails_closed(monkeypatch):
    import builtins

    original = builtins.__import__
    def without_packaging(name, *args, **kwargs):
        if name == "packaging.version":
            raise ModuleNotFoundError("packaging unavailable")
        return original(name, *args, **kwargs)
    monkeypatch.setattr(builtins, "__import__", without_packaging)
    monkeypatch.setattr(quad, "_probe", lambda _label, *_args: probe())
    assert quad.cmd_phase9(None) == 4


def test_pep440_versions_are_accepted():
    text = probe().replace("4.95.91", "5.0rc1")
    assert quad._parse_phase9_probe("LAB", text)["VER radia"] == "5.0rc1"
