"""Fail-closed cross-host probe contracts; no SSH or installed solver needed."""
import ast
import importlib.util
import subprocess
from pathlib import Path

import pytest


spec = importlib.util.spec_from_file_location(
    "quad_phase9", Path(__file__).resolve().parents[1] / "tools/release_quad.py")
quad = importlib.util.module_from_spec(spec)
spec.loader.exec_module(quad)


def probe(label):
    values = {
        "VER radia": "4.95.91", "VER cubit-mesh-export": "0.14.17",
        "VER radia-mcp": "1.4.53", "COMPAT cme -> radia": "[4.5.0, 4.999.999]",
        "COMPAT rad -> cme": "[0.5.0, 0.999.999]",
    }
    for name in ("panels/register_toolbar.py", "panels/radia_export_menu.py",
                 "simulink/application.py", "panels/calc_inductance.py",
                 "panels/calc_fem_kelvin.py", "panels/calc_fem_coilmesh.py"):
        values["SHA radia/" + name] = "0123456789ab"
    if label in ("mdx1", "mdx2"):
        for key in list(values)[1:5]:
            values[key] = "N/A"
    return "\n".join(f"{key} = {value}" for key, value in values.items())


def test_required_keys_match_both_actual_probe_file_lists():
    for script in (quad.CROSS_MACHINE_PROBE, quad.CROSS_MACHINE_PROBE_LAB):
        tree = ast.parse(script)
        paths = next(ast.literal_eval(n.iter) for n in ast.walk(tree)
                     if isinstance(n, ast.For) and isinstance(n.iter, ast.List))
        assert {"SHA radia/" + p for p in paths} == {
            k for k in quad._PHASE9_FIELDS if k.startswith("SHA ")}
    assert set(quad._parse_phase9_probe("LAB", probe("LAB"))) == set(quad._PHASE9_FIELDS)


def test_phase9_accepts_complete_shuffled_output(monkeypatch):
    monkeypatch.setattr(quad, "_probe", lambda label, *a: "\n".join(reversed(probe(label).splitlines())))
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
    lambda text: text.replace("[4.5.0, 4.999.999]", "invalid"),
    lambda text: text.replace("[4.5.0, 4.999.999]", "[1bad, 2bad]"),
    lambda text: text.replace("[4.5.0, 4.999.999]", "[5.0, 4.0]"),
    lambda text: text.replace("0123456789ab", "not-a-hash!!"),
])
def test_phase9_rejects_malformed_probe_even_when_all_hosts_agree(monkeypatch, corrupt):
    monkeypatch.setattr(quad, "_probe", lambda label, *a: corrupt(probe(label)))
    assert quad.cmd_phase9(None) == 4


@pytest.mark.parametrize("label", ["LAB", "100号機"])
def test_editable_hosts_cannot_hide_mcp_with_na(label):
    with pytest.raises(ValueError):
        quad._parse_phase9_probe(label, probe(label).replace("1.4.53", "N/A"))


def test_compute_host_must_declare_non_deployed_fields_na():
    with pytest.raises(ValueError):
        quad._parse_phase9_probe("mdx1", probe("LAB"))


def test_same_values_with_keys_swapped_cannot_hide_drift(monkeypatch):
    def run(label, *a):
        text = probe(label)
        if label == "100号機":
            text = text.replace("VER radia =", "TEMP =").replace(
                "VER radia-mcp =", "VER radia =").replace("TEMP =", "VER radia-mcp =")
        return text
    monkeypatch.setattr(quad, "_probe", run)
    assert quad.cmd_phase9(None) == 4


@pytest.mark.parametrize("failure", [None, ""])
def test_failed_host_cannot_be_omitted(monkeypatch, failure):
    monkeypatch.setattr(quad, "_probe", lambda label, *a: failure if label == "mdx2" else probe(label))
    assert quad.cmd_phase9(None) == 4


@pytest.mark.parametrize("error", [OSError("no executable"), subprocess.TimeoutExpired("probe", 120)])
def test_probe_launch_failure_is_not_acceptance(monkeypatch, error):
    def run(*args, **kwargs):
        raise error
    monkeypatch.setattr(quad.subprocess, "run", run)
    assert quad._probe("LAB", ["python", "-"]) is None


def test_failed_exit_rejects_even_complete_stdout(monkeypatch):
    monkeypatch.setattr(quad.subprocess, "run", lambda *a, **k: subprocess.CompletedProcess(
        a, 1, stdout=probe("LAB"), stderr="failed after output"))
    assert quad._probe("LAB", ["python", "-"]) is None


def test_missing_version_parser_fails_closed(monkeypatch):
    import builtins

    original = builtins.__import__
    def without_packaging(name, *args, **kwargs):
        if name == "packaging.version":
            raise ModuleNotFoundError("packaging unavailable")
        return original(name, *args, **kwargs)
    monkeypatch.setattr(builtins, "__import__", without_packaging)
    monkeypatch.setattr(quad, "_probe", lambda label, *a: probe(label))
    assert quad.cmd_phase9(None) == 4


def test_pep440_versions_and_inclusive_compatibility_bounds():
    text = probe("LAB").replace("4.95.91", "5.0rc1").replace(
        "[4.5.0, 4.999.999]", "[5.0rc1, 5.0]")
    assert quad._parse_phase9_probe("LAB", text)["VER radia"] == "5.0rc1"
    text = text.replace("[5.0rc1, 5.0]", "[5.0, 5.0]")
    assert quad._parse_phase9_probe("LAB", text)["COMPAT cme -> radia"] == "[5.0, 5.0]"


def test_mcp_preservation_does_not_waive_phase9_version_match(monkeypatch):
    monkeypatch.setenv("RADIA_RELEASE_PRESERVE_MCP_SOURCE_LAB", "C:/approved/mcp")
    monkeypatch.setenv("RADIA_RELEASE_PRESERVE_MCP_SOURCE_100", "W:/approved/mcp")
    monkeypatch.setattr(quad, "_probe", lambda label, *a: probe(label).replace(
        "1.4.53", "1.4.54") if label == "100号機" else probe(label))
    assert quad.cmd_phase9(None) == 4
