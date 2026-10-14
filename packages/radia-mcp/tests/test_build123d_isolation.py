"""build123d/CadQuery user scripts run isolated and report what they produced."""

import json

import pytest

pytest.importorskip("build123d")

from radia_mcp.build123d import server  # noqa: E402


@pytest.fixture(autouse=True)
def _no_failure_log(monkeypatch):
    monkeypatch.setattr(server._fl, "record_failure", lambda *a, **k: None)


def _inspect(script, **kwargs):
    return json.loads(server._inspect_build123d_in_process(script, **kwargs))


def test_last_reassigned_shape_is_selected_not_first_bound():
    info = _inspect("result = Box(10, 10, 10)\ntool = Cylinder(2, 20)\nresult = result - tool")
    assert info["variable"] == "result"
    assert info["volume"] < 1000.0


def test_future_import_and_triple_quoted_text_are_kept_verbatim():
    info = _inspect('from __future__ import annotations\ndoc = """\n  x\n"""\nb = Box(1, 2, 3)')
    assert info["status"] == "ok" and info["volume"] == pytest.approx(6.0)


def test_system_exit_is_a_reported_script_error():
    info = _inspect("import sys\nsys.exit(3)")
    assert info["status"] == "error" and "SystemExit" in info["error"]


def test_export_name_is_sanitized_and_file_checked(tmp_path):
    info = _inspect('b = Box(1, 1, 1)\nb.label = "yoke/pole:1"', export_dir=str(tmp_path))
    assert info["status"] == "ok"
    assert (tmp_path / "yoke_pole_1.step").is_file()
    assert info["exported"].endswith("yoke_pole_1.step")


def test_unknown_export_format_fails_instead_of_reporting_ok(tmp_path):
    info = _inspect("b = Box(1, 1, 1)", export_dir=str(tmp_path), export_format="iges")
    assert info["status"] == "error" and info["stage"] == "export"
    assert "exported" not in info


def test_process_exit_ends_only_the_isolated_child():
    info = json.loads(server._execute_build123d_sync("import os\nos._exit(7)", timeout_s=120))
    assert info["status"] == "error"
    assert "returncode 7" in info["error"]


def test_race_rejects_unknown_rule_before_spawning():
    info = json.loads(server._build123d_try_race_sync(["b = Box(1, 1, 1)"], prefer="bogus"))
    assert info["status"] == "error" and "bogus" in info["error"]
