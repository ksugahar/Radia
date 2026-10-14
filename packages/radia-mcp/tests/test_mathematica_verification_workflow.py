from __future__ import annotations

import subprocess
import sys
from types import SimpleNamespace

from radia_mcp.mathematica import tools
from radia_mcp.mathematica import helpers


def test_run_script_parses_verification_json(tmp_path, monkeypatch):
    script = tmp_path / "course_check.wls"
    script.write_text('Print["placeholder"]', encoding="utf-8")
    calls = []

    def fake_run(arguments, *, stdin, stdout, stderr, timeout):
        calls.append((arguments, timeout))
        stdout.write(b'{"ok":true,"checks":{"lesson01":true},"failures":[]}')
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(tools, "_wolframscript_path", lambda: "wolframscript")
    monkeypatch.setattr(tools.subprocess, "run", fake_run)

    result = tools.mathematica_run_script(
        str(script), result_format="json", require_ok=True
    )

    assert result["ok"] is True
    assert result["process_ok"] is True
    assert result["verification_ok"] is True
    assert result["parsed"]["checks"] == {"lesson01": True}
    assert calls == [(["wolframscript", "-file", str(script.resolve())], 300)]


def test_run_script_surfaces_failed_named_checks(tmp_path, monkeypatch):
    script = tmp_path / "course_check.wls"
    script.write_text('Print["placeholder"]', encoding="utf-8")

    def fake_run(arguments, *, stdin, stdout, stderr, timeout):
        stdout.write(b'{"ok":false,"checks":{"lesson02":false},"failures":["lesson02"]}')
        return SimpleNamespace(returncode=0)

    monkeypatch.setattr(tools, "_wolframscript_path", lambda: "wolframscript")
    monkeypatch.setattr(tools.subprocess, "run", fake_run)

    result = tools.mathematica_run_script(
        str(script), result_format="json", require_ok=True
    )

    assert result["ok"] is False
    assert result["process_ok"] is True
    assert result["verification_ok"] is False
    assert "lesson02" in result["error"]


def test_run_script_rejects_non_wolfram_file(tmp_path):
    script = tmp_path / "course_check.py"
    script.write_text("pass", encoding="utf-8")

    result = tools.mathematica_run_script(str(script))

    assert result["ok"] is False
    assert ".wls" in result["error"]


def test_batch_identity_check_uses_one_kernel(monkeypatch):
    seen = {}

    def fake_evaluate(code, timeout):
        seen["code"] = code
        seen["timeout"] = timeout
        return {
            "result": '{"ok":true,"checks":{"gradient":true,"curlGrad":true},"failures":[]}',
            "stderr": "",
            "exit_code": 0,
            "timed_out": False,
        }

    monkeypatch.setattr(helpers, "mathematica_evaluate", fake_evaluate)
    result = helpers.mathematica_check_identities(
        [
            {"name": "gradient", "lhs": "Grad[x y, {x, y}]", "rhs": "{y, x}"},
            {
                "name": "curlGrad",
                "lhs": "Curl[Grad[f[x,y,z], {x,y,z}], {x,y,z}]",
                "rhs": "{0,0,0}",
                "assumptions": "Element[{x,y,z}, Reals]",
            },
        ]
    )

    assert result["ok"] is True
    assert result["failures"] == []
    assert seen["timeout"] == 180
    assert seen["code"].count("FullSimplify") == 2
    assert "gradient" in seen["code"] and "curlGrad" in seen["code"]


def test_batch_identity_check_rejects_duplicate_names():
    result = helpers.mathematica_check_identities(
        [
            {"name": "same", "lhs": "x", "rhs": "x"},
            {"name": "same", "lhs": "y", "rhs": "y"},
        ]
    )

    assert result["ok"] is False
    assert "duplicate" in result["error"]


def test_verification_guide_exposes_course_workflow():
    result = helpers.mathematica_verification_guide("course")

    assert result["ok"] is True
    assert result["topic"] == "electromagnetics"
    assert result["small_batch"] == "mathematica_check_identities"
    assert result["tracked_suite"] == "mathematica_run_script"
    assert len(result["starter_claims"]) == 2
    assert any("Canvas" in rule for rule in result["rules"])


def test_verification_guide_rejects_unknown_topic():
    result = helpers.mathematica_verification_guide("astrology")

    assert result["ok"] is False
    assert "electromagnetics" in result["error"]


def test_helpers_can_be_imported_before_tools():
    completed = subprocess.run(
        [
            sys.executable,
            "-c",
            (
                "from radia_mcp.mathematica.helpers import "
                "mathematica_verification_guide; "
                "assert mathematica_verification_guide('course')['ok']"
            ),
        ],
        check=False,
        capture_output=True,
        text=True,
        timeout=30,
    )

    assert completed.returncode == 0, completed.stderr
