import json

import pytest

from radia_mcp.mathematica import tools
from radia_mcp.presentation import _equation_cli


@pytest.fixture
def evaluation(monkeypatch):
    calls = []

    def evaluate(code, timeout):
        calls.append((code, timeout))
        return {
            "exit_code": 0,
            "timed_out": False,
            "stderr": "",
            "result": json.dumps({"ok": True, "tex": r"\frac{1}{2}", "input_form": "1/2"}),
        }

    monkeypatch.setattr(tools, "mathematica_evaluate", evaluate)
    return calls


def test_result_evaluated_once_and_preserves_tex(evaluation):
    result = tools.mathematica_export_equation("Integrate[x, {x, 0, 1}]")
    assert result["ok"] and result["tex"] == r"\frac{1}{2}"
    assert result["input_form"] == "1/2"
    assert len(evaluation) == 1
    assert "MatrixForm[result]" in evaluation[0][0]


def test_save_creates_editor_input(evaluation, tmp_path):
    output = tmp_path / "equations" / "integral.tex"
    result = tools.mathematica_export_equation("1/2", action="save", output_path=str(output))
    assert result["ok"]
    assert output.read_text(encoding="utf-8") == "\\frac{1}{2}\n"


@pytest.mark.parametrize(
    "arguments",
    [
        {"action": "unknown"},
        {"action": "save"},
        {"action": "render", "output_path": "bad.pdf"},
        {"action": "copy", "target": "unknown"},
        {"timeout": float("nan")},
        {"timeout": True},
    ],
)
def test_bad_input_does_not_start_kernel(evaluation, arguments):
    result = tools.mathematica_export_equation("1/2", **arguments)
    assert not result["ok"] and result["stage"] == "input"
    assert evaluation == []


@pytest.mark.parametrize(
    "stdout",
    ['{"ok":false}', '"not an equation"', '{"ok":true,"tex":"","input_form":"x"}', "invalid JSON"],
)
def test_failed_evaluation_never_copies(monkeypatch, stdout):
    monkeypatch.setattr(
        tools,
        "mathematica_evaluate",
        lambda *a, **k: {"exit_code": 0, "timed_out": False, "result": stdout},
    )
    monkeypatch.setattr(
        _equation_cli,
        "presentation_copy_equation",
        lambda *a, **k: pytest.fail("editor called after failed evaluation"),
    )
    result = tools.mathematica_export_equation("x", action="copy")
    assert not result["ok"] and result["stage"] == "evaluation"


def test_clipboard_uses_calculated_result(evaluation, monkeypatch):
    calls = []

    def copy(tex, **kwargs):
        calls.append((tex, kwargs))
        return {"ok": True}

    monkeypatch.setattr(_equation_cli, "presentation_copy_equation", copy)
    result = tools.mathematica_export_equation("1/2", action="copy")
    assert result["ok"]
    assert calls == [(r"\frac{1}{2}", {"target": "office", "executable": None})]


def test_editor_failure_keeps_tex(evaluation, monkeypatch, tmp_path):
    monkeypatch.setattr(
        _equation_cli,
        "presentation_render_equation",
        lambda *a, **k: {"ok": False, "error": "unsupported notation"},
    )
    result = tools.mathematica_export_equation(
        "1/2", action="render", output_path=str(tmp_path / "eq.png")
    )
    assert not result["ok"] and result["stage"] == "editor"
    assert result["tex"] == r"\frac{1}{2}"


def test_new_tool_registered_without_kernel(monkeypatch):
    monkeypatch.setattr(
        tools, "mathematica_evaluate", lambda *a, **k: pytest.fail("discovery launched kernel")
    )
    from radia_mcp.mathematica import server

    assert "mathematica_export_equation" in server._REGISTERED


def test_wolframscript_trailing_null_preserves_matrix_lines(monkeypatch):
    tex = "\\begin{array}{cc}\n1 & 2 \\\\\n3 & 4\n\\end{array}"
    payload = json.dumps({"ok": True, "tex": tex, "input_form": "{{1,2},{3,4}}"})

    def evaluate(code, timeout):
        assert '"Compact"->True' in code
        return {"exit_code": 0, "timed_out": False, "result": payload + "\nNull"}

    monkeypatch.setattr(tools, "mathematica_evaluate", evaluate)
    result = tools.mathematica_export_equation("{{1,2},{3,4}}")
    assert result["ok"] and result["wolfram_tex"] == tex
    assert result["tex"] == tex.replace(r"\begin{array}{cc}", r"\begin{matrix}").replace(
        r"\end{array}", r"\end{matrix}"
    )


def test_noncentered_array_is_not_silently_changed(monkeypatch):
    tex = r"\begin{array}{lr}1 & 2\end{array}"
    monkeypatch.setattr(
        tools,
        "mathematica_evaluate",
        lambda *a, **k: {
            "exit_code": 0,
            "timed_out": False,
            "result": json.dumps({"ok": True, "tex": tex, "input_form": "matrix"}),
        },
    )
    result = tools.mathematica_export_equation("matrix")
    assert not result["ok"] and result["stage"] == "conversion"
    assert result["tex"] == tex
