from __future__ import annotations

from pathlib import Path
import subprocess

from radia_mcp.presentation import _equation_cli as equation_cli
from radia_mcp.presentation import tools as presentation_tools


def _fake_executable(tmp_path: Path) -> Path:
    app = tmp_path / "Eqnedit64.exe"
    app.write_bytes(b"test executable placeholder")
    return app


def test_presentation_tools_export_eqnedit64_bridge() -> None:
    assert presentation_tools.presentation_equation_policy
    assert presentation_tools.presentation_equation_backend
    assert presentation_tools.presentation_copy_equation
    assert presentation_tools.presentation_render_equation


def test_equation_policy_keeps_both_editions_in_radia() -> None:
    policy = equation_cli.presentation_equation_policy()
    assert policy["source_of_truth"]["native"] == "tools/eqnedit64/src"
    assert policy["source_of_truth"]["web"].startswith("tools/eqnedit64/web/")
    assert policy["homepage_is_source_of_truth"] is False
    assert policy["canonical_input"] == "TeX"
    assert policy["retired_formats"] == ["MTEF", ".eqn"]
    publication = policy["web_publication_contract"]
    assert publication["mode"] == "build-time import from Radia checkout"
    assert publication["checkout_environment_variable"] == "RADIA_REPOSITORY"
    assert publication["integrity"] == "SHA-256 equality after copy"
    assert publication["homepage_source_copy"] is False
    assert publication["release_gate"].endswith("run_eqnedit64_release_qa.ps1")
    assert "PowerPoint" in publication["test_scope"]


def test_copy_office_uses_utf8_file_contract(tmp_path, monkeypatch) -> None:
    app = _fake_executable(tmp_path)
    observed = {}

    def fake_invoke(command, timeout_s):
        observed["command"] = command
        observed["timeout_s"] = timeout_s
        observed["tex"] = Path(command[2]).read_text(encoding="utf-8")
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(equation_cli, "_invoke", fake_invoke)
    result = equation_cli.presentation_copy_equation(
        r"\[\frac{x}{y}\]", executable=str(app)
    )

    assert result["ok"] is True
    assert observed["command"][1] == "--copy-tex-file"
    assert observed["tex"] == r"\[\frac{x}{y}\]"
    assert "MathML" in result["formats"]
    assert "CF_ENHMETAFILE" in result["formats"]
    assert not Path(observed["command"][2]).exists()


def test_copy_target_routes_are_explicit(tmp_path, monkeypatch) -> None:
    app = _fake_executable(tmp_path)
    switches = []

    def fake_invoke(command, timeout_s):
        switches.append(command[1])
        return subprocess.CompletedProcess(command, 0, "", "")

    monkeypatch.setattr(equation_cli, "_invoke", fake_invoke)
    google = equation_cli.presentation_copy_equation(
        "x", target="google_slides", executable=str(app)
    )
    png = equation_cli.presentation_copy_equation(
        "x", target="png", executable=str(app)
    )

    assert google["ok"] and png["ok"]
    assert switches == ["--copy-google-slides-file", "--copy-png-file"]
    assert google["formats"] == ["HTML Format", "PNG"]
    assert png["formats"] == ["PNG", "CF_DIBV5"]


def test_render_png_reports_checked_artifact(tmp_path, monkeypatch) -> None:
    app = _fake_executable(tmp_path)
    output = tmp_path / "nested" / "equation.png"
    observed = {}

    def fake_invoke(command, timeout_s):
        observed["command"] = command
        Path(command[3]).write_bytes(b"\x89PNG\r\n\x1a\nfixture")
        return subprocess.CompletedProcess(command, 0, "294 243\n", "")

    monkeypatch.setattr(equation_cli, "_invoke", fake_invoke)
    result = equation_cli.presentation_render_equation(
        r"E=mc^2", str(output), executable=str(app)
    )

    assert result["ok"] is True
    assert observed["command"][1] == "--render-png-file"
    assert result["pixel_size"] == [294, 243]
    assert result["dpi"] == 300
    assert result["font_points"] == 24
    assert output.is_file()


def test_nonzero_backend_exit_fails_loudly(tmp_path, monkeypatch) -> None:
    app = _fake_executable(tmp_path)

    def fake_invoke(command, timeout_s):
        return subprocess.CompletedProcess(command, 84, "", "clipboard busy")

    monkeypatch.setattr(equation_cli, "_invoke", fake_invoke)
    result = equation_cli.presentation_copy_equation("x", executable=str(app))
    assert result["ok"] is False
    assert "exited 84" in result["error"]


def test_invalid_target_and_missing_backend_are_diagnostic(tmp_path) -> None:
    bad_target = equation_cli.presentation_copy_equation(
        "x", target="svg", executable=str(tmp_path / "Eqnedit64.exe")
    )
    missing = equation_cli.presentation_equation_backend(
        str(tmp_path / "missing.exe")
    )
    assert bad_target["ok"] is False
    assert "target must be" in bad_target["error"]
    assert missing["ok"] is False
    assert "was not found" in missing["error"]
