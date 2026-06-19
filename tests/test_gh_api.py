from __future__ import annotations

import importlib.util
from pathlib import Path
from types import SimpleNamespace


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "radia_tools_gh_api", ROOT / "tools" / "gh_api.py"
)
gh_api = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(gh_api)


def _clear_env(monkeypatch) -> None:
    monkeypatch.delenv("GH_TOKEN", raising=False)
    monkeypatch.delenv("GITHUB_TOKEN", raising=False)
    monkeypatch.delenv("RADIA_GH", raising=False)


def _redirect_home(monkeypatch, tmp_path: Path) -> None:
    def fake_expanduser(value: str) -> str:
        if value.startswith("~/"):
            return str(tmp_path / value[2:])
        return value

    monkeypatch.setattr(gh_api.os.path, "expanduser", fake_expanduser)


def test_token_prefers_gh_token_over_github_token(monkeypatch, tmp_path):
    _redirect_home(monkeypatch, tmp_path)
    monkeypatch.setenv("GH_TOKEN", "gh-token")
    monkeypatch.setenv("GITHUB_TOKEN", "github-token")

    assert gh_api.token() == "gh-token"


def test_token_reads_radia_token_file_before_gh_cli(monkeypatch, tmp_path):
    _clear_env(monkeypatch)
    _redirect_home(monkeypatch, tmp_path)
    token_file = tmp_path / ".radia" / "gh_token"
    token_file.parent.mkdir()
    token_file.write_text("file-token\n", encoding="utf-8")
    monkeypatch.setattr(gh_api.shutil, "which", lambda _name: "C:/bin/gh.exe")

    def fail_run(*_args, **_kwargs):
        raise AssertionError("gh CLI should not be queried when token file exists")

    monkeypatch.setattr(gh_api.subprocess, "run", fail_run)

    assert gh_api.token() == "file-token"


def test_token_falls_back_to_gh_cli(monkeypatch, tmp_path):
    _clear_env(monkeypatch)
    _redirect_home(monkeypatch, tmp_path)
    monkeypatch.setattr(gh_api.shutil, "which", lambda _name: "C:/bin/gh.exe")
    calls = []

    def fake_run(args, **_kwargs):
        calls.append(args)
        return SimpleNamespace(returncode=0, stdout="cli-token\n")

    monkeypatch.setattr(gh_api.subprocess, "run", fake_run)

    assert gh_api.token() == "cli-token"
    assert calls[0][1:] == ["auth", "token", "--hostname", "github.com"]


def test_token_uses_radia_gh_when_gh_is_not_on_path(monkeypatch, tmp_path):
    _clear_env(monkeypatch)
    _redirect_home(monkeypatch, tmp_path)
    monkeypatch.setenv("RADIA_GH", "C:/portable/gh.exe")
    monkeypatch.setattr(gh_api.shutil, "which", lambda _name: None)

    def fake_run(args, **_kwargs):
        return SimpleNamespace(returncode=0, stdout=f"token-from-{args[0]}\n")

    monkeypatch.setattr(gh_api.subprocess, "run", fake_run)

    assert gh_api.token() == "token-from-C:/portable/gh.exe"
