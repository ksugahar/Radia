from types import SimpleNamespace
from pathlib import Path

import pytest
from cubit_mesh_export.mcp import session


@pytest.mark.parametrize("platform", ["win32", "linux", "darwin"])
def test_headless_journal_can_select_isolated_command_plugin_directory(
    tmp_path, monkeypatch, platform
):
    monkeypatch.setenv("CUBIT_MCP_TEMP", str(tmp_path / "scratch"))
    monkeypatch.setattr(session.sys, "platform", platform)
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "coreform_cubit.com").write_bytes(b"")
    plugin_dir = tmp_path / "plugin-under-test"
    plugin_dir.mkdir()
    (plugin_dir / "exporter.ccm").write_bytes(b"candidate")
    captured = {}

    def fake_run(argv, **kwargs):
        captured["argv"] = argv
        captured["kwargs"] = kwargs
        return SimpleNamespace(returncode=0, stdout="exported", stderr="")

    monkeypatch.setattr(session, "get_cubit_bin_dir", lambda: bin_dir)
    monkeypatch.setattr(session.subprocess, "run", fake_run)

    result = session.run_headless_journal(
        ["reset", "export netgen \"mesh.vol\" overwrite"],
        working_directory=tmp_path,
        command_plugin_directory=plugin_dir,
    )

    assert result["status"] == "completed"
    assert result["gui_started"] is False
    assert result["command_plugin_directory"] == str(plugin_dir)
    assert result["user_init_loaded"] is False
    assert "-noinitfile" in result["headless_flags"]
    assert "-commandplugindir" in result["headless_flags"]
    plugin_arg = captured["argv"].index("-commandplugindir")
    assert captured["argv"][plugin_arg + 1] == str(plugin_dir)
    assert captured["argv"][-1].endswith("driver.jou")
    assert session.Path(captured["argv"][-1]).is_relative_to(tmp_path / "scratch")
    assert captured["kwargs"]["stdin"] is session.subprocess.DEVNULL


def test_headless_journal_rejects_missing_command_plugin_directory(
    tmp_path, monkeypatch
):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "coreform_cubit.com").write_bytes(b"")
    monkeypatch.setattr(session, "get_cubit_bin_dir", lambda: bin_dir)

    result = session.run_headless_journal(
        ["reset"], command_plugin_directory=tmp_path / "missing"
    )

    assert result["status"] == "error"
    assert result["stage"] == "preflight"
    assert result["kind"] == "input"
    assert result["gui_started"] is False
    assert "command plugin directory not found" in result["error"]


def test_modern_console_stages_candidate_plugins_without_modifying_sources(tmp_path, monkeypatch):
    scratch = tmp_path / "scratch"
    monkeypatch.setenv("CUBIT_MCP_TEMP", str(scratch))
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "cubitx.exe").touch()
    acis = tmp_path / "acis/code/bin"
    acis.mkdir(parents=True)
    plugins = tmp_path / "candidate"
    plugins.mkdir()
    files = {"cubit_mesh_export.ccm": b"exporter", "dependency.dll": b"dependency",
             "ignored.txt": b"not a plugin"}
    for name, data in files.items():
        (plugins / name).write_bytes(data)
    staged_paths = []

    def fake_run(argv, **kwargs):
        staged = Path(kwargs["env"]["CUBIT_PLUGIN_DIR"])
        staged_paths.append(staged)
        assert staged.is_relative_to(scratch)
        assert staged.parent == Path(argv[-1]).parent
        assert {p.name: p.read_bytes() for p in staged.iterdir()} == {
            "cubit_mesh_export.dll": b"exporter", "dependency.dll": b"dependency"}
        assert kwargs["env"]["PATH"].startswith(str(acis))
        assert "-noinitfile" in argv and "-commandplugindir" not in argv
        assert {p.name: p.read_bytes() for p in plugins.iterdir()} == files
        return SimpleNamespace(returncode=0, stdout="exported", stderr="")

    monkeypatch.setattr(session, "get_cubit_bin_dir", lambda: bin_dir)
    monkeypatch.setattr(session.subprocess, "run", fake_run)
    result = session.run_headless_journal(["reset"], command_plugin_directory=plugins)
    assert result["status"] == "completed"
    assert result["gui_started"] is False
    assert "-commandplugindir" not in result["headless_flags"]
    assert {p.name: p.read_bytes() for p in plugins.iterdir()} == files
    assert not staged_paths[0].exists()


def test_modern_console_rejects_colliding_candidate_plugin_names(tmp_path, monkeypatch):
    monkeypatch.setenv("CUBIT_MCP_TEMP", str(tmp_path / "scratch"))
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "cubitx.exe").touch()
    plugins = tmp_path / "candidate"
    plugins.mkdir()
    (plugins / "exporter.ccm").write_bytes(b"one")
    (plugins / "exporter.dll").write_bytes(b"two")
    monkeypatch.setattr(session, "get_cubit_bin_dir", lambda: bin_dir)
    monkeypatch.setattr(session.subprocess, "run", lambda *a, **k: pytest.fail("must not launch"))
    with pytest.raises(ValueError, match="Conflicting candidate plugin DLL: exporter.dll"):
        session.run_headless_journal(["reset"], command_plugin_directory=plugins)


@pytest.mark.parametrize("console", ["cubitx.exe", "coreform_cubit.com"])
@pytest.mark.parametrize("non_plugin_files", [False, True])
def test_candidate_directory_requires_plugin_files_before_launch(
        tmp_path, monkeypatch, console, non_plugin_files):
    monkeypatch.setenv("CUBIT_MCP_TEMP", str(tmp_path / "scratch"))
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / console).touch()
    plugins = tmp_path / "candidate"
    plugins.mkdir()
    if non_plugin_files:
        (plugins / "readme.txt").write_text("no plugins")
        (plugins / "not-a-file.dll").mkdir()
    monkeypatch.setattr(session, "get_cubit_bin_dir", lambda: bin_dir)
    monkeypatch.setattr(session.subprocess, "run", lambda *a, **k: pytest.fail("must not launch"))
    with pytest.raises(ValueError) as error:
        session.run_headless_journal(["reset"], command_plugin_directory=plugins)
    assert str(plugins) in str(error.value)
    assert "No .ccm/.dll candidate plugins" in str(error.value)
