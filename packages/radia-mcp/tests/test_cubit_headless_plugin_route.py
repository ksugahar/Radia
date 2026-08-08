from types import SimpleNamespace

from radia_mcp.cubit import session


def test_headless_journal_can_select_isolated_command_plugin_directory(
    tmp_path, monkeypatch
):
    bin_dir = tmp_path / "bin"
    bin_dir.mkdir()
    (bin_dir / "coreform_cubit.com").write_bytes(b"")
    plugin_dir = tmp_path / "plugin-under-test"
    plugin_dir.mkdir()
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
    assert result["persistent_gui_started"] is False
    assert result["command_plugin_directory"] == str(plugin_dir)
    assert result["user_init_loaded"] is False
    assert "-noinitfile" in result["headless_flags"]
    assert "-commandplugindir" in result["headless_flags"]
    plugin_arg = captured["argv"].index("-commandplugindir")
    assert captured["argv"][plugin_arg + 1] == str(plugin_dir)
    assert captured["argv"][-1].endswith("driver.jou")
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
    assert "command plugin directory not found" in result["error"]
