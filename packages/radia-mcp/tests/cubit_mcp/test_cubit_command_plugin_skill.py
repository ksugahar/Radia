"""The packaged cubit-command-plugin skill and its scaffold tool.

Guards: the skill and templates ship as package data, the cpp_sdk guide
serves them, and the scaffold writes exact bytes and cleans up only its own
files on failure.
Assumes: no Cubit; pure file and text checks.
"""
import re

import pytest

from radia_mcp.cubit import server
from radia_mcp.cubit.knowledge import command_plugin_skill as skill
from radia_mcp.cubit.knowledge.cpp_sdk import get_cpp_sdk_documentation


def test_skill_has_frontmatter_and_names_the_verified_release():
    text = skill.get_command_plugin_skill()
    assert text.startswith("---\nname: cubit-command-plugin\n")
    assert "Cubit 2026.8" in text
    assert "_DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR" in text


@pytest.mark.parametrize("topic", ["command_plugin", "plugin", "skill", "ccm"])
def test_cpp_sdk_guide_serves_the_skill(topic):
    assert get_cpp_sdk_documentation(topic) == skill.get_command_plugin_skill()


def test_cpp_sdk_guide_serves_every_template():
    text = get_cpp_sdk_documentation("templates")
    for name in skill.TEMPLATE_FILES:
        assert f"## templates/{name}" in text


def test_custom_commands_topic_uses_the_real_sdk_api():
    text = get_cpp_sdk_documentation("custom_commands")
    assert "CUBIT_PLUGIN(MyPlugin)" in text
    assert "execute(CubitCommandData& data)" in text
    assert "command_parser()" not in text and "CubitStatus execute" not in text


def test_template_implements_the_sdk_contract():
    plugin = skill.read_template("MyPlugin.cpp")
    assert len(re.findall(r"^CUBIT_PLUGIN\(MyPlugin\)", plugin, re.M)) == 1
    for header in ("HelloCommand.hpp", "ExportNodeCountCommand.hpp"):
        body = skill.read_template(header)
        for method in ("get_syntax()", "get_syntax_help()", "get_help()",
                       "execute(CubitCommandData& data)"):
            assert method in body, (header, method)
    cmake = skill.read_template("CMakeLists.txt")
    assert "find_package(Cubit REQUIRED CONFIG)" in cmake
    assert "_DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR" in cmake
    assert "cubit_plugin_verification;debug=no" in skill.read_template("build.ps1")


def test_scaffold_copies_bytes_and_refuses_overwrite(tmp_path):
    target = tmp_path / "plug"
    written = skill.scaffold_command_plugin(target)
    assert {p.name for p in written} == {*skill.TEMPLATE_FILES, "SKILL.md"}
    for name in skill.TEMPLATE_FILES:
        assert (target / name).read_text(encoding="utf-8") == skill.read_template(name)
    with pytest.raises(FileExistsError):
        skill.scaffold_command_plugin(target)


def test_scaffold_preflight_failure_creates_nothing(tmp_path, monkeypatch):
    real = skill._skill_dir

    class Broken:
        def joinpath(self, *parts):
            node = real().joinpath(*parts)
            if parts[-1] == "SKILL.md":
                class Missing:
                    def read_bytes(self):
                        raise OSError("gone")
                return Missing()
            return node

    monkeypatch.setattr(skill, "_skill_dir", lambda: Broken())
    with pytest.raises(FileNotFoundError):
        skill.scaffold_command_plugin(tmp_path / "plug")
    assert not (tmp_path / "plug").exists()


def _inject(monkeypatch, stage, at=3):
    """Fail the at-th exclusive create at `stage`: open, write or close."""
    import builtins

    real_open = builtins.open
    calls = {"n": 0}

    class Faulty:
        def __init__(self, inner):
            self.inner = inner

        def __enter__(self):
            return self

        def __exit__(self, *exc):
            self.inner.close()
            if stage == "close":
                raise OSError("close failed")
            return False

        def write(self, data):
            if stage == "write":
                self.inner.write(data[:5])
                raise OSError("disk full")
            return self.inner.write(data)

    def fake_open(path, mode="r", *args, **kwargs):
        if mode != "xb":
            return real_open(path, mode, *args, **kwargs)
        calls["n"] += 1
        if calls["n"] != at:
            return real_open(path, mode, *args, **kwargs)
        if stage == "open":
            raise OSError("open failed")
        return Faulty(real_open(path, mode, *args, **kwargs))

    monkeypatch.setattr(builtins, "open", fake_open)
    return real_open


@pytest.mark.parametrize("stage,message", [
    ("open", "open failed"), ("write", "disk full"), ("close", "close failed")])
def test_scaffold_failure_at_each_stage_removes_only_its_own_files(
        tmp_path, monkeypatch, stage, message):
    import builtins

    real_open = _inject(monkeypatch, stage)
    with pytest.raises(OSError, match=message):
        skill.scaffold_command_plugin(tmp_path / "plug")
    assert not (tmp_path / "plug").exists()
    monkeypatch.setattr(builtins, "open", real_open)
    assert len(skill.scaffold_command_plugin(tmp_path / "plug")) == len(skill.TEMPLATE_FILES) + 1


def test_scaffold_cleanup_failure_keeps_primary_error_and_continues(tmp_path, monkeypatch):
    from pathlib import Path

    _inject(monkeypatch, "write", at=3)
    real_unlink = Path.unlink
    stuck = []

    def unlink(self, missing_ok=False):
        if self.name == skill.TEMPLATE_FILES[0]:
            stuck.append(self)
            raise PermissionError("locked")
        return real_unlink(self, missing_ok=missing_ok)

    monkeypatch.setattr(Path, "unlink", unlink)
    with pytest.raises(OSError, match="disk full") as info:
        skill.scaffold_command_plugin(tmp_path / "plug")
    notes = "\n".join(getattr(info.value, "__notes__", []))
    assert "could not remove" in notes and skill.TEMPLATE_FILES[0] in notes
    remaining = sorted(p.name for p in (tmp_path / "plug").iterdir())
    assert remaining == [skill.TEMPLATE_FILES[0]]  # the others were still removed


def test_templates_carry_the_review_fixes():
    smoke = skill.read_template("smoke.ps1")
    assert "if ($exit -ne 0)" in smoke and "output file contents differ" in smoke
    build = skill.read_template("build.ps1")
    assert "ReparsePoint" in build and "strictly inside" in build
    export = skill.read_template("ExportNodeCountCommand.cpp")
    assert 'L"wbx"' in export and "initialize_export()" in export
    assert "id_from_handle" in export and "get_block_size(block, n)" in export


def test_scaffold_tool_validates_and_writes(tmp_path):
    assert server.cubit_scaffold_command_plugin(str(tmp_path), "../x").startswith("ERROR")
    assert server.cubit_scaffold_command_plugin(str(tmp_path / "missing")).startswith("ERROR")
    out = server.cubit_scaffold_command_plugin(str(tmp_path), "demo")
    assert "Command-plugin template written to" in out
    assert (tmp_path / "demo" / "smoke.ps1").is_file()
    assert server.cubit_scaffold_command_plugin(str(tmp_path), "demo").startswith("ERROR")


def test_scaffold_tool_is_classified_as_writing():
    assert "cubit_scaffold_command_plugin" not in server._UNCLASSIFIED_TOOLS
    assert "cubit_scaffold_command_plugin" in server._WRITING_TOOLS
