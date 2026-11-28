"""Packaged skill for native Cubit C++ command plugins.

The skill lives as files under ``radia_mcp/cubit/skills/`` so it
survives the removal of the plugin sources from the repository and ships
inside the wheel.  ``SKILL.md`` is the guide; ``templates/`` is a minimal
plugin that builds and passes ``smoke.ps1`` on Coreform Cubit 2026.8.
"""
from __future__ import annotations

from importlib import resources
from pathlib import Path

SKILL_NAME = "cubit-command-plugin"

# Every file the scaffold copies, relative to templates/.  Listed explicitly
# so a missing file in a wheel is an error, not a silently smaller scaffold.
TEMPLATE_FILES = (
    "CMakeLists.txt",
    "MyPlugin.hpp",
    "MyPlugin.cpp",
    "HelloCommand.hpp",
    "HelloCommand.cpp",
    "ExportNodeCountCommand.hpp",
    "ExportNodeCountCommand.cpp",
    "utf8_path.hpp",
    "build.ps1",
    "smoke.ps1",
)


def _skill_dir():
    return resources.files("radia_mcp.cubit").joinpath("skills", SKILL_NAME)


def get_command_plugin_skill() -> str:
    """Return SKILL.md for the command-plugin skill."""
    return _skill_dir().joinpath("SKILL.md").read_text(encoding="utf-8")


def read_template(name: str) -> str:
    if name not in TEMPLATE_FILES:
        raise KeyError(f"unknown template file: {name!r}")
    return _skill_dir().joinpath("templates", name).read_text(encoding="utf-8")


def get_command_plugin_templates() -> str:
    """Return every template file, concatenated with headings."""
    parts = []
    for name in TEMPLATE_FILES:
        fence = "powershell" if name.endswith(".ps1") else (
            "cmake" if name == "CMakeLists.txt" else "cpp")
        parts.append(f"## templates/{name}\n\n```{fence}\n{read_template(name)}```\n")
    return "\n".join(parts)


def scaffold_command_plugin(target: Path) -> list[Path]:
    """Copy the templates into a new directory ``target``.

    ``target`` must not exist yet; its parent must.  Returns the written
    paths.  Bytes are copied unchanged so the scripts keep their encoding.
    """
    target = Path(target)
    if target.exists():
        raise FileExistsError(f"refusing to overwrite an existing directory: {target}")
    if not target.parent.is_dir():
        raise FileNotFoundError(f"parent directory does not exist: {target.parent}")
    # Read every packaged resource before touching the disk, so a broken
    # install fails without leaving a half-written directory behind.
    payload = {}
    try:
        for name in TEMPLATE_FILES:
            payload[name] = _skill_dir().joinpath("templates", name).read_bytes()
        payload["SKILL.md"] = _skill_dir().joinpath("SKILL.md").read_bytes()
    except OSError as exc:
        raise FileNotFoundError(f"packaged skill resources unreadable: {exc}") from exc
    target.mkdir()  # exclusive: raises FileExistsError if it appeared meanwhile
    created = []  # every file this call created, recorded as soon as it exists
    try:
        for name, data in payload.items():
            dest = target / name
            handle = open(dest, "xb")
            created.append(dest)
            with handle:
                handle.write(data)
    except BaseException as primary:
        # Only this call created target and these files; remove exactly them.
        # Cleanup problems are attached as notes and never mask the original.
        for path in created:
            try:
                path.unlink(missing_ok=True)
            except OSError as exc:
                primary.add_note(f"cleanup: could not remove {path}: {exc}")
        try:
            target.rmdir()
        except OSError as exc:
            primary.add_note(f"cleanup: could not remove {target}: {exc}")
        raise
    return created
