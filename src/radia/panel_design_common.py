"""Small UI-neutral helpers shared by application DesignSpec contracts."""

from __future__ import annotations

from pathlib import Path


def panels_dir() -> Path:
    return Path(__file__).resolve().parent / "panels"


def calc_script(name: str, panel_dir: str | Path | None = None) -> str:
    base = Path(panel_dir) if panel_dir is not None else panels_dir()
    return str(base / name)


def _with_suffix(base_path: str | Path, suffix: str, ext: str) -> str:
    base = Path(base_path) if base_path else Path("output")
    return str(base.with_suffix("")) + suffix + ext


def msh_output(base_path: str | Path, suffix: str) -> str:
    return _with_suffix(base_path, suffix, ".msh")


def json_output(base_path: str | Path, suffix: str) -> str:
    return _with_suffix(base_path, suffix, ".json")


def append_value(cmd: list[str], flag: str, value: object) -> None:
    """Append ``flag value`` when *value* is meaningful."""

    if value is None:
        return
    text = str(value)
    if text == "":
        return
    cmd += [flag, text]


def append_switch(cmd: list[str], flag: str, enabled: bool) -> None:
    if enabled:
        cmd.append(flag)
