"""Release-wheel guards for local native backup artifacts."""

from pathlib import Path
import tomllib


ROOT = Path(__file__).resolve().parents[1]


def test_radia_wheel_excludes_native_backup_files():
    config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    excluded = set(
        config["tool"]["setuptools"]["exclude-package-data"]["radia"]
    )
    assert {
        "*.locked-old.pyd",
        "*.lockedold.pyd",
        "*.pre_*.pyd",
        "*.backup*.pyd",
    } <= excluded


def test_retired_notebook_workbenches_are_absent_from_source():
    retired = {
        "notebook_workbench.py",
        "ih_notebook.py",
        "motor_notebook.py",
        "pcb_notebook.py",
        "streamfunction_notebook.py",
    }
    present = sorted(
        name for name in retired if (ROOT / "src" / "radia" / name).exists()
    )
    assert present == []
