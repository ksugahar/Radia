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


def test_ci_native_artifact_excludes_backups_and_pre_push_does_not_upload():
    hook = (ROOT / "tools" / "git-hooks" / "pre-push").read_text(
        encoding="utf-8"
    )
    workflow = (ROOT / ".github" / "workflows" / "build-test.yml").read_text(
        encoding="utf-8"
    )
    assert "upload_release_asset.py" not in hook
    assert "!src/radia/*.locked-old.*" in workflow
    assert "!src/radia/*.lockedold.*" in workflow
    assert "!src/radia/*.pre_*.*" in workflow
    assert "!src/radia/*.backup*.*" in workflow


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
