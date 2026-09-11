"""Release-wheel guards for local native backup artifacts."""

from pathlib import Path
import tomllib


ROOT = Path(__file__).resolve().parents[1]


def test_coil_builder_dependency_is_available_without_extras():
    config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    assert "scipy>=1.11" in config["project"]["dependencies"]


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


def test_sdist_manifest_excludes_retired_and_backup_payloads():
    from setuptools._distutils.filelist import FileList

    kept = {
        str(Path("src/radia/_radia_pybind.pyd")),
        str(Path("src/radia/radia_motor_rom.dll")),
        str(Path("src/radia/panels/cubit_toolbar/icons/export.svg")),
    }
    excluded = {
        str(Path("src/radia/cubit_mesh_curver.pyd")),
        str(Path("src/radia/_radia_pybind.locked-old.pyd")),
        str(Path("src/radia/_radia_pybind.lockedold.pyd")),
        str(Path("src/radia/_radia_pybind.pre_release.pyd")),
        str(Path("src/radia/_radia_pybind.backup1.pyd")),
        str(Path("examples/retired.py")),
        str(Path("examples/cubit/retired.vol")),
    }
    files = FileList()
    files.set_allfiles(sorted(kept | excluded))
    # Also cover files found before MANIFEST.in processing (e.g. cached lists).
    files.files = sorted(kept | excluded)
    for line in (ROOT / "MANIFEST.in").read_text(encoding="utf-8").splitlines():
        if line.strip() and not line.lstrip().startswith("#"):
            files.process_template_line(line)
    assert kept <= set(files.files)
    assert not excluded.intersection(files.files)


def test_mkl_is_external_but_the_radia_owned_motor_abi_is_documented():
    config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    policies = (ROOT / "AGENTS.md").read_text(encoding="utf-8")

    assert "mkl>=2026,<2027" in config["project"]["dependencies"]
    assert "mkl>=2024.2.0" not in policies
    assert "radia_motor_rom.dll" in policies
    assert "no third-party DLLs" in policies


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


def test_wheel_does_not_package_generated_application_results():
    config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    package_data = set(config["tool"]["setuptools"]["package-data"]["radia"])

    assert "panels/samples/*.sol" not in package_data
    assert "panels/samples/*.json" not in package_data


def test_retired_desktop_viewer_is_not_packaged():
    config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    scripts = config["project"].get("scripts", {})
    gui_scripts = config["project"].get("gui-scripts", {})
    optional = config["project"]["optional-dependencies"]

    assert "radia-vol-viewer" not in scripts
    assert "radia-vol-viewer-gui" not in gui_scripts
    assert not (ROOT / "src" / "radia" / "tools" / "vol_sol_viewer.py").exists()
    assert all("pyvista" not in requirement.lower() for requirement in optional["viz"])
    assert all("pyvista" not in requirement.lower() for requirement in optional["dev"])


def test_mcp_sdk_is_optional_for_core_radia_users():
    config = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))
    project = config["project"]
    optional = project["optional-dependencies"]

    assert all(not requirement.lower().startswith("mcp") for requirement in project["dependencies"])
    assert optional["mcp"] == ["mcp>=1.0,<2"]
    assert "mcp>=1.0,<2" in optional["test"]
