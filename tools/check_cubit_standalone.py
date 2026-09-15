"""Check an installed wheel without Radia or changes to a real Cubit profile."""
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys
import tarfile
import tempfile
from unittest.mock import patch


def main():
    if importlib.util.find_spec("radia") is not None:
        raise RuntimeError("Run this gate in an isolated environment without radia")
    import cubit_mesh_export
    from cubit_mesh_export import toolbar_install
    from cubit_mesh_export.smoke_test import _find_sample_jou
    installed = Path(cubit_mesh_export.__file__).resolve().parent
    if not installed.is_relative_to(Path(sys.prefix).resolve()):
        raise RuntimeError(f"Not the installed candidate wheel: {installed}")
    fixture = _find_sample_jou()
    assert fixture.is_file() and fixture.is_relative_to(installed)
    gui = installed / "cubit_gui"
    assert (gui / "register_toolbar.py").is_file()
    temp_root = r"C:\temp" if os.name == "nt" else None
    with tempfile.TemporaryDirectory(prefix="cubit-wheel-check-", dir=temp_root) as directory:
        root = Path(directory)
        home = root / "home"
        home.mkdir()
        cubit = root / "programs/Coreform Cubit 2025.12/bin"
        (cubit / "plugins").mkdir(parents=True)
        (cubit / "cubit.py").write_text("# test profile only\n")
        environment = {
            "HOME": str(home), "USERPROFILE": str(home),
            "LOCALAPPDATA": str(root / "local"), "APPDATA": str(root / "roaming"),
            "ProgramData": str(root / "shared"), "ProgramFiles": str(root / "programs"),
            "ProgramFiles(x86)": "", "CUBIT_PATH": "", "SystemDrive": str(root / "drive"),
        }
        with patch.dict(os.environ, environment):
            assert toolbar_install.install_panels(all_users=False)
            valid, issues = toolbar_install.verify_panel_installation(all_users=False)
            assert valid, issues
        archive_path = root / "local/Radia/Cubit/radia_export_toolbar.tar.gz"
        with tarfile.open(archive_path) as archive:
            assert archive.extractfile("scripts/radia_export_menu.py").read() == (gui / "radia_export_menu.py").read_bytes()
            assert len([n for n in archive.getnames() if n.startswith("scripts/export_")]) == 6
        for module in ("install", "smoke_test", "check"):
            result = subprocess.run([sys.executable, "-m", f"cubit_mesh_export.{module}", "--help"],
                                    capture_output=True, text=True, timeout=30, cwd=root)
            assert result.returncode == 0, (module, result.stderr)
    assert importlib.util.find_spec("PySide6") is None, "Qt must remain Cubit-owned"
    print(json.dumps({"passed": True, "installed": str(installed),
                      "fixture": str(fixture), "radia_required": False,
                      "real_profile_modified": False}))


if __name__ == "__main__":
    main()
