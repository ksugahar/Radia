"""Exercise the actual generated batch guards without compiling native code."""

from pathlib import Path
import re
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
TARGETS = (
    "cln_core", "sparsesolv_ngsolve", "axifem", "_equation",
    "cubit_mesh_curver", "cubit_mesh_export_ccm",
)


@pytest.mark.parametrize("target", TARGETS)
@pytest.mark.parametrize("exit_code", [0, 7])
def test_native_build_guard_propagates_failure(tmp_path, target, exit_code):
    source = (ROOT / "Build.ps1").read_text(encoding="utf-8-sig")
    match = re.search(
        rf'"\$CMAKE_EXE" --build [^\n]*--target {target} -j\n'
        r'(?P<guard>\s*if errorlevel 1 \([^)]*\))', source,
    )
    assert match, f"missing immediate failure guard for {target}"
    guard = match.group("guard")
    assert f"ERROR: {target} build failed" in guard
    assert "exit /b 1" in guard
    if sys.platform != "win32":
        pytest.skip("cmd.exe execution requires Windows")
    batch = tmp_path / "guard.bat"
    batch.write_text(
        f'@echo off\n"{sys.executable}" -c "import sys; sys.exit({exit_code})"\n'
        f'{guard}\necho CONTINUED\nexit /b 0\n', encoding="ascii",
    )
    result = subprocess.run(
        ["cmd.exe", "/d", "/c", str(batch)], capture_output=True,
        text=True, timeout=15,
    )
    assert result.returncode == (1 if exit_code else 0)
    assert ("CONTINUED" in result.stdout) == (exit_code == 0)


def test_plugin_configuration_failure_cannot_use_stale_build_tree():
    source = (ROOT / "Build.ps1").read_text(encoding="utf-8-sig")
    commands = re.findall(
        r'"\$CMAKE_EXE" -G Ninja [^\n]*"%CUBIT_PLUGIN_SRC%"\n'
        r'(?P<guard>\s*if errorlevel 1 \([^)]*\))', source,
    )
    assert len(commands) == 2
    assert all("configuration failed" in guard and "exit /b 1" in guard
               for guard in commands)
