"""Execute the CPU discovery block with isolated installation inventories."""
from pathlib import Path
import shutil
import subprocess

import pytest


@pytest.mark.parametrize(
    "available,toolset,compiler,accepted",
    [(True, "14.44.35207", "19.44.35229.0", True),
     (False, "14.44.35207", "19.44.35229.0", False),
     (True, "14.51.10000", "19.51.10000.0", False),
     (True, "14.44.35207", "19.51.10000.0", False)],
)
def test_cpu_build_rejects_unsupported_installations(tmp_path, available, toolset, compiler, accepted):
    pwsh = shutil.which("pwsh")
    if pwsh is None:
        pytest.skip("Windows PowerShell 7 discovery contract")
    root = Path(__file__).resolve().parents[1]
    text = (root / "Build.ps1").read_text(encoding="utf-8")
    start = text.index("# CPU toolchain follows")
    block = text[start:text.index("if (-not $CMAKE_EXE)", start)]
    install = tmp_path / "VS2022"
    identity = install / "VC/Auxiliary/Build/Microsoft.VCToolsVersion.default.txt"
    identity.parent.mkdir(parents=True)
    identity.write_text(toolset, encoding="ascii")
    binary = install / f"VC/Tools/MSVC/{toolset}/bin/Hostx64/x64/cl.exe"
    binary.parent.mkdir(parents=True)
    binary.touch()
    probe = tmp_path / "inventory.ps1"
    probe.write_text(
        "if ($args -notcontains '-version' -or $args -notcontains '[17.0,18.0)') { throw 'Unrestricted discovery' }\n"
        + (f"'{install.as_posix()}'\n" if available else "# VS2026 is the only installed product\n"),
        encoding="utf-8",
    )
    block = block.replace(
        '$VSWHERE = "${env:ProgramFiles(x86)}\\Microsoft Visual Studio\\Installer\\vswhere.exe"',
        f"$VSWHERE = '{probe.as_posix()}'",
    )
    script = tmp_path / "check.ps1"
    script.write_text(
        "$ErrorActionPreference = 'Stop'\n"
        # File-version querying is the sole mock: no fixture compiler is executed.
        + f"function Get-Item {{ param($LiteralPath) [pscustomobject]@{{VersionInfo=[pscustomobject]@{{FileVersion='{compiler}'}}}} }}\n"
        + block,
        encoding="utf-8",
    )
    result = subprocess.run([pwsh, "-NoProfile", "-File", str(script)], capture_output=True, text=True)
    assert (result.returncode == 0) is accepted, result.stdout + result.stderr
