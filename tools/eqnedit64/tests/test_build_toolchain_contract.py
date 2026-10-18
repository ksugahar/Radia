"""Guard all Eqnedit64 native entry points against an older toolchain bypass."""
from pathlib import Path
import re

ROOT = Path(__file__).resolve().parents[3]
BUILD = ROOT / "tools/eqnedit64/build"


def test_every_batch_build_uses_the_same_mandatory_selector():
    for name in ("build_eqnedit64.bat", "build_pymodule.bat", "build_asan.bat"):
        text = (BUILD / name).read_text()
        assert 'call "%~dp0select_vs2022.bat"\nif errorlevel 1 exit /b 1' in text
        assert "goto have_cl" not in text
        assert "vswhere" not in text.lower()


def test_selector_rejects_inherited_or_older_compilers():
    text = (BUILD / "select_vs2022.bat").read_text()
    assert "-version '[17.0,18.0)'" in text
    assert 'if not "%VisualStudioVersion%"=="17.0"' in text
    assert 'if /I not "%EQN_FIRST_CL%"=="%VCToolsInstallDir%bin\\Hostx64\\x64\\cl.exe"' in text
    assert "-prerelease" not in text
    assert "-vcvars_ver=14.44.35207" in text
    assert '"%VCToolsVersion%"=="14.44.35207"' in text
    assert 'assert_msvc1944.ps1" -Compiler "%EQN_FIRST_CL%"' in text
    assert "vcvars64.bat is missing" in text


def test_each_ci_cmake_configure_selects_pinned_vs2022_explicitly():
    workflow = (ROOT / ".github/workflows/eqnedit64.yml").read_text()
    commands = re.findall(r"^.*cmake -S .*", workflow, re.M)
    assert commands
    assert all('-G "Visual Studio 17 2022"' in line and '-T "v143,version=14.44.35207"' in line for line in commands)



def test_every_cmake_configure_checks_actual_compiler():
    workflow = (ROOT / ".github/workflows/eqnedit64.yml").read_text()
    configurations = workflow.count("cmake -S . -B build_cmake")
    assert configurations == 4
    assert workflow.count("assert_msvc1944.ps1 -BuildDirectory build_cmake") == configurations
