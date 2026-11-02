"""Guards: a Ninja/MSVC build tree whose msvc_deps_prefix does not match the
compiler's /showIncludes bytes, or whose _radia_pybind objects carry no
recorded header dependencies, fails Build.ps1 before the provenance receipt.

Assumes: the byte strings are the observed cl.exe 19.44 Japanese prefix in
cp932 and UTF-8 and its UTF-8-read-as-cp932 mojibake; no compilation runs.
"""

from pathlib import Path
import importlib.util
import re

import pytest


ROOT = Path(__file__).resolve().parents[1]
_SPEC = importlib.util.spec_from_file_location("check_ninja_msvc_deps", ROOT / "tools" / "check_ninja_msvc_deps.py")
deps = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(deps)

NOTE = "\u30e1\u30e2: \u30a4\u30f3\u30af\u30eb\u30fc\u30c9 \u30d5\u30a1\u30a4\u30eb:  "
UTF8 = NOTE.encode("utf-8")
CP932 = NOTE.encode("cp932")
MOJIBAKE = UTF8.decode("cp932", "replace").encode("utf-8")


def _tree(tmp_path, prefix):
    rules = tmp_path / "CMakeFiles" / "rules.ninja"
    rules.parent.mkdir(parents=True)
    rules.write_bytes(b"ninja_required_version = 1.5\nmsvc_deps_prefix = " + prefix + b"\r\n\nrule CXX\n  deps = msvc\n")
    return tmp_path


@pytest.mark.parametrize("prefix", [UTF8, CP932])
def test_prefix_matching_compiler_bytes_passes(tmp_path, prefix):
    configured = deps.configured_prefix(_tree(tmp_path, prefix))
    assert configured == prefix, (configured, prefix)
    deps.check_prefix(configured, prefix + rb"C:\src\probe_header.h")


@pytest.mark.parametrize("configured,written", [(MOJIBAKE, UTF8), (MOJIBAKE, CP932), (UTF8, CP932), (CP932, UTF8)])
def test_prefix_in_another_encoding_fails(tmp_path, configured, written):
    with pytest.raises(deps.DepsError, match="does not match"):
        deps.check_prefix(deps.configured_prefix(_tree(tmp_path, configured)), written + rb"C:\src\probe_header.h")


@pytest.mark.parametrize("body", [b"rule CXX\n  deps = msvc\n", b"msvc_deps_prefix = \n",
                                  b"msvc_deps_prefix = A: \nmsvc_deps_prefix = B: \n"])
def test_missing_empty_or_ambiguous_prefix_fails(tmp_path, body):
    rules = tmp_path / "CMakeFiles" / "rules.ninja"
    rules.parent.mkdir(parents=True)
    rules.write_bytes(body)
    with pytest.raises(deps.DepsError, match="one non-empty msvc_deps_prefix"):
        deps.configured_prefix(tmp_path)


def _objects(tmp_path, names):
    for name in names:
        path = tmp_path / "CMakeFiles" / "_radia_pybind.dir" / "src" / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_bytes(b"")
    return [f"CMakeFiles/_radia_pybind.dir/src/{name}" for name in names]


def _deps_output(counts):
    return "".join(f"{path}: #deps {count}, deps mtime 8153543202164041 (VALID)\n    src/h.h\n\n"
                   for path, count in counts.items())


def test_target_with_recorded_deps_passes(tmp_path):
    objects = _objects(tmp_path, ["a.cpp.obj", "core/b.cpp.obj"])
    records = deps.recorded_deps(_deps_output(dict.fromkeys(objects, 3)))
    assert records == dict.fromkeys(objects, 3)
    assert deps.check_target_deps(tmp_path, "_radia_pybind", records) == 2


@pytest.mark.parametrize("record", ["zero", "missing"])
def test_target_object_without_deps_fails(tmp_path, record):
    objects = _objects(tmp_path, ["a.cpp.obj", "core/b.cpp.obj"])
    counts = {objects[0]: 3}
    if record == "zero":
        counts[objects[1]] = 0
    with pytest.raises(deps.DepsError, match=r"no header dependencies for 1/2 _radia_pybind objects"):
        deps.check_target_deps(tmp_path, "_radia_pybind", deps.recorded_deps(_deps_output(counts)))


def test_target_without_objects_fails(tmp_path):
    with pytest.raises(deps.DepsError, match="holds no compiled objects"):
        deps.check_target_deps(tmp_path, "_radia_pybind", {})


def test_build_checks_header_deps_around_radia_pybind():
    source = (ROOT / "Build.ps1").read_text(encoding="utf-8-sig")
    guard = r'\nif errorlevel 1 \(\n\s*echo ERROR: [^\n]*(?:configuration|build) failed\n\s*exit /b 1\n\)'
    order = re.search(
        r'-DCMAKE_BUILD_TYPE=Release\$PythonCMakeArgs\$NGSolveCMakeArgs\n'
        r'(?:.|\n)*?echo ERROR: CMake configuration failed(?:.|\n)*?\n\)\n'
        r'"\$PythonExecutable" "\$PROJECT_DIR\\tools\\check_ninja_msvc_deps\.py" prefix [^\n]*'
        r'--compiler "\$CompilerPath"' + guard +
        r'(?:.|\n)*?--target _radia_pybind -j \$BuildJobs\n(?:.|\n)*?\n\)\n'
        r'"\$PythonExecutable" "\$PROJECT_DIR\\tools\\check_ninja_msvc_deps\.py" deps [^\n]*'
        r'--ninja "\$NINJA_EXE" --target _radia_pybind' + guard, source)
    assert order, "Build.ps1 must check the deps prefix after configure and the recorded deps after _radia_pybind"
    batch_run = source.index('& cmd.exe /c "$BatchFile')
    alignment = source.index("[Console]::InputEncoding = [Console]::OutputEncoding")
    assert alignment < batch_run < source.index("Write-NativeBuildProvenance `"), (alignment, batch_run)
