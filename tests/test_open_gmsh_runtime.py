"""Guards: the Radia panel Gmsh launcher accepts only a Gmsh 5 runtime.

Assumes: no native Gmsh is needed; fake modules stand in for the runtime, and
the launcher must return before any GUI call when the major version differs.
"""
import importlib.util
import sys
import types
from pathlib import Path

import pytest

LAUNCHER = Path(__file__).resolve().parents[1] / "src" / "radia" / "panels" / "open_gmsh.py"


def _load():
    spec = importlib.util.spec_from_file_location("open_gmsh_under_test", LAUNCHER)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _fake_gmsh(version):
    calls = []
    fake = types.SimpleNamespace(__version__=version, __file__="gmsh.py", libpath="gmsh.dll",
                                 initialize=lambda *a, **k: calls.append("initialize"),
                                 open=lambda *a: calls.append("open"))
    return fake, calls


def test_identity_reports_runtime_and_rejects_other_majors(monkeypatch):
    launcher = _load()
    monkeypatch.setenv("RADIA_GMSH_PROVENANCE", "record.json")
    good, _ = _fake_gmsh("5.0.0")
    identity = launcher._gmsh5_identity(good)
    assert identity["gmsh_version"] == "5.0.0", identity
    assert identity["provenance_record"] == "record.json", identity
    for version in ("4.15.2", "", "6.0.0"):
        with pytest.raises(RuntimeError, match="major 5"):
            launcher._gmsh5_identity(_fake_gmsh(version)[0])


def test_launcher_refuses_gmsh4_before_opening(monkeypatch, tmp_path):
    launcher = _load()
    msh = tmp_path / "case.msh"
    msh.write_text("placeholder", encoding="utf-8")
    fake, calls = _fake_gmsh("4.15.2")
    monkeypatch.setitem(sys.modules, "gmsh", fake)
    native = types.ModuleType("radia._gmsh_native_fields")
    native.validate_native_msh_fields = lambda path: None
    monkeypatch.setitem(sys.modules, "radia._gmsh_native_fields", native)
    status = launcher.main(["open_gmsh.py", str(msh)])
    assert (status, calls) == (6, []), f"status={status}, calls={calls}"


VIEWER = Path(__file__).resolve().parents[1] / "docs" / "stream_function" / "view_sf_coil_gx_gmsh.py"


def _viewer(monkeypatch, exe_file, returncode, stdout):
    spec = importlib.util.spec_from_file_location("view_sf_gmsh_under_test", VIEWER)
    viewer = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(viewer)
    seen = []

    def run(cmd, **kw):
        seen.append(cmd)
        return types.SimpleNamespace(returncode=returncode, stdout=stdout, stderr="")

    monkeypatch.setattr(viewer.subprocess, "run", run)
    monkeypatch.delenv("GMSH_EXECUTABLE", raising=False)
    viewer._GMSH_EXE = str(exe_file)
    return viewer, seen


def test_viewer_resolves_and_reports_the_executable_path(monkeypatch, tmp_path):
    """Guards: the viewer reports the resolved absolute executable it probed."""
    exe = tmp_path / "gmsh.exe"
    exe.write_bytes(b"")
    exe.chmod(0o755)
    viewer, seen = _viewer(monkeypatch, exe, 0, "5.0.0-git-798153893\n")
    viewer._GMSH_EXE = "gmsh"  # bare name: must be resolved through PATH
    monkeypatch.setenv("PATH", str(tmp_path))
    path, version = viewer._gmsh_executable()
    assert path == str(exe.resolve()), path
    assert version == "5.0.0-git-798153893", version
    assert seen == [[path, "-version"]], seen


@pytest.mark.parametrize("rc, out", [(1, "5.0.0\n"), (0, "4.15.2\n"), (0, "")])
def test_viewer_rejects_nonzero_exit_or_other_major(monkeypatch, tmp_path, rc, out):
    """Guards: a failing -version is refused even when its text says 5."""
    exe = tmp_path / "gmsh.exe"
    exe.write_bytes(b"")
    exe.chmod(0o755)
    viewer, _ = _viewer(monkeypatch, exe, rc, out)
    with pytest.raises(SystemExit) as err:
        viewer._gmsh_executable()
    want = "exited with" if rc else "major 5"
    assert want in str(err.value), str(err.value)


def test_viewer_explicit_missing_executable_fails(monkeypatch, tmp_path):
    """Guards: an explicit --gmsh-exe that does not exist is an error, not a skip."""
    viewer, seen = _viewer(monkeypatch, tmp_path / "absent.exe", 0, "5.0.0")
    with pytest.raises(SystemExit, match="not found"):
        viewer._gmsh_executable()
    assert seen == []
