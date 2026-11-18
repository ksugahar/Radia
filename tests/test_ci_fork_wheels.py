import hashlib
import importlib.util
import json
from pathlib import Path
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("ci_fork_wheels", ROOT / "tools/ci_fork_wheels.py")
fork = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fork)


def test_modified_wheel_cannot_reach_pip(tmp_path, monkeypatch):
    filename = "netgen_mesher-test.whl"
    (tmp_path / filename).write_bytes(b"modified")
    lock = {"wheels": {"netgen-mesher": {
        "filename": filename, "sha256": hashlib.sha256(b"verified").hexdigest()
    }}}
    monkeypatch.setattr(sys, "prefix", str(tmp_path))
    monkeypatch.setattr(fork.subprocess, "run", lambda *a, **k: pytest.fail("pip received an unverified wheel"))
    with pytest.raises(RuntimeError, match="Incorrect wheel SHA-256"):
        fork.install(lock, tmp_path)


def test_lock_identifies_both_cpu_wheels_and_native_dlls():
    from packaging.specifiers import SpecifierSet
    import tomllib

    lock = json.loads(fork.LOCK.read_text(encoding="utf-8"))
    deps = tomllib.loads((ROOT / "pyproject.toml").read_text(encoding="utf-8"))["project"]["dependencies"]
    assert set(lock["wheels"]) == {"netgen-mesher", "ngsolve"}
    assert set(lock["dlls"]) == {"nglib.dll", "ngcore.dll", "libngsolve.dll"}
    for package, wheel in lock["wheels"].items():
        pin, = [dep[len(package):] for dep in deps if dep.startswith(package + "==")]
        assert SpecifierSet(pin).contains(wheel["version"])
        expected = lock["sources"]["netgen" if package == "netgen-mesher" else "ngsolve"]
        assert "+radia." + expected[:8] in wheel["filename"]
    assert lock["wheels"]["ngsolve"]["version"].endswith(".cpu")
