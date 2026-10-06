"""Admission failures must prevent a mismatched runtime/job from executing."""
import importlib.util
import json
from pathlib import Path
import pytest

path = Path(__file__).resolve().parents[1] / "tools/compute_runtime.py"
spec = importlib.util.spec_from_file_location("compute_runtime", path)
runtime = importlib.util.module_from_spec(spec)
spec.loader.exec_module(runtime)


def test_changed_lock_is_rejected(tmp_path):
    value = {"schema": runtime.SCHEMA, "python": "3.12.10"}
    import hashlib
    value["lock_sha256"] = hashlib.sha256(runtime.canonical(value)).hexdigest()
    p = tmp_path / "lock.json"
    p.write_text(json.dumps(value))
    assert runtime.read_lock(p) == value
    value["python"] = "3.12.11"
    p.write_text(json.dumps(value))
    with pytest.raises(ValueError, match="modified"):
        runtime.read_lock(p)


def test_mismatched_or_changed_job_is_rejected(tmp_path):
    p = tmp_path / "entry.py"
    p.write_text("raise SystemExit(0)")
    lock = {"lock_sha256": "a"}
    job = dict(schema="radia.compute-job.v1", lock_sha256="a", source_commit="commit",
               physics_contract={"units": "SI"}, files={"entry.py": runtime.digest(p)}, entry="entry.py")
    runtime.verify_job(tmp_path, job, lock)
    with pytest.raises(ValueError, match="lock mismatch"):
        runtime.verify_job(tmp_path, job, {"lock_sha256": "b"})
    p.write_text("raise SystemExit(1)")
    with pytest.raises(ValueError, match="changed"):
        runtime.verify_job(tmp_path, job, lock)


def test_bundle_cannot_escape_root(tmp_path):
    with pytest.raises(ValueError, match="Escaping"):
        runtime.safe_member(tmp_path, "../private.py")


def test_changed_wheel_bytes_are_rejected(tmp_path):
    p = tmp_path / "solver.whl"
    p.write_bytes(b"original wheel")
    lock = {"wheels": [{"file": p.name, "sha256": runtime.digest(p)}]}
    runtime.verify_wheels(lock, tmp_path)
    p.write_bytes(b"different wheel")
    with pytest.raises(ValueError, match="Wheel hash mismatch"):
        runtime.verify_wheels(lock, tmp_path)


def test_recovered_file_must_match_manifest(tmp_path):
    import zipfile
    archive = tmp_path / "recovery.zip"
    with zipfile.ZipFile(archive, "w") as z:
        z.writestr("result.json", "changed result")
        z.writestr("recovery_manifest.json", json.dumps(dict(root="owned", host="host", files={"result.json": "wrong hash"})))
    with pytest.raises(ValueError, match="Recovered file hash mismatch"):
        runtime.verify_recovery(archive, tmp_path / "recovered", runtime.digest(archive))


def test_environment_cannot_inherit_external_mkl_or_source_override(tmp_path, monkeypatch):
    monkeypatch.setenv("MKLROOT", "C:/external/mkl")
    monkeypatch.setenv("MKL_THREADING_LAYER", "INTEL")
    monkeypatch.setenv("PYTHONPATH", "C:/another/source")
    env = runtime.environment({"threads": runtime.THREADS}, tmp_path)
    assert Path(env["MKLROOT"]) == tmp_path.resolve() / "Library"
    assert env["MKL_THREADING_LAYER"] == "TBB"
    assert "PYTHONPATH" not in env


def test_memory_includes_windows_interpreter_child():
    from types import SimpleNamespace
    class Child:
        def memory_info(self):
            return SimpleNamespace(private=1000, rss=800)
    class Launcher(Child):
        def children(self, recursive):
            assert recursive
            return [Child(), Child()]
    assert runtime.process_tree_private_bytes(Launcher()) == 3000
