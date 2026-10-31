"""Guards: tag-only approval and immutable-wheel publication cannot be bypassed.

Assumes: service metadata/ZIPs are controlled in memory; the real package wheel
verifier, GitHub environment approval and actual upload are separate checks.
No network, native runtime, package installation or publication is executed.
"""
import base64
import copy
import csv
from email.message import Message
import importlib.util
import io
import json
from pathlib import Path
import zipfile
from urllib import error

import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location("optuna_publication_contract", ROOT / "tools/verify_publication.py")
assert SPEC is not None and SPEC.loader is not None
gate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gate)
SOURCE = "1" * 40
ANNOTATED = "2" * 40
VERSION = "0.5.0"
FILENAME = "radia_optuna-0.5.0-cp312-cp312-win_amd64.whl"


def context():
    return {"repository": "owner/repo", "ref": "refs/tags/radia-optuna-v0.5.0",
            "ref_type": "tag", "sha": SOURCE}


def run():
    return {"id": 7, "status": "completed", "conclusion": "success", "head_sha": SOURCE,
            "head_branch": "main", "path": ".github/workflows/radia-optuna.yml",
            "event": "push", "head_repository": {"full_name": "owner/repo"}, "run_attempt": 1}


def jobs():
    return [{"id": i+10, "name": name, "status": "completed", "conclusion": "success"}
            for i, name in enumerate(gate.REQUIRED_JOBS)]


def wheel(source=SOURCE):
    files = {"radia_optuna/native_build.json": json.dumps({"source_sha": source,
             "build_route": "nvcc-shared-library", "wrapper": "plain-msvc-mex"}).encode(),
             "radia_optuna-0.5.0.dist-info/METADATA": b"Version: 0.5.0\n"}
    output = io.StringIO(); writer = csv.writer(output, lineterminator="\n")
    record = "radia_optuna-0.5.0.dist-info/RECORD"
    for name, body in files.items():
        value = base64.urlsafe_b64encode(bytes.fromhex(gate.digest(body))).decode().rstrip("=")
        writer.writerow([name, "sha256="+value, len(body)])
    writer.writerow([record, "", ""]); files[record] = output.getvalue().encode()
    data = io.BytesIO()
    with zipfile.ZipFile(data, "w") as archive:
        for name, body in files.items(): archive.writestr(name, body)
    return data.getvalue()


def artifact(data, name=FILENAME):
    output = io.BytesIO()
    with zipfile.ZipFile(output, "w") as archive: archive.writestr(name, data)
    return output.getvalue()


class Services:
    def __init__(self):
        self.run = run(); self.jobs = jobs(); self.wheel = wheel(); self.archive = artifact(self.wheel)
        self.artifacts = [{"id": 31, "name": "radia-optuna-wheel", "expired": False, "size_in_bytes": len(self.archive)}]
        self.reference = {"ref": context()["ref"], "object": {"type": "tag", "sha": ANNOTATED}}
        self.pypi = None; self.release = None; self.asset = self.wheel; self.calls = []

    def json(self, path, **kwargs):
        self.calls.append(("json", path))
        if "/git/ref/tags/" in path: return self.reference
        if "/git/tags/" in path: return {"sha": ANNOTATED, "object": {"type": "commit", "sha": SOURCE}}
        if path.endswith("/actions/runs/7"): return self.run
        if path.startswith("https://pypi.org/"): return self.pypi
        if "/releases/tags/" in path: return self.release
        raise AssertionError(path)

    def pages(self, path, key):
        self.calls.append(("pages", path))
        return self.jobs if key == "jobs" else self.artifacts

    def read(self, path, **kwargs):
        self.calls.append(("read", path))
        return self.asset if "/releases/assets/" in path else self.archive


@pytest.fixture
def package(monkeypatch, tmp_path):
    package = tmp_path / "package"; package.mkdir()
    (package / "pyproject.toml").write_text('[project]\nversion="0.5.0"\n', encoding="utf8")
    # This sentinel proves delegation/flags, not that synthetic bytes are real PE.
    (package / "verify_wheel.py").write_text(
        'def verify(wheel, *, release_candidate=False):\n'
        '    assert release_candidate\n'
        '    wheel.with_suffix(".checked").write_text("checked")\n'
        '    return {"ok": True}\n', encoding="utf8")
    monkeypatch.setattr(gate, "PACKAGE", package)
    return package


def test_tag_ci_artifact_and_approval_are_rechecked(package, tmp_path):
    service = Services()
    proof = gate.verify(service, context(), 7, gate.digest(service.wheel), tmp_path/"before")
    assert proof["source_sha"] == SOURCE and proof["artifact_id"] == 31
    assert (tmp_path/"before"/FILENAME).with_suffix(".checked").is_file()
    after = gate.verify(service, context(), 7, proof["wheel_sha256"], tmp_path/"after", proof)
    gate.compare_proofs(after, proof)
    assert not after["pypi_exists"] and not after["github_asset_exists"]


@pytest.mark.parametrize("field,value", [("ref_type", "branch"), ("ref", "refs/heads/main"),
    ("ref", "refs/tags/radia-optuna-v0.4.0"), ("sha", "--upload"),
    ("repository", "owner/repo;echo bad"), ("repository", "../repo")])
def test_wrong_ref_and_argument_injection_stop_before_service_reads(package, tmp_path, field, value):
    service = Services(); c = context(); c[field] = value
    with pytest.raises(gate.GateError): gate.verify(service, c, 7, gate.digest(service.wheel), tmp_path/"out")
    assert service.calls == []


def test_annotated_tags_must_peel_to_actual_triggering_commit(package, tmp_path):
    service = Services(); c = context(); c["sha"] = "3"*40
    with pytest.raises(gate.GateError, match="peeled tag commit"):
        gate.verify(service, c, 7, gate.digest(service.wheel), tmp_path/"out")
    with pytest.raises(gate.GateError, match="cyclic"):
        gate.peel(service.reference, lambda _: {"sha": ANNOTATED, "object": {"type": "tag", "sha": ANNOTATED}})


@pytest.mark.parametrize("field,value", [("head_sha", "3"*40), ("head_branch", "feature"),
    ("status", "in_progress"), ("conclusion", "failure"), ("event", "pull_request"),
    ("path", ".github/workflows/other.yml"), ("head_repository", {"full_name": "fork/repo"})])
def test_only_the_exact_successful_main_ci_is_accepted(field, value):
    candidate = run(); candidate[field] = value
    with pytest.raises(gate.GateError): gate.run_identity(candidate, jobs(), "owner/repo", SOURCE, 7)


@pytest.mark.parametrize("mode", ["missing", "duplicate", "skipped", "failed"])
def test_both_required_job_names_have_exactly_one_success(mode):
    values = jobs()
    if mode == "missing": values.pop()
    elif mode == "duplicate": values.append(copy.deepcopy(values[0]))
    else: values[0]["conclusion"] = "skipped" if mode == "skipped" else "failure"
    with pytest.raises(gate.GateError): gate.run_identity(run(), values, "owner/repo", SOURCE, 7)


@pytest.mark.parametrize("mode", ["missing", "duplicate", "expired", "oversized"])
def test_artifact_must_be_unique_live_and_bounded(mode):
    values = Services().artifacts
    if mode == "missing": values = []
    elif mode == "duplicate": values *= 2
    elif mode == "expired": values[0]["expired"] = True
    else: values[0]["size_in_bytes"] = gate.MAX_ARCHIVE+1
    with pytest.raises(gate.GateError): gate.choose_artifact(values)


@pytest.mark.parametrize("name", ["../"+FILENAME, "nested/"+FILENAME, "C:\\"+FILENAME,
    "radia_optuna-0.5.0-py3-none-win_amd64.whl", "radia_optuna-0.5.0-cp311-cp312-win_amd64.whl"])
def test_archive_paths_and_abi_cannot_be_substituted(name):
    with pytest.raises(gate.GateError): gate.unpack_wheel(artifact(wheel(), name), VERSION)


def test_duplicate_and_zip_bomb_members_fail():
    data = io.BytesIO()
    with zipfile.ZipFile(data, "w") as archive:
        archive.writestr(FILENAME, b"a")
        with pytest.warns(UserWarning, match="Duplicate name"):
            archive.writestr(FILENAME, b"b")
    with pytest.raises(gate.GateError, match="duplicate"): gate.unpack_wheel(data.getvalue(), VERSION)
    data = io.BytesIO()
    with zipfile.ZipFile(data, "w", compression=zipfile.ZIP_DEFLATED) as archive:
        archive.writestr(FILENAME, b"0"*100000)
    with pytest.raises(gate.GateError, match="size"): gate.unpack_wheel(data.getvalue(), VERSION)


def test_native_source_and_record_hash_are_checked():
    gate.wheel_source(wheel(), SOURCE)
    with pytest.raises(gate.GateError, match="source SHA"): gate.wheel_source(wheel("3"*40), SOURCE)
    original = zipfile.ZipFile(io.BytesIO(wheel())); output = io.BytesIO()
    with zipfile.ZipFile(output, "w") as archive:
        for name in original.namelist():
            body = original.read(name)
            if name.endswith("METADATA"): body += b"tampered\n"
            archive.writestr(name, body)
    with pytest.raises(gate.GateError, match="RECORD"): gate.wheel_source(output.getvalue(), SOURCE)


@pytest.mark.parametrize("field", ["artifact_id", "run_attempt", "source_sha", "wheel_sha256", "tag_object_sha"])
def test_post_approval_identity_changes_cannot_resume(package, tmp_path, field):
    service = Services(); before = gate.verify(service, context(), 7, gate.digest(service.wheel), tmp_path/"before")
    altered = copy.deepcopy(before); altered[field] = "changed"
    with pytest.raises(gate.GateError, match="Approved proof changed"):
        gate.verify(service, context(), 7, before["wheel_sha256"], tmp_path/"after", altered)
    assert not (tmp_path/"after").exists()


def test_same_hash_can_resume_but_existing_different_bytes_cannot(package, tmp_path):
    service = Services(); expected = gate.digest(service.wheel)
    service.pypi = {"info": {"version": VERSION}, "urls": [{"filename": FILENAME, "digests": {"sha256": expected}}]}
    service.release = {"tag_name": "radia-optuna-v0.5.0", "draft": False, "prerelease": False,
                       "assets": [{"id": 55, "name": FILENAME, "state": "uploaded", "size": len(service.wheel)}]}
    proof = gate.verify(service, context(), 7, expected, tmp_path/"same")
    assert proof["pypi_exists"] and proof["github_asset_exists"]
    service.asset = service.wheel[::-1]
    with pytest.raises(gate.GateError, match="overwrite prohibited"):
        gate.verify(service, context(), 7, expected, tmp_path/"wrong-asset")
    service.pypi["urls"][0]["digests"]["sha256"] = "0"*64
    with pytest.raises(gate.GateError, match="overwrite prohibited"):
        gate.verify(service, context(), 7, expected, tmp_path/"wrong-pypi")


class Response:
    def __init__(self, body): self.body = body; self.headers = {"Content-Length": str(len(body))}
    def read(self, limit): return self.body[:limit]
    def __enter__(self): return self
    def __exit__(self, *unused): return False


class Opener:
    def __init__(self, results): self.results = iter(results); self.requests = []
    def open(self, req, timeout):
        self.requests.append(req); value = next(self.results)
        if isinstance(value, Exception): raise value
        return value


def redirect(target):
    headers = Message(); headers["Location"] = target
    return error.HTTPError("https://api.github.com/repos/owner/repo/actions/artifacts/31/zip", 302, "redirect", headers, None)


def test_github_authorization_is_removed_on_signed_download_redirect():
    opener = Opener([redirect("https://production.blob.core.windows.net/path?signature=secret"), Response(b"zip")])
    assert gate.Remote("not-a-real-token", opener).read("https://api.github.com/repos/owner/repo/actions/artifacts/31/zip", binary=True) == b"zip"
    assert opener.requests[0].get_header("Authorization") == "Bearer not-a-real-token"
    assert opener.requests[1].get_header("Authorization") is None
    assert all(req.get_method() == "GET" for req in opener.requests)


@pytest.mark.parametrize("url", ["http://api.github.com/x", "https://api.github.com.evil.test/x",
    "https://token@api.github.com/x", "https://api.github.com:444/x", "file:///tmp/x"])
def test_api_host_and_credentials_are_validated_before_io(url):
    opener = Opener([])
    with pytest.raises(gate.GateError): gate.Remote("token", opener).read(url)
    assert opener.requests == []


def test_untrusted_redirect_and_unknown_service_status_fail_closed():
    opener = Opener([redirect("https://evil.test/steal")])
    with pytest.raises(gate.GateError, match="host"):
        gate.Remote("token", opener).read("https://api.github.com/download", binary=True)
    assert len(opener.requests) == 1
    failure = error.URLError("offline")
    with pytest.raises(gate.GateError, match="absence is not established"):
        gate.Remote("token", Opener([failure])).json("https://pypi.org/pypi/radia-optuna/0.5.0/json", public=True, absent=True)


def test_download_size_bound_is_checked_before_consuming_body(monkeypatch):
    monkeypatch.setattr(gate, "MAX_ARCHIVE", 2)
    with pytest.raises(gate.GateError, match="too large"):
        gate.Remote("token", Opener([Response(b"big")])).read("https://api.github.com/download", binary=True)
