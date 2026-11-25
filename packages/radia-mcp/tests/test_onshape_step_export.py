"""Guards: immutable Onshape export, credential scope and local STEP handoff."""

import asyncio
import hashlib
import io
import json
from pathlib import Path
from types import SimpleNamespace
from urllib.error import HTTPError

import pytest

from radia_mcp.onshape import server


DOCUMENT, ELEMENT, VERSION = "a" * 24, "b" * 24, "c" * 24
TRANSLATION, RESULT = "d" * 24, "e" * 24


@pytest.mark.parametrize("element_type, collection", [
    ("part_studio", "partstudios"), ("assembly", "assemblies"),
])
def test_pinned_step_export_and_provenance(tmp_path, monkeypatch, element_type, collection):
    cad = pytest.importorskip("build123d")
    from radia_mcp.build123d import server as build123d

    reference = tmp_path / "reference.step"
    assert cad.export_step(cad.Box(10, 8, 6), reference)
    payload = reference.read_bytes()
    requests = []
    storage_url = "https://example-bucket.s3.amazonaws.com/signed-object?signature=dummy"
    prefix = "https://cad.onshape.com/api/v11"
    paths = [
        f"{prefix}/{collection}/d/{DOCUMENT}/v/{VERSION}/e/{ELEMENT}/export/step",
        f"{prefix}/translations/{TRANSLATION}",
        f"{prefix}/documents/d/{DOCUMENT}/externaldata/{RESULT}",
        storage_url,
    ]
    requested_unit = "MILLIMETER" if element_type == "assembly" else None

    class MockOpener:
        def open(self, request, timeout):
            index = len(requests)
            requests.append(request)
            assert request.full_url == paths[index]
            assert request.get_method() == ("POST" if index == 0 else "GET")
            assert bool(request.get_header("Authorization")) == (index != 3)
            if index == 0:
                expected = {
                    "storeInDocument": False, "grouping": True, "isYAxisUp": False,
                    "notifyUser": False, "triggerAutoDownload": False,
                    "advancedParams": {"configuration": "size%3Dlarge"},
                }
                if requested_unit is not None:
                    expected["stepUnit"] = requested_unit
                assert json.loads(request.data) == expected
                data = {"id": TRANSLATION, "requestState": "ACTIVE"}
            elif index == 1:
                data = {"id": TRANSLATION, "requestState": "DONE",
                        "resultExternalDataIds": [RESULT]}
            elif index == 2:
                raise HTTPError(request.full_url, 307, "Redirect",
                                {"Location": storage_url}, io.BytesIO())
            else:
                return io.BytesIO(payload)
            return io.BytesIO(json.dumps(data).encode())

    monkeypatch.setenv("ONSHAPE_ACCESS_KEY", "test-access")
    monkeypatch.setenv("ONSHAPE_SECRET_KEY", "test-secret")
    monkeypatch.setattr(server, "build_opener", lambda *args: MockOpener())
    monkeypatch.setattr(server.time, "sleep", lambda _: None)
    result = asyncio.run(server.onshape_export_cad(
        DOCUMENT, ELEMENT, VERSION, element_type, str(tmp_path), "step", "size%3Dlarge",
        step_unit=requested_unit,
    ))
    assert result["status"] == "ok", result
    assert len(requests) == 4
    provenance = Path(result["provenance"]).read_text(encoding="utf-8")
    record = json.loads(provenance)
    assert record["document_id"] == DOCUMENT and record["element_id"] == ELEMENT
    assert record["version_id"] == VERSION and record["configuration"] == "size%3Dlarge"
    assert record["step_sha256"] == hashlib.sha256(payload).hexdigest()
    assert record["requested_step_unit"] == requested_unit
    assert record["declared_length_units"] == ["millimetre"]
    assert record["measurement_length_unit"] == "mm"
    assert "test-secret" not in provenance and "signature=dummy" not in provenance
    assert Path(result["step"]).read_bytes() == payload

    inspection = json.loads(build123d._inspect_step_readback(result["step"]))
    assert inspection["solids"][0]["volume_mm3"] == pytest.approx(480)
    diff = json.loads(build123d._diff_step_files(str(reference), result["step"]))
    assert diff["assembly"]["symmetric_difference_volume_mm3"] == pytest.approx(0)
    render = json.loads(build123d._render_step_views(
        result["step"], str(tmp_path), ["svg"], "xy", None
    ))
    assert render["status"] == "ok", render
    assert len(render["views"]) == 5


def test_untrusted_download_redirect_is_rejected(tmp_path):
    class RedirectOpener:
        def open(self, request, timeout):
            assert request.full_url == "https://cad.onshape.com/export"
            raise HTTPError(request.full_url, 307, "Redirect",
                            {"Location": "https://other.example/export"}, io.BytesIO())

    target = tmp_path / "model.part"
    with pytest.raises(ValueError, match="redirect destination"):
        server._download(RedirectOpener(), "https://cad.onshape.com/export", "Basic dummy", target)
    assert not target.exists()


def test_invalid_version_returns_before_auth(tmp_path, monkeypatch):
    monkeypatch.delenv("ONSHAPE_ACCESS_KEY", raising=False)
    monkeypatch.setattr(server, "_authorization", lambda: pytest.fail("unexpected authentication"))
    result = asyncio.run(server.onshape_export_cad(
        DOCUMENT, ELEMENT, "workspace", "assembly", str(tmp_path), "step"
    ))
    assert result["stage"] == "input"


@pytest.mark.parametrize("element_type,collection", [
    ("part_studio", "partstudios"), ("assembly", "assemblies"),
])
@pytest.mark.parametrize("suffix,api_format,payload", [
    ("x_t", "PARASOLID", b"**PARASOLID !\n**TEST HEADER\n"),
    ("sat", "ACIS", b"500 0 1 0\nTest header\n1 0 0\nEnd-of-ACIS-data\n"),
])
def test_native_export_pins_version_and_configuration(
    tmp_path, monkeypatch, element_type, collection, suffix, api_format, payload,
):
    configuration = "Size=large;Style=round"
    requests = []

    class MockOpener:
        def open(self, request, timeout):
            requests.append((request.get_method(), request.full_url))
            if request.get_method() == "POST":
                options = json.loads(request.data)
                assert options["formatName"] == api_format
                assert options["configuration"] == configuration
                assert options["storeInDocument"] is False
                assert options["grouping"] is True
                assert options.get("parasolidMode") == ("TEXT" if suffix == "x_t" else None)
                return io.BytesIO(json.dumps({
                    "id": TRANSLATION, "requestState": "DONE", "resultExternalDataIds": [RESULT]
                }).encode())
            return io.BytesIO(payload)

    monkeypatch.setenv("ONSHAPE_ACCESS_KEY", "test-access")
    monkeypatch.setenv("ONSHAPE_SECRET_KEY", "test-secret")
    monkeypatch.setattr(server, "build_opener", lambda *args: MockOpener())
    result = asyncio.run(server.onshape_export_cad(
        DOCUMENT, ELEMENT, VERSION, element_type, str(tmp_path), suffix, configuration
    ))
    assert result["status"] == "ok", result
    base = "https://cad.onshape.com/api/v11"
    assert requests == [
        ("POST", f"{base}/{collection}/d/{DOCUMENT}/v/{VERSION}/e/{ELEMENT}/translations"),
        ("GET", f"{base}/documents/d/{DOCUMENT}/externaldata/{RESULT}"),
    ]
    artifact = Path(result["file"])
    assert artifact.suffix == f".{suffix}" and artifact.read_bytes() == payload
    provenance = json.loads(Path(result["provenance"]).read_bytes())
    assert {key: provenance[key] for key in ("document_id", "element_id", "version_id",
                                             "configuration", "format", "file", "file_sha256")} == {
        "document_id": DOCUMENT, "element_id": ELEMENT, "version_id": VERSION,
        "configuration": configuration, "format": suffix, "file": artifact.name,
        "file_sha256": hashlib.sha256(payload).hexdigest(),
    }
    assert "step" not in result and "step_sha256" not in result


@pytest.mark.parametrize("suffix", ["x_t", "sat"])
def test_native_export_rejects_zip_payload(suffix):
    payload = b"PK\x03\x04 archive"
    with pytest.raises(ValueError, match=f"expected .{suffix}"):
        server._validate_download(payload, payload, suffix)


@pytest.fixture
def run_export(tmp_path, monkeypatch):
    monkeypatch.setenv("ONSHAPE_ACCESS_KEY", "test-access")
    monkeypatch.setenv("ONSHAPE_SECRET_KEY", "test-secret")

    def run(opener, **options):
        monkeypatch.setattr(server, "build_opener", lambda *args: opener)
        return asyncio.run(server.onshape_export_cad(
            DOCUMENT, ELEMENT, VERSION, "part_studio", str(tmp_path), "step", **options
        ))

    return run


@pytest.mark.parametrize("content_length,payload,limit,message", [
    ("100", b"", 16, "16-byte limit"),
    (None, b"a" * 17, 16, "16-byte limit"),
    ("8", b"a" * 17, 16, "16-byte limit"),
    (None, b"<html>not a STEP</html>", 64, "expected .step"),
    (None, b"ISO-10303-21;\nEND-ISO-10303-21;\n", 64, "read STEP length units"),
])
def test_failed_download_or_provenance_cleans_owned_staging(
    tmp_path, run_export, content_length, payload, limit, message,
):
    existing = tmp_path / "onshape-step-existing"
    existing.mkdir()
    (existing / "model.step").write_bytes(b"keep")

    class MockOpener:
        def open(self, request, timeout):
            if request.get_method() == "POST":
                return io.BytesIO(json.dumps({
                    "id": TRANSLATION, "requestState": "DONE", "resultExternalDataIds": [RESULT]
                }).encode())
            response = io.BytesIO(payload)
            response.headers = {} if content_length is None else {"Content-Length": content_length}
            return response

    if message == "read STEP length units":
        pytest.importorskip("OCP")
    result = run_export(MockOpener(), max_download_bytes=limit)
    assert result["status"] == "error" and message in result["error"], result
    assert set(tmp_path.iterdir()) == {existing}
    assert (existing / "model.step").read_bytes() == b"keep"


def test_poll_timeout_keeps_id_and_never_resubmits(run_export, monkeypatch):
    now = 0
    requests = []

    def sleep(seconds):
        nonlocal now
        now += seconds

    monkeypatch.setattr(server, "time", SimpleNamespace(monotonic=lambda: now, sleep=sleep))

    class MockOpener:
        def open(self, request, timeout):
            requests.append((request.get_method(), timeout))
            return io.BytesIO(json.dumps({"id": TRANSLATION, "requestState": "ACTIVE"}).encode())

    result = run_export(MockOpener(), poll_timeout_s=3)
    assert result == {"status": "error", "stage": "poll", "error": "Onshape polling deadline exceeded.",
                      "translation_id": TRANSLATION}
    assert requests == [("POST", 30), ("GET", 3), ("GET", 1)]


def test_http_diagnostics_never_echo_upstream_body_or_url(run_export):
    class MockOpener:
        def open(self, request, timeout):
            raise HTTPError("https://cad.onshape.com/?signature=secret",
                            403, "test-secret", {},
                            io.BytesIO(b'{"message":"test-access:test-secret"}'))

    result = run_export(MockOpener())
    assert result == {"status": "error", "stage": "export", "http_status": 403,
                      "error": "Check export permissions for the fixed version."}


def test_download_forbidden_reports_download_context(run_export):
    class MockOpener:
        def open(self, request, timeout):
            if request.get_method() == "POST":
                return io.BytesIO(json.dumps({
                    "id": TRANSLATION, "requestState": "DONE", "resultExternalDataIds": [RESULT]
                }).encode())
            raise HTTPError(request.full_url, 403, "Forbidden", {},
                            io.BytesIO(b'{"message":"test-secret"}'))

    result = run_export(MockOpener())
    assert (result["stage"], result["http_status"]) == ("download", 403)
    assert "signed download link may have expired" in result["error"]
    assert "test-secret" not in result["error"]


def test_export_tool_declares_positive_integer_limits():
    [tool] = [t for t in asyncio.run(server.mcp.list_tools()) if t.name == "onshape_export_cad"]
    for name in ("max_download_bytes", "poll_timeout_s"):
        argument = tool.inputSchema["properties"][name]
        assert argument["type"] == "integer"
        assert argument["exclusiveMinimum"] == 0


def test_fixed_step_cache_reuses_verified_bytes_and_isolates_options(tmp_path, run_export, monkeypatch):
    cad = pytest.importorskip("build123d")
    reference = tmp_path / "reference.step"
    assert cad.export_step(cad.Box(10, 8, 6), reference)
    payload = reference.read_bytes()
    cache_root = tmp_path / "cache"
    monkeypatch.setenv("ONSHAPE_STEP_CACHE_DIR", str(cache_root))
    requests = []

    class MockOpener:
        def open(self, request, timeout):
            requests.append(request)
            if request.get_method() == "POST":
                return io.BytesIO(json.dumps({
                    "id": TRANSLATION, "requestState": "DONE",
                    "resultExternalDataIds": [RESULT],
                }).encode())
            return io.BytesIO(payload)

    first = run_export(MockOpener(), step_unit="MILLIMETER", configuration="size=large")
    changed = run_export(MockOpener(), step_unit="INCH", configuration="size=large")
    different_config = run_export(MockOpener(), step_unit="INCH", configuration="size=small")
    assert [r.get("cache") for r in (first, changed, different_config)] == ["miss"] * 3
    assert len(requests) == 6
    assert len(list(cache_root.iterdir())) == 3
    assert str(cache_root) not in json.dumps(first)

    monkeypatch.delenv("ONSHAPE_ACCESS_KEY")
    monkeypatch.delenv("ONSHAPE_SECRET_KEY")
    hit = run_export(MockOpener(), step_unit="MILLIMETER", configuration="size=large")
    assert hit["cache"] == "hit" and hit["status"] == "ok", hit
    assert hit["file"] != first["file"]
    assert Path(hit["file"]).read_bytes() == payload
    assert hit["file_sha256"] == first["file_sha256"]
    assert len(requests) == 6

    (cache_root / first["cache_key"] / "model.step").write_bytes(b"corrupted")
    broken = run_export(MockOpener(), step_unit="MILLIMETER", configuration="size=large")
    assert (broken["status"], broken["stage"], broken["cache"]) == ("error", "cache", "error")
    assert len(requests) == 6


def test_workspace_export_is_never_cached(tmp_path, monkeypatch):
    cache_root = tmp_path / "cache"
    monkeypatch.setenv("ONSHAPE_STEP_CACHE_DIR", str(cache_root))
    monkeypatch.setattr(server, "_authorization", lambda: pytest.fail("unexpected authentication"))
    result = asyncio.run(server.onshape_export_cad(
        DOCUMENT, ELEMENT, "workspace", "part_studio", str(tmp_path), "step",
    ))
    assert (result["status"], result["stage"], result["cache"]) == ("error", "input", "not_eligible")
    assert not cache_root.exists()
