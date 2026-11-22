"""Guards: Onshape read-only discovery and pinned export handoff."""

import asyncio
import io
import json
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

import pytest

from radia_mcp.onshape import server


HOST = "cad.onshape.com"
DID, WID, VID, EID, MID = (digit * 24 for digit in "abcde")
BASE = f"https://cad.onshape.com/documents/{DID}"


@pytest.fixture
def discover(monkeypatch):
    monkeypatch.setenv("ONSHAPE_ACCESS_KEY", "test-access")
    monkeypatch.setenv("ONSHAPE_SECRET_KEY", "test-secret")
    return lambda operation, **kwargs: asyncio.run(server.onshape_discover(operation, **kwargs))


@pytest.mark.parametrize("kind,ref,immutable", [("w", WID, False), ("v", VID, True), ("m", MID, True)])
def test_parse_tab_url_without_credentials_or_network(monkeypatch, kind, ref, immutable):
    monkeypatch.delenv("ONSHAPE_ACCESS_KEY", raising=False)
    monkeypatch.setattr(server, "build_opener", lambda *args: pytest.fail("URL parse opened network"))
    result = asyncio.run(server.onshape_discover(
        "parse_url", url=f"{BASE}/{kind}/{ref}/e/{EID}?configuration=Size%3Dlarge%3BType%3DA"
    ))
    assert result == {"status": "ok", "reference": {
        "document_id": DID, "reference_kind": kind, "reference_id": ref, "element_id": EID,
        "immutable": immutable, "configuration": "Size=large;Type=A",
    }}


@pytest.mark.parametrize("url", [
    f"https://cad.onshape.com.evil.test/documents/{DID}/v/{VID}",
    f"https://user@{HOST}/documents/{DID}/v/{VID}",  # userinfo is rejected
    f"http://cad.onshape.com/documents/{DID}/v/{VID}",
    f"{BASE}/v/{VID}?configuration=A&configuration=B",
    f"{BASE}/v/{VID}/e/invalid",
])
def test_invalid_url_fails_before_network(discover, monkeypatch, url):
    monkeypatch.setattr(server, "build_opener", lambda *args: pytest.fail("invalid URL opened network"))
    assert discover("elements", url=url)["stage"] == "input"


def test_document_search_pages_without_following_server_next(discover, monkeypatch):
    calls = []

    class MockOpener:
        def open(self, request, timeout):
            assert request.get_method() == "GET"
            assert request.get_header("Authorization").startswith("Basic ")
            url = urlsplit(request.full_url)
            assert (url.netloc, url.path) == ("cad.onshape.com", "/api/v11/documents")
            query = parse_qs(url.query)
            assert query | {"offset": ["0"]} == {
                "filter": ["0"], "q": ["Teaching CAD & fixture"], "offset": ["0"], "limit": ["2"]
            }
            calls.append(query["offset"][0])
            response = {"items": [{"id": DID, "name": "Teaching CAD", "resourceType": "document"},
                                  {"id": WID, "name": "Folder", "resourceType": "folder"}],
                        "next": "https://not-onshape.example/next"}
            return io.BytesIO(json.dumps(response if len(calls) == 1 else {"items": []}).encode())

    monkeypatch.setattr(server, "build_opener", lambda *args: MockOpener())
    first = discover("documents", search="Teaching CAD & fixture", limit=2)
    assert first == {"status": "ok", "documents": [
        {"document_id": DID, "name": "Teaching CAD", "url": BASE}], "next_offset": 2}
    last = discover("documents", search="Teaching CAD & fixture", limit=2, offset=first["next_offset"])
    assert last == {"status": "ok", "documents": [], "next_offset": None}
    assert calls == ["0", "2"]


def test_discovery_to_fixed_step_export(discover, monkeypatch, tmp_path):
    cad = pytest.importorskip("build123d")
    fixture = tmp_path / "lesson.step"
    assert cad.export_step(cad.Box(10, 8, 6), fixture)
    calls = []
    translation, artifact = "1" * 24, "2" * 24
    workspace_path = f"/api/v11/documents/d/{DID}/w/{WID}/elements"
    version_path = f"/api/v11/documents/d/{DID}/v/{VID}/elements"
    versions_path = f"/api/v11/documents/d/{DID}/versions"
    export_path = f"/api/v11/partstudios/d/{DID}/v/{VID}/e/{EID}/export/step"
    download_path = f"/api/v11/documents/d/{DID}/externaldata/{artifact}"

    class MockOpener:
        def open(self, request, timeout):
            url = urlsplit(request.full_url)
            assert url.netloc == "cad.onshape.com"
            calls.append((request.get_method(), url.path))
            if url.path == versions_path:
                assert parse_qs(url.query) == {"offset": ["0"], "limit": ["20"]}
                body = [{"id": VID, "name": "Review", "microversion": MID}]
            elif url.path in {workspace_path, version_path}:
                body = [{"id": EID, "elementType": "PARTSTUDIO", "name": "Bracket", "microversionId": MID},
                        {"id": WID, "elementType": "ASSEMBLY", "name": "Fixture"},
                        {"id": MID, "elementType": "DRAWING", "name": "Sheet"}]
            elif url.path == export_path:
                options = json.loads(request.data)
                assert options["storeInDocument"] is False
                assert options["advancedParams"] == {"configuration": "Size=large"}
                body = {"id": translation, "requestState": "DONE", "resultExternalDataIds": [artifact]}
            elif url.path == download_path:
                return io.BytesIO(fixture.read_bytes())
            else:
                pytest.fail(f"unexpected mock request: {url.path}")
            return io.BytesIO(json.dumps(body).encode())

    monkeypatch.setattr(server, "build_opener", lambda *args: MockOpener())
    from_workspace = discover("elements", url=f"{BASE}/w/{WID}")
    assert len(from_workspace["elements"]) == 2
    assert all("export_args" not in row for row in from_workspace["elements"])
    versions = discover("versions", url=BASE)
    assert versions["versions"] == [{"version_id": VID, "name": "Review",
                                     "microversion_id": MID, "url": f"{BASE}/v/{VID}"}]
    version_url = versions["versions"][0]["url"] + f"/e/{EID}?configuration=Size%3Dlarge"
    fixed = discover("elements", url=version_url)
    [tab] = fixed["elements"]
    assert fixed["reference"]["element_id"] == tab["element_id"] == EID
    assert tab["microversion_id"] == MID
    output = asyncio.run(server.onshape_export_cad(
        **tab["export_args"], export_dir=str(tmp_path), file_format="step"
    ))
    assert output["status"] == "ok", output
    assert output["version_id"] == VID and output["configuration"] == "Size=large"
    assert Path(output["file"]).read_bytes() == fixture.read_bytes()
    assert calls == [("GET", workspace_path), ("GET", versions_path), ("GET", version_path),
                     ("POST", export_path), ("GET", download_path)]
