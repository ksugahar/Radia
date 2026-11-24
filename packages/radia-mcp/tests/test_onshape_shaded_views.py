"""Guards: pinned GET-only shaded previews, PNG identity and batch cleanup."""

import asyncio
import base64
import hashlib
import io
import json
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

import pytest

from radia_mcp.onshape import server


DID, VID, EID = (letter * 24 for letter in "abc")
CONFIGURATION = "Size=large;Gap=2%20mm"
URL = f"https://cad.onshape.com/documents/{DID}/v/{VID}/e/{EID}?configuration=Size%3Dlarge%3BGap%3D2%2520mm"


@pytest.fixture(autouse=True)
def credentials(monkeypatch):
    monkeypatch.setenv("ONSHAPE_ACCESS_KEY", "test-access")
    monkeypatch.setenv("ONSHAPE_SECRET_KEY", "test-secret")


@pytest.fixture
def png():
    image = pytest.importorskip("PIL.Image")
    output = io.BytesIO()
    image.new("RGB", (16, 16), (20, 40, 70)).save(output, format="PNG")
    return output.getvalue()


@pytest.mark.parametrize("api_type,collection", [("PARTSTUDIO", "partstudios"), ("ASSEMBLY", "assemblies")])
def test_discovered_version_and_configuration_flow_to_shaded_png(
    tmp_path, monkeypatch, png, api_type, collection,
):
    calls = []

    class MockOpener:
        def open(self, request, timeout):
            assert request.get_method() == "GET"
            assert request.get_header("Authorization").startswith("Basic ")
            url = urlsplit(request.full_url)
            assert url.netloc == "cad.onshape.com"
            if url.path == f"/api/v11/documents/d/{DID}/v/{VID}/elements":
                body = [{"id": EID, "elementType": api_type}]
            else:
                assert url.path == f"/api/v11/{collection}/d/{DID}/v/{VID}/e/{EID}/shadedviews"
                options = parse_qs(url.query)
                assert options.pop("configuration") == [CONFIGURATION]
                assert options.pop("pixelSize") == ["0"]
                assert options.pop("outputWidth") == options.pop("outputHeight") == ["512"]
                assert options.pop("useAntiAliasing") == ["true"]
                [view] = options.pop("viewMatrix")
                assert not options
                calls.append(view)
                encoded = base64.b64encode(png).decode()
                if api_type == "ASSEMBLY":
                    encoded = [base64.b64encode(chunk).decode() for chunk in (png[:13], png[13:])]
                body = {"images": [encoded]}
            return io.BytesIO(json.dumps(body).encode())

    monkeypatch.setattr(server, "build_opener", lambda *args: MockOpener())
    selected = asyncio.run(server.onshape_discover("elements", url=URL))
    result = asyncio.run(server.onshape_render_shaded_views(
        **selected["elements"][0]["export_args"], export_dir=str(tmp_path)
    ))
    assert result["status"] == "ok", result
    assert result["version_id"] == VID and result["configuration"] == CONFIGURATION
    assert set(result["views"]) == {"isometric", "front"}
    [iso, front] = calls
    matrix = [float(value) for value in iso.split(",")]
    assert len(matrix) == 12 and matrix[3::4] == [0, 0, 0]
    directions = [matrix[i:i + 3] for i in (0, 4, 8)]
    for i, a in enumerate(directions):
        for j, b in enumerate(directions):
            assert sum(x * y for x, y in zip(a, b)) == pytest.approx(float(i == j))
    assert directions[2][0] > 0 and directions[2][1] < 0 and directions[2][2] > 0
    assert front == "front"
    for name, artifact in result["views"].items():
        assert Path(artifact["path"]).name == f"{name}.png"
        assert Path(artifact["path"]).read_bytes() == png
        assert artifact["sha256"] == hashlib.sha256(png).hexdigest()


@pytest.mark.parametrize("invalid", ["not-base64!", "ZW1wdHk=", "oversized", "truncated"])
def test_failed_batch_leaves_no_published_artifacts(tmp_path, monkeypatch, png, invalid):
    image = pytest.importorskip("PIL.Image")
    buffer = io.BytesIO()
    image.new("RGB", (513, 1)).save(buffer, format="PNG")
    invalid = {
        "oversized": base64.b64encode(buffer.getvalue()).decode(),
        "truncated": base64.b64encode(png[:-24]).decode(),
    }.get(invalid, invalid)
    requests = []

    class MockOpener:
        def open(self, request, timeout):
            requests.append(request)
            payload = base64.b64encode(png).decode() if len(requests) == 1 else invalid
            return io.BytesIO(json.dumps({"images": [payload]}).encode())

    monkeypatch.setattr(server, "build_opener", lambda *args: MockOpener())
    result = asyncio.run(server.onshape_render_shaded_views(
        DID, EID, VID, "assembly", str(tmp_path), views=["right", "top"]
    ))
    assert (result["status"], result["stage"]) == ("error", "shaded_views")
    assert [parse_qs(urlsplit(r.full_url).query)["viewMatrix"][0] for r in requests] == ["right", "top"]
    assert not list(tmp_path.iterdir())


def test_duplicate_view_names_fail_before_network(tmp_path, monkeypatch):
    monkeypatch.setattr(server, "build_opener", lambda *args: pytest.fail("duplicate views opened network"))
    result = asyncio.run(server.onshape_render_shaded_views(
        DID, EID, VID, "assembly", str(tmp_path), views=["front", "front"]
    ))
    assert result["stage"] == "input"
