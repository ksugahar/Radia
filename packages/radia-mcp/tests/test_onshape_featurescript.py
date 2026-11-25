"""Guards: fixed-version FeatureScript evaluation accepts only a face query."""

import asyncio
import hashlib
import io
import json
from urllib.parse import parse_qs, urlsplit

import pytest

from radia_mcp.onshape import server


DID, VID, EID = "a" * 24, "b" * 24, "c" * 24


@pytest.mark.parametrize("face_id,version", [
    ('opDeleteBodies(context, id, {});', VID),
    ("JHD", "workspace"),
])
def test_rejects_non_query_input_before_auth(monkeypatch, face_id, version):
    monkeypatch.setattr(server, "_authorization", lambda: pytest.fail("unexpected authentication"))
    result = asyncio.run(server.onshape_measure_planar_face(DID, EID, version, face_id))
    assert (result["status"], result["stage"]) == ("error", "input")


def test_planar_face_measurement_posts_only_a_fixed_query(monkeypatch):
    monkeypatch.setenv("ONSHAPE_ACCESS_KEY", "dummy")
    monkeypatch.setenv("ONSHAPE_SECRET_KEY", "dummy")
    sent_scripts = []

    class MockOpener:
        def open(self, request, timeout):
            url = urlsplit(request.full_url)
            assert url.netloc == "cad.onshape.com"
            assert url.path == f"/api/v11/partstudios/d/{DID}/v/{VID}/e/{EID}/featurescript"
            assert parse_qs(url.query) == {"rollbackBarIndex": ["-1"], "configuration": ["gap=2 mm"]}
            assert request.get_method() == "POST"
            body = json.loads(request.data)
            assert body["queries"] == {"face": ["JHD"]}
            assert "evArea(" in body["script"] and "evApproximateCentroid(" in body["script"]
            assert "opDeleteBodies" not in body["script"]
            sent_scripts.append(body["script"])
            return io.BytesIO(json.dumps({
                "result": {"value": [{"value": n} for n in (0.00012, 0.002, -0.003, 0.0)]},
            }).encode())

    monkeypatch.setattr(server, "build_opener", lambda *args: MockOpener())
    result = asyncio.run(server.onshape_measure_planar_face(DID, EID, VID, "JHD", "gap=2 mm"))
    assert result["status"] == "ok", result
    assert result["version_id"] == VID and result["configuration"] == "gap=2 mm"
    assert result["coordinate_frame"] == "part_studio"
    assert result["area_m2"] == pytest.approx(0.00012)
    assert result["approximate_centroid_m"] == pytest.approx([0.002, -0.003, 0.0])
    assert len(sent_scripts) == 1
    assert result["script_sha256"] == hashlib.sha256(sent_scripts[0].encode()).hexdigest()
