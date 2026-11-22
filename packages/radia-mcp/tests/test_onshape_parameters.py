"""Guards: read-only fixed-version parameter discovery and encoded STEP handoff."""

import asyncio
import io
import json
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

import pytest

from radia_mcp.onshape import server


DID, VID, EID, MID = (char * 24 for char in "abcd")
REFERENCE = f"/d/{DID}/v/{VID}/e/{EID}"
ENCODED = "Size=large;Gap=2%20mm"


@pytest.fixture(autouse=True)
def credentials(monkeypatch):
    monkeypatch.setenv("ONSHAPE_ACCESS_KEY", "test-access")
    monkeypatch.setenv("ONSHAPE_SECRET_KEY", "test-secret")


def test_discovery_parameters_and_export_use_one_pinned_version(monkeypatch, tmp_path):
    cad = pytest.importorskip("build123d")
    sample = tmp_path / "lesson.step"
    assert cad.export_step(cad.Box(10, 8, 6), sample)
    calls = []
    translation, artifact = "1" * 24, "2" * 24
    definitions = [{"parameterId": "Size", "options": [{"option": "large"}]},
                   {"parameterId": "Gap", "defaultValue": "2 mm"}]
    table = {"variables": [
        {"name": "gap", "expression": "#plate / 2", "value": "2 mm", "type": "LENGTH"},
        {"name": "angle", "expression": "60 deg", "value": "60 deg", "type": "ANGLE"},
        {"name": "count", "expression": "2 + 1", "value": "3", "type": "NUMBER"},
        {"name": "shape", "expression": "vector(1,2)", "value": "[1, 2]", "type": "ANY"},
        {"name": "waiting", "expression": "#missing", "type": "ANY"},
    ]}

    class MockOpener:
        def open(self, request, timeout):
            url = urlsplit(request.full_url)
            method = request.get_method()
            path = url.path
            calls.append((method, path))
            assert url.netloc == "cad.onshape.com"
            assert request.get_header("Authorization").startswith("Basic ")
            if path == f"/api/v11/documents/d/{DID}/v/{VID}/elements":
                body = [{"id": EID, "elementType": "PARTSTUDIO", "microversionId": MID}]
            elif path == f"/api/v11/elements{REFERENCE}/configuration":
                body = {"configurationParameters": definitions, "sourceMicroversion": MID,
                        "currentConfiguration": []}
            elif path == f"/api/v11/elements/d/{DID}/e/{EID}/configurationencodings":
                assert parse_qs(url.query) == {"versionId": [VID]}
                assert json.loads(request.data) == {"parameters": [
                    {"parameterId": "Gap", "parameterValue": "2 mm"},
                    {"parameterId": "Size", "parameterValue": "large"},
                ]}
                body = {"encodedId": ENCODED}
            elif path == f"/api/v11/variables{REFERENCE}/variables":
                assert parse_qs(url.query) == {
                    "configuration": [ENCODED], "includeValuesAndReferencedVariables": ["true"]
                }
                body = [table]
            elif path == f"/api/v11/partstudios{REFERENCE}/export/step":
                options = json.loads(request.data)
                assert options["advancedParams"] == {"configuration": ENCODED}
                assert options["storeInDocument"] is False
                body = {"id": translation, "requestState": "DONE", "resultExternalDataIds": [artifact]}
            elif path == f"/api/v11/documents/d/{DID}/externaldata/{artifact}":
                return io.BytesIO(sample.read_bytes())
            else:
                pytest.fail(f"unexpected mock request: {path}")
            return io.BytesIO(json.dumps(body).encode())

    monkeypatch.setattr(server, "build_opener", lambda *args: MockOpener())
    discovery = asyncio.run(server.onshape_discover(
        "elements", url=f"https://cad.onshape.com/documents/{DID}/v/{VID}/e/{EID}"
    ))
    selection = discovery["elements"][0]["export_args"]
    result = asyncio.run(server.onshape_read_parameters(
        **selection, parameter_values={"Size": "large", "Gap": "2 mm"}
    ))
    assert result["status"] == "ok", result
    assert result["source_microversion_id"] == MID
    assert result["configuration"] == {
        "parameters": definitions, "encodedId": ENCODED, "current_configuration": [],
    }
    [returned] = result["variable_tables"]
    gap, angle, count, shape, missing = returned["variables"]
    assert gap["expression"] == "#plate / 2"
    assert [(row["value"], row["unit"]) for row in (gap, angle, count, shape, missing)] == [
        (2.0, "mm"), (60.0, "deg"), (3.0, None), (None, None), (None, None)
    ]
    assert shape["value_text"] == "[1, 2]" and missing["value_text"] is None
    exported = asyncio.run(server.onshape_export_cad(
        **result["export_args"], file_format="step", export_dir=str(tmp_path)
    ))
    assert exported["status"] == "ok", exported
    assert Path(exported["file"]).read_bytes() == sample.read_bytes()
    assert [method for method, _ in calls] == ["GET", "GET", "POST", "GET", "POST", "GET"]


def test_assembly_reads_configuration_without_encoding_when_already_selected(monkeypatch):
    class MockOpener:
        def open(self, request, timeout):
            assert request.get_method() == "GET"
            assert urlsplit(request.full_url).path == f"/api/v11/elements{REFERENCE}/configuration"
            return io.BytesIO(b'{"configurationParameters":[]}')

    monkeypatch.setattr(server, "build_opener", lambda *args: MockOpener())
    result = asyncio.run(server.onshape_read_parameters(DID, EID, VID, "assembly", ENCODED))
    assert result["status"] == "ok", result
    assert result["configuration"]["encodedId"] == ENCODED
    assert result["variable_tables"] is None


def test_unknown_configuration_parameter_does_not_call_encoder(monkeypatch):
    class MockOpener:
        def open(self, request, timeout):
            assert request.get_method() == "GET"
            return io.BytesIO(b'{"configurationParameters":[]}')

    monkeypatch.setattr(server, "build_opener", lambda *args: MockOpener())
    result = asyncio.run(server.onshape_read_parameters(
        DID, EID, VID, "assembly", parameter_values={"unknown": "value", "bogus": "1"}
    ))
    assert (result["status"], result["stage"]) == ("error", "input")
    assert result["unknown_parameter_ids"] == ["bogus", "unknown"]
