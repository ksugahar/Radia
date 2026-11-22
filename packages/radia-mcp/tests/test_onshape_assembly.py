"""Guards: fixed-version assembly reads and a STEP centroid/volume cross-check."""

import asyncio
import io
import json
from urllib.parse import parse_qs, urlsplit

import pytest

from radia_mcp.onshape import server


DID, VID, EID, ROOT_MID, SUB_EID, SUB_MID, PART_EID, PART_MID = (
    char * 24 for char in "abcdef12"
)
CONFIG = "Size=large;Gap=2%20mm"
PART_CONFIG = "Width=10%20mm"
PART_ID = "JHD"
ROOT = f"/api/v11/assemblies/d/{DID}/v/{VID}/e/{EID}"
PART = f"/api/v11/parts/d/{DID}/m/{PART_MID}/e/{PART_EID}/partid/{PART_ID}/massproperties"


def _source(element, microversion, configuration):
    return {"documentId": DID, "elementId": element, "documentMicroversion": microversion,
            "configuration": configuration, "fullConfiguration": configuration}


def _definition():
    part = _source(PART_EID, PART_MID, PART_CONFIG) | {"type": "Part", "partId": PART_ID}
    sub = _source(SUB_EID, SUB_MID, "default") | {"type": "Assembly", "id": "nested"}
    return {
        "rootAssembly": _source(EID, ROOT_MID, CONFIG) | {
            "instances": [
                part | {"id": "first", "name": "Fixed"},
                sub,
                part | {"id": "suppressed", "suppressed": True},
            ],
            "occurrences": [
                {"path": ["first"], "transform": [1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1]},
                {"path": ["nested", "second"], "transform": [
                    0, -1, 0, 0.070, 1, 0, 0, -0.003, 0, 0, 1, 0.005, 0, 0, 0, 1
                ]},
            ],
            "features": [
                {"id": "mate-1", "featureType": "mate", "suppressed": False,
                 "featureData": {"mateType": "FASTENED", "matedEntities": [
                     {"matedOccurrence": ["first"]}, {"matedOccurrence": ["nested", "second"]},
                 ]}},
                {"id": "skip", "featureType": "pattern"},
                {"id": "connector", "featureType": "mateConnector"},
                {"id": "group", "featureType": "mateGroup"},
                {"id": "relation", "featureType": "mateRelation"},
            ],
        },
        "subAssemblies": [sub | {"instances": [part | {"id": "second", "name": "Rotated"}]}],
    }


@pytest.fixture(autouse=True)
def credentials(monkeypatch):
    monkeypatch.setenv("ONSHAPE_ACCESS_KEY", "test-access")
    monkeypatch.setenv("ONSHAPE_SECRET_KEY", "test-secret")


def test_pinned_assembly_top_level_placements_match_step_readback(tmp_path, monkeypatch):
    cad = pytest.importorskip("build123d")
    from radia_mcp.build123d import server as build123d

    free_part = cad.Pos(10, 2, 3) * cad.Box(10, 4, 3)
    # Construct the expected world-space solid independently of the API matrix.
    second = cad.Pos(68, 7, 8) * cad.Box(4, 10, 3)
    free_part.label, second.label = "fixed_block", "rotated_block"
    step = tmp_path / "assembly.step"
    assert cad.export_step(cad.Compound(children=[free_part, second]), step)
    calls = []

    class MockOpener:
        def open(self, request, timeout):
            assert request.get_method() == "GET"
            assert request.get_header("Authorization").startswith("Basic ")
            url = urlsplit(request.full_url)
            assert url.netloc == "cad.onshape.com"
            calls.append(url.path)
            query = parse_qs(url.query)
            if url.path == f"/api/v11/documents/d/{DID}/v/{VID}/elements":
                body = [{"id": EID, "elementType": "ASSEMBLY"}]
            elif url.path == ROOT:
                assert query == {
                    "configuration": [CONFIG], "includeMateFeatures": ["true"],
                    "includeNonSolids": ["false"], "excludeSuppressed": ["false"],
                }
                body = _definition()
            elif url.path == PART:
                assert query == {"configuration": [PART_CONFIG], "useMassPropertyOverrides": ["false"]}
                body = {"microversionId": PART_MID, "bodies": {PART_ID: {
                    "volume": [120e-9, 119e-9, 121e-9], "centroid": [
                        .010, .002, .003, .009, .001, .002, .011, .003, .004
                    ], "mass": [.12], "hasMass": True,
                }}}
            elif url.path == f"{ROOT}/massproperties":
                assert query == {"configuration": [CONFIG]}
                body = {"volume": [240e-9], "centroid": [.039, .0045, .0055], "hasMass": True, "mass": [.24]}
            else:
                pytest.fail(f"unexpected mock request: {url.path}")
            return io.BytesIO(json.dumps(body).encode())

    monkeypatch.setattr(server, "build_opener", lambda *args: MockOpener())
    url = (f"https://cad.onshape.com/documents/{DID}/v/{VID}/e/{EID}"
           "?configuration=Size%3Dlarge%3BGap%3D2%2520mm")
    discovery = asyncio.run(server.onshape_discover("elements", url=url))
    args = discovery["elements"][0]["export_args"]
    result = asyncio.run(server.onshape_inspect_assembly(**args))
    assert result["status"] == "ok", result
    assert result["source_microversion_id"] == ROOT_MID
    assert result["version_id"] == VID and result["configuration"] == CONFIG
    assert result["unresolved"] == []
    assert calls == [f"/api/v11/documents/d/{DID}/v/{VID}/elements",
                     ROOT, PART, f"{ROOT}/massproperties"]
    mate, relation = result["definitions"][0]["mates"]
    assert relation["id"] == "relation"
    assert mate["featureData"]["mateType"] == "FASTENED"
    assert result["definitions"][0]["instances"][-1]["suppressed"]
    assert result["assembly_properties"]["volume_mm3"] == pytest.approx(240)

    readback = json.loads(build123d._inspect_step_readback(str(step)))
    assert readback["status"] == "ok", readback
    expected_labels = {("first",): "fixed_block", ("nested", "second"): "rotated_block"}
    geometry = {body["label"]: body for body in readback["solids"]}
    assert len(result["part_occurrences"]) == 2
    for row in result["part_occurrences"]:
        body = geometry[expected_labels[tuple(row["path"])]]
        assert row["properties"]["volume_mm3"] == pytest.approx(body["volume_mm3"])
        assert row["properties"]["centroid_mm"] == pytest.approx(body["centroid_mm"])
        assert row["properties"]["mass_kg"] == pytest.approx(.12)


def test_request_limit_and_unresolved_source_are_visible(monkeypatch):
    source = _definition()
    root = source["rootAssembly"]
    root["instances"].append(root["instances"][0] | {"id": "over_budget", "partId": "OTHER"})
    root["occurrences"].append(root["occurrences"][0] | {"path": ["over_budget"]})
    root["occurrences"].append(root["occurrences"][0] | {"path": ["missing"]})

    class MockOpener:
        def open(self, request, timeout):
            assert request.get_method() == "GET"
            path = urlsplit(request.full_url).path
            if path == ROOT:
                body = source
            elif path == PART:
                body = {"bodies": {}}
            elif path == f"{ROOT}/massproperties":
                body = {"volume": [], "centroid": []}
            else:
                pytest.fail(f"unexpected mock request: {path}")
            return io.BytesIO(json.dumps(body).encode())

    monkeypatch.setattr(server, "build_opener", lambda *args: MockOpener())
    result = asyncio.run(server.onshape_inspect_assembly(
        DID, EID, VID, "assembly", max_part_requests=1
    ))
    assert result["status"] == "ok", result
    assert result["assembly_properties"] is None
    assert result["unresolved"] == [
        {"path": ["first"], "reason": "part_properties_unavailable"},
        {"path": ["nested", "second"], "reason": "part_properties_unavailable"},
        {"path": ["over_budget"], "reason": "part_request_limit"},
        {"path": ["missing"], "reason": "unknown_occurrence"},
    ]
