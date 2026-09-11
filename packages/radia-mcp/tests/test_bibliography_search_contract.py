"""Offline response and candidate contracts for Crossref search."""
import io
import json
from urllib.parse import parse_qs, urlsplit

import pytest

from radia_mcp.bibliography.plans import T3_search_crossref as adapter


def run(monkeypatch, payload):
    stream = io.BytesIO(json.dumps(payload).encode())
    monkeypatch.setattr(adapter.urllib.request, "urlopen", lambda *a, **k: stream)
    result = adapter.bibliography_search_crossref("test")
    assert stream.closed
    return result


@pytest.mark.parametrize("payload", [None, [], {}, {"message": None}, {"message": {}},
    {"message": {"items": None}}, {"message": {"items": {}}},
    {"message": {"items": [None]}}, {"status": "failed", "message": {"items": []}},
    {"error": "unavailable", "message": {"items": []}}])
def test_invalid_search_response_is_not_no_results(monkeypatch, payload):
    assert run(monkeypatch, payload).startswith("Error:")


def test_explicit_empty_results_are_valid(monkeypatch):
    assert run(monkeypatch, {"message": {"items": []}}).startswith("No results")


@pytest.mark.parametrize("change", [{"title": None}, {"title": []}, {"DOI": None},
                                     {"DOI": "not-a-doi"}, {"title": [42]}])
def test_invalid_candidate_blocks_partial_list(monkeypatch, change):
    valid = {"title": ["A test paper"], "DOI": "10.1234/test"}
    broken = {**valid, **change}
    result = run(monkeypatch, {"message": {"items": [valid, broken]}})
    assert result.startswith("Error:") and "candidate 2" in result


@pytest.mark.parametrize("optional", [
    {}, {"author": None, "issued": None, "container-title": None},
    {"author": [None], "issued": {"date-parts": [[True]]}},
    {"author": [{"family": 42}], "issued": "invalid"},
    {"issued": {"date-parts": [[None]]}, "created": {"date-parts": [[2026]]}},
])
def test_missing_optional_metadata_stays_unknown(monkeypatch, optional):
    item = {"title": ["A test paper"], "DOI": "10.1234/test", **optional}
    result = run(monkeypatch, {"message": {"items": [item]}})
    assert "[?] ? — A test paper" in result and "DOI: 10.1234/test" in result


def test_query_encoding_and_success_fields(monkeypatch):
    query = "C++ & 電磁界 #1 50%"
    def open_url(req, **kwargs):
        assert parse_qs(urlsplit(req.full_url).query)["query"] == [query]
        return io.BytesIO(json.dumps({"message": {"items": [{
            "title": ["A test paper"], "DOI": "https://doi.org/10.1234/test",
            "author": [{"name": "Research Group"}], "issued": {"date-parts": [[2024]]},
            "container-title": ["Journal"]}]}}).encode())
    monkeypatch.setattr(adapter.urllib.request, "urlopen", open_url)
    result = adapter.bibliography_search_crossref(query)
    assert "[2024] Research Group" in result and "Journal" in result
    assert "DOI: 10.1234/test" in result


@pytest.mark.parametrize("query,limit", [("", 5), (" ", 5), ("test", True), ("test", "5")])
def test_invalid_inputs_do_not_fetch(monkeypatch, query, limit):
    def unexpected(*args, **kwargs):
        pytest.fail("invalid input must not fetch")
    monkeypatch.setattr(adapter.urllib.request, "urlopen", unexpected)
    assert adapter.bibliography_search_crossref(query, limit).startswith("Error:")


@pytest.mark.parametrize("raw", [b"invalid JSON", b"\xff"])
def test_decode_failure_closes_response(monkeypatch, raw):
    stream = io.BytesIO(raw)
    monkeypatch.setattr(adapter.urllib.request, "urlopen", lambda *a, **k: stream)
    assert adapter.bibliography_search_crossref("test").startswith("Error:")
    assert stream.closed
