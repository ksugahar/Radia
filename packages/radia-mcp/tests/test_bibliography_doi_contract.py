"""Offline contracts for the standalone bibliography DOI adapter."""
import json

import pytest

from radia_mcp.bibliography.plans import T1_doi_to_bibtex as adapter


def metadata():
    return {"DOI": "10.1234/test", "type": "journal-article", "title": ["R&D 20%"],
            "author": [{"name": "Research and Development Group"}],
            "issued": {"date-parts": [[2024]]}, "page": "1--5"}


def test_complete_bibtex_protects_corporate_name_and_text(monkeypatch):
    monkeypatch.setattr(adapter, "_fetch_crossref", lambda doi: metadata())
    result = adapter.bibliography_doi_to_bibtex("https://doi.org/10.1234/test")
    assert not result.startswith("Error:")
    assert r"title = {R\&D 20\%}" in result
    assert "author = {{Research and Development Group}}" in result
    assert "pages = {1--5}" in result


@pytest.mark.parametrize("change", [
    {"DOI": "10.1234/another"}, {"DOI": None}, {"title": []}, {"title": ["$x$"]},
    {"author": []}, {"author": [{"name": "Group"}, {}]}, {"author": [None]},
    {"author": [{"family": "Doe", "given": 42}]},
    {"issued": None}, {"issued": {"date-parts": [[None]]}},
    {"issued": {"date-parts": [[True]]}},
    {"issued": {}, "created": {"date-parts": [[2024]]}},
])
def test_incomplete_or_wrong_metadata_is_not_an_entry(monkeypatch, change):
    data = metadata()
    data.update(change)
    monkeypatch.setattr(adapter, "_fetch_crossref", lambda doi: data)
    assert adapter.bibliography_doi_to_bibtex("10.1234/test").startswith("Error:")


def test_published_date_and_null_given_name_are_supported(monkeypatch):
    data = metadata()
    data.update(issued=None, published={"date-parts": [[2025]]},
                author=[{"family": "Doe", "given": None}])
    monkeypatch.setattr(adapter, "_fetch_crossref", lambda doi: data)
    result = adapter.bibliography_doi_to_bibtex("10.1234/test")
    assert "year = {2025}" in result and "author = {Doe}" in result


@pytest.mark.parametrize("value", ["DOI:10.1234/A#B?C%D",
                                  "HTTPS://DX.DOI.ORG/10.1234/A%23B%3FC%25D"])
def test_crossref_request_encodes_doi_and_closes_response(monkeypatch, value):
    closed = []
    class Response:
        def __enter__(self):
            return self
        def __exit__(self, *args):
            closed.append(True)
        def read(self):
            return json.dumps({"message": metadata()}).encode()
    def open_url(req, **kwargs):
        assert req.full_url == "https://api.crossref.org/works/10.1234/A%23B%3FC%25D"
        return Response()
    monkeypatch.setattr(adapter.urllib.request, "urlopen", open_url)
    assert adapter._fetch_crossref(value) == metadata()
    assert closed == [True]


@pytest.mark.parametrize("payload", [[], None, {"message": []}, {"message": None}])
def test_invalid_response_shape_is_lookup_failure(monkeypatch, payload):
    from io import BytesIO
    monkeypatch.setattr(adapter.urllib.request, "urlopen", lambda *a, **k:
                        BytesIO(json.dumps(payload).encode()))
    assert adapter._fetch_crossref("10.1234/test") is None
