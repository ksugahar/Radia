"""Offline external-response contracts; no live API calls or downloads."""
from types import SimpleNamespace
from urllib.parse import parse_qs, urlsplit

import pytest

from radia_mcp.paper_writing import _arxiv_source as arxiv
from radia_mcp.paper_writing._citation_verify import paper_writing_verify_citation


EMPTY_FEED = '<feed xmlns="http://www.w3.org/2005/Atom"/>'


def mock_response(monkeypatch, text, status=200):
    def check_status():
        if status >= 400:
            raise RuntimeError(f"HTTP {status}")
    response = SimpleNamespace(text=text, raise_for_status=check_status)
    monkeypatch.setattr(arxiv, "_require_requests", lambda: SimpleNamespace(get=lambda *a, **kw: response))


@pytest.mark.parametrize("query", ["A & B", "C++", "field #1", "電磁界 + 最適化", "50% accuracy"])
def test_arxiv_query_special_characters_roundtrip(monkeypatch, query):
    requests = pytest.importorskip("requests")
    def get(url, **kwargs):
        prepared = requests.Request("GET", url, params=kwargs["params"]).prepare()
        parsed = urlsplit(prepared.url)
        assert parsed.fragment == ""
        assert parse_qs(parsed.query)["search_query"] == [f"all:{query}"]
        assert kwargs["timeout"] == 30
        return SimpleNamespace(text=EMPTY_FEED, raise_for_status=lambda: None)
    monkeypatch.setattr(arxiv, "_require_requests", lambda: SimpleNamespace(get=get))
    result = arxiv.paper_writing_arxiv_search(query, categories="all")
    assert result["n_results"] == 0
    assert result["papers"] == []


@pytest.mark.parametrize("body", [
    "not xml", "<html><body>Service unavailable</body></html>", "<feed/>",
    '<feed xmlns="http://www.w3.org/2005/Atom"><entry><title>Error</title><summary>Bad query</summary></entry></feed>',
    '<feed xmlns="http://www.w3.org/2005/Atom"><entry><id>http://arxiv.org/api/errors#bad_query</id><title>Invalid query</title></entry></feed>',
    '<feed xmlns="http://www.w3.org/2005/Atom"><entry><id>http://arxiv.org/abs/2603.17339</id></entry></feed>',
])
def test_arxiv_error_documents_are_not_papers_or_empty_results(monkeypatch, body):
    mock_response(monkeypatch, body)
    result = arxiv.paper_writing_arxiv_search("query")
    assert result["error"]
    assert "papers" not in result
    assert "n_results" not in result


@pytest.mark.parametrize("status", [400, 429, 503])
def test_arxiv_http_error_not_empty_success(monkeypatch, status):
    mock_response(monkeypatch, EMPTY_FEED, status)
    assert "error" in arxiv.paper_writing_arxiv_search("query")


def test_arxiv_old_style_identifier_remains_supported(monkeypatch):
    mock_response(monkeypatch, '<feed xmlns="http://www.w3.org/2005/Atom"><entry>'
                  '<id>http://arxiv.org/abs/physics/0501123v2</id><title>Sample</title>'
                  '</entry></feed>')
    result = arxiv.paper_writing_arxiv_search("query")
    assert result["papers"][0]["arxiv_id"] == "physics/0501123"


@pytest.mark.parametrize("mode", ["exception", "error", "none", "empty", "inconsistent", "bad_paper", "bad_title"])
def test_citation_search_failure_is_not_no_candidate(monkeypatch, tmp_path, mode):
    def search(*a, **kw):
        if mode == "exception":
            raise RuntimeError("injected transport failure")
        return {"error": {"error": "HTTP 503"}, "none": None, "empty": {},
                "inconsistent": {"n_results": 0, "papers": [{}]},
                "bad_paper": {"n_results": 1, "papers": [{}]},
                "bad_title": {"n_results": 1, "papers": [{"title": 42, "arxiv_id": "2603.17339"}]}}[mode]
    monkeypatch.setattr(arxiv, "paper_writing_arxiv_search", search)
    bib = tmp_path / "references.bib"
    bib.write_text("", encoding="utf-8")
    result = paper_writing_verify_citation(claim="Electromagnetic field analysis", bib_path=str(bib))
    assert result["verdict"] == "error"
    assert result["suggested_bibtex"] is None
    assert result["verification_method"] == "arxiv-search"
    assert "Do not infer" in result["advice"]


def test_real_search_error_propagates_to_citation_verification(monkeypatch, tmp_path):
    mock_response(monkeypatch, "<html>Unavailable</html>")
    bib = tmp_path / "references.bib"
    bib.write_text("", encoding="utf-8")
    result = paper_writing_verify_citation(claim="Electromagnetic analysis", bib_path=str(bib))
    assert result["verdict"] == "error"


def test_valid_empty_feed_means_no_candidate(monkeypatch, tmp_path):
    mock_response(monkeypatch, EMPTY_FEED)
    bib = tmp_path / "references.bib"
    bib.write_text("", encoding="utf-8")
    result = paper_writing_verify_citation(claim="Electromagnetic analysis", bib_path=str(bib))
    assert result["verdict"] == "no_candidate_found"
