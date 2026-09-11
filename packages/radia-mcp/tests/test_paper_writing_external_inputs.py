"""Offline external-response contracts; no live API calls or downloads."""
from types import SimpleNamespace
import gzip
import io
import tarfile
import sys
from pathlib import Path
from urllib.parse import parse_qs, urlsplit

import pytest

from radia_mcp.paper_writing import _arxiv_source as arxiv
from radia_mcp.paper_writing._citation_verify import paper_writing_verify_citation
from radia_mcp.paper_writing import paper_download as download


def _crossref_stub(monkeypatch, message):
    seen = []
    def get(url, **kwargs):
        seen.append(url)
        return SimpleNamespace(status_code=200, json=lambda: {"message": message})
    monkeypatch.setattr(download, "_require_requests", lambda: SimpleNamespace(get=get))
    return seen


@pytest.mark.parametrize("value", [
    "10.1109/a#b?c%d", "DOI:10.1109/a#b?c%d",
    "HTTPS://DX.DOI.ORG/10.1109/a%23b%3Fc%25d",
    "doi.org/10.1109/a%23b%3Fc%25d",
])
def test_doi_consumers_share_lossless_encoding(monkeypatch, value):
    seen = _crossref_stub(monkeypatch, {})
    result = download.paper_writing_resolve_doi(value)
    assert result["doi"] == "10.1109/a#b?c%d"
    assert seen == ["https://api.crossref.org/works/10.1109/a%23b%3Fc%25d"]
    assert result["url"] == "https://doi.org/10.1109/a%23b%3Fc%25d"
    def get(url, **kwargs):
        seen.append(url)
        return SimpleNamespace(status_code=200,
                               url="https://ieeexplore.ieee.org/document/123/")
    monkeypatch.setattr(download, "_require_requests", lambda: SimpleNamespace(
        Session=lambda: SimpleNamespace(get=get)))
    ieee = download.paper_writing_ieee_doi_to_arnumber(value)
    assert ieee["ok"] and ieee["doi"] == result["doi"]
    assert seen[-1] == result["url"]


@pytest.mark.parametrize("message", [None, [], "bad", 42])
def test_crossref_invalid_message_is_not_success(monkeypatch, message):
    _crossref_stub(monkeypatch, message)
    result = download.paper_writing_resolve_doi("10.1234/test")
    assert not result["ok"] and result["temporary_failure"]


@pytest.mark.parametrize("date", [None, {}, {"date-parts": []},
                                   {"date-parts": [[None]]},
                                   {"date-parts": [[True]]}, "bad"])
def test_crossref_created_date_is_never_publication_year(monkeypatch, date):
    _crossref_stub(monkeypatch, {
        "title": ["Test"], "author": [{"name": "Research Group"}],
        "published": date, "created": {"date-parts": [[2026]]},
    })
    result = download.paper_writing_doi_to_bibtex("10.1234/test")
    assert not result["ok"]
    assert result["error_kind"] == "incomplete_metadata"
    assert result["missing_fields"] == ["year"]
    assert result["metadata"]["year"] is None
    assert "bibtex" not in result


def test_crossref_issued_fallback_and_corporate_author(monkeypatch):
    _crossref_stub(monkeypatch, {
        "title": ["Test"], "type": "journal-article",
        "author": [{"name": "Research and Development Group"},
                   {"family": "Doe", "given": None}],
        "published": {"date-parts": [[]]}, "issued": {"date-parts": [[2024]]},
    })
    result = download.paper_writing_doi_to_bibtex("10.1234/test", "test2024")
    assert result["ok"]
    assert result["metadata"]["year"] == 2024
    assert "author  = {{Research and Development Group} and Doe}" in result["bibtex"]


@pytest.mark.parametrize("overrides,missing", [
    ({"title": "Not a title list"}, "title"),
    ({"author": None}, "authors"),
    ({"author": [None, {}, {"name": 42}]}, "authors"),
])
def test_partial_crossref_metadata_cannot_generate_bibtex(monkeypatch, overrides, missing):
    message = {"title": ["Test"], "author": [{"name": "Group"}],
               "issued": {"date-parts": [[2024]]}}
    message.update(overrides)
    _crossref_stub(monkeypatch, message)
    result = download.paper_writing_doi_to_bibtex("10.1234/test")
    assert not result["ok"] and result["missing_fields"] == [missing]
    assert "bibtex" not in result


@pytest.mark.parametrize("bad", [None, {}, {"name": 12},
                                 {"family": "Doe", "given": 42}])
def test_crossref_partial_author_list_blocks_bibtex(monkeypatch, bad):
    _crossref_stub(monkeypatch, {"title": ["Test"],
        "author": [{"name": "Group"}, bad], "issued": {"date-parts": [[2024]]}})
    result = download.paper_writing_doi_to_bibtex("10.1234/test")
    assert not result["ok"] and result["missing_fields"] == ["authors"]
    assert result["metadata"]["invalid_author_indices"] == [1]
    assert "bibtex" not in result


def test_crossref_escapes_text_after_decoding_entities(monkeypatch):
    _crossref_stub(monkeypatch, {
        "title": ["&lt;i&gt;R&amp;amp;D&lt;/i&gt; at 20%: A_B #1"],
        "type": "journal-article", "container-title": ["Science &amp; Design"],
        "author": [{"name": "R&D Group"}], "issued": {"date-parts": [[2024]]}})
    result = download.paper_writing_doi_to_bibtex("10.1234/test", "test2024")
    assert result["ok"]
    assert r"title   = {R\&D at 20\%: A\_B \#1}" in result["bibtex"]
    assert r"journal = {Science \& Design}" in result["bibtex"]
    assert r"author  = {{R\&D Group}}" in result["bibtex"]


@pytest.mark.parametrize("title", ["<math>x</math>", r"A $x$ field", r"A \alpha field",
                                  "<i></i>"])
def test_crossref_unsupported_or_empty_markup_requires_review(monkeypatch, title):
    _crossref_stub(monkeypatch, {"title": [title], "author": [{"name": "Group"}],
                               "issued": {"date-parts": [[2024]]}})
    result = download.paper_writing_doi_to_bibtex("10.1234/test")
    assert not result["ok"] and result["error_kind"] == "unsupported_metadata"
    assert "bibtex" not in result


def _arxiv_citation(monkeypatch, tmp_path, lookup):
    monkeypatch.setattr(arxiv, "paper_writing_semantic_scholar_lookup", lambda *_: lookup)
    bib = tmp_path / "references.bib"
    bib.write_text("", encoding="utf-8")
    return paper_writing_verify_citation(claim="", bib_path=str(bib),
                                        candidate_arxiv_id="2603.17339")


def test_arxiv_bibtex_keeps_all_authors_and_escapes_text(monkeypatch, tmp_path):
    names = [f"Author {n}" for n in range(7)]
    result = _arxiv_citation(monkeypatch, tmp_path, {"metadata": {
        "title": "R&amp;D at 20%", "authors": [{"name": name} for name in names],
        "year": 2024, "externalIds": None}})
    assert result["verdict"] == "ready_to_insert"
    assert "author       = {" + " and ".join(names) + "}" in result["suggested_bibtex"]
    assert r"title        = {R\&D at 20\%}" in result["suggested_bibtex"]
    assert "Semantic Scholar" in result["advice"]


@pytest.mark.parametrize("override", [
    {"authors": None}, {"authors": []}, {"authors": [{"name": "A"}, {}]},
    {"authors": [{"name": "A"}, {"name": 12}]}, {"year": None},
    {"authors": [{"name": "<i></i>"}]},
    {"year": True}, {"year": "2024"}, {"title": None},
    {"title": "<i></i>"}, {"title": "$x$"}, {"externalIds": []},
    {"externalIds": {"DOI": 42}},
])
def test_arxiv_incomplete_metadata_is_not_insertion_ready(monkeypatch, tmp_path, override):
    metadata = {"title": "Test", "authors": [{"name": "A"}], "year": 2024}
    metadata.update(override)
    result = _arxiv_citation(monkeypatch, tmp_path, {"metadata": metadata})
    assert result["verdict"] == "error"
    assert result["suggested_bibtex"] is None


@pytest.mark.parametrize("lookup", [None, [], {"error": "HTTP 429"},
                                     {"metadata": None}, {"metadata": []}])
def test_arxiv_failed_lookup_does_not_mean_no_candidate(monkeypatch, tmp_path, lookup):
    result = _arxiv_citation(monkeypatch, tmp_path, lookup)
    assert result["verdict"] == "error"
    assert result["suggested_bibtex"] is None


def test_arxiv_lookup_exception_is_reported(monkeypatch, tmp_path):
    def fail(*_):
        raise ValueError("invalid JSON")
    monkeypatch.setattr(arxiv, "paper_writing_semantic_scholar_lookup", fail)
    bib = tmp_path / "references.bib"
    bib.write_text("", encoding="utf-8")
    result = paper_writing_verify_citation(claim="", bib_path=str(bib),
                                         candidate_arxiv_id="2603.17339")
    assert result["verdict"] == "error" and "invalid JSON" in result["advice"]
    assert result["suggested_bibtex"] is None


def test_arxiv_doi_route_still_delegates(monkeypatch, tmp_path):
    from radia_mcp.paper_writing import _citation_verify as verify
    original = verify.paper_writing_verify_citation
    seen = []
    def delegated(**kwargs):
        seen.append(kwargs)
        return {"delegated": True}
    monkeypatch.setattr(verify, "paper_writing_verify_citation", delegated)
    monkeypatch.setattr(arxiv, "paper_writing_semantic_scholar_lookup", lambda *_: {
        "metadata": {"externalIds": {"DOI": "10.1234/test"}}})
    bib = tmp_path / "references.bib"
    bib.write_text("", encoding="utf-8")
    result = original(claim="claim", bib_path=str(bib), candidate_arxiv_id="2603.17339")
    assert result == {"delegated": True}
    assert seen[0]["candidate_doi"] == "10.1234/test"


@pytest.mark.parametrize("response", [None, [], {}, {"ok": "yes"},
    {"ok": False, "error": "connection failed"}, {"ok": False, "error_kind": "temporary"}])
def test_doi_unknown_resolver_failure_is_error(monkeypatch, tmp_path, response):
    monkeypatch.setattr(download, "paper_writing_resolve_doi", lambda *_: response)
    bib = tmp_path / "references.bib"
    bib.write_text("", encoding="utf-8")
    result = paper_writing_verify_citation("", str(bib), candidate_doi="10.1234/test")
    assert result["verdict"] == "error"
    assert result["suggested_bibtex"] is None


@pytest.mark.parametrize("response", [None, [], {}, {"ok": True},
    {"ok": True, "bibtex": "", "citation_key": "test"},
    {"ok": True, "bibtex": "entry", "citation_key": None},
    {"ok": "yes", "bibtex": "entry", "citation_key": "test"}])
def test_doi_incomplete_generator_response_is_error(monkeypatch, tmp_path, response):
    monkeypatch.setattr(download, "paper_writing_resolve_doi", lambda *_: {"ok": True})
    monkeypatch.setattr(download, "paper_writing_doi_to_bibtex", lambda *_: response)
    bib = tmp_path / "references.bib"
    bib.write_text("", encoding="utf-8")
    result = paper_writing_verify_citation("", str(bib), candidate_doi="10.1234/test")
    assert result["verdict"] == "error"
    assert result["suggested_bibtex"] is None


@pytest.mark.parametrize("stage", ["resolver", "generator"])
def test_doi_dependency_exception_is_error(monkeypatch, tmp_path, stage):
    def fail(*_):
        raise RuntimeError("dependency unavailable")
    monkeypatch.setattr(download, "paper_writing_resolve_doi",
                        fail if stage == "resolver" else lambda *_: {"ok": True})
    monkeypatch.setattr(download, "paper_writing_doi_to_bibtex", fail)
    bib = tmp_path / "references.bib"
    bib.write_text("", encoding="utf-8")
    result = paper_writing_verify_citation("", str(bib), candidate_doi="10.1234/test")
    assert result["verdict"] == "error" and "dependency unavailable" in result["advice"]
    assert result["suggested_bibtex"] is None


@pytest.mark.parametrize("suffix", [".", ",", ";", "#?%"])
def test_doi_equality_preserves_identifier_punctuation(suffix):
    from radia_mcp.paper_writing._citation_verify import _doi_already_cited
    entries = {"target": {"doi": "10.1234/test" + suffix}}
    assert _doi_already_cited(entries, "10.1234/test") is None
    from urllib.parse import quote
    url = "HTTPS://DX.DOI.ORG/" + quote("10.1234/TEST" + suffix, safe="/.")
    assert _doi_already_cited(entries, url) == "target"


@pytest.mark.parametrize("mode", ["missing", "directory", "encoding", "denied"])
def test_unreadable_bibliography_blocks_external_lookup(monkeypatch, tmp_path, mode):
    from radia_mcp.paper_writing import _citation_verify as verify
    bib = tmp_path / "references.bib"
    if mode == "directory":
        bib.mkdir()
    elif mode == "encoding":
        bib.write_bytes(b"\xff")
    elif mode == "denied":
        bib.write_text("", encoding="utf-8")
        def denied(*_):
            raise PermissionError("access denied")
        monkeypatch.setattr(verify, "_parse_bib_lightweight", denied)
    def unexpected(*_):
        pytest.fail("Unreadable bibliography must stop before network lookup")
    monkeypatch.setattr(download, "paper_writing_resolve_doi", unexpected)
    result = verify.paper_writing_verify_citation("", str(bib), candidate_doi="10.1234/test")
    assert result["verdict"] == "error"
    assert result["suggested_bibtex"] is None


class StreamResponse:
    status_code = 200

    def __init__(self, content, fail=False):
        self.content = content
        self.fail = fail
        self.closed = False

    def raise_for_status(self):
        pass

    def iter_content(self, chunk_size=65536):
        yield self.content
        if self.fail:
            raise OSError("injected interrupted download")

    def close(self):
        self.closed = True


@pytest.mark.parametrize("mode", ["corrupt", "interrupted", "oversize", "success"])
def test_verified_pdf_save_preserves_destination_until_valid(monkeypatch, tmp_path, mode):
    dest = tmp_path / "paper.pdf"
    dest.write_bytes(b"previous user file")
    response = StreamResponse(b"downloaded bytes", fail=mode == "interrupted")
    monkeypatch.setattr(download, "_verify_pdf", lambda path: {"ok": mode != "corrupt", "error": "injected"})
    if mode == "oversize":
        monkeypatch.setattr(download, "_MAX_PDF_DOWNLOAD_BYTES", 4)
    result = download._save_verified_pdf(response, str(dest))
    assert result["ok"] is (mode == "success")
    assert dest.read_bytes() == (b"downloaded bytes" if mode == "success" else b"previous user file")
    assert list(tmp_path.glob(".radia-download-*")) == []
    assert response.closed


def test_pdf_magic_alone_does_not_validate_corruption(tmp_path):
    pytest.importorskip("fitz")
    path = tmp_path / "corrupt.pdf"
    path.write_bytes(b"%PDF-1.7\nnot a document\n")
    assert download._verify_pdf(str(path))["ok"] is False


def test_pdf_verification_requires_parser(monkeypatch, tmp_path):
    path = tmp_path / "unverified.pdf"
    path.write_bytes(b"%PDF-1.7\n")
    monkeypatch.setitem(sys.modules, "fitz", None)
    assert download._verify_pdf(str(path))["ok"] is False


@pytest.mark.parametrize("fetch", [download.paper_writing_ieee_download_pdf,
    download.paper_writing_sciencedirect_download_pdf, download.paper_writing_emerald_download_pdf])
def test_publishers_preserve_previous_file_then_allow_valid_retry(monkeypatch, tmp_path, fetch):
    response = StreamResponse(b"corrupt")
    session = SimpleNamespace(get=lambda *a, **kw: response)
    monkeypatch.setattr(download, "_require_requests", lambda: SimpleNamespace(Session=lambda: session))
    monkeypatch.setattr(download, "_verify_pdf", lambda path: {"ok": Path(path).read_bytes() == b"verified"})
    dest = tmp_path / "paper.pdf"
    dest.write_bytes(b"previous")
    result = fetch("1234567", str(dest), overwrite=True)
    assert result["ok"] is False
    assert dest.read_bytes() == b"previous"
    response = StreamResponse(b"verified")
    result = fetch("1234567", str(dest), overwrite=True)
    assert result["ok"] is True
    assert dest.read_bytes() == b"verified"
    assert not list(tmp_path.glob(".radia-download-*"))


@pytest.mark.parametrize("kind", ["single", "tar", "expanded", "member_count", "download", "interrupted", "member_size", "single_size"])
def test_arxiv_source_resource_limits(monkeypatch, kind):
    if kind in {"single", "expanded", "single_size"}:
        blob = gzip.compress(b"\\documentclass{article}" + b"x" * 100)
    else:
        buffer = io.BytesIO()
        with tarfile.open(fileobj=buffer, mode="w:gz") as archive:
            for name in ["main.tex", "appendix.tex"]:
                body = b"\\documentclass{article}"
                info = tarfile.TarInfo(name)
                info.size = len(body)
                archive.addfile(info, io.BytesIO(body))
        blob = buffer.getvalue()
    response = StreamResponse(blob, fail=kind == "interrupted")
    def get(*args, **kwargs):
        assert kwargs["stream"] is True
        return response
    monkeypatch.setattr(arxiv, "_require_requests", lambda: SimpleNamespace(get=get))
    if kind == "expanded":
        monkeypatch.setattr(arxiv, "_MAX_EXPANDED_BYTES", 32)
    if kind == "member_count":
        monkeypatch.setattr(arxiv, "_MAX_SOURCE_MEMBERS", 1)
    if kind == "download":
        monkeypatch.setattr(arxiv, "_MAX_SOURCE_BYTES", 4)
    result = arxiv.paper_writing_arxiv_fetch_latex_source(
        "2603.17339", max_chars_per_file=1 if kind in {"member_size", "single_size"} else 200_000)
    assert response.closed
    if kind in {"single", "tar"}:
        assert "error" not in result
        assert result["files"]
    else:
        assert result["error"]
        assert "files" not in result


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
