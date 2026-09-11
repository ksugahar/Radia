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
