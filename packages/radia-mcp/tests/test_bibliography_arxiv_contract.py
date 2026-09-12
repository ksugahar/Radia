"""Offline identity/completeness contracts for the bibliography arXiv adapter."""
import pytest

from radia_mcp.bibliography.plans import T2_arxiv_to_bibtex as adapter


def feed(aid="physics/0501123v2", title="A paper", published="2024-01-02T00:00:00Z",
         authors="<author><name>A Author</name></author>"):
    return f'''<feed xmlns="http://www.w3.org/2005/Atom"
      xmlns:arxiv="http://arxiv.org/schemas/atom"><entry>
      <id>http://arxiv.org/abs/{aid}</id><title>{title}</title>
      <published>{published}</published>{authors}
      <arxiv:primary_category term="physics.class-ph"/></entry></feed>'''


@pytest.mark.parametrize("value", ["physics/0501123", "arXiv:physics/0501123",
    "https://arxiv.org/abs/physics/0501123", "HTTP://ARXIV.ORG/pdf/physics/0501123.pdf"])
def test_old_style_identifier_and_primary_category_retained(monkeypatch, value):
    monkeypatch.setattr(adapter, "_fetch_arxiv", lambda aid: feed())
    result = adapter.bibliography_arxiv_to_bibtex(value)
    assert not result.startswith("Error:")
    assert "physics/0501123v2" in result and "physics.class-ph" in result


@pytest.mark.parametrize("changes", [{"title": ""}, {"title": "Error"},
    {"published": ""}, {"published": "2024-02-30"}, {"authors": ""},
    {"authors": "<author><name>A Author</name></author><author/>"}, {"aid": "invalid"}])
def test_incomplete_response_is_not_bibtex(monkeypatch, changes):
    monkeypatch.setattr(adapter, "_fetch_arxiv", lambda aid: feed(**changes))
    assert adapter.bibliography_arxiv_to_bibtex("physics/0501123").startswith("Error:")


@pytest.mark.parametrize("requested,returned,success", [
    ("2603.17339", "2603.17339v2", True), ("2603.17339v2", "2603.17339v2", True),
    ("2603.17339", "2603.17340", False), ("2603.17339v2", "2603.17339v3", False)])
def test_requested_paper_and_version_must_match(monkeypatch, requested, returned, success):
    monkeypatch.setattr(adapter, "_fetch_arxiv", lambda aid: feed(aid=returned))
    result = adapter.bibliography_arxiv_to_bibtex(requested)
    assert result.startswith("Error:") is not success


@pytest.mark.parametrize("value", ["", "not-an-id", "2603.17339,2603.17340"])
def test_invalid_input_does_not_fetch(monkeypatch, value):
    def unexpected(*_):
        pytest.fail("invalid ID must stop before fetching")
    monkeypatch.setattr(adapter, "_fetch_arxiv", unexpected)
    assert adapter.bibliography_arxiv_to_bibtex(value).startswith("Error:")


@pytest.mark.parametrize("body", ["bad XML", "<feed/>",
    feed().replace("<feed ", "<wrong ").replace("</feed>", "</wrong>"),
    feed().replace("</entry>", "</entry><entry/>"), None])
def test_malformed_or_multiple_entries_rejected(body):
    assert adapter._atom_to_bibentry(body) is None
