"""Source-preserving surname convention checks; no canonical data mutations."""
import pytest

from radia_mcp.bibliography.plans.T13_check_surname_braces import (
    _wrap_surname, bibliography_check_surname_braces,
)


@pytest.mark.parametrize("before,after", [
    ("Doe, Jane", "{Doe}, Jane"),
    ("  Doe , Jane  ", "  {Doe} , Jane  "),
    ("Jane Doe", "Jane {Doe}"),
    ("H. A. van  der Vorst", "H. A. {van  der Vorst}"),
    ("H. A. {van der Vorst}", "H. A. {van der Vorst}"),
    ("{Research and Development Group}", "{Research and Development Group}"),
    ("{Doe}, Jr., Jane", "{Doe}, Jr., Jane"),
    ("Doe, Jr., Jane", "{Doe}, Jr., Jane"),
    ("others", "others"),
    ("", ""),
])
def test_surname_grouping_is_brace_aware_and_idempotent(before, after):
    assert _wrap_surname(before) == after
    assert _wrap_surname(after) == after


@pytest.mark.parametrize("newline", ["\n", "\r\n"])
def test_fix_preserves_comments_macros_bom_and_author_whitespace(tmp_path, newline):
    source = ("\ufeff% author={Comment Name}\n"
              "@string{note = {author={Macro Name}}}\n"
              "@misc(one, author = {Jane Doe AND\n  {Research and Development Group} and others},\n"
              "title={Nested {Title}}, note=note # { suffix})\n").replace("\n", newline)
    path = tmp_path / "fixture.bib"
    path.write_bytes(source.encode("utf-8"))
    assert "FAIL" in bibliography_check_surname_braces(str(path))
    assert path.read_bytes() == source.encode("utf-8")
    assert "FIXED" in bibliography_check_surname_braces(str(path), fix=True)
    expected = source.replace("Jane Doe", "Jane {Doe}").encode("utf-8")
    assert path.read_bytes() == expected
    assert "OK" in bibliography_check_surname_braces(str(path), fix=True)
    assert path.read_bytes() == expected


@pytest.mark.parametrize("expression", ['names', '{Jane} # { Doe}', '"Jane" # " Doe"'])
def test_macro_author_is_unavailable_and_blocks_entire_fix(tmp_path, expression):
    data = ('@string{names={Jane Doe}}\n@misc{a,author={John Roe}}\n'
            '@misc{b,author=' + expression + '}').encode()
    path = tmp_path / "fixture.bib"
    path.write_bytes(data)
    result = bibliography_check_surname_braces(str(path), fix=True)
    assert "UNAVAILABLE" in result and "FIXED" not in result
    assert path.read_bytes() == data


def test_quoted_nested_author_literal_is_fixed(tmp_path):
    path = tmp_path / "fixture.bib"
    path.write_text('@misc{a, author="Jane Doe and {Alpha and Beta}"}', encoding="utf-8")
    assert "FIXED" in bibliography_check_surname_braces(str(path), fix=True)
    assert path.read_text(encoding="utf-8") == '@misc{a, author="Jane {Doe} and {Alpha and Beta}"}'


@pytest.mark.parametrize("data", [b"\xff", b"not a bibliography", b"@misc{a, author={Jane Doe}"])
def test_invalid_input_is_never_reencoded_or_rewritten(tmp_path, data):
    path = tmp_path / "fixture.bib"
    path.write_bytes(data)
    assert bibliography_check_surname_braces(str(path), fix=True).startswith("Error:")
    assert path.read_bytes() == data


def test_failed_surname_publication_preserves_source(tmp_path, monkeypatch):
    from radia_mcp.bibliography import _source_edit
    path = tmp_path / "fixture.bib"
    data = b"@misc{a,author={Jane Doe}}"
    path.write_bytes(data)
    def fail(*args):
        raise PermissionError("locked")
    monkeypatch.setattr(_source_edit.os, "replace", fail)
    assert "Error:" in bibliography_check_surname_braces(str(path), fix=True)
    assert path.read_bytes() == data
    assert not list(tmp_path.glob(".radia-bib-*.tmp"))
