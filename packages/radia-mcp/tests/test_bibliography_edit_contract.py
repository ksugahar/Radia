"""Bibliography key and non-destructive edit contracts; synthetic inputs only."""
import pytest

from radia_mcp.bibliography._bibparse import first_author_lastname, first_title_word
from radia_mcp.bibliography.plans.T6_canonicalize_keys import bibliography_canonicalize_keys
from radia_mcp.bibliography.plans.T10_normalize_journals import bibliography_normalize_journal_names
from radia_mcp.bibliography._bibparse import parse_bib
from radia_mcp.bibliography import _source_edit


@pytest.mark.parametrize("authors,expected", [
    ("{Research and Development Group} and Doe, Jane", "researchanddevelopmentgroup"),
    ("{Institute, Division A} AND Doe, Jane", "institutedivisiona"),
    ("Ludwig van Beethoven and Doe, Jane", "vanbeethoven"),
    ("Charles de la Vallee Poussin", "delavalleepoussin"),
    ("Smith, Jr., John and Doe, Jane", "smith"),
    ("Doe, Jane AND Roe, John", "doe"),
    ("Müller, Hans", "muller"), (r'M{\"u}ller, Hans', "muller"),
    ("John {van der Waals}", "vanderwaals"), ("", "unknown"),
])
def test_first_author_respects_bibtex_name_boundaries(authors, expected):
    assert first_author_lastname(authors) == expected


@pytest.mark.parametrize("title,expected", [
    (r"\textit{Magnetic} analysis", "magnetic"),
    ("The Étude of fields", "etude"), ("On the field", "field")])
def test_title_key_ignores_format_commands(title, expected):
    assert first_title_word(title) == expected


@pytest.mark.parametrize("dry_run", [True, False])
@pytest.mark.parametrize("keys", [("one", "two"), ("same", "same")])
def test_key_collision_never_writes(tmp_path, dry_run, keys):
    path = tmp_path / "references.bib"
    original = "\n".join(f"@article{{{key},author={{Doe, Jane}},year={{2024}},title={{Magnetic fields}}}}"
                         for key in keys).encode()
    path.write_bytes(original)
    result = bibliography_canonicalize_keys(str(path), dry_run=dry_run)
    assert result.startswith("Error:") and path.read_bytes() == original


def test_invalid_utf8_never_becomes_replacement_text(tmp_path):
    path = tmp_path / "references.bib"
    path.write_bytes(b"@article{old,title={\xff}}")
    result = bibliography_canonicalize_keys(str(path), dry_run=False)
    assert result.startswith("Error:")
    assert path.read_bytes() == b"@article{old,title={\xff}}"


@pytest.mark.parametrize("newline", ["\n", "\r\n"])
@pytest.mark.parametrize("dry_run", [True, False])
def test_key_edit_preserves_source_and_updates_direct_crossref(tmp_path, newline, dry_run):
    text = ('% Keep comments and macro syntax exactly.\n'
            '@string{prefix = "Shared "}\n'
            '@article{old, author={Doe, Jane}, year=2024, title={Magnetic fields},\n'
            ' note=prefix # {detail}} % trailing comment\n'
            '@inproceedings{roe2025electric,author={Roe, John},year=2025,\n'
            ' title={Electric fields}, crossref="old"}\n').replace("\n", newline)
    original = text.encode("utf-8")
    path = tmp_path / "references.bib"
    path.write_bytes(original)
    result = bibliography_canonicalize_keys(str(path), dry_run=dry_run)
    assert not result.startswith("Error:"), result
    expected = text if dry_run else text.replace("@article{old,", "@article{doe2024magnetic,").replace(
        'crossref="old"', 'crossref={doe2024magnetic}')
    assert path.read_bytes() == expected.encode("utf-8")


@pytest.mark.parametrize("newline", ["\n", "\r\n"])
@pytest.mark.parametrize("dry_run", [True, False])
def test_journal_edit_preserves_non_target_expression(tmp_path, newline, dry_run):
    text = ('\ufeff% journal maintenance\n@string(prefix = "Shared ) text")\n'
            '@article(old, title=prefix # {rest}, journal="IEEE Transactions on Magnetics",\n'
            ' author={Doe, Jane},year=2024) % retain\n').replace("\n", newline)
    path = tmp_path / "references.bib"
    path.write_bytes(text.encode())
    result = bibliography_normalize_journal_names(str(path), dry_run=dry_run)
    assert not result.startswith("Error:"), result
    expected = text if dry_run else text.replace('"IEEE Transactions on Magnetics"', '{IEEE Trans. Magn.}')
    assert path.read_bytes() == expected.encode()


@pytest.mark.parametrize("tool", [bibliography_canonicalize_keys, bibliography_normalize_journal_names])
@pytest.mark.parametrize("body", [
    '@article{old,title={unterminated', '@article{old,title="unterminated',
    '@article{old,title={A},title={B}}', '@article{old,title=}',
    '@article{old,title={A} #}', '@article{old,title={A}',
])
def test_malformed_source_is_not_rewritten(tmp_path, tool, body):
    path = tmp_path / "references.bib"
    path.write_bytes(body.encode())
    result = tool(str(path), dry_run=False)
    assert result.startswith("Error:")
    assert path.read_bytes() == body.encode()


@pytest.mark.parametrize("mode", ["replace_failure", "source_changed"])
def test_staged_edit_does_not_overwrite_on_failure(tmp_path, monkeypatch, mode):
    path = tmp_path / "references.bib"
    path.write_bytes(b"@misc{old,title={A}}")
    original, text, entries = _source_edit.read_source(path)
    edits = [(*entries[0].key_span, "new")]
    if mode == "replace_failure":
        def fail(*args):
            raise OSError("locked")
        monkeypatch.setattr(_source_edit.os, "replace", fail)
        expected = original
    else:
        expected = b"another session's content"
        path.write_bytes(expected)
    with pytest.raises((OSError, ValueError)):
        _source_edit.write_source_edits(path, original, text, edits)
    assert path.read_bytes() == expected
    assert not list(tmp_path.glob(".radia-bib-*"))


def test_parser_source_spans_preserve_escaped_braces_and_quotes():
    text = r'@article(key, title="A \"quote\" and \{literal\}", year=2024)'
    entry = parse_bib(text)[0]
    assert text[slice(*entry.key_span)] == "key"
    assert text[slice(*entry.field_spans["title"])] == '"' + entry.fields["title"] + '"'
