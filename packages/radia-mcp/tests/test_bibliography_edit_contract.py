"""Bibliography key and non-destructive edit contracts; synthetic inputs only."""
import pytest

from radia_mcp.bibliography._bibparse import first_author_lastname, first_title_word
from radia_mcp.bibliography.plans.T6_canonicalize_keys import bibliography_canonicalize_keys


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
