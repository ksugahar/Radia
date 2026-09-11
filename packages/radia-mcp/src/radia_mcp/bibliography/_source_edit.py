"""Narrow UTF-8 source edits: preserve unrelated BibTeX text and line endings."""
import os
from pathlib import Path
import tempfile

from ._bibparse import parse_bib


def read_source(path: Path):
    original = path.read_bytes()
    text = original.decode("utf-8", errors="strict")
    return original, text, parse_bib(text)


def write_source_edits(path: Path, original: bytes, text: str, edits: list) -> None:
    """Stage validated non-overlapping ranges; never serialize the entire database."""
    edits = sorted(edits)
    previous_end = 0
    for start, end, replacement in edits:
        if not 0 <= previous_end <= start <= end <= len(text):
            raise ValueError("overlapping or invalid bibliography edit ranges")
        previous_end = end
    for start, end, replacement in reversed(edits):
        text = text[:start] + replacement + text[end:]
    # Reject a malformed result before touching the existing file.
    parse_bib(text)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(dir=path.parent, prefix=".radia-bib-",
                                         suffix=".tmp", delete=False) as output:
            temporary = Path(output.name)
            output.write(text.encode("utf-8", errors="strict"))
        if path.read_bytes() != original:
            raise ValueError("bibliography changed during operation; retry from current source")
        os.replace(temporary, path)
        temporary = None
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
