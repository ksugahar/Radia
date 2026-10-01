"""Narrow UTF-8 source edits: preserve unrelated BibTeX text and line endings."""
import os
from pathlib import Path
import tempfile

from ._bibparse import parse_bib
from ._write_lock import target_lock


def literal_value(expression: str) -> str | None:
    """Return one complete braced/quoted literal; never expand a macro or #."""
    if not expression or expression[0] not in '{"':
        return None
    quoted = expression[0] == '"'
    depth = 0 if quoted else 1
    i = 1
    while i < len(expression):
        char = expression[i]
        if char == "\\":
            i += 2
            continue
        if char == "{":
            depth += 1
        elif char == "}":
            depth -= 1
            if not quoted and depth == 0:
                return expression[1:i] if i == len(expression) - 1 else None
        elif quoted and char == '"' and depth == 0:
            return expression[1:i] if i == len(expression) - 1 else None
        i += 1
    return None


def read_source(path: Path):
    original = path.read_bytes()
    text = original.decode("utf-8", errors="strict")
    return original, text, parse_bib(text)


def require_canonical_target(path: Path) -> None:
    """Only the canonical parent bibliography is ever rewritten.

    Manuscript folders cite canonical keys and ship a generated ``.bbl``; a
    local ``.bib`` copy is not maintained, so editing one would fork the data.
    """
    from .plans import T14_canonical

    canonical = Path(T14_canonical.CANONICAL).resolve()
    if Path(path).resolve() != canonical:
        raise ValueError(
            f"refusing to rewrite {path}: only the canonical bibliography "
            f"({canonical}) is edited; manuscripts use bibliography_make_bbl")


def write_source_edits(path: Path, original: bytes, text: str, edits: list) -> None:
    require_canonical_target(path)
    with target_lock(path):
        _write_source_edits_unlocked(path, original, text, edits)


def _write_source_edits_unlocked(path: Path, original: bytes, text: str, edits: list) -> None:
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
