r"""Tier 3 -- bibliography_check_surname_braces.

Lab style Rule 1 (paper_writing skill.md, 2026-06-14, from 2023 Compumag
review feedback): every author surname in a ``.bib`` ``author = {...}``
field must be wrapped in ``{}`` to defend against ``.bst`` case-mangling.

The bst files of IEEEtran (and many others) transform author names
during BibTeX output: forcing case, removing accents, ALL-CAPS surnames,
etc.  Text inside a brace group is preserved verbatim, so wrapping the
surname in ``{}`` (or the whole name in ``{}``) prevents:

    McDonald      -> Mcdonald     (default IEEEtran style)
    Bíró          -> bíró         (when accent macro is not protected)
    Köster        -> köster       (same problem)
    van der Vorst -> van der vorst

Acceptable forms (both produce safe output):

    author = {{Hollaus}, K.}                    # Lastname-first + brace
    author = {K. {Hollaus}}                     # First-then-Lastname + brace
    author = {{B{\'\i}r{\'o}}, Oszk{\'a}r}      # accent macro inside the brace
    author = {{van der Vorst}, H. A.}           # compound surname grouped

Violations:

    author = {Hollaus, K.}                      # surname NOT braced
    author = {K. Hollaus}                       # same
    author = {B{\'\i}r{\'o}, Oszk{\'a}r}        # surname accent unprotected
"""

from __future__ import annotations

import pathlib
import re

from .._bibparse import read_bib_file


# Find the FIRST top-level comma inside an author atom (after skipping any
# leading whitespace).  Returns -1 if not found.
def _first_top_comma(s: str) -> int:
    depth = 0
    for i, ch in enumerate(s):
        if ch == "{":
            depth += 1
        elif ch == "}":
            depth -= 1
        elif ch == "," and depth == 0:
            return i
    return -1


def _is_surname_braced(name: str) -> bool:
    """Return True if the surname is already wrapped in ``{}``."""
    stripped = name.lstrip()
    if not stripped:
        return True  # empty -> nothing to flag

    comma_idx = _first_top_comma(stripped)
    if comma_idx >= 0:
        # Lastname, Firstname form -- check if the part BEFORE the comma is braced.
        surname = stripped[:comma_idx].strip()
        if surname.startswith("{") and surname.endswith("}"):
            return True
        # Also accept whole-name brace ``{Lastname, Firstname}`` as a fallback
        # (some authors wrap the whole entry instead).
        if stripped.startswith("{") and stripped.endswith("}"):
            return True
        return False

    # No comma -- assume ``Firstname Lastname`` form; the last whitespace-separated
    # token is the surname, which should be ``{Lastname}``.
    tokens = stripped.split()
    if not tokens:
        return True
    last = tokens[-1]
    return last.startswith("{") and last.endswith("}")


def _wrap_surname(name: str) -> str:
    """Return the same name with the surname wrapped in ``{}``.  If already
    wrapped, return unchanged.
    """
    if _is_surname_braced(name):
        return name

    leading_ws = len(name) - len(name.lstrip())
    body = name[leading_ws:]

    comma_idx = _first_top_comma(body)
    if comma_idx >= 0:
        surname = body[:comma_idx]
        rest = body[comma_idx:]
        return name[:leading_ws] + "{" + surname + "}" + rest

    # No comma: wrap the LAST whitespace-separated token.
    tokens = body.split()
    if not tokens:
        return name
    last = tokens[-1]
    # Replace last occurrence of `last` in body with `{last}`.
    idx = body.rfind(last)
    return name[:leading_ws] + body[:idx] + "{" + last + "}" + body[idx + len(last):]


def _split_authors(value: str) -> list[str]:
    """Split a BibTeX author value by ``and`` (whitespace-bounded)."""
    return re.split(r"\s+and\s+", value)


def _join_authors(authors: list[str]) -> str:
    return " and ".join(authors)


# Match a full `author = {...}` field, value possibly multi-line.  Brace
# depth-1 only -- BibTeX values use simple nested braces.
_AUTHOR_RE = re.compile(
    r"(?P<key>author\s*=\s*)\{(?P<value>[^{}]*(?:\{[^{}]*\}[^{}]*)*)\}",
    re.IGNORECASE,
)


def _check_text(text: str) -> tuple[list[tuple[str, str]], str, int]:
    """Scan and (optionally) rewrite the .bib body.

    Returns:
        (violations, rewritten_text, n_fields_fixed)
        violations: list of (cite_key_or_index, original_author_value) for
        every field that has at least one un-braced surname.
        rewritten_text: text after wrapping all un-braced surnames in ``{}``.
        n_fields_fixed: number of author fields modified.
    """
    violations: list[tuple[str, str]] = []
    n_fixed = 0

    # Per-entry cite keys for nicer reporting.
    # First pass: walk through `@kind{key,` matches to know which key each
    # author field belongs to.  Cheap enough to redo in this scan.
    entry_keys: list[tuple[int, str]] = []
    for m in re.finditer(r"@\w+\s*\{\s*([^,\s]+)\s*,", text):
        entry_keys.append((m.start(), m.group(1)))

    def _key_at(pos: int) -> str:
        # Last entry starting at or before pos.
        key = "<unknown>"
        for start, k in entry_keys:
            if start <= pos:
                key = k
            else:
                break
        return key

    def _repl(match: re.Match) -> str:
        nonlocal n_fixed
        original = match.group("value")
        single = re.sub(r"\s+", " ", original).strip()

        authors = _split_authors(single)
        has_violation = any(not _is_surname_braced(a) for a in authors)
        if has_violation:
            violations.append((_key_at(match.start()), single))
            new_authors = [_wrap_surname(a) for a in authors]
            new_value = _join_authors(new_authors)
            if new_value != single:
                n_fixed += 1
                return match.group("key") + "{" + new_value + "}"
        return match.group(0)

    new_text = _AUTHOR_RE.sub(_repl, text)
    return violations, new_text, n_fixed


def bibliography_check_surname_braces(bib_path: str, fix: bool = False) -> str:
    """Check (and optionally fix) Rule 1: surname must be brace-protected.

    Args:
        bib_path: Path to the ``.bib`` source.
        fix: When False (default), report violations without modifying the
            file.  When True, rewrite the file in place with every
            un-braced surname wrapped in ``{}``.

    Returns:
        Multi-line lint / fix report.
    """
    p = pathlib.Path(bib_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"

    # Validate the file actually parses as .bib (gives an error if the path
    # points to something unrelated).
    try:
        read_bib_file(p)
    except Exception as exc:
        return f"Error: cannot parse {p} as .bib: {exc}"

    # Encoding fallback: many lab .bib files predate the UTF-8 convention
    # and contain cp932 mojibake (broken Japanese journal names from
    # Word -> BibTeX copy-paste).  Try UTF-8 first, then fall back to
    # cp932 / latin-1 so we can at least lint the surname format.
    text: str | None = None
    used_encoding = "utf-8"
    for enc in ("utf-8", "cp932", "latin-1"):
        try:
            text = p.read_text(encoding=enc)
            used_encoding = enc
            break
        except UnicodeDecodeError:
            continue
    if text is None:
        return f"Error: cannot decode {p} (tried utf-8, cp932, latin-1)"
    violations, rewritten, n_fixed = _check_text(text)

    n_entries_with_author = sum(1 for _ in _AUTHOR_RE.finditer(text))
    n_violations = len(violations)

    lines: list[str] = [
        f"bibliography_check_surname_braces: {p}",
        f"  author fields total : {n_entries_with_author}",
        f"  un-braced surnames  : {n_violations}",
    ]

    if n_violations == 0:
        lines.append("OK -- every surname is brace-protected (Rule 1 compliant).")
        return "\n".join(lines)

    # Report up to 50 violations to keep output reasonable.
    SHOW = 50
    lines.append(f"FAIL -- {n_violations} entry/entries violate Rule 1"
                 f" (showing up to {SHOW}):")
    for key, value in violations[:SHOW]:
        lines.append(f"  [HIGH] {key}: {value}")
    if n_violations > SHOW:
        lines.append(f"  ... and {n_violations - SHOW} more")

    if fix:
        p.write_text(rewritten, encoding=used_encoding)
        lines.append("")
        lines.append(f"FIXED -- rewrote {n_fixed} author field(s) in place.")
        lines.append("Re-run with fix=False to confirm 0 violations remain.")
    else:
        lines.append("")
        lines.append("To rewrite the file in place, call with fix=True.")

    return "\n".join(lines)
