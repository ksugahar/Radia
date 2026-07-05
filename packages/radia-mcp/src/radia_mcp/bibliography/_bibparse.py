"""Minimal BibTeX parser / serializer for the bibliography sub-skill.

We deliberately don't depend on ``bibtexparser`` — its installation is
heavyweight and its API churn has caused issues. The grammar handled
here covers:

    - ``@type{key, field = value, …}``
    - balanced-brace values, ``"..."`` values, and bare values
    - string concatenation with ``#``
    - line comments starting with ``%`` (TeX-style)

What it does NOT handle (deliberate omissions):

    - ``@string`` / ``@preamble`` definitions (preserved verbatim as
      ``BibEntry(kind='@string'|'@preamble', body=…)``)
    - cross-references (kept as a normal field; the user can opt in)
    - URL escape semantics — we treat fields opaquely
"""
from __future__ import annotations

import dataclasses
import pathlib
import re


@dataclasses.dataclass
class BibEntry:
    """A single ``@type{key, ...}`` entry."""
    kind: str             # e.g. "article", "inproceedings", "@string"
    key: str              # citation key; "" for @string/@preamble
    fields: dict[str, str]
    raw_body: str = ""    # original body for @string/@preamble


_ENTRY_HEAD = re.compile(r"@(?P<kind>\w+)\s*[\{(]\s*", re.IGNORECASE)


def _skip_ws_and_comments(s: str, i: int) -> int:
    while i < len(s):
        if s[i].isspace():
            i += 1
        elif s[i] == "%":
            # Line comment until newline.
            nl = s.find("\n", i)
            i = len(s) if nl == -1 else nl + 1
        else:
            break
    return i


def _read_value(s: str, i: int) -> tuple[str, int]:
    """Read a single value (braced / quoted / bare) starting at ``s[i]``."""
    i = _skip_ws_and_comments(s, i)
    if i >= len(s):
        return "", i
    parts: list[str] = []
    while True:
        if i >= len(s):
            break
        c = s[i]
        if c == "{":
            depth = 1
            j = i + 1
            while j < len(s) and depth > 0:
                if s[j] == "{":
                    depth += 1
                elif s[j] == "}":
                    depth -= 1
                j += 1
            parts.append(s[i + 1:j - 1])
            i = j
        elif c == '"':
            j = i + 1
            depth = 0
            while j < len(s):
                if s[j] == "{":
                    depth += 1
                elif s[j] == "}":
                    depth -= 1
                elif s[j] == '"' and depth == 0:
                    break
                j += 1
            parts.append(s[i + 1:j])
            i = j + 1
        else:
            # Bare token (number, single word, or @string reference).
            j = i
            while j < len(s) and s[j] not in ",}#)\n":
                j += 1
            parts.append(s[i:j].strip())
            i = j
        i = _skip_ws_and_comments(s, i)
        if i < len(s) and s[i] == "#":
            i = _skip_ws_and_comments(s, i + 1)
            continue
        break
    return "".join(parts), i


def parse_bib(text: str) -> list[BibEntry]:
    """Parse a BibTeX source into a list of ``BibEntry``."""
    entries: list[BibEntry] = []
    i = 0
    while i < len(text):
        i = _skip_ws_and_comments(text, i)
        if i >= len(text):
            break
        m = _ENTRY_HEAD.match(text, i)
        if not m:
            # Lone text outside an entry — skip a char.
            i += 1
            continue
        kind = m.group("kind").lower()
        i = m.end()
        if kind in ("string", "preamble", "comment"):
            # Treat as opaque body; consume until the matching closing brace.
            depth = 1
            j = i
            while j < len(text) and depth > 0:
                if text[j] == "{":
                    depth += 1
                elif text[j] == "}":
                    depth -= 1
                j += 1
            entries.append(BibEntry(kind=f"@{kind}", key="",
                                    fields={}, raw_body=text[i:j - 1]))
            i = j
            continue
        # Normal entry: key , field = value , ...
        # Read key (up to ',').
        j = i
        while j < len(text) and text[j] not in ",}":
            j += 1
        key = text[i:j].strip()
        fields: dict[str, str] = {}
        if j < len(text) and text[j] == ",":
            j += 1
        while True:
            j = _skip_ws_and_comments(text, j)
            if j >= len(text) or text[j] == "}":
                break
            # Read field name until '='.
            k = j
            while k < len(text) and text[k] != "=":
                k += 1
            field = text[j:k].strip().lower()
            if not field:
                break
            j = k + 1
            value, j = _read_value(text, j)
            fields[field] = value
            j = _skip_ws_and_comments(text, j)
            if j < len(text) and text[j] == ",":
                j += 1
        # Consume closing brace.
        if j < len(text) and text[j] == "}":
            j += 1
        entries.append(BibEntry(kind=kind, key=key, fields=fields))
        i = j
    return entries


def write_bib(entries: list[BibEntry], indent: str = "  ") -> str:
    """Serialize entries back to .bib text. Pretty but deterministic."""
    out: list[str] = []
    for e in entries:
        if e.kind.startswith("@"):
            out.append(f"{e.kind}{{{e.raw_body}}}\n")
            continue
        out.append(f"@{e.kind}{{{e.key},\n")
        for field, value in e.fields.items():
            # Wrap value in braces; escape inner braces is the user's problem.
            out.append(f"{indent}{field} = {{{value}}},\n")
        if out[-1].endswith(",\n"):
            out[-1] = out[-1][:-2] + "\n"  # drop trailing comma on last field
        out.append("}\n\n")
    return "".join(out)


def read_bib_file(path: str | pathlib.Path) -> list[BibEntry]:
    p = pathlib.Path(path)
    text = p.read_text(encoding="utf-8", errors="replace")
    return parse_bib(text)


# Citation-key handling -----------------------------------------------------

def first_author_lastname(authors_field: str) -> str:
    """Return the first author's lastname from a BibTeX ``author`` field.

    Handles "Last, First and Last2, First2" and "First Last and First2 Last2".
    Returns lowercase ASCII when possible.
    """
    if not authors_field:
        return "unknown"
    first = authors_field.split(" and ")[0].strip()
    if "," in first:
        last = first.split(",", 1)[0].strip()
    else:
        # Take the last whitespace-separated token.
        last = first.split()[-1] if first.split() else "unknown"
    # Strip non-letters for the cite-key alphabet.
    return re.sub(r"[^A-Za-z]", "", last).lower() or "unknown"


def first_title_word(title_field: str) -> str:
    """First content word of a title, stop-word-stripped."""
    stop = {
        "a", "an", "the", "on", "in", "of", "and", "or", "for", "with",
        "to", "by", "from", "into", "via", "using", "based",
    }
    for raw in re.split(r"[^A-Za-z]+", title_field):
        if not raw:
            continue
        low = raw.lower()
        if low in stop or len(low) < 3:
            continue
        return low
    return "untitled"


_LAB_KEY_RE = re.compile(r"^[a-z]+\d{4}[a-z]+$")


def is_lab_style_key(key: str) -> bool:
    """Return True if ``key`` matches the lab convention ``<author><year><word>``.

    Lab convention (from public-safe curated corpus): all lowercase, format
    ``lastname + 4-digit year + semantic keyword``, e.g. ``bobbio1997play``,
    ``sugahara2026binput``. The 3rd word is hand-picked, not the title's
    first word — so this checker validates the SHAPE only.
    """
    return bool(_LAB_KEY_RE.match(key or ""))


def make_cite_key(entry: BibEntry, style: str = "author_year_word",
                  keyword_override: str | None = None) -> str:
    """Generate a deterministic cite-key from entry fields.

    Styles:
        - ``author_year_word``: ``sugahara2025kelvin`` (lab default)
        - ``Author_year``: ``Sugahara2025``

    The 3rd word defaults to the title's first content word
    (algorithmically derived). The lab convention prefers semantic
    keywords like ``play`` or ``preisach``; pass ``keyword_override`` to
    use a hand-picked keyword instead.
    """
    f = entry.fields
    author = first_author_lastname(f.get("author", "") or f.get("editor", ""))
    year = (f.get("year", "") or "").strip()[:4] or "nodate"
    word = keyword_override or first_title_word(f.get("title", ""))
    word = re.sub(r"[^a-z]", "", word.lower()) or first_title_word(f.get("title", ""))
    if style == "author_year_word":
        return f"{author}{year}{word}"
    if style == "Author_year":
        return author.capitalize() + year
    return f"{author}{year}{word}"
