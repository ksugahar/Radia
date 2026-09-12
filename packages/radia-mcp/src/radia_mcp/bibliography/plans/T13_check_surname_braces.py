"""Check local surname protection without rewriting unrelated BibTeX source.

This is syntax, not author identity verification or a guarantee about every
bibliography style. Protected literal names are preserved; author macros and
concatenations require explicit resolution before automatic editing.
"""
from __future__ import annotations

import pathlib

from .._bibparse import _split_name_parts, parse_bib
from .._source_edit import literal_value, read_source, write_source_edits


def _protected(value: str) -> bool:
    value = value.strip()
    if not value.startswith("{"):
        return False
    depth = 0
    i = 0
    while i < len(value):
        if value[i] == "\\":
            i += 2
            continue
        if value[i] == "{":
            depth += 1
        elif value[i] == "}":
            depth -= 1
            if depth == 0:
                return i == len(value) - 1
        i += 1
    return False


def _surname_span(name: str) -> tuple[int, int]:
    start = len(name) - len(name.lstrip())
    end = len(name.rstrip())
    body = name[start:end]
    commas = _split_name_parts(body, r",")
    if len(commas) > 1:
        return start, start + len(commas[0])
    tokens = _split_name_parts(body, r"\s+")
    particle = next((i for i, token in enumerate(tokens[:-1])
                     if token and token[0].islower()), len(tokens) - 1)
    position = 0
    for token in tokens[:particle]:
        position = body.index(token, position) + len(token)
    position = body.index(tokens[particle], position)
    return start + position, end


def _is_surname_braced(name: str) -> bool:
    if not name.strip() or name.strip().casefold() == "others" or _protected(name):
        return True
    start, end = _surname_span(name)
    return _protected(name[start:end])


def _wrap_surname(name: str) -> str:
    if _is_surname_braced(name):
        return name
    start, end = _surname_span(name)
    return name[:start] + "{" + name[start:end] + "}" + name[end:]


def _split_authors(value: str) -> list[str]:
    return _split_name_parts(value, r"\s+and\s+")


def _source_changes(text: str):
    violations, edits, unresolved = [], [], []
    total = 0
    for entry in parse_bib(text):
        if "author" not in entry.fields:
            continue
        total += 1
        start, end = entry.field_spans["author"]
        expression = text[start:end]
        value = literal_value(expression)
        if value is None:
            unresolved.append(entry.key)
            continue
        position = 0
        changed = False
        for author in _split_authors(value):
            offset = value.index(author, position)
            position = offset + len(author)
            replacement = _wrap_surname(author)
            if replacement != author:
                edits.append((start + 1 + offset, start + 1 + position, replacement))
                changed = True
        if changed:
            violations.append((entry.key, value))
    return violations, edits, unresolved, total


def _check_text(text: str) -> tuple[list[tuple[str, str]], str, int]:
    violations, edits, unresolved, _ = _source_changes(text)
    if unresolved:
        raise ValueError("author expressions require explicit resolution: " + ", ".join(unresolved))
    for start, end, replacement in reversed(sorted(edits)):
        text = text[:start] + replacement + text[end:]
    return violations, text, len(violations)


def bibliography_check_surname_braces(bib_path: str, fix: bool = False) -> str:
    """Check local surname protection; fix literal author fields only on request."""
    path = pathlib.Path(bib_path).absolute()
    try:
        original, text, entries = read_source(path)
        if not entries:
            raise ValueError("bibliography contains no entries")
        violations, edits, unresolved, total = _source_changes(text)
    except (OSError, UnicodeError, ValueError) as exc:
        return f"Error: cannot inspect bibliography: {exc}"
    lines = [f"bibliography_check_surname_braces: {path}",
             f"  author fields total : {total}",
             f"  un-braced surnames  : {len(violations)}"]
    if unresolved:
        lines.append("UNAVAILABLE -- author macros/concatenations require explicit resolution: "
                     + ", ".join(unresolved))
        if fix:
            return "\n".join(lines + ["No changes written."])
    elif not violations:
        return "\n".join(lines + ["OK -- every inspected surname is brace-protected (local Rule 1)."])
    for key, value in violations[:50]:
        lines.append(f"  [HIGH] {key}: {value}")
    if fix:
        try:
            write_source_edits(path, original, text, edits)
        except (OSError, UnicodeError, ValueError) as exc:
            return "\n".join(lines + [f"Error: cannot apply surname edits: {exc}"])
        lines.append(f"FIXED -- rewrote {len(violations)} author field(s) in place.")
    elif violations:
        lines.append("FAIL -- local surname protection convention; call with fix=True to apply.")
    return "\n".join(lines)
