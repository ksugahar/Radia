"""Terminology normalization helpers for paper writing.

The existing notation-variant lints report suspicious spelling drift.
This module adds a conservative, TeX-friendly normalizer for cases where
the preferred term is already known, such as using ``立方体`` instead of
``cube`` in Japanese prose while preserving English captions and file paths.
"""
from __future__ import annotations

import json
import pathlib
import re
from dataclasses import dataclass


@dataclass(frozen=True)
class _TermRule:
    variant: str
    preferred: str
    note: str = ""


_DEFAULT_RULES: tuple[_TermRule, ...] = (
    _TermRule(
        "cube",
        "立方体",
        "Japanese technical prose should use 立方体; keep English captions as cube.",
    ),
)


_CJK_RE = re.compile(r"[\u3040-\u30ff\u3400-\u9fff\uff00-\uffef]")


def _has_cjk_near(text: str, start: int, end: int, window: int) -> bool:
    left = max(0, start - window)
    right = min(len(text), end + window)
    return bool(_CJK_RE.search(text[left:start]) or _CJK_RE.search(text[end:right]))


def _parse_rules(rules: str) -> list[_TermRule]:
    raw = (rules or "").strip()
    if not raw:
        return list(_DEFAULT_RULES)

    # JSON forms accepted by MCP callers:
    #   {"cube": "立方体"}
    #   [{"variant": "cube", "preferred": "立方体", "note": "..."}]
    if raw[0] in "[{":
        data = json.loads(raw)
        if isinstance(data, dict):
            return [_TermRule(str(k), str(v)) for k, v in data.items()]
        if isinstance(data, list):
            parsed: list[_TermRule] = []
            for item in data:
                if isinstance(item, dict):
                    parsed.append(
                        _TermRule(
                            str(item["variant"]),
                            str(item["preferred"]),
                            str(item.get("note", "")),
                        )
                    )
                elif isinstance(item, (list, tuple)) and len(item) >= 2:
                    parsed.append(_TermRule(str(item[0]), str(item[1])))
                else:
                    raise ValueError(f"unsupported terminology rule: {item!r}")
            return parsed
        raise ValueError("rules JSON must be an object or a list")

    parsed = []
    for line in raw.splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        if "=>" in line:
            variant, preferred = line.split("=>", 1)
        elif "=" in line:
            variant, preferred = line.split("=", 1)
        else:
            raise ValueError(
                "terminology rules must use 'variant=>preferred' or JSON"
            )
        parsed.append(_TermRule(variant.strip(), preferred.strip()))
    return parsed


def _term_pattern(variant: str) -> re.Pattern[str]:
    escaped = re.escape(variant)
    if re.fullmatch(r"[A-Za-z0-9_ -]+", variant):
        return re.compile(r"(?<![A-Za-z0-9_\\])" + escaped + r"(?![A-Za-z0-9_])")
    return re.compile(escaped)


def _read_text_with_fallback(path: pathlib.Path, encoding: str) -> tuple[str, str]:
    candidates = [encoding, "utf-8-sig", "utf-8", "cp932"]
    seen: set[str] = set()
    for enc in candidates:
        if enc in seen:
            continue
        seen.add(enc)
        try:
            return path.read_text(encoding=enc), enc
        except UnicodeDecodeError:
            continue
    return path.read_text(encoding=encoding, errors="replace"), encoding


def paper_writing_normalize_terminology(
    text: str,
    rules: str = "",
    japanese_context_only: bool = True,
    context_window: int = 8,
) -> dict:
    """Normalize known terminology variants in paper text.

    Args:
        text: TeX or plain text to normalize.
        rules: Replacement rules. Empty uses the lab default
            ``cube=>立方体``. Accepted forms are line-based
            ``variant=>preferred`` or JSON dict/list.
        japanese_context_only: When true, replace Latin variants only
            if CJK text appears nearby. This preserves English captions,
            bibliography entries, labels, and graphics file paths.
        context_window: Number of characters used for the CJK-nearby guard.

    Returns:
        dict containing the normalized text and replacement report.
    """
    parsed_rules = _parse_rules(rules)
    normalized = text
    reports: list[dict] = []
    total = 0

    for rule in parsed_rules:
        pattern = _term_pattern(rule.variant)
        examples: list[str] = []
        count = 0

        def repl(match: re.Match[str]) -> str:
            nonlocal count
            if japanese_context_only and re.fullmatch(r"[A-Za-z0-9_ -]+", rule.variant):
                if not _has_cjk_near(normalized, match.start(), match.end(), context_window):
                    return match.group(0)
            count += 1
            if len(examples) < 5:
                cs = max(0, match.start() - 20)
                ce = min(len(normalized), match.end() + 20)
                examples.append(normalized[cs:ce].replace("\n", " "))
            return rule.preferred

        normalized = pattern.sub(repl, normalized)
        total += count
        reports.append(
            {
                "variant": rule.variant,
                "preferred": rule.preferred,
                "count": count,
                "note": rule.note,
                "examples": examples,
            }
        )

    return {
        "changed": normalized != text,
        "n_replacements": total,
        "rules": [
            {
                "variant": r.variant,
                "preferred": r.preferred,
                "note": r.note,
            }
            for r in parsed_rules
        ],
        "replacements": reports,
        "normalized_text": normalized,
        "policy": (
            "Latin variants are replaced only in nearby Japanese context by default; "
            "English captions, paths, labels, and references are preserved."
        ),
    }


def paper_writing_normalize_terminology_file(
    tex_path: str,
    rules: str = "",
    dry_run: bool = True,
    encoding: str = "utf-8",
    japanese_context_only: bool = True,
    context_window: int = 8,
) -> dict:
    """Normalize known terminology variants in a TeX/text file.

    ``dry_run`` defaults to true. Set ``dry_run=False`` only after
    reviewing the returned replacement counts/examples.
    """
    path = pathlib.Path(tex_path)
    if not path.exists():
        return {"error": f"file not found: {tex_path}"}

    text, used_encoding = _read_text_with_fallback(path, encoding)
    result = paper_writing_normalize_terminology(
        text,
        rules=rules,
        japanese_context_only=japanese_context_only,
        context_window=context_window,
    )
    if result["changed"] and not dry_run:
        path.write_text(result["normalized_text"], encoding=used_encoding)

    return {
        "file": str(path),
        "encoding": used_encoding,
        "dry_run": dry_run,
        "changed": result["changed"],
        "n_replacements": result["n_replacements"],
        "rules": result["rules"],
        "replacements": result["replacements"],
        "policy": result["policy"],
    }
