"""Tier 3 — bibliography_cite_validation.

Cross-check ``\\cite{key}`` references in .tex against keys defined in
.bib. The 2 failure modes:
    - cite-key used in .tex but not in .bib (compilation error)
    - bib key defined but never cited (bloat / forgotten reference)
"""
from __future__ import annotations

import pathlib
import re

from .._bibparse import read_bib_file


_CITE_CMD_RE = re.compile(r"\\(?:no)?cite[a-zA-Z]*\*?(?:\[[^\]]*\])?\{([^}]+)\}")


def _collect_cited_keys(tex_path: pathlib.Path) -> set[str]:
    """Walk a .tex tree (single file or via \\input/\\include) for \\cite keys."""
    seen_files: set[str] = set()
    keys: set[str] = set()

    def visit(path: pathlib.Path) -> None:
        if not path.exists() or str(path) in seen_files:
            return
        seen_files.add(str(path))
        try:
            text = path.read_text(encoding="utf-8", errors="replace")
        except OSError:
            return
        for m in _CITE_CMD_RE.finditer(text):
            for k in m.group(1).split(","):
                k = k.strip()
                if k:
                    keys.add(k)
        for inc_match in re.finditer(r"\\(?:input|include|InputIfFileExists)\{([^}]+)\}", text):
            target = inc_match.group(1)
            cand = path.parent / target
            if not cand.suffix:
                cand = cand.with_suffix(".tex")
            visit(cand)

    visit(tex_path)
    return keys


def bibliography_cite_validation(tex_path: str, bib_path: str) -> str:
    """Cross-check ``\\cite{}`` keys in .tex against entries in .bib.

    Args:
        tex_path: Root .tex file (follows ``\\input{}`` / ``\\include{}``).
        bib_path: Path to the .bib source.

    Returns:
        - Keys cited but missing from .bib (HIGH — compile error)
        - Keys in .bib but never cited (LOW — bloat)
    """
    tp = pathlib.Path(tex_path)
    if not tp.is_absolute():
        tp = pathlib.Path.cwd() / tp
    bp = pathlib.Path(bib_path)
    if not bp.is_absolute():
        bp = pathlib.Path.cwd() / bp
    if not tp.exists():
        return f"Error: tex file not found: {tp}"
    if not bp.exists():
        return f"Error: bib file not found: {bp}"

    cited = _collect_cited_keys(tp)
    bib_keys = {e.key for e in read_bib_file(bp) if not e.kind.startswith("@")}

    missing = cited - bib_keys
    uncited = bib_keys - cited

    lines = [f"bibliography_cite_validation: tex={tp}, bib={bp}",
             f"  cited keys: {len(cited)}  bib keys: {len(bib_keys)}"]
    if missing:
        lines.append(f"[HIGH] {len(missing)} key(s) cited in .tex but missing from .bib:")
        lines.extend(f"  - {k}" for k in sorted(missing))
    if uncited:
        lines.append(f"[LOW] {len(uncited)} key(s) defined in .bib but never cited:")
        lines.extend(f"  - {k}" for k in sorted(uncited))
    if not missing and not uncited:
        lines.append("PASS — cite-key set matches bib-key set")
    return "\n".join(lines)
