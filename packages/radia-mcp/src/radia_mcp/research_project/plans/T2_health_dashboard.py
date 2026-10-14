"""Tier 1 — research_project_health_dashboard.

For each detected artifact in a project, call the matching sub-skill's
``*_health_report`` and aggregate into a single dashboard.
"""
from __future__ import annotations

import pathlib

from .T1_scan import _HEURISTICS


def research_project_health_dashboard(directory: str,
                                      recursive: bool = True) -> str:
    """Run each artifact's health report in one pass.

    Args:
        directory: Project root.
        recursive: descend or not.

    Returns:
        A combined dashboard with per-artifact verdict.
    """
    root = pathlib.Path(directory)
    if not root.is_absolute():
        root = pathlib.Path.cwd() / root
    if not root.is_dir():
        return f"Error: not a directory: {root}"

    # Lazy-import the health reports to avoid cyclic deps.
    from ..tools_handlers import handlers  # type: ignore

    iterator = root.rglob("*") if recursive else root.glob("*")
    files = [p for p in iterator if p.is_file()]

    found: dict[str, list[pathlib.Path]] = {}
    for p in files:
        name = p.name.lower()
        ext = p.suffix.lower()
        kind = next((k for k, fn in _HEURISTICS if fn(name, ext)), None)
        if kind:
            found.setdefault(kind, []).append(p)

    lines = [f"research_project_health_dashboard: {root}"]
    for kind, paths in found.items():
        handler = handlers.get(kind)
        if not handler:
            continue
        lines.append(f"\n=== {kind} ({len(paths)}) ===")
        for path in paths:
            try:
                rep = handler(str(path))
            except Exception as e:
                lines.append(f"  - {path.name}: ERROR ({type(e).__name__}: {e})")
                continue
            # Extract score / verdict lines.
            for line in rep.splitlines():
                line = line.strip()
                if "score:" in line.lower() or "verdict:" in line.lower():
                    lines.append(f"  - {path.relative_to(root)}: {line}")
                    break
    return "\n".join(lines)
