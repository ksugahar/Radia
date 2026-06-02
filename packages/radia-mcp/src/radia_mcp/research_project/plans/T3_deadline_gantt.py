"""Tier 2 — research_project_deadline_gantt.

Read deadlines from a YAML / JSON / .md plan and render an ASCII Gantt
chart showing days remaining per task.
"""
from __future__ import annotations

import datetime
import json
import pathlib
import re


_DATE_RE = re.compile(r"(20\d{2})[-/](\d{1,2})[-/](\d{1,2})")


def _parse_deadlines(path: pathlib.Path) -> list[tuple[str, datetime.date]]:
    text = path.read_text(encoding="utf-8", errors="replace")
    out: list[tuple[str, datetime.date]] = []
    # Try JSON / YAML-ish dict-of-deadlines
    if path.suffix.lower() == ".json":
        try:
            data = json.loads(text)
            for k, v in data.items():
                if isinstance(v, str):
                    m = _DATE_RE.search(v)
                    if m:
                        out.append((k, datetime.date(int(m.group(1)),
                                                       int(m.group(2)),
                                                       int(m.group(3)))))
        except (ValueError, AttributeError):
            pass
        return out
    # Markdown / plain text: look for "task: 2026-05-15" lines
    for line in text.splitlines():
        m = _DATE_RE.search(line)
        if not m:
            continue
        head = line[:m.start()].strip(" \t-•*:#|")
        if not head:
            continue
        try:
            d = datetime.date(int(m.group(1)),
                               int(m.group(2)),
                               int(m.group(3)))
        except ValueError:
            continue
        out.append((head[:50], d))
    return out


def research_project_deadline_gantt(plan_path: str,
                                    reference_date: str | None = None) -> str:
    """Render an ASCII Gantt-ish chart of upcoming deadlines.

    Args:
        plan_path: A .json / .md / .txt file with lines like
                   ``task: 2026-05-15`` or a JSON ``{"task": "2026-05-15"}``.
        reference_date: ISO date string ``YYYY-MM-DD`` for "today".
                        Defaults to actual today.

    Returns:
        ASCII chart with relative-to-today positions.
    """
    p = pathlib.Path(plan_path)
    if not p.is_absolute():
        p = pathlib.Path.cwd() / p
    if not p.exists():
        return f"Error: file not found: {p}"

    deadlines = _parse_deadlines(p)
    if not deadlines:
        return f"research_project_deadline_gantt: {p}\n  no deadline entries detected"

    today = datetime.date.today()
    if reference_date:
        try:
            today = datetime.date.fromisoformat(reference_date)
        except ValueError:
            pass

    # Compute span and bin to columns.
    days = [(t, (d - today).days) for t, d in deadlines]
    span = max(d for _, d in days) - min(0, min(d for _, d in days))
    span = max(span, 30)
    cols = 50  # width of the chart

    lines = [f"research_project_deadline_gantt: {p}",
             f"  today: {today.isoformat()}, deadlines: {len(deadlines)}",
             "",
             "  " + ("." * cols) + " (span: {} days)".format(span)]
    for label, dd in sorted(days, key=lambda x: x[1]):
        pos = max(0, min(cols - 1, int(dd / span * cols)))
        bar = " " * pos + "|"
        suffix = "TODAY" if dd == 0 else (f"in {dd}d" if dd > 0 else f"{-dd}d ago")
        lines.append(f"  {bar:50s} {label} ({suffix})")
    return "\n".join(lines)
