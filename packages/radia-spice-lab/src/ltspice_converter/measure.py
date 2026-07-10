"""Public LTspice .measure log parsing helpers.

The functions in this module intentionally parse text logs only.  They turn
LTspice scalar `.measure` output into a small, public-safe table schema that
MCP tools and validation artifacts can reuse without copying private RAW data.
"""
from __future__ import annotations

import math
import re
from typing import Any


_NUMBER = r"[-+]?(?:(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?|inf|nan)"
_SCALAR = rf"{_NUMBER}[A-Za-zµμ]*"
_WINDOW_RE = re.compile(rf"\bFROM\s+(?P<from>{_SCALAR})\s+TO\s+(?P<to>{_SCALAR})\b", re.IGNORECASE)


def parse_spice_scalar(value: str) -> float | None:
    """Parse a SPICE engineering scalar such as ``1k`` or ``22u``.

    Unknown suffixes return ``None`` instead of guessing.
    """
    text = str(value or "").strip().strip("{}")
    match = re.match(
        r"^\s*(?P<number>[-+]?(?:\d+(?:\.\d*)?|\.\d+)(?:[eE][-+]?\d+)?|inf|nan)(?P<suffix>[A-Za-zµμ]*)",
        text,
        re.IGNORECASE,
    )
    if not match:
        return None
    number_text = match.group("number")
    base = float(number_text)
    suffix = match.group("suffix").strip()
    if not suffix:
        return base
    suffix_l = suffix.lower().replace("μ", "u").replace("µ", "u")
    if suffix_l.startswith("meg"):
        return base * 1e6
    scale = {
        "f": 1e-15,
        "p": 1e-12,
        "n": 1e-9,
        "u": 1e-6,
        "m": 1e-3,
        "k": 1e3,
        "g": 1e9,
        "t": 1e12,
    }.get(suffix_l[:1])
    return base * scale if scale is not None else None


def _as_float(value: str) -> float:
    parsed = parse_spice_scalar(value)
    if parsed is None:
        return float(value)
    return parsed


def parse_ltspice_measure_lines(lines: list[str]) -> list[dict[str, Any]]:
    """Parse scalar LTspice `.measure` result lines from log text lines.

    Supported public schema rows include:

    - scalar value rows: ``gain: mag(v(out))=0.707``
    - event rows: ``trise: v(out)=0.5 AT 0.001``
    - AC rows: ``gain: mag(v(out))=(-3.01dB,-45deg) at 1k``
    - window annotations: ``FROM ... TO ...``
    """
    measures: list[dict[str, Any]] = []
    at_re = re.compile(
        rf"^(?P<name>[A-Za-z_]\w*)\s*:\s*(?P<expr>.*?)\s*=\s*(?P<target>{_NUMBER})\s+AT\s+(?P<at>{_NUMBER})(?P<rest>.*)$",
        re.IGNORECASE,
    )
    ac_db_re = re.compile(
        rf"^(?P<name>[A-Za-z_]\w*)\s*:\s*(?P<expr>.*?)\s*=\s*"
        rf"\(\s*(?P<value>{_NUMBER})\s*dB\s*,\s*(?P<phase>{_NUMBER})[^)]*\)"
        rf"(?P<rest>.*)$",
        re.IGNORECASE,
    )
    value_re = re.compile(
        rf"^(?P<name>[A-Za-z_]\w*)\s*:\s*(?P<expr>.*?)\s*=\s*(?P<value>{_NUMBER})(?P<rest>.*)$",
        re.IGNORECASE,
    )
    scalar_re = re.compile(
        rf"^(?P<name>[A-Za-z_]\w*)\s*:\s*(?P<kind>time|freq|frequency)\s*=\s*(?P<value>{_NUMBER})(?P<rest>.*)$",
        re.IGNORECASE,
    )
    for raw_line in lines:
        line = str(raw_line).strip()
        if not line:
            continue
        at_match = at_re.match(line)
        if at_match:
            expr = at_match.group("expr").strip()
            expr_lower = expr.lower()
            measures.append({
                "name": at_match.group("name"),
                "kind": "time_at" if expr_lower in {"time", "freq", "frequency"} else "value_at",
                "value": float(at_match.group("target")),
                "at": float(at_match.group("at")),
                "target_value": float(at_match.group("target")),
                "expression": expr,
                "line": line,
            })
            continue
        ac_match = ac_db_re.match(line)
        if ac_match:
            rest = ac_match.group("rest") or ""
            at_in_rest = re.search(rf"\bat\s+(?P<at>{_SCALAR})\b", rest, re.IGNORECASE)
            window = _WINDOW_RE.search(rest)
            item = {
                "name": ac_match.group("name"),
                "kind": "ac_db_phase_at" if at_in_rest else "ac_db_phase",
                "value": float(ac_match.group("value")),
                "unit": "dB",
                "phase_deg": float(ac_match.group("phase")),
                "expression": ac_match.group("expr").strip(),
                "line": line,
            }
            if at_in_rest:
                item["at"] = _as_float(at_in_rest.group("at"))
            if window:
                item["from"] = _as_float(window.group("from"))
                item["to"] = _as_float(window.group("to"))
            measures.append(item)
            continue
        match = value_re.match(line) or scalar_re.match(line)
        if not match:
            continue
        item: dict[str, Any] = {
            "name": match.group("name"),
            "value": float(match.group("value")),
            "line": line,
        }
        expr = match.groupdict().get("expr")
        kind = match.groupdict().get("kind")
        if expr is not None:
            item["expression"] = expr.strip()
        if kind is not None:
            item["kind"] = kind.lower()
        rest = match.groupdict().get("rest") or ""
        window = _WINDOW_RE.search(rest)
        if window:
            item["from"] = _as_float(window.group("from"))
            item["to"] = _as_float(window.group("to"))
        measures.append(item)
    return measures


def parse_ltspice_step_lines(lines: list[str]) -> list[dict[str, Any]]:
    """Parse LTspice `.step` assignment lines from log text lines."""
    steps: list[dict[str, Any]] = []
    step_re = re.compile(r"^\s*\.step\s+(?P<body>.+)$", re.IGNORECASE)
    assign_re = re.compile(r"(?P<name>[A-Za-z_]\w*)\s*=\s*(?P<value>[^\s]+)")
    for raw_line in lines:
        line = str(raw_line).strip()
        match = step_re.match(line)
        if not match:
            continue
        body = match.group("body").strip()
        assignments = {
            item.group("name"): item.group("value")
            for item in assign_re.finditer(body)
        }
        numeric_assignments: dict[str, float] = {}
        for key, value in assignments.items():
            parsed = parse_spice_scalar(value)
            if parsed is not None and math.isfinite(parsed):
                numeric_assignments[key] = parsed
        steps.append({
            "line": line,
            "body": body,
            "assignments": assignments,
            "numeric_assignments": numeric_assignments,
        })
    return steps


def summarize_measure_log(log_text: str) -> dict[str, Any]:
    """Summarize LTspice log text as a public-safe measurement evidence table."""
    lines = [line.strip() for line in str(log_text or "").splitlines() if line.strip()]
    measures = parse_ltspice_measure_lines(lines)
    steps = parse_ltspice_step_lines(lines)
    duplicate_names = sorted({
        item["name"]
        for item in measures
        if sum(1 for other in measures if other.get("name") == item.get("name")) > 1
    })
    warnings: list[str] = []
    if not measures:
        warnings.append("no LTspice .measure result rows were parsed")
    if duplicate_names:
        warnings.append(f"duplicate .measure result names: {', '.join(duplicate_names)}")
    return {
        "schema": "radia-spice-lab.measure-log.v1",
        "ok": bool(measures) and not duplicate_names,
        "measure_count": len(measures),
        "step_count": len(steps),
        "measure_names": [str(item.get("name")) for item in measures],
        "duplicate_measure_names": duplicate_names,
        "measures": measures,
        "steps": steps,
        "warnings": warnings,
    }
