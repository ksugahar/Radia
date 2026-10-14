#!/usr/bin/env python3
"""Regression tests for the LLM usability-trace evidence pipeline."""

import json
import sys
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from tools.analyze_usability_trace import (
    build_review_bundle,
    detect_candidates,
    parse_operation_log_text,
)


def row(
    seq: int,
    millis: int,
    event: str,
    details: str = "",
    latex: str = "x",
    extra: str = "",
) -> str:
    stamp = f"2026-08-21 12:00:{millis // 1000:02d}.{millis % 1000:03d}"
    fields = [str(seq), stamp, event]
    if details:
        fields.append(details)
    fields.extend([
        "caret=root:1", "selection=<none>", "dirty=1",
        f"elapsed_ms={millis}", f"delta_ms={millis if seq == 1 else 100}",
        "focus=canvas", "input_style=M", "alignment=-1",
        "zoom_percent=175", "equation_mode=equation",
        "shortcut_prefix=-", f"latex={latex}",
    ])
    if extra:
        fields.insert(-1, extra)
    return "\t".join(fields)


def main() -> None:
    text = "\n".join([
        "\ufeff# Eqnedit64 operation log v2",
        row(1, 0, "debug.start", "Eqnedit64 3.0.7 path=C:\\\\redacted"),
        row(2, 100, "template.frac", latex="\\\\frac{}{}"),
        row(3, 500, "edit.undo", latex="x"),
        row(4, 600, "caret.right", "shift=0 changed=1"),
        row(5, 750, "caret.left", "shift=0 changed=1"),
        row(6, 900, "shortcut.prefix.invalid", "prefix=T key=90 shift=0"),
        row(7, 1000, "edit.delete.no_change"),
        row(8, 1100, "edit.backspace", latex=""),
        row(9, 1200, "edit.delete", latex=""),
        row(10, 1300, "user.marker", "F12: problem noticed here"),
        row(11, 1400, "debug.stop"),
    ]) + "\n"
    version, events = parse_operation_log_text(text)
    assert version == 2
    assert len(events) == 11
    assert events[1].state["latex"] == "\\frac{}{}"
    assert events[2].elapsed_ms == 500

    detectors = {candidate.detector for candidate in detect_candidates(events)}
    required = {
        "immediate_undo", "navigation_reversal", "invalid_shortcut",
        "ineffective_command", "correction_burst", "explicit_user_marker",
    }
    assert required <= detectors, (required, detectors)

    structure = build_review_bundle(Path("operation-test.log"), version, events, "structure")
    assert structure["schema"] == "eqnedit64.usability-review-bundle.v1"
    assert structure["summary"]["event_count"] == 11
    assert structure["summary"]["explicit_marker_count"] == 1
    contexts = [event for item in structure["candidates"] for event in item["context"]]
    latex_values = [event["state"].get("latex") for event in contexts]
    assert any(isinstance(value, dict) and value.get("redacted") for value in latex_values)
    assert not any(value == "\\frac{}{}" for value in latex_values)

    full = build_review_bundle(Path("operation-test.log"), version, events, "full")
    full_contexts = [event for item in full["candidates"] for event in item["context"]]
    assert any(event["state"].get("latex") == "\\frac{}{}" for event in full_contexts)

    v1 = "\n".join([
        "# Eqnedit64 operation log v1",
        "1\t2026-08-21 23:59:59.900\ttext.insert\tx\tcaret=root:1"
        "\tselection=<none>\tdirty=1\tlatex=x",
        "2\t2026-08-22 00:00:00.100\tedit.undo\tcaret=root:0"
        "\tselection=<none>\tdirty=1\tlatex=",
    ])
    old_version, old_events = parse_operation_log_text(v1)
    assert old_version == 1
    assert old_events[1].elapsed_ms == 200
    assert "immediate_undo" in {
        candidate.detector for candidate in detect_candidates(old_events)
    }

    preference_path = PROJECT_ROOT / "docs" / "USABILITY_PREFERENCES.jsonl"
    preferences = [
        json.loads(line) for line in preference_path.read_text(encoding="utf-8").splitlines()
        if line.strip()
    ]
    ids = [item["decision_id"] for item in preferences]
    assert len(ids) == len(set(ids)) and len(ids) >= 4
    for item in preferences:
        assert item["schema"] == "eqnedit64.usability-preference.v1"
        assert item["status"] in {"accepted", "rejected", "superseded"}
        assert item["decision"] and item["authority"] and item["regressions"]

    print(
        "ok    usability trace v1/v2 parsing, privacy, 6 friction detectors, "
        f"and {len(preferences)} preference decisions"
    )


if __name__ == "__main__":
    main()
