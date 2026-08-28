#!/usr/bin/env python3
"""Build a privacy-aware LLM review bundle from an Eqnedit64 operation log.

The deterministic detectors do not claim to decide what feels natural.  They
locate evidence windows where an LLM can compare the event sequence, semantic
editor state, Eqnedit32 reference, and an explicit user marker.  A human still
accepts or rejects every behavioural change.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from dataclasses import asdict, dataclass
from datetime import datetime
from pathlib import Path
from typing import Iterable, Sequence


LOG_VERSION_RE = re.compile(r"Eqnedit64 operation log v(\d+)")
COMMAND_RE = re.compile(r"\\([A-Za-z]+)")
ENVIRONMENT_RE = re.compile(r"\\begin\{([^}]+)\}")


@dataclass(frozen=True)
class TraceEvent:
    seq: int
    timestamp: str
    event: str
    details: str
    state: dict[str, str]
    elapsed_ms: int
    delta_ms: int


@dataclass(frozen=True)
class Candidate:
    detector: str
    confidence: str
    anchor_seq: int
    evidence_seq: list[int]
    reason: str
    llm_question: str


def debug_unescape(text: str) -> str:
    """Reverse Eqnedit64's tab-safe debug_escape representation."""
    output: list[str] = []
    index = 0
    escapes = {"t": "\t", "r": "\r", "n": "\n", "\\": "\\"}
    while index < len(text):
        if text[index] == "\\" and index + 1 < len(text):
            following = text[index + 1]
            if following in escapes:
                output.append(escapes[following])
                index += 2
                continue
        output.append(text[index])
        index += 1
    return "".join(output)


def _timestamp_ms(value: str) -> int:
    parsed = datetime.strptime(value, "%Y-%m-%d %H:%M:%S.%f")
    midnight = parsed.replace(hour=0, minute=0, second=0, microsecond=0)
    return int((parsed - midnight).total_seconds() * 1000)


def parse_operation_log_text(text: str) -> tuple[int, list[TraceEvent]]:
    version = 0
    raw_events: list[tuple[int, str, str, str, dict[str, str]]] = []
    for raw_line in text.lstrip("\ufeff").splitlines():
        if raw_line.startswith("#"):
            match = LOG_VERSION_RE.search(raw_line)
            if match:
                version = int(match.group(1))
            continue
        if not raw_line.strip():
            continue
        fields = raw_line.split("\t")
        if len(fields) < 4:
            raise ValueError(f"malformed operation-log row: {raw_line!r}")
        try:
            seq = int(fields[0])
        except ValueError as exc:
            raise ValueError(f"invalid sequence number: {fields[0]!r}") from exc
        state_start = next(
            (i for i in range(3, len(fields)) if fields[i].startswith("caret=")),
            None,
        )
        if state_start is None:
            raise ValueError(f"operation-log row has no caret state: seq={seq}")
        details = debug_unescape("\t".join(fields[3:state_start]))
        state: dict[str, str] = {}
        for field in fields[state_start:]:
            if "=" not in field:
                continue
            key, value = field.split("=", 1)
            state[key] = debug_unescape(value)
        raw_events.append((seq, fields[1], fields[2], details, state))

    if version not in (1, 2):
        raise ValueError(f"unsupported Eqnedit64 operation-log version: {version}")
    if not raw_events:
        return version, []

    first_timestamp = _timestamp_ms(raw_events[0][1])
    previous_elapsed = 0
    events: list[TraceEvent] = []
    day_offset = 0
    previous_clock = first_timestamp
    for seq, timestamp, event, details, state in raw_events:
        clock = _timestamp_ms(timestamp)
        if clock < previous_clock:
            day_offset += 24 * 60 * 60 * 1000
        previous_clock = clock
        derived_elapsed = day_offset + clock - first_timestamp
        elapsed = int(state.get("elapsed_ms", derived_elapsed))
        delta = int(state.get("delta_ms", max(0, elapsed - previous_elapsed)))
        previous_elapsed = elapsed
        events.append(TraceEvent(
            seq=seq,
            timestamp=timestamp,
            event=event,
            details=details,
            state=state,
            elapsed_ms=elapsed,
            delta_ms=delta,
        ))
    return version, events


def parse_operation_log(path: Path) -> tuple[int, list[TraceEvent]]:
    return parse_operation_log_text(path.read_text(encoding="utf-8-sig"))


def _is_noise(event: TraceEvent) -> bool:
    return (
        event.event.startswith("debug.")
        or event.event.startswith("focus.")
        or event.event == "shortcut.coach"
        or event.event == "shortcut.sequence"
    )


def _previous_meaningful(events: Sequence[TraceEvent], index: int) -> int | None:
    for candidate in range(index - 1, -1, -1):
        if not _is_noise(events[candidate]):
            return candidate
    return None


def detect_candidates(events: Sequence[TraceEvent]) -> list[Candidate]:
    found: list[Candidate] = []
    seen: set[tuple[str, int]] = set()

    def add(candidate: Candidate) -> None:
        key = (candidate.detector, candidate.anchor_seq)
        if key not in seen:
            seen.add(key)
            found.append(candidate)

    opposites = {
        "caret.left": "caret.right", "caret.right": "caret.left",
        "caret.up": "caret.down", "caret.down": "caret.up",
        "caret.next_slot": "caret.prev_slot",
        "caret.prev_slot": "caret.next_slot",
    }
    mutation_prefixes = (
        "text.insert", "template.", "symbol.insert", "edit.paste",
        "edit.cut", "edit.delete", "edit.backspace", "edit.new_line",
        "edit.alignment", "shortcut.prime", "shortcut.double_prime",
        "source.edit",
    )

    for index, event in enumerate(events):
        if event.event == "user.marker":
            start = max(0, index - 8)
            end = min(len(events), index + 3)
            add(Candidate(
                detector="explicit_user_marker",
                confidence="high",
                anchor_seq=event.seq,
                evidence_seq=[item.seq for item in events[start:end]],
                reason="The user explicitly marked that the interaction felt wrong here.",
                llm_question=(
                    "What expectation was violated immediately before the marker, and "
                    "what is the smallest behavioural change that would satisfy it?"
                ),
            ))

        if event.event == "shortcut.prefix.invalid":
            add(Candidate(
                detector="invalid_shortcut",
                confidence="high",
                anchor_seq=event.seq,
                evidence_seq=[item.seq for item in events[max(0, index - 3):index + 2]],
                reason="A remembered or guessed shortcut chord was rejected.",
                llm_question=(
                    "Is this a missing learned chord, an ambiguous prefix, or a discoverability "
                    "problem? Compare the Eqnedit32 shortcut reference before proposing a key."
                ),
            ))

        if event.event.endswith(".no_change"):
            add(Candidate(
                detector="ineffective_command",
                confidence="medium",
                anchor_seq=event.seq,
                evidence_seq=[item.seq for item in events[max(0, index - 3):index + 2]],
                reason="A command completed without changing the structural editor state.",
                llm_question=(
                    "Was no-op behaviour expected at this structural boundary, or should the "
                    "command move, explain why it is disabled, or do something else?"
                ),
            ))

        if event.event == "edit.undo":
            previous = _previous_meaningful(events, index)
            if previous is not None:
                prior = events[previous]
                elapsed = event.elapsed_ms - prior.elapsed_ms
                if prior.event.startswith(mutation_prefixes) and 0 <= elapsed <= 2500:
                    add(Candidate(
                        detector="immediate_undo",
                        confidence="medium",
                        anchor_seq=event.seq,
                        evidence_seq=[item.seq for item in events[max(0, previous - 2):index + 2]],
                        reason=f"A structural change was undone after only {elapsed} ms.",
                        llm_question=(
                            "Did insertion produce an unexpected structure/caret destination, or "
                            "was this an intentional correction? Use surrounding repeats and markers."
                        ),
                    ))

        if event.event in opposites:
            previous = _previous_meaningful(events, index)
            if previous is not None:
                prior = events[previous]
                elapsed = event.elapsed_ms - prior.elapsed_ms
                if prior.event == opposites[event.event] and 0 <= elapsed <= 1500:
                    add(Candidate(
                        detector="navigation_reversal",
                        confidence="low",
                        anchor_seq=event.seq,
                        evidence_seq=[item.seq for item in events[max(0, previous - 2):index + 2]],
                        reason=f"Caret direction was reversed after {elapsed} ms.",
                        llm_question=(
                            "Was the first move an overshoot caused by an unintuitive structural "
                            "boundary, or normal visual inspection? Do not change behaviour from this alone."
                        ),
                    ))

    correction_events = {"edit.undo", "edit.backspace", "edit.delete"}
    for index, event in enumerate(events):
        if event.event not in correction_events:
            continue
        window = [
            item for item in events[:index + 1]
            if item.event in correction_events
            and 0 <= event.elapsed_ms - item.elapsed_ms <= 4000
        ]
        if len(window) >= 3:
            add(Candidate(
                detector="correction_burst",
                confidence="medium",
                anchor_seq=event.seq,
                evidence_seq=[item.seq for item in events[max(0, index - 6):index + 2]],
                reason=f"{len(window)} correction commands occurred within four seconds.",
                llm_question=(
                    "Does this burst reveal an incorrect insertion unit, selection unit, or caret "
                    "destination? Separate ordinary text correction from structural friction."
                ),
            ))

    return sorted(found, key=lambda candidate: (candidate.anchor_seq, candidate.detector))


def _content_summary(value: str) -> dict[str, object]:
    commands = sorted(set(COMMAND_RE.findall(value)))
    environments = sorted(set(ENVIRONMENT_RE.findall(value)))
    return {
        "redacted": True,
        "utf8_bytes": len(value.encode("utf-8")),
        "sha256_12": hashlib.sha256(value.encode("utf-8")).hexdigest()[:12],
        "commands": commands,
        "environments": environments,
    }


def _safe_details(event: TraceEvent) -> bool:
    safe_events = (
        "caret.", "view.", "focus.", "debug.", "mouse.",
        "input.style", "tex.numbering", "shortcut.prefix",
    )
    return event.event.startswith(safe_events) and "path=" not in event.details


def event_for_review(event: TraceEvent, privacy: str) -> dict[str, object]:
    result = asdict(event)
    if privacy == "full":
        return result
    state = dict(event.state)
    for key in ("latex", "selection"):
        if key in state and state[key] not in ("", "<none>"):
            state[key] = _content_summary(state[key])
    result["state"] = state
    if event.details and not _safe_details(event):
        result["details"] = _content_summary(event.details)
    return result


def build_review_bundle(
    source: Path,
    version: int,
    events: Sequence[TraceEvent],
    privacy: str,
) -> dict[str, object]:
    candidates = detect_candidates(events)
    by_seq = {event.seq: event for event in events}
    review_candidates: list[dict[str, object]] = []
    for number, candidate in enumerate(candidates, 1):
        item = asdict(candidate)
        item["candidate_id"] = f"C{number:03d}"
        item["context"] = [
            event_for_review(by_seq[seq], privacy)
            for seq in candidate.evidence_seq if seq in by_seq
        ]
        review_candidates.append(item)

    event_counts: dict[str, int] = {}
    for event in events:
        event_counts[event.event] = event_counts.get(event.event, 0) + 1
    marker_count = event_counts.get("user.marker", 0)
    duration = events[-1].elapsed_ms if events else 0
    return {
        "schema": "eqnedit64.usability-review-bundle.v1",
        "source": source.name,
        "operation_log_version": version,
        "privacy": privacy,
        "summary": {
            "event_count": len(events),
            "duration_ms": duration,
            "explicit_marker_count": marker_count,
            "candidate_count": len(candidates),
            "event_counts": dict(sorted(event_counts.items())),
        },
        "reference_files": [
            "docs/SHORTCUTS.md",
            "docs/GUI_SPEC.md",
            "docs/USABILITY_PREFERENCES.jsonl",
        ],
        "review_contract": {
            "principle": (
                "A detector identifies evidence, not a usability defect. The LLM proposes; "
                "a human accepts; an accepted behaviour becomes an executable regression test."
            ),
            "required_finding_fields": [
                "candidate_id", "verdict", "hypothesis", "evidence_seq",
                "legacy_evidence", "proposed_minimal_change", "regression_test",
                "human_question",
            ],
            "verdict_values": ["friction", "expected", "uncertain"],
        },
        "llm_instruction": (
            "Review each candidate against the accepted compatibility behaviour and current GUI "
            "contract. Cite exact event sequence numbers. Never infer equation content hidden by "
            "structure privacy. Prefer the smallest reversible interaction change. Do not alter "
            "shortcuts or behaviour without a human preference decision. For every accepted "
            "finding, specify a background replay regression that does not control the desktop."
        ),
        "candidates": review_candidates,
    }


def write_bundle(input_path: Path, output_path: Path, privacy: str) -> dict[str, object]:
    version, events = parse_operation_log(input_path)
    bundle = build_review_bundle(input_path, version, events, privacy)
    output_path.write_text(
        json.dumps(bundle, ensure_ascii=False, indent=2) + "\n",
        encoding="utf-8",
    )
    return bundle


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Build an LLM usability-review bundle from an Eqnedit64 operation log."
    )
    parser.add_argument("operation_log", type=Path)
    parser.add_argument("--output", "-o", type=Path, required=True)
    parser.add_argument(
        "--privacy", choices=("structure", "full"), default="structure",
        help="structure redacts equation text/details; full keeps them for local review",
    )
    return parser


def main(argv: Iterable[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    bundle = write_bundle(args.operation_log, args.output, args.privacy)
    summary = bundle["summary"]
    print(
        "PASS: usability bundle "
        f"events={summary['event_count']} candidates={summary['candidate_count']} "
        f"markers={summary['explicit_marker_count']} privacy={args.privacy}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
