#!/usr/bin/env python3
"""Create a lossless, static inventory of an ACIS SAT file.

This utility intentionally does not classify materials or infer electrical
windings.  SAT name attributes are emitted as best-effort display labels plus
hex-encoded source bytes, so the audit remains useful for legacy encodings.
"""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


BODY_RE = re.compile(rb"^-(\d+)\s+body\s+\$(-?\d+)\b")
NAME_RE = re.compile(
    rb"^-(\d+)\s+name_attrib-gen-attrib\b.*?@(\d+)\s+(.*?)\s+#\s*$"
)


def _display_name(value: bytes) -> str:
    for encoding in ("utf-8", "cp932"):
        try:
            return value.decode(encoding)
        except UnicodeDecodeError:
            pass
    return value.decode("latin-1")


def _parse_header(lines: list[bytes]) -> dict[str, Any]:
    if not lines:
        raise ValueError("SAT file is empty")
    tokens = lines[0].split()
    declared_body_count = None
    if len(tokens) >= 3:
        try:
            declared_body_count = int(tokens[2])
        except ValueError:
            pass
    return {
        "header_line_1": lines[0].decode("latin-1"),
        "header_line_2": lines[1].decode("latin-1") if len(lines) > 1 else None,
        "declared_body_count": declared_body_count,
    }


def audit_sat(path: Path) -> dict[str, Any]:
    data = path.read_bytes()
    lines = data.splitlines()
    header = _parse_header(lines)
    names: dict[int, dict[str, str]] = {}
    bodies: list[dict[str, Any]] = []

    for line in lines:
        name = NAME_RE.match(line)
        if name:
            attribute_id = int(name.group(1))
            raw_name = name.group(3)
            names[attribute_id] = {
                "display_name": _display_name(raw_name),
                "source_bytes_hex": raw_name.hex(),
            }
            continue
        body = BODY_RE.match(line)
        if body:
            attribute_id = abs(int(body.group(2)))
            bodies.append(
                {
                    "sat_body_id": int(body.group(1)),
                    "name_attribute_id": attribute_id,
                }
            )

    for body in bodies:
        body["name"] = names.get(body["name_attribute_id"])

    return {
        "schema": "radia.sat-audit/v1",
        "source_sat": str(path.resolve()),
        "sha256": hashlib.sha256(data).hexdigest(),
        "size_bytes": len(data),
        "audited_at_utc": datetime.now(timezone.utc).isoformat(),
        "header": header,
        "parsed_body_count": len(bodies),
        "body_count_matches_header": (
            header["declared_body_count"] is None
            or header["declared_body_count"] == len(bodies)
        ),
        "bodies": bodies,
        "warning": (
            "This is a structural SAT inventory only. It does not identify "
            "materials, coil centerlines, turns, current direction, or "
            "electrical connectivity."
        ),
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sat", type=Path, help="input ACIS SAT file")
    parser.add_argument("--output", type=Path, required=True,
                        help="JSON audit output path")
    args = parser.parse_args()

    if not args.sat.is_file():
        parser.error(f"SAT file does not exist: {args.sat}")
    if args.output.resolve() == args.sat.resolve():
        parser.error("output must not overwrite the SAT file")

    report = audit_sat(args.sat)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(
        json.dumps(report, ensure_ascii=False, indent=2) + "\n", encoding="utf-8"
    )
    print(
        "SAT audit: "
        f"declared={report['header']['declared_body_count']} "
        f"parsed={report['parsed_body_count']} "
        f"match={report['body_count_matches_header']}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
