"""Field checks shared by the versioned public-identity validators.

Each validator module binds these under its historical private names
(``_digest``, ``_result``, ``_generations`` ...); the definitions used to be
copied verbatim into every module.
"""
from __future__ import annotations

import math
from collections.abc import Mapping, Sequence


def is_sha256(value: object) -> bool:
    """A 64-character hexadecimal digest (either case)."""
    if not isinstance(value, str):
        return False
    text = value.lower()
    return len(text) == 64 and all(character in "0123456789abcdef" for character in text)


def result_digest_accepted(row: Mapping[str, object]) -> bool:
    """The row's result digest is well formed and equals the accepted one."""
    return is_sha256(row.get("result_sha256")) and row.get("accepted_result_sha256") == row.get("result_sha256")


def is_finite_number(value: object) -> bool:
    return isinstance(value, (int, float)) and not isinstance(value, bool) and math.isfinite(float(value))


def generation_closed(row: Mapping[str, object], *names: str) -> bool:
    """Every named field carries the row's non-empty ``generation``."""
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def generation_closed_fields(row: Mapping[str, object], names: Sequence[str]) -> bool:
    """:func:`generation_closed` with the field names as one sequence."""
    generation = str(row.get("generation") or "")
    return bool(generation) and all(row.get(name) == generation for name in names)


def generation_closed_stripped(row: Mapping[str, object], fields: Sequence[str]) -> bool:
    """:func:`generation_closed_fields` that ignores surrounding whitespace in ``generation``."""
    generation = str(row.get("generation", "")).strip()
    return bool(generation) and all(row.get(field) == generation for field in fields)


def finite_vector(value: object, length: int | None = None) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and (length is None or len(value) == length)
        and bool(value)
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def finite_sequence(value: object, *, minimum: int = 1) -> bool:
    return (
        isinstance(value, Sequence)
        and not isinstance(value, (str, bytes))
        and len(value) >= minimum
        and all(isinstance(item, (int, float)) and math.isfinite(float(item)) for item in value)
    )


def is_integer(value: object) -> bool:
    return isinstance(value, int) and not isinstance(value, bool)


def prefixed_equal(value: Mapping[str, object], left: str, right: str, prefix: str) -> bool:
    return str(value.get(left, "")).startswith(prefix) and value.get(left) == value.get(right)


def public_rows(payload: Mapping[str, object]) -> list[Mapping[str, object]]:
    """Reference rows plus every measured family's rows (mappings only)."""
    rows: list[Mapping[str, object]] = []
    reference = payload.get("reference")
    if isinstance(reference, Sequence) and not isinstance(reference, (str, bytes)):
        rows.extend(item for item in reference if isinstance(item, Mapping))
    measured = payload.get("measured")
    if isinstance(measured, Mapping):
        for family in measured.values():
            if isinstance(family, Sequence) and not isinstance(family, (str, bytes)):
                rows.extend(item for item in family if isinstance(item, Mapping))
    return rows


def gate_report(policy: str, checks: dict[str, bool]) -> dict[str, object]:
    return {
        "policy": policy,
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, accepted in checks.items() if not accepted],
    }
