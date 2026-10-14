"""Pure-Python convection-boundary schema shared by IH front ends."""
from __future__ import annotations

import json
import math
from dataclasses import dataclass
from pathlib import Path
from typing import Mapping

KELVIN_OFFSET_C = 273.15


@dataclass(frozen=True)
class ConvectionBoundary:
    label: str
    h_W_m2K: float
    ambient_C: float


def _pairs_no_duplicates(pairs):
    result = {}
    for key, value in pairs:
        if key in result:
            raise ValueError(f"duplicate convection boundary label {key!r}")
        result[key] = value
    return result


def _finite(value, name, *, nonnegative=False):
    if isinstance(value, bool):
        raise ValueError(f"{name} must be a finite number")
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{name} must be a finite number") from exc
    if not math.isfinite(result) or (nonnegative and result < 0.0):
        suffix = " and nonnegative" if nonnegative else ""
        raise ValueError(f"{name} must be finite{suffix}")
    return result


def normalize_convection_map(mapping) -> tuple[ConvectionBoundary, ...]:
    """Validate the exact-label JSON/Python convection map."""
    if mapping is None:
        return ()
    if not isinstance(mapping, Mapping):
        raise ValueError("convection_map must be a mapping")
    result = []
    for label, values in mapping.items():
        if not isinstance(label, str) or not label:
            raise ValueError("convection_map boundary labels must be nonempty strings")
        if not isinstance(values, Mapping):
            raise ValueError(f"convection_map[{label!r}] must be an object")
        expected = {"h_W_m2K", "ambient_C"}
        extra = set(values) - expected
        missing = expected - set(values)
        if missing or extra:
            raise ValueError(
                f"convection_map[{label!r}] requires exactly {sorted(expected)!r}; "
                f"missing={sorted(missing)!r}, extra={sorted(extra)!r}")
        ambient = _finite(values["ambient_C"],
                          f"convection_map[{label!r}].ambient_C")
        if ambient < -KELVIN_OFFSET_C:
            raise ValueError(
                f"convection_map[{label!r}].ambient_C is below absolute zero")
        result.append(ConvectionBoundary(
            label=label,
            h_W_m2K=_finite(values["h_W_m2K"],
                            f"convection_map[{label!r}].h_W_m2K",
                            nonnegative=True),
            ambient_C=ambient))
    return tuple(result)


def load_convection_map(path_or_mapping) -> dict[str, dict[str, float]]:
    """Load and canonicalize a JSON/Python map, rejecting duplicate JSON keys."""
    if isinstance(path_or_mapping, (str, Path)):
        try:
            raw = json.loads(Path(path_or_mapping).read_text(encoding="utf-8-sig"),
                             object_pairs_hook=_pairs_no_duplicates)
        except json.JSONDecodeError as exc:
            raise ValueError(f"invalid convection-map JSON: {exc}") from exc
        records = normalize_convection_map(raw)
    else:
        records = normalize_convection_map(path_or_mapping)
    return {record.label: {"h_W_m2K": record.h_W_m2K,
                           "ambient_C": record.ambient_C}
            for record in records}
