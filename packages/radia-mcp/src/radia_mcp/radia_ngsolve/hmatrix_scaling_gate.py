"""Solver-neutral quality gate for educational H-matrix scaling studies."""

from __future__ import annotations

import math
from typing import Any


def _number(row: dict[str, Any], snake: str, camel: str) -> float:
    value = row[snake] if snake in row else row.get(camel)
    if value is None:
        raise ValueError(f"row is missing {snake}")
    result = float(value)
    if not math.isfinite(result):
        raise ValueError(f"{snake} must be finite")
    return result


def hmatrix_compression_scaling_gate(
    rows: list[dict[str, Any]],
    *,
    max_matvec_relative_error: float = 1.0e-8,
    max_rank: int = 20,
    max_storage_growth_exponent: float = 1.25,
    min_dense_growth_exponent: float = 1.9,
    max_compression_ratio_consistency_error: float = 1.0e-12,
) -> dict[str, Any]:
    """Require accurate low-rank matvecs and subquadratic stored-entry growth."""

    if len(rows) < 3:
        raise ValueError("at least three scaling rows are required")
    tolerances = [
        float(max_matvec_relative_error),
        float(max_storage_growth_exponent),
        float(min_dense_growth_exponent),
        float(max_compression_ratio_consistency_error),
    ]
    if int(max_rank) < 1:
        raise ValueError("max_rank must be positive")
    if any(not math.isfinite(value) or value < 0.0 for value in tolerances):
        raise ValueError("tolerances must be finite and nonnegative")

    normalized = []
    for row in rows:
        if not isinstance(row, dict):
            raise ValueError("every scaling row must be an object")
        normalized.append(
            {
                "point_count": int(_number(row, "point_count", "pointCount")),
                "max_rank": int(_number(row, "max_rank", "maxRank")),
                "low_rank_blocks": int(_number(row, "low_rank_blocks", "lowRankBlocks")),
                "stored_entries": int(_number(row, "stored_entries", "storedEntries")),
                "dense_entries": int(_number(row, "dense_entries", "denseEntries")),
                "compression_ratio": _number(row, "compression_ratio", "compressionRatio"),
                "matvec_relative_error": _number(
                    row, "matvec_relative_error", "matvecRelativeError"
                ),
                "build_seconds": _number(row, "build_seconds", "buildSeconds"),
                "matvec_seconds": _number(row, "matvec_seconds", "matvecSeconds"),
                "dense_reference_seconds": _number(
                    row, "dense_reference_seconds", "denseReferenceSeconds"
                ),
            }
        )

    counts = [row["point_count"] for row in normalized]
    ranks = [row["max_rank"] for row in normalized]
    low_rank_blocks = [row["low_rank_blocks"] for row in normalized]
    stored = [row["stored_entries"] for row in normalized]
    dense = [row["dense_entries"] for row in normalized]
    compression = [row["compression_ratio"] for row in normalized]
    errors = [row["matvec_relative_error"] for row in normalized]
    timings = [
        row[key]
        for row in normalized
        for key in ("build_seconds", "matvec_seconds", "dense_reference_seconds")
    ]
    storage_exponents = [
        math.log(stored[i] / stored[i - 1]) / math.log(counts[i] / counts[i - 1])
        for i in range(1, len(rows))
        if counts[i] > counts[i - 1] and stored[i] > 0 and stored[i - 1] > 0
    ]
    dense_exponents = [
        math.log(dense[i] / dense[i - 1]) / math.log(counts[i] / counts[i - 1])
        for i in range(1, len(rows))
        if counts[i] > counts[i - 1] and dense[i] > 0 and dense[i - 1] > 0
    ]
    ratio_consistency = [
        abs(ratio - entries / reference)
        for ratio, entries, reference in zip(compression, stored, dense)
        if reference > 0
    ]

    checks = {
        "point_counts_strictly_increase": all(a < b for a, b in zip(counts, counts[1:])),
        "point_counts_positive": all(value > 0 for value in counts),
        "low_rank_blocks_present": all(value > 0 for value in low_rank_blocks),
        "rank_is_bounded": all(1 <= value <= int(max_rank) for value in ranks),
        "stored_entries_are_compressed": all(
            0 < compressed < reference for compressed, reference in zip(stored, dense)
        ),
        "compression_ratio_matches_storage": (
            len(ratio_consistency) == len(rows)
            and max(ratio_consistency) <= float(max_compression_ratio_consistency_error)
        ),
        "compression_improves_with_size": all(
            a > b for a, b in zip(compression, compression[1:])
        ),
        "matvec_matches_dense_reference": all(
            0.0 <= value <= float(max_matvec_relative_error) for value in errors
        ),
        "stored_entry_growth_is_subquadratic": (
            len(storage_exponents) == len(rows) - 1
            and max(storage_exponents) <= float(max_storage_growth_exponent)
        ),
        "dense_reference_growth_is_quadratic": (
            len(dense_exponents) == len(rows) - 1
            and min(dense_exponents) >= float(min_dense_growth_exponent)
        ),
        "timings_are_nonnegative": all(value >= 0.0 for value in timings),
    }
    return {
        "policy": "hmatrix_compression_scaling_gate_v1",
        "status": "ok" if all(checks.values()) else "needs_attention",
        "checks": checks,
        "issues": [name for name, ok in checks.items() if not ok],
        "metrics": {
            "point_counts": counts,
            "max_ranks": ranks,
            "stored_entries": stored,
            "dense_entries": dense,
            "compression_ratios": compression,
            "max_matvec_relative_error": max(errors),
            "storage_growth_exponents": storage_exponents,
            "dense_growth_exponents": dense_exponents,
            "max_compression_ratio_consistency_error": max(ratio_consistency, default=math.inf),
        },
        "notes": [
            "Accuracy alone is insufficient: an H-matrix study must also show bounded far-field rank and subquadratic storage growth.",
            "Wall-clock monotonicity is deliberately not required because first-call JIT and cache effects dominate small educational cases.",
        ],
    }
