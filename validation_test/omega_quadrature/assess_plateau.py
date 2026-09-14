"""Assess the predeclared ESRF6 order-2 high-rule quadrature plateau."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np


def assess(payload: dict, low_bonus: int = 8, high_bonus: int = 12) -> dict:
    rows = {(int(row["bonus"]), int(row["order"])): row for row in payload["rows"]}
    low = rows[(low_bonus, 2)]
    high = rows[(high_bonus, 2)]
    fields = []
    checks = {
        "completed": payload.get("completed") is True,
        "source_unchanged": payload.get("source_unchanged") is True,
        "mesh_unchanged": payload.get("mesh_unchanged") is True,
        "wheel_mode": payload.get("runtime", {}).get("mode") == "wheel",
        "wheel_archive_hash_present": bool(
            payload.get("runtime", {}).get("direct_url", {}).get("archive_info", {}).get("hashes", {}).get("sha256")
        ),
    }
    harmonic = []
    for bonus, row in ((low_bonus, low), (high_bonus, high)):
        source = row.get("source_hodge", {})
        checks[f"hodge_bonus_{bonus}"] = source.get("bonus_intorder") == bonus
        harmonic.append(float(source["relative_harmonic_norm"]))
        residual = row["linear_residual"]["free_dofs"]["relative"]
        checks[f"free_residual_{bonus}"] = (
            residual is not None and np.isfinite(residual) and 0 <= residual <= 1e-8
        )
        action = row["block_action_residual"]["blocks"]
        checks[f"block_action_{bonus}"] = all(
            value.get("action_relative") is not None
            and np.isfinite(value["action_relative"])
            and 0 <= value["action_relative"] <= 1e-8
            for value in action.values()
        )
        field = np.asarray(row["field_observations"]["B_samples_T"], dtype=float)
        checks[f"finite_field_{bonus}"] = field.ndim == 2 and field.shape[1] == 3 and np.isfinite(field).all()
        fields.append(field)
    checks["same_samples"] = (
        low["field_observations"]["samples_m"] == high["field_observations"]["samples_m"]
        and fields[0].shape == fields[1].shape
    )
    delta = fields[0] - fields[1]
    high_norm = float(np.linalg.norm(fields[1]))
    field_relative_rms = float(np.linalg.norm(delta) / high_norm) if high_norm else float("inf")
    field_scale = float(np.sqrt(np.mean(np.sum(fields[1] ** 2, axis=1))))
    maximum_scaled_change = (
        float(np.max(np.linalg.norm(delta, axis=1)) / field_scale) if field_scale else float("inf")
    )
    harmonic_relative_change = abs(harmonic[0] - harmonic[1]) / max(abs(harmonic[1]), 1e-30)
    checks.update(
        field_relative_rms=field_relative_rms <= 1e-3,
        maximum_scaled_change=maximum_scaled_change <= 5e-3,
        harmonic_relative_change=harmonic_relative_change <= 1e-4,
    )
    checks = {name: bool(value) for name, value in checks.items()}
    return {
        "schema": "radia.validation.omega-quadrature-plateau.v1",
        "passed": all(checks.values()),
        "scope": "ESRF6 nominal order-2 mixed Omega; not BDM2/IMA/general acceptance",
        "bonuses": [low_bonus, high_bonus],
        "thresholds": {
            "field_relative_rms": 1e-3,
            "maximum_scaled_change": 5e-3,
            "harmonic_relative_change": 1e-4,
            "residual": 1e-8,
        },
        "observed": {
            "field_relative_rms": field_relative_rms,
            "maximum_scaled_change": maximum_scaled_change,
            "harmonic_relative_change": harmonic_relative_change,
        },
        "checks": checks,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--low-bonus", type=int, default=8)
    parser.add_argument("--high-bonus", type=int, default=12)
    args = parser.parse_args()
    if not 0 <= args.low_bonus < args.high_bonus:
        parser.error("bonuses must be nonnegative and strictly increasing")
    result = assess(json.loads(args.input.read_text(encoding="utf-8")),
                    args.low_bonus, args.high_bonus)
    args.output.write_text(json.dumps(result, indent=2, allow_nan=False) + "\n", encoding="utf-8")
    return 0 if result["passed"] else 2


if __name__ == "__main__":
    raise SystemExit(main())
