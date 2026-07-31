"""P1 electrostatic BEM convergence for two finite square electrodes."""

from __future__ import annotations

import hashlib
import json
import platform
import subprocess
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

import numpy as np


ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT / "src"))

from radia.bem.electrostatic_p1 import solve_prescribed_potential_p1


OUT_JSON = Path(__file__).with_name(
    "validation_parallel_square_plates_p1_bem_summary.json"
)
SIZE_M = 10e-3
GAP_M = 5e-3
VOLTAGE_V = 1000.0


def _plate(n: int, z: float, offset: int) -> tuple[np.ndarray, list[list[int]]]:
    grid = np.linspace(-0.5 * SIZE_M, 0.5 * SIZE_M, n + 1)
    vertices = np.array([[x, y, z] for y in grid for x in grid], dtype=float)
    triangles: list[list[int]] = []
    for j in range(n):
        for i in range(n):
            a = offset + j * (n + 1) + i
            b = a + 1
            c = a + n + 1
            d = c + 1
            triangles.extend(([a, b, d], [a, d, c]))
    return vertices, triangles


def _identity_digest() -> str:
    identity = {
        "geometry": "two_aligned_finite_square_electrodes",
        "size_m": SIZE_M,
        "gap_m": GAP_M,
        "potentials_v": [VOLTAGE_V, -VOLTAGE_V],
        "medium": "vacuum",
    }
    payload = json.dumps(identity, sort_keys=True, separators=(",", ":"))
    return hashlib.sha256(payload.encode("ascii")).hexdigest()


def _git_head() -> str:
    completed = subprocess.run(
        ["git", "-C", str(ROOT), "rev-parse", "HEAD"],
        check=True,
        capture_output=True,
        text=True,
    )
    return completed.stdout.strip()


def solve_level(n: int) -> dict[str, float | int]:
    positive, positive_tris = _plate(n, +0.5 * GAP_M, 0)
    negative, negative_tris = _plate(n, -0.5 * GAP_M, len(positive))
    vertices = np.vstack((positive, negative))
    triangles = np.asarray(positive_tris + negative_tris, dtype=np.int64)
    potential = np.r_[
        np.full(len(positive), VOLTAGE_V),
        np.full(len(negative), -VOLTAGE_V),
    ]
    started = time.perf_counter()
    result = solve_prescribed_potential_p1(vertices, triangles, potential)
    duration = time.perf_counter() - started
    positive_ids = np.arange(len(positive), dtype=np.int64)
    negative_ids = np.arange(len(positive), len(vertices), dtype=np.int64)
    positive_charge = result.charge_on_vertices(positive_ids)
    negative_charge = result.charge_on_vertices(negative_ids)
    return {
        "subdivisions_per_side": n,
        "vertices": len(vertices),
        "triangles": len(triangles),
        "positive_charge_c": positive_charge,
        "negative_charge_c": negative_charge,
        "charge_balance_c": positive_charge + negative_charge,
        "capacitance_f": positive_charge / (2.0 * VOLTAGE_V),
        "duration_s": duration,
    }


def main() -> int:
    levels = [solve_level(n) for n in (4, 8, 12)]
    relative_last_change = abs(
        levels[-1]["positive_charge_c"] - levels[-2]["positive_charge_c"]
    ) / abs(levels[-1]["positive_charge_c"])
    max_balance = max(abs(row["charge_balance_c"]) for row in levels)
    checks = {
        "opposite_conductor_charges": all(
            row["positive_charge_c"] > 0.0 and row["negative_charge_c"] < 0.0
            for row in levels
        ),
        "charge_conservation": max_balance < 1e-20,
        "positive_capacitance": all(row["capacitance_f"] > 0.0 for row in levels),
        "last_refinement_change_below_1_percent": relative_last_change < 0.01,
    }
    artifact = {
        "schema": "radia.validation.electrostatic.parallel-square-plates-p1-bem.v1",
        "executed_at_utc": datetime.now(timezone.utc).isoformat(),
        "execution_version": {
            "producer": Path(__file__).name,
            "python": platform.python_version(),
            "radia_source_head": _git_head(),
        },
        "physics_identity_digest": _identity_digest(),
        "discretization": "continuous P1 surface charge on flat triangles",
        "levels": levels,
        "metrics": {
            "relative_last_refinement_change": relative_last_change,
            "maximum_charge_imbalance_c": max_balance,
        },
        "timing_breakdown_s": {
            f"solve_n{row['subdivisions_per_side']}": row["duration_s"]
            for row in sorted(levels, key=lambda item: item["duration_s"], reverse=True)
        },
        "checks": checks,
        "pass": all(checks.values()),
    }
    OUT_JSON.write_text(json.dumps(artifact, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(artifact, indent=2))
    return 0 if artifact["pass"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
