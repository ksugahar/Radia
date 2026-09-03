"""High-order quadrature reference for mapped HEX BDM2 HDiv-VIM.

Run this expensive companion on hibino first, or on mdx only when hibino is
unavailable and the mdx CI queue is idle. The ordinary production gate
uses q9/q12 against q10/q16; this script checks that q10/q16 itself approaches
the still richer q11/q20 rule.
"""

from __future__ import annotations

import json
import platform
from datetime import UTC, datetime
from pathlib import Path

import ngsolve as ng
from validate_mapped_hex_bdm2_production import _mesh, _operator_sweep

import radia

DEFAULT_OUTPUT = Path(__file__).with_name(
    "mapped_hex_bdm2_quadrature_reference_summary.json"
)


def run() -> dict:
    with ng.TaskManager():
        operator = _operator_sweep(_mesh(half=False), [(10, 16), (11, 20)])
    comparison = operator["comparison"]
    checks = {
        "both_rules_have_physical_spectrum": all(
            row["eigenvalues_outside_physical_interval"] == 0
            for row in operator["rules"]
        ),
        "q10_q16_has_sub_half_per_mille_material_response_error": (
            comparison["material_solution_mass_relative"] <= 5.0e-4
        ),
    }
    return {
        "schema": "radia.validation.mapped-hex-bdm2-quadrature-reference.v1",
        "created_at_utc": datetime.now(UTC).isoformat(),
        "machine": platform.node(),
        "versions": {
            "radia": getattr(radia, "__version__", "unknown"),
            "ngsolve": getattr(ng, "__version__", "unknown"),
            "python": platform.python_version(),
        },
        "operator_quadrature": operator,
        "checks": checks,
        "pass": all(checks.values()),
    }


def main() -> None:
    summary = run()
    DEFAULT_OUTPUT.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(
        json.dumps(
            {
                "output": str(DEFAULT_OUTPUT),
                "pass": summary["pass"],
                "checks": summary["checks"],
            },
            indent=2,
        ),
        flush=True,
    )
    if not summary["pass"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
