from __future__ import annotations

import json
import math
from pathlib import Path


RESULT = Path(__file__).with_name("cogging_skew_demo_results.json")


def test_cogging_skew_demo_evidence_preserves_physical_checks():
    result = json.loads(RESULT.read_text(encoding="utf-8"))
    assert result["schema"] == "radia.validation.electric_machine.cogging_skew_demo.v1"
    assert all(result["checks"].values())
    assert result["dominant_order"] == 2
    assert result["mean_fraction"] < 0.15
    assert result["peak_torque_Nm"] > 0.0
    for row in result["skew_rows"]:
        assert math.isclose(
            row["amplitude_ratio"],
            row["predicted_ratio"],
            rel_tol=2.0e-2,
            abs_tol=1.0e-12,
        )
    assert abs(result["skew_rows"][-1]["amplitude_ratio"]) < 1.0e-12
