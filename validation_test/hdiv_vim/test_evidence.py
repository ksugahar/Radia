import json
import math
from pathlib import Path


HERE = Path(__file__).resolve().parent
REPORT = HERE / "c_type_three_formulation_tosca_mixed_iron.vol-check.json"


def test_c_type_iron_mesh_preflight_evidence():
    report = json.loads(REPORT.read_text(encoding="utf-8"))

    assert report["schema"] == "cubit-mesh-export.vol-check.v1"
    assert report["passed"] is True
    assert report["labels"]["passed"] is True
    assert report["labels"]["materials"] == ["iron"]
    assert report["labels"]["boundaries"] == ["iron_boundary"]
    assert report["quality"]["passed"] is True
    assert report["quality"]["minimum_scaled_jacobian"] > 0.1
    assert report["boundary_domain_ownership"]["passed"] is True
    assert report["mesh"]["n_elements"] == 1688
    assert report["mesh"]["n_points"] == 500

    errors = [
        abs(float(entry["error_pct"]))
        for family in ("materials", "boundaries")
        for entry in report[family]
    ]
    assert errors and all(math.isfinite(error) and error < 1.0e-8 for error in errors)
