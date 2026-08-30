import json
from pathlib import Path

PRODUCTION = Path(__file__).with_name("mapped_hex_bdm2_production_summary.json")
REFERENCE = Path(__file__).with_name(
    "mapped_hex_bdm2_quadrature_reference_summary.json"
)


def _load(path: Path) -> dict:
    return json.loads(path.read_text(encoding="utf-8"))


def test_mapped_hex_bdm2_production_summary_is_green():
    result = _load(PRODUCTION)
    operator = result["operator_quadrature"]

    assert result["schema"] == "radia.validation.mapped-hex-bdm2-production.v1"
    assert result["machine"] == "mdx"
    assert result["geometry"]["geometry_order"] == 1
    assert "trilinear" in result["geometry"]["map"]
    assert [(row["outer_order"], row["inner_order"]) for row in operator["rules"]] == [
        (9, 12),
        (10, 16),
    ]
    assert all(
        row["eigenvalues_outside_physical_interval"] == 0 for row in operator["rules"]
    )
    assert operator["comparison"]["material_solution_mass_relative"] < 1.0e-3
    assert result["prescribed_source_roundoff"]["field_error_in_machine_eps"] < 10.0
    assert all(result["checks"].values())
    assert result["pass"] is True


def test_mapped_hex_bdm2_quadrature_reference_is_green():
    result = _load(REFERENCE)
    operator = result["operator_quadrature"]

    assert result["schema"] == (
        "radia.validation.mapped-hex-bdm2-quadrature-reference.v1"
    )
    assert result["machine"] == "mdx"
    assert [(row["outer_order"], row["inner_order"]) for row in operator["rules"]] == [
        (10, 16),
        (11, 20),
    ]
    assert all(
        row["eigenvalues_outside_physical_interval"] == 0 for row in operator["rules"]
    )
    assert operator["comparison"]["material_solution_mass_relative"] < 5.0e-4
    assert all(result["checks"].values())
    assert result["pass"] is True
