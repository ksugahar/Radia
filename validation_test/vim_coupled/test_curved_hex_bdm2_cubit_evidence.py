import json
from pathlib import Path

RESULT = Path(__file__).with_name("curved_hex_bdm2_cubit_summary.json")


def _result() -> dict:
    return json.loads(RESULT.read_text(encoding="utf-8"))


def test_cubit_curve_two_mesh_is_checked_and_nonaffine():
    result = _result()
    mesh = result["mesh"]

    assert result["schema"] == "radia.validation.curved-hex-bdm2-cubit.v1"
    assert result["machine"] == "mdx"
    assert result["versions"]["radia"] == "4.95.71"
    assert result["versions"]["cubit_mesh_export"] == "0.14.14"
    assert mesh["curve_order"] == 2
    assert mesh["elements"] == 4
    assert mesh["affinity"]["nonaffine_cell_count"] == mesh["elements"]
    assert mesh["minimum_scaled_jacobian"] > 0.0
    assert len(mesh["sha256"]) == 64
    assert len(mesh["sidecar_sha256"]) == 64


def test_cubit_curve_two_bdm2_solve_and_field_are_production_green():
    result = _result()

    assert result["linear"]["converged"] is True
    assert result["linear"]["linear_solver"] == "mass-riesz-cg"
    assert result["nonlinear"]["converged"] is True
    assert result["nonlinear"]["linear_solver"] == "energy-newton-cpp"
    assert result["equivalent_material_relative_difference"] < 1.0e-7
    assert result["prescribed_source"]["maximum_relative_error"] < 2.0e-10
    assert result["prescribed_source"]["source_stats"]["curve_order"] == 2
    assert result["shape_derivative"]["rejected"] is True
    assert all(result["checks"].values())
    assert result["pass"] is True
