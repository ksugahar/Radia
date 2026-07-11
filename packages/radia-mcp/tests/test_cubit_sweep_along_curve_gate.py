import json

from radia_mcp.cubit.server import cubit_sweep_along_curve_gate as mcp_gate
from radia_mcp.cubit.vol_inventory import cubit_sweep_along_curve_gate


def _summary():
    return {
        "source_journal": "sweep_along-curve.jou",
        "execution_mode": "python_api_headless",
        "persistent_gui_started": False,
        "version": "2025.12",
        "process_exit_code": 2,
        "startup_diagnostics": [
            "Could not open file: <install>/plugins",
            "Could not open file: -commandplugindir",
        ],
        "script_error_lines": [],
        "result_artifact_fresh": True,
        "source_quad_count": 50,
        "sweep_interval_count": 20,
        "element_counts": {"hex": 1000, "pyramid": 0, "wedge": 0, "tet": 0},
        "quality": {
            "scaled_jacobian": {"min": 0.9999999999999999},
            "shape": {"min": 0.5772367461702522},
        },
        "cad_volume_by_body": {"1": 2000.0},
        "total_cad_volume": 2000.0,
        "analytic_volume": 2000.0,
    }


def test_accepts_source_quad_times_interval_hex_sweep():
    result = cubit_sweep_along_curve_gate(_summary())
    assert result["status"] == "ok"
    assert result["launcher_classification"] == "allowlisted_startup_diagnostic_with_clean_script"
    assert json.loads(mcp_gate(_summary()))["status"] == "ok"


def test_rejects_layer_count_or_quality_regression():
    bad = _summary()
    bad["element_counts"] = {**bad["element_counts"], "hex": 999}
    bad["quality"] = {**bad["quality"], "scaled_jacobian": {"min": -0.1}}
    result = cubit_sweep_along_curve_gate(bad)
    assert result["checks"]["sweep_layer_count_conserved"] is False
    assert result["checks"]["scaled_jacobian_above_threshold"] is False


def test_rejects_nonallowlisted_exit_or_script_error():
    bad = {**_summary(), "startup_diagnostics": ["unrelated failure"], "script_error_lines": ["mesh failed"]}
    result = cubit_sweep_along_curve_gate(bad)
    assert result["launcher_classification"] == "execution_error"
    assert result["checks"]["script_error_lines_empty"] is False
