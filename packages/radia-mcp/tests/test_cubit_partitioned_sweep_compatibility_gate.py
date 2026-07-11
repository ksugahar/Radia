import copy
import json

from radia_mcp.cubit.server import cubit_partitioned_sweep_compatibility_gate as mcp_gate
from radia_mcp.cubit.vol_inventory import cubit_partitioned_sweep_compatibility_gate


def _summary():
    first = {str(index): 1.0 + index / 100.0 for index in range(1, 13)}
    copied = {str(index + 12): value for index, value in enumerate(first.values(), 1)}
    volumes = {**first, **copied}
    return {
        "source_journal": "partitioned_sweep.jou",
        "source_sha256": "a" * 64,
        "execution_mode": "python_api_headless",
        "persistent_gui_started": False,
        "process_exit_code": 2,
        "startup_diagnostics": [
            "Could not open file: <install>/plugins",
            "Could not open file: -commandplugindir",
        ],
        "script_error_lines": [],
        "result_artifact_fresh": True,
        "source_command_count": 98,
        "executed_command_count": 98,
        "compatibility_transforms": {
            "surface 22 scheme quad_dominant": "surface 22 scheme pave"
        },
        "webcut_count": 11,
        "base_volume_count": 12,
        "volume_count": 24,
        "base_hex_count": 126,
        "mesh_copy_factor": 2,
        "element_counts": {"hex": 252, "pyramid": 0, "wedge": 0, "tet": 0},
        "quality": {
            "scaled_jacobian": {"min": 0.2438},
            "shape": {"min": 0.3060},
        },
        "cad_volume_by_body": volumes,
        "total_cad_volume": sum(volumes.values()),
    }


def test_accepts_legacy_partitioned_sweep_and_mesh_copy():
    result = cubit_partitioned_sweep_compatibility_gate(_summary())
    assert result["status"] == "ok"
    assert result["launcher_classification"] == "allowlisted_startup_diagnostic_with_clean_script"
    assert json.loads(mcp_gate(_summary()))["status"] == "ok"


def test_rejects_silent_legacy_skip_and_copy_count_drift():
    bad = copy.deepcopy(_summary())
    bad["compatibility_transforms"] = {}
    bad["element_counts"]["hex"] = 251
    result = cubit_partitioned_sweep_compatibility_gate(bad)
    assert result["status"] == "needs_attention"
    assert set(result["issues"]) >= {
        "legacy_quad_dominant_promoted_to_pave",
        "hex_copy_count_conserved",
    }


def test_rejects_mixed_volume_elements_and_unrelated_launcher_error():
    bad = copy.deepcopy(_summary())
    bad["element_counts"]["tet"] = 1
    bad["startup_diagnostics"] = ["unrelated startup failure"]
    result = cubit_partitioned_sweep_compatibility_gate(bad)
    assert result["checks"]["hex_only_volume_mesh"] is False
    assert result["launcher_classification"] == "execution_error"
