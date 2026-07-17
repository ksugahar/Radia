from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent))

from ltspice_converter.ideal_transformer_gate import ideal_transformer_identity_gate
from test_generalization_v23 import _v23


_PROMOTED_CASE_IDS = (
    "v24_public_hierarchical_step_parameter_scope_model_bin_temperature_sample_mismatch",
    "v24_public_ac_noise_source_normalization_node_alias_complex_axis_generation_mismatch",
)


def _v24():
    summary = _v23()
    positive = summary["metrics"]["positive"]
    positive["hierarchical_step_parameter_scope_model_bin_temperature_sample_generation_identity"] = {
        "step_generation_id": "hier-step-101",
        "scope_step_generation_id": "hier-step-101",
        "model_bin_step_generation_id": "hier-step-101",
        "temperature_step_generation_id": "hier-step-101",
        "sample_row_step_generation_id": "hier-step-101",
        "result_step_generation_id": "hier-step-101",
        "hierarchy_path": "XTOP.XAMP",
        "result_hierarchy_path": "XTOP.XAMP",
        "parameter_scope": "XTOP.XAMP:RLOAD",
        "result_parameter_scope": "XTOP.XAMP:RLOAD",
        "model_bin": "NMOS_BIN_2",
        "result_model_bin": "NMOS_BIN_2",
        "step_sample_ids": [1, 2, 3],
        "result_step_sample_ids": [1, 2, 3],
        "parameter_values_ohm": [1000.0, 2000.0, 3000.0],
        "result_parameter_values_ohm": [1000.0, 2000.0, 3000.0],
        "temperatures_c": [25.0, 50.0, 75.0],
        "result_temperatures_c": [25.0, 50.0, 75.0],
        "sample_values_v": [1.2, 1.0, 0.8],
        "result_sample_values_v": [1.2, 1.0, 0.8],
        "step_table_sha256": "1" * 64,
        "result_step_table_sha256": "1" * 64,
    }
    positive["ac_noise_source_normalization_node_alias_complex_axis_generation_identity"] = {
        "analysis_generation_id": "ac-noise-101",
        "source_analysis_generation_id": "ac-noise-101",
        "node_alias_analysis_generation_id": "ac-noise-101",
        "complex_axis_analysis_generation_id": "ac-noise-101",
        "frequency_grid_analysis_generation_id": "ac-noise-101",
        "result_analysis_generation_id": "ac-noise-101",
        "source_id": "V1",
        "result_source_id": "V1",
        "source_normalization": "1_V_ac",
        "result_source_normalization": "1_V_ac",
        "node_aliases": [["out", "V(n003)"], ["in", "V(n001)"]],
        "result_node_aliases": [["out", "V(n003)"], ["in", "V(n001)"]],
        "complex_axis_convention": "real_imaginary",
        "result_complex_axis_convention": "real_imaginary",
        "frequency_grid_hz": [10.0, 100.0, 1000.0],
        "result_frequency_grid_hz": [10.0, 100.0, 1000.0],
        "transfer_function_ri": [[1.0, -0.1], [0.7, -0.3], [0.2, -0.4]],
        "result_transfer_function_ri": [[1.0, -0.1], [0.7, -0.3], [0.2, -0.4]],
        "output_noise_v_per_sqrt_hz": [1.0e-9, 1.5e-9, 2.0e-9],
        "result_output_noise_v_per_sqrt_hz": [1.0e-9, 1.5e-9, 2.0e-9],
        "ac_noise_table_sha256": "2" * 64,
        "result_ac_noise_table_sha256": "2" * 64,
    }
    return summary


def test_v24_positive():
    assert ideal_transformer_identity_gate(_v24())["status"] == "ok"


def test_v24_public_hierarchical_step_parameter_scope_model_bin_temperature_sample_mismatch():
    summary = _v24()
    contract = summary["metrics"]["positive"][
        "hierarchical_step_parameter_scope_model_bin_temperature_sample_generation_identity"
    ]
    contract.update(
        {
            "scope_step_generation_id": "hier-step-100",
            "result_hierarchy_path": "XTOP.XOLD",
            "result_step_sample_ids": [1, 3, 4],
            "result_temperatures_c": [25.0, 27.0, 75.0],
            "result_step_table_sha256": "a" * 64,
        }
    )
    result = ideal_transformer_identity_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "hierarchical_steps_use_current_scope_model_bin_temperature_and_sample_rows"
    ]


def test_v24_public_ac_noise_source_normalization_node_alias_complex_axis_generation_mismatch():
    summary = _v24()
    contract = summary["metrics"]["positive"][
        "ac_noise_source_normalization_node_alias_complex_axis_generation_identity"
    ]
    contract.update(
        {
            "source_analysis_generation_id": "ac-noise-100",
            "result_source_id": "I1",
            "result_source_normalization": "1_A_ac",
            "result_complex_axis_convention": "magnitude_phase_deg",
            "result_frequency_grid_hz": [10.0, 200.0, 1000.0],
            "result_ac_noise_table_sha256": "b" * 64,
        }
    )
    result = ideal_transformer_identity_gate(summary)
    assert result["status"] == "needs_attention"
    assert not result["checks"][
        "ac_noise_uses_current_source_normalization_aliases_complex_axis_and_grid"
    ]
