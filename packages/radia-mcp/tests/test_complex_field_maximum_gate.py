import copy
import json

from radia_mcp.radia_ngsolve.server import complex_vector_field_maximum_gate


def good() -> dict:
    return {
        "cases": [
            {
                "case_id": "field-a",
                "frequency_hz": 1000.0,
                "ampere_turns": 100.0,
                "field_unit": "T",
                "rows": [
                    {"part": "real", "element": 1, "material_id": 2, "bx_t": 3.0, "by_t": 4.0, "bz_t": 0.0, "bmag_t": 5.0},
                    {"part": "real", "element": 2, "material_id": 2, "bx_t": 0.0, "by_t": 0.0, "bz_t": 2.0, "bmag_t": 2.0},
                    {"part": "imaginary", "element": 1, "material_id": 2, "bx_t": 0.0, "by_t": 0.0, "bz_t": 7.0, "bmag_t": 7.0},
                ],
                "maxima": [
                    {"part": "real", "element": 1, "material_id": 2, "bmax_t": 5.0},
                    {"part": "imaginary", "element": 1, "material_id": 2, "bmax_t": 7.0},
                ],
            }
        ]
    }


def test_accepts_component_norms_and_material_maxima():
    result = json.loads(complex_vector_field_maximum_gate(json.dumps(good())))
    assert result["status"] == "ok"
    assert result["metrics"]["maximum_component_magnitude_relative_error"] == 0.0


def test_rejects_wrong_norm_and_stale_maximum_element():
    payload = copy.deepcopy(good())
    payload["cases"][0]["rows"][0]["bmag_t"] = 4.0
    payload["cases"][0]["maxima"][1]["element"] = 99
    result = json.loads(complex_vector_field_maximum_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["case_1_component_magnitudes_close"] is False
    assert result["checks"]["case_1_maximum_elements_match"] is False


def test_rejects_missing_imaginary_part_and_wrong_unit():
    payload = good()
    payload["cases"][0]["rows"] = [
        row for row in payload["cases"][0]["rows"] if row["part"] == "real"
    ]
    payload["cases"][0]["maxima"] = [payload["cases"][0]["maxima"][0]]
    payload["cases"][0]["field_unit"] = "mT"
    result = json.loads(complex_vector_field_maximum_gate(json.dumps(payload)))
    assert result["status"] == "needs_attention"
    assert result["checks"]["case_1_real_and_imaginary_parts_present"] is False
    assert result["checks"]["case_1_tesla_frequency_excitation_recorded"] is False
