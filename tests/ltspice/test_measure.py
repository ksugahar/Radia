from __future__ import annotations

import math

from radia.ltspice.measure import (
    parse_ltspice_measure_lines,
    parse_ltspice_stepped_measure_tables,
    parse_ltspice_step_lines,
    parse_spice_scalar,
    summarize_measure_log,
    summarize_stepped_measure_log,
)
from radia.ltspice.mcp_server import balanced_learning_profile


def test_balanced_learning_profile_has_ten_unique_controlled_stages():
    profile = balanced_learning_profile()
    assert profile["policy"] == "equal_capability_gain_v1"
    assert profile["stage_count"] == 10
    assert len({row["capability_id"] for row in profile["stages"]}) == 10
    assert all(row["positive_control"] and row["negative_control"] for row in profile["stages"])
    assert profile["self_check"]["status"] == "ok"
    from radia.ltspice.learning_quality import validate_balanced_learning_profile
    bad = {**profile, "stages": [dict(row) for row in profile["stages"]]}
    bad["stages"][7]["negative_control"] = ""
    assert validate_balanced_learning_profile(bad)["status"] == "needs_attention"


def test_parse_spice_scalar_engineering_suffixes():
    assert parse_spice_scalar("1k") == 1000.0
    assert parse_spice_scalar("2.2u") == 2.2e-6
    assert parse_spice_scalar("3Meg") == 3e6
    assert math.isclose(parse_spice_scalar("{47n}"), 47e-9)
    assert parse_spice_scalar("bad") is None


def test_parse_ltspice_measure_value_event_ac_and_window_rows():
    lines = [
        "gain_mid: mag(v(out))=0.707106 from 1k to 2k",
        "t_rise1: V(in)=0.5 AT 0.000201",
        "low_gain: mag(V(out))=(-3.0103dB,-45.0deg) at 1591.55",
    ]
    rows = parse_ltspice_measure_lines(lines)
    by_name = {row["name"]: row for row in rows}

    assert by_name["gain_mid"]["value"] == 0.707106
    assert by_name["gain_mid"]["from"] == 1000.0
    assert by_name["gain_mid"]["to"] == 2000.0
    assert by_name["t_rise1"]["kind"] == "value_at"
    assert by_name["t_rise1"]["target_value"] == 0.5
    assert by_name["t_rise1"]["at"] == 0.000201
    assert by_name["low_gain"]["unit"] == "dB"
    assert math.isclose(by_name["low_gain"]["phase_deg"], -45.0)
    assert math.isclose(by_name["low_gain"]["at"], 1591.55)
    assert math.isclose(by_name["low_gain"]["linear_value"], 10 ** (-3.0103 / 20.0))
    assert by_name["low_gain"]["ac_measure_semantics_safe"] is True


def test_ac_scalar_phase_measure_is_rejected_because_sign_is_lost():
    log_text = "\n".join([
        "mag_vout_1k: mag(V(out)) =(44.0364026328dB,0deg) at 1000",
        "phase_vout_1k: ph(V(out)) =(39.0848501888dB,0deg) at 1000",
    ])
    summary = summarize_measure_log(log_text)
    rows = {row["name"]: row for row in summary["measures"]}

    assert summary["ok"] is False
    assert summary["unsafe_ac_measure_names"] == ["phase_vout_1k"]
    assert math.isclose(rows["mag_vout_1k"]["linear_value"], 159.1549430919, rel_tol=1e-10)
    assert math.isclose(rows["phase_vout_1k"]["encoded_abs_phase_deg"], 90.0, rel_tol=1e-10)
    assert rows["phase_vout_1k"]["ac_measure_semantics_safe"] is False
    assert any("complex RAW" in warning for warning in summary["warnings"])


def test_parse_ltspice_step_lines_with_numeric_assignments():
    steps = parse_ltspice_step_lines([
        ".step RLOAD=1k COUT=22u label=fast",
        "not a step",
    ])
    assert len(steps) == 1
    assert steps[0]["assignments"]["label"] == "fast"
    assert steps[0]["numeric_assignments"]["RLOAD"] == 1000.0
    assert steps[0]["numeric_assignments"]["COUT"] == 22e-6
    assert "label" not in steps[0]["numeric_assignments"]


def test_summarize_measure_log_schema_and_duplicate_gate():
    summary = summarize_measure_log(
        "\n".join([
            "gain: mag(v(out))=1",
            "gain: mag(v(out))=2",
            ".step RLOAD=1k",
        ])
    )
    assert summary["schema"] == "radia-spice-lab.measure-log.v1"
    assert summary["measure_count"] == 2
    assert summary["step_count"] == 1
    assert summary["ok"] is False
    assert summary["duplicate_measure_names"] == ["gain"]
    assert any("duplicate" in warning for warning in summary["warnings"])


def test_summarize_measure_log_requires_measure_rows():
    summary = summarize_measure_log("Elapsed time: 0.1 seconds")
    assert summary["ok"] is False
    assert summary["measure_count"] == 0
    assert summary["warnings"] == ["no LTspice .measure result rows were parsed"]


def test_parse_and_summarize_stepped_ac_measure_table():
    log_text = "\n".join([
        ".step rval=1000",
        ".step rval=2000",
        ".step rval=4000",
        "Measurement: gain",
        "  step MAX(mag(V(out)) ) FROM TO",
        "     1 (-1.2060607353dB,0deg) 90 110",
        "     2 (-3.57679299092dB,0deg) 90 110",
        "     3 (-7.85991318292dB,0deg) 90 110",
    ])
    tables = parse_ltspice_stepped_measure_tables(log_text.splitlines())
    assert len(tables) == 1
    assert tables[0]["row_count"] == 3
    assert tables[0]["rows"][1]["numeric_step_assignments"]["rval"] == 2000.0
    assert math.isclose(tables[0]["rows"][2]["value"], -7.85991318292)
    assert tables[0]["rows"][2]["unit"] == "dB"

    summary = summarize_stepped_measure_log(log_text)
    assert summary["schema"] == "radia-spice-lab.stepped-measure-log.v1"
    assert summary["ok"] is True
    assert summary["checks"] == {
        "row_counts_match_steps": True,
        "values_finite": True,
        "step_assignments_complete": True,
    }


def test_stepped_measure_summary_rejects_missing_result_row():
    summary = summarize_stepped_measure_log("\n".join([
        ".step rval=1000",
        ".step rval=2000",
        "Measurement: gain",
        " step MAX(mag(V(out)))",
        " 1 -1.0",
    ]))
    assert summary["ok"] is False
    assert summary["checks"]["row_counts_match_steps"] is False
    assert any("row counts" in warning for warning in summary["warnings"])
