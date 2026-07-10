from __future__ import annotations

import math

from ltspice_converter.measure import (
    parse_ltspice_measure_lines,
    parse_ltspice_step_lines,
    parse_spice_scalar,
    summarize_measure_log,
)


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
