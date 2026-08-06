"""Fail-fast contracts for Gmsh field-line tools (no Gmsh required)."""

import math

import pytest
from radia_mcp.gmsh import post_process


@pytest.fixture
def input_file(tmp_path):
    path = tmp_path / "field.msh"
    path.write_text("placeholder", encoding="ascii")
    return path


@pytest.fixture(autouse=True)
def forbid_subprocess(monkeypatch):
    def unexpected(*_args, **_kwargs):
        raise AssertionError("invalid input reached the Gmsh subprocess")

    monkeypatch.setattr(post_process, "_run_post", unexpected)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"step_size": 0.0}, "step_size must be positive"),
        ({"n_seeds": 0}, "n_seeds must be positive"),
        ({"max_steps": 0}, "max_steps must be positive"),
        ({"max_turn_deg": 0.0}, "max_turn_deg must be in"),
        ({"arrows_every": -1}, "arrows_every must be non-negative"),
        ({"time_step": -1}, "time_step must be non-negative"),
    ],
)
def test_streamlines_rejects_invalid_controls(input_file, kwargs, message):
    result = post_process.streamlines(
        input_file, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0], **kwargs)
    assert result["ok"] is False
    assert message in result["error"]


def test_streamlines_rejects_nonfinite_seed(input_file):
    result = post_process.streamlines(
        input_file, [math.nan, 0.0, 0.0], [1.0, 0.0, 0.0])
    assert result == {"ok": False,
                      "error": "seed_start must be a finite number"}


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"n_levels": 0}, "n_levels must be positive"),
        ({"levels": []}, "levels must not be empty"),
        ({"levels": [0.1, math.inf]}, "levels must be a finite number"),
        ({"levels": [0.1, 0.1]}, "levels must be unique"),
        ({"result_name": "  "}, "result_name must not be empty"),
    ],
)
def test_flux_lines_rejects_invalid_levels(input_file, kwargs, message):
    result = post_process.flux_lines(input_file, **kwargs)
    assert result["ok"] is False
    assert message in result["error"]


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"d_sep": 0.0}, "d_sep must be positive"),
        ({"d_sep": 1.0, "d_test": 1.1},
         "d_test must not exceed d_sep"),
        ({"d_sep": 1.0, "d_test": 0.5, "step_size": 0.6},
         "step_size must not exceed d_test"),
        ({"first_seed": [1.1, 0.5]}, "first_seed must lie inside"),
        ({"max_lines": 0}, "max_lines must be positive"),
        ({"max_steps": 0}, "max_steps must be positive"),
        ({"max_total_steps": 0}, "max_total_steps must be positive"),
    ],
)
def test_streamlines_2d_rejects_invalid_controls(input_file, kwargs,
                                                  message):
    result = post_process.streamlines_2d(
        input_file, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
        [0.0, 1.0, 0.0], **kwargs)
    assert result["ok"] is False
    assert message in result["error"]


def test_streamlines_2d_rejects_degenerate_plane(input_file):
    result = post_process.streamlines_2d(
        input_file, [0.0, 0.0, 0.0], [1.0, 0.0, 0.0],
        [2.0, 0.0, 0.0])
    assert result["ok"] is False
    assert "must not be collinear" in result["error"]
