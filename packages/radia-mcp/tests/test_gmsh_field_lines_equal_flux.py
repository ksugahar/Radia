"""Field-line seeding contracts of the gmsh adapter that need no Radia core.

Guards: seed_line stays uniform and is labelled not flux-weighted; ambiguous
or malformed seeding (schema-accepted JSON such as null, bool, nan, wrong
lengths, huge integers, over-budget tables or line counts) fails with a
field error before probing.
Assumes: the gmsh package; B = (0, 0, 1 + x/2 + y/4) as P1 node data on the
Kuhn cube. The cases that place seeds with radia.flux_seeding live in
validation_test/radia_mcp/test_gmsh_seeded_field_lines.py.
"""

import importlib.util

import numpy as np
import pytest
from msh42_fixtures import field_section
from test_gmsh_flux_streamlines import _CUBE, _P

from radia_mcp.gmsh.equal_flux import trace_field_lines

pytestmark = pytest.mark.skipif(importlib.util.find_spec("gmsh") is None,
                                reason="gmsh package not installed")

PATCH = {"origin": [-0.4, -0.4, 0.0], "e_xi": [0.8, 0.0, 0.0], "e_eta": [0.0, 0.8, 0.0],
         "n_u": 3, "n_v": 2, "table_nodes": 5}
PLANE_Y0 = {"origin": [-0.4, 0.0, -0.2], "e_xi": [0.8, 0.0, 0.0], "e_eta": [0.0, 0.0, 0.4],
            "n_lines": 4, "table_nodes": 5}
TRACE = {"step_size": 0.05, "max_steps": 40, "return_points": True}


def _bz(p):
    return 1 + p[0] / 2 + p[1] / 4


def _write(tmp_path, name, rows, scale=1.0):
    cube = _CUBE
    if scale != 1.0:
        lines = cube.split("\n")
        start = lines.index("$Nodes") + 3 + len(_P)
        for k in range(start, start + len(_P)):
            lines[k] = " ".join(repr(float(v) * scale) for v in lines[k].split())
        cube = "\n".join(lines)
    path = tmp_path / name
    path.write_text(cube + field_section("NodeData", "B", 3, rows), encoding="utf-8")
    return path


def _msh(tmp_path, scale=1.0, name="linear.msh"):
    """The unit cube in file units of 1/scale metres; B depends on metres."""
    return _write(tmp_path, name, [(i, [0.0, 0.0, float(_bz(p))]) for i, p in _P.items()], scale)


def test_seed_line_is_uniform_and_labelled_not_flux_weighted(tmp_path):
    result = trace_field_lines(_msh(tmp_path), seed_line={"start": [-0.4, 0, 0], "end": [0.4, 0, 0],
                                                          "n": 3}, **TRACE)
    assert result["ok"] is True, result.get("error")
    assert (result["seeding"]["method"], result["seeding"]["flux_weighted"]) == ("uniform_line", False)
    starts = [np.asarray(line["points"])[0, 0] for line in result["polylines"]]
    np.testing.assert_allclose(sorted(starts), [-0.4, 0.0, 0.4], atol=1e-9)


@pytest.mark.parametrize("kwargs, message", [
    ({}, "exactly one of seed_patch"),
    ({"seed_patch": PATCH, "seed_line": {"start": [0, 0, 0], "end": [1, 0, 0], "n": 2}},
     "exactly one of seed_patch"),
    ({"seed_patch": {**PATCH, "seeds": []}}, "unknown ['seeds']"),
    ({"seed_patch": {"origin": [0, 0, 0], "n_lines": 2}}, "missing ['e_eta', 'e_xi']"),
    ({"seed_patch": {**PATCH, "n_lines": 4}}, "either n_lines"),
    ({"seed_patch": {**PATCH, "sign": "positive"}}, "sign applies to n_lines only"),
    ({"seed_patch": {**PATCH, "shape": "sphere"}}, "shape must be one of"),
    ({"seed_patch": {"shape": "cylinder", "axis": [0, 0, 1], "radius": 1, "n_lines": 2}},
     "missing ['axis_origin']"),
    ({"seed_line": {"start": [0, 0, 0], "end": [1, 0, 0]}}, "seed_line needs exactly"),
    ({"seed_patch": PATCH, "view": 0, "seeds": [[0, 0, 0]]}, "unsupported tracer options ['seeds']"),
])
def test_ambiguous_or_malformed_seeding_fails_before_probing(tmp_path, kwargs, message):
    absent = tmp_path / "absent.msh"
    absent.write_text(_CUBE, encoding="utf-8")  # no view: probing would fail differently
    result = trace_field_lines(absent, **kwargs)
    assert result["ok"] is False
    assert message in result["error"], result["error"]


DISK = {"shape": "annulus", "center": [0, 0, 0], "normal": [0, 0, 1], "r_outer": 0.3,
        "n_lines": 3}


@pytest.mark.parametrize("seeding, message", [
    ({"seed_patch": {**PATCH, "n_u": None}}, "seed_patch n_u must be an integer in [1, 1000], got None"),
    ({"seed_patch": {**PATCH, "table_nodes": None}}, "seed_patch table_nodes must be an integer in [2, 316]"),
    ({"seed_patch": {**PATCH, "quadrature_order": -1}}, "quadrature_order must be an integer in [0, 64]"),
    ({"seed_patch": {**PLANE_Y0, "n_lines": True}}, "n_lines must be an integer in [1, 1000], got True"),
    ({"seed_patch": {**PLANE_Y0, "n_lines": 2.5}}, "n_lines must be an integer"),
    ({"seed_patch": {**PLANE_Y0, "sign": None}}, "seed_patch sign must be"),
    ({"seed_patch": {**PATCH, "origin": [0, 0]}}, "seed_patch origin must be a list of 3"),
    ({"seed_patch": {**PATCH, "e_xi": [0, float("nan"), 0]}}, "seed_patch e_xi must be a list of 3"),
    ({"seed_patch": {**PATCH, "origin": None}}, "seed_patch origin must be a list of 3"),
    ({"seed_patch": {**DISK, "r_outer": "0.3"}}, "seed_patch r_outer must be a finite number"),
    ({"seed_patch": {**DISK, "theta_range": [0]}}, "seed_patch theta_range must be"),
    ({"seed_patch": {**DISK, "shape": ["annulus"]}}, "shape must be one of"),
    ({"seed_line": {"start": [0, 0, 0], "end": [1, 0, 0], "n": None}},
     "seed_line n must be an integer in [1, 1000], got None"),
    ({"seed_line": {"start": None, "end": [1, 0, 0], "n": 2}}, "seed_line start must be a list"),
    ({"seed_line": None, "seed_patch": None}, "exactly one of seed_patch"),
    ({"seed_patch": {**DISK, "r_outer": 2**64 * 10**300}}, "seed_patch r_outer must be a finite number"),
    ({"seed_patch": {**PATCH, "origin": [10**400, 0, 0]}}, "seed_patch origin must be a list of 3"),
    ({"seed_patch": {**PLANE_Y0, "n_lines": 10**400}},
     "n_lines must be an integer in [1, 1000], got 'a huge integer'"),
    ({"seed_patch": {**PATCH, "table_nodes": 2**64}}, "table_nodes must be an integer in [2, 316]"),
    ({"seed_patch": {**PATCH, "table_nodes": 300, "quadrature_order": 8}},
     "needs about 5,811,664 field samples"),
    ({"seed_patch": {**PATCH, "n_u": 40, "n_v": 30}}, "n_u x n_v = 1200 exceeds 1000 lines"),
])
def test_malformed_values_return_field_errors_before_probing(tmp_path, monkeypatch,
                                                             seeding, message):
    """Schema-accepted JSON values (null, bool, nan, wrong length, huge ints) never reach int()."""
    def unexpected(*args, **kwargs):
        pytest.fail("a malformed payload reached the probe or the tracer")

    monkeypatch.setattr("radia_mcp.gmsh.equal_flux.probe_field", unexpected)
    monkeypatch.setattr("radia_mcp.gmsh.equal_flux.streamlines", unexpected)
    from radia_mcp.gmsh.server import gmsh_equal_flux_seeds, gmsh_streamlines
    trace = getattr(gmsh_streamlines, "fn", gmsh_streamlines)
    result = trace(str(_msh(tmp_path)), **seeding)
    assert result["ok"] is False and message in result["error"], result
    if seeding.get("seed_patch") is not None:
        place = getattr(gmsh_equal_flux_seeds, "fn", gmsh_equal_flux_seeds)
        placed = place(str(_msh(tmp_path)), seeding["seed_patch"])
        assert placed["ok"] is False and message in placed["error"], placed
