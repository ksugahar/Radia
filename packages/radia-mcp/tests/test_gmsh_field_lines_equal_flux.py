"""Field-line tracing seeds by equal flux unless uniform seeds are asked for.

Guards: `trace_field_lines` and the `gmsh_streamlines` tool place seed_patch
seeds with radia.flux_seeding (grid or balanced bisection, parallelogram or
annulus) and trace exactly those seeds; seed_line stays uniform and is
labelled not flux-weighted; ambiguous or malformed seeding fails before
probing; the stated unit scales only fluxes; the flux account reports how
lines end, repeated crossings and the drawn share of the opposite flux.
Assumes, with P1 node data exact for linear fields on the Kuhn cube:
- B = (0, 0, 1 + x/2 + y/4): the square [-0.4, 0.4]^2 at z = 0 carries
  0.64 Wb and the disk r < 0.3 at z = 0 carries pi 0.09 Wb;
- B = (-y, x, 0): lines are circles about z, so a line seeded at x < 0 on
  the plane y = 0 returns through x > 0 exactly once, in the opposite sense;
  the patch [-0.4, 0.4] x [-0.2, 0.2] (x, z) carries 0.032 Wb each way;
- B = (-y, x, 0.02): helices that cross that plane again in the same sense.
"""

import importlib
import importlib.metadata
import importlib.util

import numpy as np
import pytest
from msh42_fixtures import field_section
from test_gmsh_flux_streamlines import _CUBE, _P

pytestmark = pytest.mark.skipif(importlib.util.find_spec("gmsh") is None,
                                reason="gmsh package not installed")


def _flux_seeding():
    try:
        module = importlib.import_module("radia.flux_seeding")
    except ImportError as exc:
        pytest.skip(f"radia {importlib.metadata.version('radia')} lacks radia.flux_seeding ({exc})")
    if not hasattr(module, "balanced_equal_flux_seeds"):
        pytest.skip(f"radia {importlib.metadata.version('radia')} lacks balanced seeding")
    return module


from radia_mcp.gmsh.equal_flux import place_equal_flux_seeds, trace_field_lines  # noqa: E402

PATCH = {"origin": [-0.4, -0.4, 0.0], "e_xi": [0.8, 0.0, 0.0], "e_eta": [0.0, 0.8, 0.0],
         "n_u": 3, "n_v": 2, "table_nodes": 5}
PLANE_Y0 = {"origin": [-0.4, 0.0, -0.2], "e_xi": [0.8, 0.0, 0.0], "e_eta": [0.0, 0.0, 0.4],
            "n_lines": 4, "table_nodes": 5}
TRACE = {"step_size": 0.05, "max_steps": 40, "return_points": True}
LOOP = {"step_size": 0.02, "max_steps": 400, "return_points": True}


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


def _rotation(tmp_path, pitch=0.0):
    return _write(tmp_path, f"rot_{pitch}.msh",
                  [(i, [-p[1], p[0], pitch]) for i, p in _P.items()])


def _by_seed(result):
    return {poly["seed_index"]: np.asarray(poly["points"]) for poly in result["polylines"]}


def test_seed_patch_traces_equal_flux_lines(tmp_path):
    _flux_seeding()
    path = _msh(tmp_path)
    result = trace_field_lines(path, seed_patch=PATCH, **TRACE)
    assert result["ok"] is True, result.get("error")
    seeding = result["seeding"]
    assert (seeding["method"], seeding["flux_weighted"], seeding["layout"]) == \
        ("equal_flux", True, "grid")
    assert seeding["total_flux_Wb"] == pytest.approx(0.64)
    assert seeding["flux_per_line_Wb"] == pytest.approx(0.64 / 6)
    assert seeding["equal_flux_quadrature"]["max_relative_cell_deviation"] < 1e-12
    # The traced lines start from exactly the placed equal-flux seeds.
    placed = place_equal_flux_seeds(path, {**PATCH, "quadrature_order": 0})
    lines = _by_seed(result)
    assert sorted(lines) == list(range(6))
    for index, seed in enumerate(placed["seeds"]):
        xy = lines[index][:, :2]
        np.testing.assert_allclose(xy, np.broadcast_to(seed["point"][:2], xy.shape), atol=1e-9)
    account = result["flux_account"]
    assert account["lines_by_end"] == {"closed": 0, "open": 6, "truncated": 0, "stagnation": 0}
    assert account["opposite_crossing_flux_ratio"] is None and result["warnings"] == []


def test_balanced_annulus_seeds_carry_the_disk_flux(tmp_path):
    _flux_seeding()
    disk = {"shape": "annulus", "center": [0, 0, 0], "normal": [0, 0, 1], "r_outer": 0.3,
            "n_lines": 7, "table_nodes": 9}
    result = trace_field_lines(_msh(tmp_path), seed_patch=disk, **TRACE)
    assert result["ok"] is True, result.get("error")
    seeding = result["seeding"]
    assert (seeding["layout"], seeding["n_lines"]) == ("balanced", 7)
    assert seeding["total_flux_Wb"] == pytest.approx(np.pi * 0.09, rel=1e-9)
    points = np.asarray([s["point"] for s in seeding["seeds"]])
    assert np.all(np.hypot(points[:, 0], points[:, 1]) < 0.3) and np.allclose(points[:, 2], 0)
    assert result["n_polylines"] == 7


def test_closed_lines_draw_the_opposite_flux_once(tmp_path):
    """Each seeded loop returns through the patch once with the other sign."""
    _flux_seeding()
    path = _rotation(tmp_path)
    mixed = trace_field_lines(path, seed_patch={k: v for k, v in PLANE_Y0.items()
                                                if k != "n_lines"} | {"n_u": 2, "n_v": 2},
                              **LOOP)
    assert mixed["ok"] is False and "changes sign" in mixed["error"]
    auto = trace_field_lines(path, seed_patch=PLANE_Y0, **LOOP)
    assert auto["ok"] is False and "sign='positive'" in auto["error"]

    result = trace_field_lines(path, seed_patch={**PLANE_Y0, "sign": "positive"}, **LOOP)
    assert result["ok"] is True, result.get("error")
    seeding, account = result["seeding"], result["flux_account"]
    assert (seeding["positive_flux_Wb"], seeding["negative_flux_Wb"]) == \
        pytest.approx((0.032, -0.032), abs=1e-14)
    assert all(s["point"][0] < 0 for s in seeding["seeds"])
    assert account["lines_by_end"]["closed"] == 4
    assert (account["repeat_crossings"], account["opposite_crossings"]) == (0, 4)
    assert account["opposite_crossing_flux_ratio"] == pytest.approx(1.0, abs=1e-12)
    assert result["warnings"] == []


def test_extent_stops_lines_and_loses_the_return_flux(tmp_path):
    _flux_seeding()
    extent = {"min": [-0.5, -0.5, -0.5], "max": [0.1, 0.5, 0.5]}
    result = trace_field_lines(_rotation(tmp_path), seed_patch={**PLANE_Y0, "sign": "positive"},
                               extent=extent, **LOOP)
    assert result["ok"] is True, result.get("error")
    account = result["flux_account"]
    assert account["lines_by_end"]["open"] >= 1, account
    for points in _by_seed(result).values():
        assert points[:, 0].max() <= 0.1
    stopped = [line for line in result["lines"] if "left_extent" in line["reasons"].values()]
    assert len(stopped) == account["lines_by_end"]["open"]
    assert account["opposite_crossing_flux_ratio"] < 1


def test_truncated_and_repeating_lines_are_warned(tmp_path):
    _flux_seeding()
    short = trace_field_lines(_rotation(tmp_path), seed_patch={**PLANE_Y0, "sign": "positive"},
                              step_size=0.02, max_steps=3)
    assert short["flux_account"]["lines_by_end"]["truncated"] == 4
    assert any("max_steps" in w for w in short["warnings"])
    assert short["polylines"] is None  # points were used internally, not returned
    forward = trace_field_lines(_msh(tmp_path), seed_patch=PATCH, both_directions=False, **TRACE)
    assert any("forward only" in w for w in forward["warnings"]), forward["warnings"]

    helix = trace_field_lines(_rotation(tmp_path, pitch=0.02),
                              seed_patch={**PLANE_Y0, "e_eta": [0, 0, 0.8],
                                          "origin": [-0.4, 0, -0.4], "sign": "positive"},
                              step_size=0.02, max_steps=300)
    assert helix["ok"] is True, helix.get("error")
    assert helix["flux_account"]["repeat_crossings"] > 0
    assert any("over-represents" in w for w in helix["warnings"])


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


def test_millimetre_coordinates_give_the_same_flux_in_wb(tmp_path):
    _flux_seeding()
    metres = trace_field_lines(_msh(tmp_path), seed_patch=PATCH, **TRACE)
    mm_patch = {**PATCH, **{k: [1000 * v for v in PATCH[k]] for k in ("origin", "e_xi", "e_eta")}}
    mm_file = _msh(tmp_path, scale=1000.0, name="linear_mm.msh")
    mm = trace_field_lines(mm_file, seed_patch=mm_patch, length_unit="mm",
                           step_size=50.0, max_steps=40, return_points=True)
    assert mm["ok"] is True, mm.get("error")
    assert mm["seeding"]["total_flux_Wb"] == pytest.approx(metres["seeding"]["total_flux_Wb"])
    assert mm["seeding"]["provenance"]["length_unit"] == "mm"
    lines = _by_seed(mm)
    for index, (a, b) in enumerate(zip(metres["seeding"]["seeds"], mm["seeding"]["seeds"])):
        np.testing.assert_allclose(b["point"], 1000 * np.asarray(a["point"]), atol=1e-9)
        np.testing.assert_allclose(b["point_m"], a["point_m"], atol=1e-12)
        # Lines are traced from the seeds in the file's own (mm) coordinates.
        xy = lines[index][:, :2]
        np.testing.assert_allclose(xy, np.broadcast_to(b["point"][:2], xy.shape), atol=1e-6)
    # The unit changes only the reported flux: the same seeds, a million times the Wb.
    unstated = trace_field_lines(mm_file, seed_patch=mm_patch, step_size=50.0, max_steps=40)
    assert unstated["seeding"]["total_flux_Wb"] == pytest.approx(1e6 * metres["seeding"]["total_flux_Wb"])
    np.testing.assert_allclose([s["point"] for s in unstated["seeding"]["seeds"]],
                               [s["point"] for s in mm["seeding"]["seeds"]], atol=1e-9)
    bad = trace_field_lines(mm_file, seed_patch=mm_patch, length_unit="inch")
    assert bad["ok"] is False and "length_unit" in bad["error"]


def test_mcp_tool_defaults_to_equal_flux_seeding(tmp_path):
    _flux_seeding()
    from radia_mcp.gmsh.server import gmsh_equal_flux_seeds, gmsh_streamlines
    tool = getattr(gmsh_streamlines, "fn", gmsh_streamlines)
    result = tool(str(_msh(tmp_path)), seed_patch=PATCH, step_size=0.05, max_steps=40)
    assert result["ok"] is True, result.get("error")
    assert result["seeding"]["method"] == "equal_flux" and result["n_polylines"] == 6
    assert tool(str(_msh(tmp_path)))["ok"] is False
    place = getattr(gmsh_equal_flux_seeds, "fn", gmsh_equal_flux_seeds)
    balanced = {k: v for k, v in PATCH.items() if k not in ("n_u", "n_v")} | {"n_lines": 5}
    placed = place(str(_msh(tmp_path)), balanced)
    assert placed["ok"] is True and placed["layout"] == "balanced" and len(placed["seeds"]) == 5


