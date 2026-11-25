"""Field-line tracing seeds by equal flux unless uniform seeds are asked for.

Guards: `trace_field_lines` and the `gmsh_streamlines` tool place seed_patch
seeds with radia.flux_seeding (each line carries total/(n_u*n_v) Wb) and
trace exactly those seeds; seed_line stays uniform and is labelled not
flux-weighted; ambiguous or malformed seeding fails before probing; the
stated coordinate unit scales fluxes to Wb.
Assumes: P1 node data reproduce B = (0, 0, 1 + x/2 + y/4) with x, y in
metres exactly, so the patch [-0.4, 0.4]^2 at z = 0 carries 0.64 Wb.
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
        return importlib.import_module("radia.flux_seeding")
    except ImportError as exc:
        pytest.skip(f"radia {importlib.metadata.version('radia')} lacks radia.flux_seeding ({exc})")


from radia_mcp.gmsh.equal_flux import equal_flux_seeds_on_plane, trace_field_lines  # noqa: E402

PATCH = {"origin": [-0.4, -0.4, 0.0], "e_xi": [0.8, 0.0, 0.0], "e_eta": [0.0, 0.8, 0.0],
         "n_u": 3, "n_v": 2, "table_nodes": 5}
TRACE = {"step_size": 0.05, "max_steps": 40, "return_points": True}


def _bz(p):
    return 1 + p[0] / 2 + p[1] / 4


def _msh(tmp_path, scale=1.0, name="linear.msh"):
    """The unit cube in file units of 1/scale metres; B depends on metres."""
    cube = _CUBE
    if scale != 1.0:
        lines = cube.split("\n")
        start = lines.index("$Nodes") + 3 + len(_P)
        for k in range(start, start + len(_P)):
            lines[k] = " ".join(repr(float(v) * scale) for v in lines[k].split())
        cube = "\n".join(lines)
    rows = [(i, [0.0, 0.0, float(_bz(p))]) for i, p in _P.items()]
    path = tmp_path / name
    path.write_text(cube + field_section("NodeData", "B", 3, rows), encoding="utf-8")
    return path


def test_seed_patch_traces_equal_flux_lines(tmp_path):
    _flux_seeding()
    path = _msh(tmp_path)
    result = trace_field_lines(path, seed_patch=PATCH, **TRACE)
    assert result["ok"] is True, result.get("error")
    seeding = result["seeding"]
    assert (seeding["method"], seeding["flux_weighted"]) == ("equal_flux", True)
    assert seeding["total_flux_Wb"] == pytest.approx(0.64)
    assert seeding["flux_per_line_Wb"] == pytest.approx(0.64 / 6)
    assert seeding["equal_flux_quadrature"]["max_relative_cell_deviation"] < 1e-12
    # The traced lines start from exactly the placed equal-flux seeds.
    placed = equal_flux_seeds_on_plane(path, PATCH["origin"], PATCH["e_xi"], PATCH["e_eta"],
                                       n_u=3, n_v=2, table_nodes=5, quadrature_order=0)
    expected = np.asarray([s["point"] for s in placed["seeds"]])
    assert len(result["polylines"]) == 6
    for seed, line in zip(expected, result["polylines"]):
        xy = np.asarray(line["points"])[:, :2]
        np.testing.assert_allclose(xy, np.broadcast_to(seed[:2], xy.shape), atol=1e-9)


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
    ({"seed_patch": {"origin": [0, 0, 0]}}, "missing ['e_eta', 'e_xi', 'n_u', 'n_v']"),
    ({"seed_line": {"start": [0, 0, 0], "end": [1, 0, 0]}}, "seed_line needs exactly"),
    ({"seed_patch": PATCH, "view": 0, "seeds": [[0, 0, 0]]}, "unsupported tracer options ['seeds']"),
])
def test_ambiguous_or_malformed_seeding_fails_before_probing(tmp_path, kwargs, message):
    result = trace_field_lines(tmp_path / "absent.msh", **kwargs)
    assert result["ok"] is False
    # Validation runs first, so a missing file is not what is reported.
    assert message in result["error"], result["error"]


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
    for a, b, line in zip(metres["seeding"]["seeds"], mm["seeding"]["seeds"], mm["polylines"]):
        np.testing.assert_allclose(b["point"], 1000 * np.asarray(a["point"]), atol=1e-9)
        np.testing.assert_allclose(b["point_m"], a["point_m"], atol=1e-12)
        # Lines are traced from the seeds in the file's own (mm) coordinates.
        xy = np.asarray(line["points"])[:, :2]
        np.testing.assert_allclose(xy, np.broadcast_to(b["point"][:2], xy.shape), atol=1e-6)
    # Reading the same mm file as metres would inflate the flux a million times.
    unstated = trace_field_lines(mm_file, seed_patch=mm_patch, step_size=50.0, max_steps=40)
    assert unstated["seeding"]["total_flux_Wb"] == pytest.approx(1e6 * metres["seeding"]["total_flux_Wb"])
    bad = trace_field_lines(mm_file, seed_patch=mm_patch, length_unit="inch")
    assert bad["ok"] is False and "length_unit" in bad["error"]


def test_mcp_tool_defaults_to_equal_flux_seeding(tmp_path):
    _flux_seeding()
    from radia_mcp.gmsh.server import gmsh_streamlines
    tool = getattr(gmsh_streamlines, "fn", gmsh_streamlines)
    result = tool(str(_msh(tmp_path)), seed_patch=PATCH, step_size=0.05, max_steps=40)
    assert result["ok"] is True, result.get("error")
    assert result["seeding"]["method"] == "equal_flux" and result["n_polylines"] == 6
    assert tool(str(_msh(tmp_path)))["ok"] is False
