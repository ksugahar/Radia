"""Tests for the streamline upgrades and flux-line extraction.

Fixtures: a Kuhn 6-tet cube [-0.5, 0.5]^3 carrying analytic linear
fields -- the rotational field (-y, x, 0) has EXACT circular field
lines (P1 interpolation is exact for linear fields), so closure,
radius preservation, evenly spaced ring placement, and termination
reasons all have exact expected answers.
"""

import importlib.util
import math
from itertools import pairwise

import pytest
from radia_mcp.gmsh.post_process import (
    flux_lines,
    probe_field,
    streamlines,
    streamlines_2d,
)

_GMSH_AVAILABLE = importlib.util.find_spec("gmsh") is not None

pytestmark = pytest.mark.skipif(not _GMSH_AVAILABLE,
                                reason="gmsh package not installed")

# Kuhn decomposition of the cube [-0.5, 0.5]^3: 6 positively oriented
# tets sharing the diagonal node1 -> node8.
_CUBE = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 8 1 8
3 1 0 8
1
2
3
4
5
6
7
8
-0.5 -0.5 -0.5
0.5 -0.5 -0.5
-0.5 0.5 -0.5
0.5 0.5 -0.5
-0.5 -0.5 0.5
0.5 -0.5 0.5
-0.5 0.5 0.5
0.5 0.5 0.5
$EndNodes
$Elements
1 6 1 6
3 1 4 6
1 1 2 4 8
2 1 2 8 6
3 1 3 8 4
4 1 3 7 8
5 1 5 6 8
6 1 5 8 7
$EndElements
"""

_P = {1: (-0.5, -0.5, -0.5), 2: (0.5, -0.5, -0.5), 3: (-0.5, 0.5, -0.5),
      4: (0.5, 0.5, -0.5), 5: (-0.5, -0.5, 0.5), 6: (0.5, -0.5, 0.5),
      7: (-0.5, 0.5, 0.5), 8: (0.5, 0.5, 0.5)}


def _nodedata(name, ncomp, rows, time=0.0, step=0):
    lines = ["$NodeData", "1", f'"{name}"', "1", str(time), "3",
             str(step), str(ncomp), str(len(rows))]
    for tag, vals in rows:
        lines.append(str(tag) + " " + " ".join(f"{v:.16g}" for v in vals))
    lines.append("$EndNodeData")
    return "\n".join(lines) + "\n"


_ROT = _CUBE + _nodedata(  # (-y, x, 0): field lines are circles
    "B", 3, [(i, [-p[1], p[0], 0.0]) for i, p in _P.items()])

_RAD = _CUBE + _nodedata(  # (x, y, z): zero at the origin
    "R", 3, [(i, list(map(float, p))) for i, p in _P.items()])

_UNI = _CUBE + _nodedata(  # uniform (1, 0, 0)
    "U", 3, [(i, [1.0, 0.0, 0.0]) for i in _P])

# unit square in z=0 (2 tris) with A = x: isolines are vertical lines
_SQUARE_AX = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Entities
0 0 1 0
1 0 0 0 1 1 0 0 0
$EndEntities
$Nodes
1 4 1 4
2 1 0 4
1
2
3
4
0 0 0
1 0 0
1 1 0
0 1 0
$EndNodes
$Elements
1 2 1 2
2 1 2 2
1 1 2 3
2 1 3 4
$EndElements
""" + _nodedata("A", 1, [(1, [0.0]), (2, [1.0]), (3, [1.0]), (4, [0.0])])


def _write(tmp_path, text, name="case.msh"):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


# ----------------------------------------------------------------------
# 3D tracer: closure, reasons, arrows
# ----------------------------------------------------------------------

def test_closed_field_line_detected_and_exact(tmp_path):
    msh = _write(tmp_path, _ROT, "rot.msh")
    out = tmp_path / "loop.pos"
    result = streamlines(msh, [0.3, 0.0, 0.0], [0.3, 0.0, 0.0],
                         n_seeds=1, step_size=0.02, max_steps=400,
                         return_points=True, out_file=out)
    assert result["ok"] is True, result.get("error")
    assert result["n_closed"] == 1
    line = result["lines"][0]
    assert line["closed"] is True
    assert line["reasons"] == {"forward": "closed"}
    # circumference 2*pi*0.3 within 2%
    assert line["arc_length"] == pytest.approx(2 * math.pi * 0.3,
                                               rel=0.02)
    pts = result["polylines"][0]["points"]
    assert pts[0] == pytest.approx(pts[-1])  # exactly closed
    for p in pts:
        r = (p[0] ** 2 + p[1] ** 2) ** 0.5
        assert r == pytest.approx(0.3, abs=2e-3)  # circles stay circles
        assert abs(p[2]) < 1e-9


def test_termination_reasons_left_data_and_stagnation(tmp_path):
    msh = _write(tmp_path, _RAD, "rad.msh")
    result = streamlines(msh, [0.3, 0.0, 0.0], [0.3, 0.0, 0.0],
                         n_seeds=1, step_size=0.02, max_steps=200,
                         return_points=True,
                         out_file=tmp_path / "rad.pos")
    assert result["ok"] is True, result.get("error")
    line = result["lines"][0]
    assert line["reasons"]["forward"] == "left_data"
    assert line["reasons"]["backward"] == "stagnation"
    # the backward march converges onto the field zero at the origin
    first = result["polylines"][0]["points"][0]
    assert (first[0] ** 2 + first[1] ** 2 + first[2] ** 2) ** 0.5 < 0.05


def test_streamline_arrows_companion_view(tmp_path):
    msh = _write(tmp_path, _ROT, "rot.msh")
    out = tmp_path / "arrows.pos"
    result = streamlines(msh, [0.25, 0.0, 0.0], [0.4, 0.0, 0.0],
                         n_seeds=2, step_size=0.02, arrows_every=10,
                         out_file=out)
    assert result["ok"] is True, result.get("error")
    assert result["n_arrows"] > 0
    probe = probe_field(out, [[0.0, 0.0, 0.0]])
    names = [v["name"] for v in probe["views"]]
    assert "streamlines" in names
    assert "streamline_arrows" in names


# ----------------------------------------------------------------------
# evenly spaced streamlines on a plane (Jobard-Lefer)
# ----------------------------------------------------------------------

def test_streamlines_2d_concentric_circles(tmp_path):
    msh = _write(tmp_path, _ROT, "rot.msh")
    out = tmp_path / "rings.pos"
    result = streamlines_2d(msh, [-0.45, -0.45, 0.0],
                            [0.45, -0.45, 0.0], [-0.45, 0.45, 0.0],
                            d_sep=0.12, return_points=True,
                            out_file=out)
    assert result["ok"] is True, result.get("error")
    assert result["n_lines"] >= 3
    assert result["n_closed"] >= 2
    assert result["budget_exceeded"] is False
    for line in result["lines"]:
        pts = line["points"]
        radii = [(p[0] ** 2 + p[1] ** 2) ** 0.5 for p in pts]
        if line["closed"]:
            # closed rings are true circles: radius spread stays tiny
            assert max(radii) - min(radii) < 0.02
        for p in pts:
            assert abs(p[2]) < 1e-9  # stays in the z=0 plane


def test_streamlines_2d_uniform_field_even_spacing(tmp_path):
    msh = _write(tmp_path, _UNI, "uni.msh")
    result = streamlines_2d(msh, [-0.45, -0.45, 0.0],
                            [0.45, -0.45, 0.0], [-0.45, 0.45, 0.0],
                            d_sep=0.12, return_points=True,
                            out_file=tmp_path / "uni.pos")
    assert result["ok"] is True, result.get("error")
    # extent 0.9 / d_sep 0.12 -> about 7-8 parallel lines
    assert 5 <= result["n_lines"] <= 10
    ys = []
    for line in result["lines"]:
        pts = line["points"]
        y0 = pts[0][1]
        for p in pts:
            assert p[1] == pytest.approx(y0, abs=1e-6)  # straight
        ys.append(y0)
    ys.sort()
    gaps = [b - a for a, b in pairwise(ys)]
    # Jobard-Lefer spacing: every gap within [0.5, 1.6] * d_sep
    assert all(0.06 <= g <= 0.20 for g in gaps), gaps


def test_streamlines_2d_honors_total_step_budget(tmp_path):
    msh = _write(tmp_path, _UNI, "budget.msh")
    result = streamlines_2d(
        msh, [-0.45, -0.45, 0.0], [0.45, -0.45, 0.0],
        [-0.45, 0.45, 0.0], d_sep=0.12, max_steps=100,
        max_total_steps=3, out_file=tmp_path / "budget.pos")
    assert result["ok"] is True, result.get("error")
    assert result["budget_exceeded"] is True
    assert result["steps_used"] == 3


# ----------------------------------------------------------------------
# flux lines (multi-level isolines)
# ----------------------------------------------------------------------

def test_flux_lines_equal_levels_on_az(tmp_path):
    msh = _write(tmp_path, _SQUARE_AX, "square.msh")
    out = tmp_path / "flux.pos"
    result = flux_lines(msh, n_levels=3, out_file=out)
    assert result["ok"] is True, result.get("error")
    assert result["levels"] == pytest.approx([0.25, 0.5, 0.75])
    assert result["pieces"].get("SL", 0) >= 3
    # every level produced at least one line piece
    assert all(per.get("SL", 0) >= 1
               for per in result["pieces_per_level"])
    # a point ON the mid line carries the level value
    probe = probe_field(out, [[0.5, 0.3, 0.0]])
    entry = probe["views"][0]["points"][0]
    assert entry["found"] is True
    row = entry.get("values") or entry["steps"][0]
    assert row == [pytest.approx(0.5)]


def test_flux_lines_explicit_levels_3d_isosurfaces(tmp_path):
    # 3D scalar: the same verb stacks isosurfaces (levels echoed)
    cube_t = _CUBE + _nodedata(
        "T", 1, [(i, [1.0 + p[0] + 2 * p[1] + 3 * p[2]])
                 for i, p in _P.items()])
    msh = _write(tmp_path, cube_t, "cube_t.msh")
    out = tmp_path / "iso3.pos"
    result = flux_lines(msh, levels=[0.5, 1.0, 1.5], out_file=out)
    assert result["ok"] is True, result.get("error")
    assert result["levels"] == pytest.approx([0.5, 1.0, 1.5])
    assert result["pieces"].get("ST", 0) >= 3


# ----------------------------------------------------------------------
# adaptive isosurface extraction on high-order (TET10) data
# ----------------------------------------------------------------------

def _tet10_r2_msh():
    """Unit TET10 carrying the EXACT quadratic T = x^2 + y^2 + z^2."""
    corners = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 1, 0), 4: (0, 0, 1)}
    edges = {5: (1, 2), 6: (2, 3), 7: (1, 3), 8: (1, 4), 9: (3, 4),
             10: (2, 4)}
    coords = dict(corners)
    for nid, (a, b) in edges.items():
        pa, pb = corners[a], corners[b]
        coords[nid] = tuple((pa[k] + pb[k]) / 2 for k in range(3))
    lines = ["$MeshFormat", "4.1 0 8", "$EndMeshFormat", "$Nodes",
             "1 10 1 10", "3 1 0 10"]
    lines += [str(i) for i in range(1, 11)]
    lines += [f"{coords[i][0]} {coords[i][1]} {coords[i][2]}"
              for i in range(1, 11)]
    lines += ["$EndNodes", "$Elements", "1 1 1 1", "3 1 11 1",
              "1 " + " ".join(str(i) for i in range(1, 11)),
              "$EndElements"]
    rows = [(i, [sum(c * c for c in coords[i])]) for i in range(1, 11)]
    return "\n".join(lines) + "\n" + _nodedata("T", 1, rows)


def _st_radial_errors(pos_path, r0):
    import re
    errs = []
    text = pos_path.read_text(encoding="utf-8", errors="replace")
    for m in re.finditer(r"ST\(([^)]*)\)", text):
        nums = [float(t) for t in m.group(1).split(",")]
        for k in range(3):
            r = (nums[3 * k] ** 2 + nums[3 * k + 1] ** 2
                 + nums[3 * k + 2] ** 2) ** 0.5
            errs.append(abs(r - r0))
    return errs


def test_isosurface_adaptive_follows_quadratic_field(tmp_path):
    from radia_mcp.gmsh.post_process import isosurface

    msh = _write(tmp_path, _tet10_r2_msh(), "tet10.msh")
    # plain P1 cut: one flat triangle, chord error ~0.21 on r = 0.3
    flat_out = tmp_path / "flat.pos"
    flat = isosurface(msh, 0.09, out_file=flat_out)
    assert flat["ok"] is True, flat.get("error")
    assert flat["pieces"].get("ST") == 1
    flat_errs = _st_radial_errors(flat_out, 0.3)
    assert max(flat_errs) > 0.15

    # adaptive: the extraction follows the quadratic interpolant
    fine_out = tmp_path / "fine.pos"
    fine = isosurface(msh, 0.09, recur_level=4, target_error=1e-6,
                      out_file=fine_out)
    assert fine["ok"] is True, fine.get("error")
    assert fine["pieces"].get("ST", 0) > 50
    fine_errs = _st_radial_errors(fine_out, 0.3)
    assert max(fine_errs) < 0.01


def test_flux_lines_pass_recur_level_through(tmp_path):
    msh = _write(tmp_path, _tet10_r2_msh(), "tet10b.msh")
    out = tmp_path / "adaptive_stack.pos"
    result = flux_lines(msh, levels=[0.09], recur_level=3, out_file=out)
    assert result["ok"] is True, result.get("error")
    assert result["pieces"].get("ST", 0) > 10  # 1 without adaptivity
