"""Tests for radia_mcp.gmsh.post_process against analytic fields.

The base fixture carries T = 1 + x + 2y + 3z on a unit tet (nodes
1..4 hold 1,2,3,4) with a second step at T+0.5, so every probe,
integral, and transform has an exact expected value.
"""

import importlib.util

import pytest

from radia_mcp.gmsh.post_process import (
    cut_plane_extract,
    harmonic_to_time,
    integrate_view,
    isosurface,
    line_profile,
    math_eval,
    probe_field,
    streamlines,
)

_GMSH_AVAILABLE = importlib.util.find_spec("gmsh") is not None
_MPL_AVAILABLE = importlib.util.find_spec("matplotlib") is not None

pytestmark = pytest.mark.skipif(not _GMSH_AVAILABLE,
                                reason="gmsh package not installed")

_TET_HEADER = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 4 1 4
3 1 0 4
1
2
3
4
0 0 0
1 0 0
0 1 0
0 0 1
$EndNodes
$Elements
1 1 1 1
3 1 4 1
1 1 2 3 4
$EndElements
"""

_SCALAR_2STEP = _TET_HEADER + """$NodeData
1
"T"
1
0.0
3
0
1
4
1 1.0
2 2.0
3 3.0
4 4.0
$EndNodeData
$NodeData
1
"T"
1
0.5
3
1
1
4
1 1.5
2 2.5
3 3.5
4 4.5
$EndNodeData
"""

_VECTOR_UNIFORM_Z = _TET_HEADER + """$NodeData
1
"D"
1
0.0
3
0
3
4
1 0.0 0.0 1.0
2 0.0 0.0 1.0
3 0.0 0.0 1.0
4 0.0 0.0 1.0
$EndNodeData
"""

_VECTOR_MIXED = _TET_HEADER + """$NodeData
1
"D"
1
0.0
3
0
3
4
1 3.0 4.0 0.0
2 0.0 0.0 1.0
3 1.0 0.0 0.0
4 0.0 2.0 0.0
$EndNodeData
"""


def _write(tmp_path, text, name="case.msh"):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def test_probe_interpolates_all_steps_and_reports_outside(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    result = probe_field(msh, [[0.25, 0.25, 0.25], [5.0, 5.0, 5.0]])

    assert result["ok"] is True, result.get("error")
    view = result["views"][0]
    assert view["name"] == "T"
    inside = view["points"][0]
    assert inside["found"] is True
    assert inside["steps"][0] == [pytest.approx(2.5)]
    assert inside["steps"][1] == [pytest.approx(3.0)]
    outside = view["points"][1]
    assert outside["found"] is False
    assert outside["distance"] > 1.0


def test_line_profile_recovers_linear_field(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    # interior-leaning segment: T = 1 + x + 2y + 3z
    result = line_profile(msh, [0.1, 0.1, 0.1], [0.1, 0.1, 0.6], n=6,
                          view="T", step=0)

    assert result["ok"] is True, result.get("error")
    assert result["length"] == pytest.approx(0.5)
    pts = result["views"][0]["points"]
    assert all(p["found"] for p in pts)
    # T along the line: 1 + 0.1 + 0.2 + 3z = 1.3 + 3z, z = 0.1..0.6
    expected = [1.6 + 1.5 * i / 5 for i in range(6)]
    got = [p["values"][0] for p in pts]
    assert got == pytest.approx(expected)


@pytest.mark.skipif(not _MPL_AVAILABLE, reason="matplotlib not installed")
def test_line_profile_writes_plot(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    png = tmp_path / "profile.png"
    result = line_profile(msh, [0.1, 0.1, 0.1], [0.1, 0.1, 0.6], n=6,
                          plot_png=png)
    assert result["ok"] is True, result.get("error")
    assert png.is_file()
    assert png.stat().st_size > 5000


# Mixed-dimension variant needs $Entities so gmsh can resolve the
# surface entity the boundary tri lives on.
_SCALAR_MIXED_DIM = _SCALAR_2STEP.replace(
    "$EndMeshFormat\n",
    "$EndMeshFormat\n"
    "$Entities\n0 0 1 1\n"
    "1 0 0 0 1 1 0 0 0\n"
    "1 0 0 0 1 1 1 0 1 1\n"
    "$EndEntities\n").replace(
    "$Elements\n1 1 1 1\n3 1 4 1\n1 1 2 3 4\n$EndElements",
    "$Elements\n2 2 1 2\n2 1 2 1\n1 1 2 3\n3 1 4 1\n2 1 2 3 4\n$EndElements")


def test_integrate_volume_matches_analytic(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)

    vol = integrate_view(msh, dimension=3)
    assert vol["ok"] is True, vol.get("error")
    # linear field on a simplex: integral = vertex mean * volume
    assert vol["integral_per_step"][0] == pytest.approx(2.5 / 6.0, rel=1e-9)
    assert vol["integral_per_step"][1] == pytest.approx(3.0 / 6.0, rel=1e-9)

    # dimension=-1 SUMS all element dimensions (locked, documented trap):
    # with a boundary tri present, the volume AND surface integrals add up.
    mixed = _write(tmp_path, _SCALAR_MIXED_DIM, "mixed.msh")
    both = integrate_view(mixed, dimension=-1)
    assert both["integral_per_step"][0] == pytest.approx(
        2.5 / 6.0 + 2.0 * 0.5, rel=1e-9)  # + surface tri mean 2 * area 0.5
    only_vol = integrate_view(mixed, dimension=3)
    assert only_vol["integral_per_step"][0] == pytest.approx(2.5 / 6.0,
                                                             rel=1e-9)


def test_math_eval_magnitude_is_nodal(tmp_path):
    msh = _write(tmp_path, _VECTOR_MIXED)
    out = tmp_path / "mag.pos"
    result = math_eval(msh, ["Sqrt(v0^2+v1^2+v2^2)"], out_file=out,
                       result_name="magD")

    assert result["ok"] is True, result.get("error")
    assert out.is_file()

    # nodal semantics: |D| at node 1 = |(3,4,0)| = 5 exactly
    probe = probe_field(out, [[0.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    assert probe["ok"] is True
    values = [p["values"][0] for p in probe["views"][0]["points"]]
    assert values[0] == pytest.approx(5.0)
    assert values[1] == pytest.approx(1.0)


def test_isosurface_extracts_midplane(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    out = tmp_path / "iso.pos"
    result = isosurface(msh, 2.5, view="T", out_file=out)

    assert result["ok"] is True, result.get("error")
    assert out.is_file()
    assert result["pieces"].get("ST", 0) >= 1  # triangles inside the tet


def test_cut_plane_extracts_section(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    out = tmp_path / "cut.pos"
    result = cut_plane_extract(msh, [0, 0, 1], -0.25, out_file=out)

    assert result["ok"] is True, result.get("error")
    assert out.is_file()
    assert result["pieces"].get("ST", 0) >= 1

    section = probe_field(out, [[0.2, 0.2, 0.25]])
    assert section["ok"] is True
    entry = section["views"][0]["points"][0]
    assert entry["found"] is True
    # T on the z=0.25 plane at (0.2, 0.2): 1 + 0.2 + 0.4 + 0.75 = 2.35
    assert entry["steps"][0] == [pytest.approx(2.35)]


def test_harmonic_to_time_expands_cosine(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)  # step0 = re, step1 = im
    out = tmp_path / "time.pos"
    result = harmonic_to_time(msh, n_steps=8, out_file=out)

    assert result["ok"] is True, result.get("error")
    assert result["n_steps"] == 8

    probe = probe_field(out, [[0.25, 0.25, 0.25]])
    steps = probe["views"][0]["points"][0]["steps"]
    assert len(steps) == 8
    # t=0: re*cos(0) - im*sin(0) = re = 2.5
    assert steps[0] == [pytest.approx(2.5)]
    # quarter period: -im = -3.0
    assert steps[2] == [pytest.approx(-3.0)]


def test_streamlines_trace_straight_in_uniform_field(tmp_path):
    msh = _write(tmp_path, _VECTOR_UNIFORM_Z)
    out = tmp_path / "stream.pos"
    result = streamlines(msh, [0.15, 0.15, 0.05], [0.25, 0.25, 0.05],
                         n_seeds=2, step_size=0.05, max_steps=50,
                         return_points=True, out_file=out)

    assert result["ok"] is True, result.get("error")
    assert out.is_file()
    assert result["n_polylines"] >= 2
    assert result["n_segments"] > 5
    # uniform (0,0,1) field: every polyline is a vertical line through
    # its seed, and the carried |v| is exactly 1
    for line in result["polylines"]:
        pts = line["points"]
        assert len(pts) >= 2
        x0, y0 = pts[0][0], pts[0][1]
        for p in pts:
            assert p[0] == pytest.approx(x0, abs=1e-9)
            assert p[1] == pytest.approx(y0, abs=1e-9)
        zs = [p[2] for p in pts]
        assert zs == sorted(zs)  # forward+backward merged, ascending
        # the tracer STOPS at the data boundary (vector views report
        # nearest values at any distance; the distance gate ends the
        # march): the tet spans 0 <= z <= 1 - x - y
        assert zs[0] >= -0.051
        assert zs[-1] <= 1.0 - x0 - y0 + 0.051
        assert all(n == pytest.approx(1.0) for n in line["field_norms"])
    longest = max(len(line["points"]) for line in result["polylines"])
    assert longest > 5  # the forward march crosses the tet interior


def test_probe_missing_view_lists_available(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    result = probe_field(msh, [[0.2, 0.2, 0.2]], view="nope")
    assert result["ok"] is False
    assert "'T'" in result["error"] or "T" in result["error"]
