"""Tests for the ParaView-parity post verbs against analytic fields.

Fixture: 2-tet mesh (unit tet + its face-neighbor reaching (1,1,1)).
T = 1 + x + 2y + 3z carries two time steps (T, T + 0.5), so gradients,
thresholds, mirror copies, histograms, and CSV dumps all have exact
expected values.  Plugin semantics were measured on gmsh 4.15.2
(Transform / Warp / Smooth / ModulusPhase run IN PLACE; ExtractElements
selects on the ELEMENT MEAN).
"""

import csv
import importlib.util
import math

import pytest

from radia_mcp.gmsh.post_process import (
    curve_profile,
    derived_field,
    export_view_csv,
    extract_skin,
    field_histogram,
    mirror_expand,
    modulus_phase,
    point_history,
    probe_field,
    resample_grid,
    smooth_to_nodes,
    threshold,
    transform_view,
    view_min_max,
    warp_view,
)

_GMSH_AVAILABLE = importlib.util.find_spec("gmsh") is not None
_MPL_AVAILABLE = importlib.util.find_spec("matplotlib") is not None
_PIL_AVAILABLE = importlib.util.find_spec("PIL") is not None

pytestmark = pytest.mark.skipif(not _GMSH_AVAILABLE,
                                reason="gmsh package not installed")

# 2 tets: (1,2,3,4) and (2,3,5,4); p5 = (1,1,1).
_MESH = """$MeshFormat
4.1 0 8
$EndMeshFormat
$Nodes
1 5 1 5
3 1 0 5
1
2
3
4
5
0 0 0
1 0 0
0 1 0
0 0 1
1 1 1
$EndNodes
$Elements
1 2 1 2
3 1 4 2
1 1 2 3 4
2 2 3 5 4
$EndElements
"""

_P = {1: (0, 0, 0), 2: (1, 0, 0), 3: (0, 1, 0), 4: (0, 0, 1), 5: (1, 1, 1)}


def _nodedata(name, ncomp, rows, time, step):
    lines = ["$NodeData", "1", f'"{name}"', "1", str(time), "3",
             str(step), str(ncomp), str(len(rows))]
    for tag, vals in rows:
        lines.append(str(tag) + " " + " ".join(f"{v:.16g}" for v in vals))
    lines.append("$EndNodeData")
    return "\n".join(lines) + "\n"


def _elementdata(name, ncomp, rows, time=0.0, step=0):
    lines = ["$ElementData", "1", f'"{name}"', "1", str(time), "3",
             str(step), str(ncomp), str(len(rows))]
    for tag, vals in rows:
        lines.append(str(tag) + " " + " ".join(f"{v:.16g}" for v in vals))
    lines.append("$EndElementData")
    return "\n".join(lines) + "\n"


def _t(x, y, z):
    return 1.0 + x + 2.0 * y + 3.0 * z


_SCALAR_2STEP = _MESH + _nodedata(
    "T", 1, [(i, [_t(*p)]) for i, p in _P.items()], 0.0, 0) + _nodedata(
    "T", 1, [(i, [_t(*p) + 0.5]) for i, p in _P.items()], 0.5, 1)

_VROT = _MESH + _nodedata(  # (-y, x, 0): curl = (0,0,2)
    "Vrot", 3, [(i, [-p[1], p[0], 0.0]) for i, p in _P.items()], 0.0, 0)

_VDIV = _MESH + _nodedata(  # (x, y, z): div = 3
    "Vdiv", 3, [(i, list(map(float, p))) for i, p in _P.items()], 0.0, 0)

_TENSOR = _MESH + _nodedata(
    "S", 9, [(i, [1, 0, 0, 0, 2, 0, 0, 0, 3]) for i in _P], 0.0, 0)

_COMPLEX = _MESH + _nodedata(
    "Z", 1, [(i, [3.0]) for i in _P], 0.0, 0) + _nodedata(
    "Z", 1, [(i, [4.0]) for i in _P], 1.0, 1)

_EDATA = _MESH + _elementdata("E", 1, [(1, [10.0]), (2, [20.0])])

_VZ = _MESH + _nodedata(
    "W", 3, [(i, [0.0, 0.0, 1.0]) for i in _P], 0.0, 0)


def _write(tmp_path, text, name="case.msh"):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def _probe_one(path, point, **kw):
    result = probe_field(path, [point], **kw)
    assert result["ok"] is True, result.get("error")
    return result["views"][0]["points"][0]


# ----------------------------------------------------------------------
# derived fields
# ----------------------------------------------------------------------

def test_gradient_of_linear_scalar(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    out = tmp_path / "grad.pos"
    result = derived_field(msh, "gradient", view="T", out_file=out,
                           check_point=[0.2, 0.2, 0.2])
    assert result["ok"] is True, result.get("error")
    assert result["probe"][:3] == pytest.approx([1.0, 2.0, 3.0])
    entry = _probe_one(out, [0.6, 0.6, 0.6])
    assert entry["found"] is True
    row = entry.get("values") or entry["steps"][0]
    assert row[:3] == pytest.approx([1.0, 2.0, 3.0])


def test_curl_and_divergence(tmp_path):
    rot = _write(tmp_path, _VROT, "rot.msh")
    result = derived_field(rot, "curl", check_point=[0.2, 0.2, 0.2],
                           out_file=tmp_path / "curl.pos")
    assert result["ok"] is True, result.get("error")
    assert result["probe"] == pytest.approx([0.0, 0.0, 2.0], abs=1e-12)

    div = _write(tmp_path, _VDIV, "div.msh")
    result = derived_field(div, "divergence", check_point=[0.2, 0.2, 0.2],
                           out_file=tmp_path / "div.pos")
    assert result["ok"] is True, result.get("error")
    assert result["probe"] == pytest.approx([3.0])


def test_eigenvalues_writes_three_views(tmp_path):
    msh = _write(tmp_path, _TENSOR)
    out = tmp_path / "eig.pos"
    result = derived_field(msh, "eigenvalues", out_file=out)
    assert result["ok"] is True, result.get("error")
    assert len(result["views"]) == 3

    probe = probe_field(out, [[0.2, 0.2, 0.2]])
    assert probe["ok"] is True
    by_name = {v["name"]: v["points"][0] for v in probe["views"]}
    assert len(by_name) == 3
    values = sorted(entry["values"][0] for entry in by_name.values())
    assert values == pytest.approx([1.0, 2.0, 3.0])


# ----------------------------------------------------------------------
# threshold / skin
# ----------------------------------------------------------------------

def test_threshold_selects_on_element_mean(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    out = tmp_path / "thresh.pos"
    # tet1 nodal T: 1,2,3,4 (mean 2.5); tet2: 2,3,7,4 (mean 4.0)
    result = threshold(msh, 2.6, 99.0, view="T", out_file=out)
    assert result["ok"] is True, result.get("error")
    assert result["n_kept"] == 1

    inside_kept = _probe_one(out, [0.6, 0.6, 0.6])
    assert inside_kept["found"] is True
    dropped = _probe_one(out, [0.1, 0.1, 0.1])
    assert dropped["found"] is False


def test_skin_extracts_boundary_with_field(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    out = tmp_path / "skin.pos"
    result = extract_skin(msh, view="T", out_file=out)
    assert result["ok"] is True, result.get("error")
    # 2 tets x 4 faces - 2 copies of the shared face = 6 boundary tris
    assert result["pieces"].get("ST") == 6
    face = _probe_one(out, [0.2, 0.2, 0.0])
    assert face["found"] is True
    assert face["steps"][0] == [pytest.approx(1.6)]


# ----------------------------------------------------------------------
# mirror expansion / affine transform
# ----------------------------------------------------------------------

def test_mirror_expand_scalar_keeps_both_sides_and_steps(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    out = tmp_path / "full.pos"
    result = mirror_expand(msh, ["x"], parity="scalar", view="T",
                           out_file=out)
    assert result["ok"] is True, result.get("error")
    assert result["n_copies"] == 1
    assert result["n_steps"] == 2

    probe = probe_field(out, [[0.2, 0.2, 0.2], [-0.2, 0.2, 0.2]])
    assert probe["ok"] is True
    assert len(probe["views"]) == 1  # combined into ONE view
    orig, mirrored = probe["views"][0]["points"]
    assert orig["steps"][0] == [pytest.approx(2.2)]
    assert orig["steps"][1] == [pytest.approx(2.7)]
    assert mirrored["steps"][0] == [pytest.approx(2.2)]
    assert mirrored["steps"][1] == [pytest.approx(2.7)]


def test_mirror_expand_vector_parities(tmp_path):
    rot = _write(tmp_path, _VROT, "rot.msh")
    # source value at (0.2, 0.2, 0.2): (-0.2, 0.2, 0)
    polar = mirror_expand(rot, ["x"], parity="vector",
                          out_file=tmp_path / "polar.pos")
    assert polar["ok"] is True, polar.get("error")
    entry = _probe_one(tmp_path / "polar.pos", [-0.2, 0.2, 0.2])
    row = entry.get("values") or entry["steps"][0]
    assert row == pytest.approx([0.2, 0.2, 0.0], abs=1e-12)

    axial = mirror_expand(rot, ["x"], parity="pseudovector",
                          out_file=tmp_path / "axial.pos")
    assert axial["ok"] is True, axial.get("error")
    entry = _probe_one(tmp_path / "axial.pos", [-0.2, 0.2, 0.2])
    row = entry.get("values") or entry["steps"][0]
    assert row == pytest.approx([-0.2, -0.2, 0.0], abs=1e-12)


def test_mirror_expand_two_planes_covers_four_quadrants(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    out = tmp_path / "quad.pos"
    result = mirror_expand(msh, ["x", "y"], view="T", out_file=out)
    assert result["ok"] is True, result.get("error")
    assert result["n_copies"] == 3
    pts = [[0.2, 0.2, 0.2], [-0.2, 0.2, 0.2],
           [0.2, -0.2, 0.2], [-0.2, -0.2, 0.2]]
    probe = probe_field(out, pts)
    for entry in probe["views"][0]["points"]:
        assert entry["found"] is True
        assert entry["steps"][0] == [pytest.approx(2.2)]


def test_transform_view_translates_and_rewrites_values(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    out = tmp_path / "moved.pos"
    result = transform_view(msh, [1, 0, 0, 0, 1, 0, 0, 0, 1],
                            translation=[10.0, 0.0, 0.0], view="T",
                            value_expressions=["2*v0"], out_file=out)
    assert result["ok"] is True, result.get("error")
    entry = _probe_one(out, [10.2, 0.2, 0.2])
    assert entry["found"] is True
    assert entry["steps"][0] == [pytest.approx(4.4)]
    # source location no longer carries data in the output view
    gone = _probe_one(out, [0.2, 0.2, 0.2])
    assert gone["found"] is False


# ----------------------------------------------------------------------
# warp / smooth / modulus-phase / min-max
# ----------------------------------------------------------------------

def test_warp_displaces_geometry(tmp_path):
    msh = _write(tmp_path, _VZ, "vz.msh")
    out = tmp_path / "warp.pos"
    result = warp_view(msh, factor=0.5, out_file=out)
    assert result["ok"] is True, result.get("error")
    shifted = _probe_one(out, [0.2, 0.2, 0.55])
    assert shifted["found"] is True
    old = _probe_one(out, [0.2, 0.2, 0.05])
    assert old["found"] is False


def test_smooth_averages_element_data_to_nodes(tmp_path):
    msh = _write(tmp_path, _EDATA, "edata.msh")
    out = tmp_path / "nodal.pos"
    result = smooth_to_nodes(msh, out_file=out)
    assert result["ok"] is True, result.get("error")
    corner = _probe_one(out, [0.0, 0.0, 0.0])
    row = corner.get("values") or corner["steps"][0]
    assert row == [pytest.approx(10.0)]
    shared = _probe_one(out, [1.0, 0.0, 0.0])
    row = shared.get("values") or shared["steps"][0]
    assert row == [pytest.approx(15.0)]


def test_modulus_phase_two_steps(tmp_path):
    msh = _write(tmp_path, _COMPLEX, "cplx.msh")
    out = tmp_path / "mp.pos"
    result = modulus_phase(msh, out_file=out)
    assert result["ok"] is True, result.get("error")
    assert result["n_steps"] == 2
    entry = _probe_one(out, [0.2, 0.2, 0.2])
    assert entry["steps"][0] == [pytest.approx(5.0)]
    assert entry["steps"][1] == [pytest.approx(math.atan2(4.0, 3.0))]


def test_view_min_max_locates_extrema(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    result = view_min_max(msh, view="T")
    assert result["ok"] is True, result.get("error")
    assert result["min"]["point"] == pytest.approx([0.0, 0.0, 0.0])
    assert result["max"]["point"] == pytest.approx([1.0, 1.0, 1.0])
    assert result["min"]["values"][0] == pytest.approx(1.0)
    assert result["max"]["values"][0] == pytest.approx(7.0)


# ----------------------------------------------------------------------
# curve profile / grid resample
# ----------------------------------------------------------------------

def test_curve_profile_recovers_linear_field(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    csv_out = tmp_path / "curve.csv"
    result = curve_profile(msh, "u", "0.05", "0.05", 0.0, 0.8, n=5,
                           view="T", csv_out=csv_out,
                           out_file=tmp_path / "curve.pos")
    assert result["ok"] is True, result.get("error")
    assert result["u"] == pytest.approx([0.0, 0.2, 0.4, 0.6, 0.8])
    # T = 1 + u + 0.1 + 0.15 = 1.25 + u
    got = [vals[0][0] for vals in result["values"]]
    assert got == pytest.approx([1.25 + u for u in result["u"]])

    with open(csv_out, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert rows[0][:4] == ["u", "x", "y", "z"]
    assert len(rows) == 6  # header + 5 samples
    assert float(rows[1][4]) == pytest.approx(1.25)


def test_resample_grid_exact_on_linear_field(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    csv_out = tmp_path / "grid.csv"
    result = resample_grid(msh, [0.05, 0.05, 0.05],
                           [0.35, 0.05, 0.05], [0.05, 0.35, 0.05],
                           [0.05, 0.05, 0.35], 2, 2, 2, view="T",
                           csv_out=csv_out)
    assert result["ok"] is True, result.get("error")
    assert len(result["points"]) == 8
    for pt, vals in zip(result["points"], result["values"]):
        assert vals[0][0] == pytest.approx(_t(*pt))
    with open(csv_out, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert len(rows) == 9  # header + 8 samples


# ----------------------------------------------------------------------
# pure-Python verbs: CSV, histogram, point history
# ----------------------------------------------------------------------

def test_export_view_csv_nodes(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    out = tmp_path / "nodes.csv"
    result = export_view_csv(msh, out)
    assert result["ok"] is True, result.get("error")
    assert result["kind"] == "nodes"
    assert result["n_rows"] == 5
    with open(out, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    assert rows[0] == ["tag", "x", "y", "z", "T_s0_c0", "T_s1_c0"]
    node5 = next(r for r in rows[1:] if r[0] == "5")
    assert float(node5[4]) == pytest.approx(7.0)
    assert float(node5[5]) == pytest.approx(7.5)


def test_export_view_csv_element_centroids(tmp_path):
    msh = _write(tmp_path, _EDATA, "edata.msh")
    out = tmp_path / "cells.csv"
    result = export_view_csv(msh, out, kind="elements")
    assert result["ok"] is True, result.get("error")
    assert result["n_rows"] == 2
    with open(out, newline="", encoding="utf-8") as f:
        rows = list(csv.reader(f))
    tet1 = next(r for r in rows[1:] if r[0] == "1")
    assert float(tet1[1]) == pytest.approx(0.25)  # centroid of unit tet
    assert float(tet1[4]) == pytest.approx(10.0)


def test_export_view_csv_reports_available_views(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    result = export_view_csv(msh, tmp_path / "x.csv", view="nope")
    assert result["ok"] is False
    assert "T" in result["error"]


def test_field_histogram_scalar_counts(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    result = field_histogram(msh, view="T", step=0, bins=2,
                             value_range=[1.0, 7.0])
    assert result["ok"] is True, result.get("error")
    # node values 1, 2, 3, 4, 7 -> bins [1,4): 3 samples, [4,7]: 2
    assert result["n_samples"] == 5
    assert result["counts"] == [3, 2]
    assert result["stats"]["min"] == pytest.approx(1.0)
    assert result["stats"]["max"] == pytest.approx(7.0)


def test_field_histogram_vector_magnitude(tmp_path):
    msh = _write(tmp_path, _VDIV, "div.msh")
    result = field_histogram(msh, bins=4)
    assert result["ok"] is True, result.get("error")
    assert result["n_samples"] == 5
    assert result["stats"]["max"] == pytest.approx(math.sqrt(3.0))
    assert result["stats"]["min"] == pytest.approx(0.0)


def test_point_history_returns_steps_and_times(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    result = point_history(msh, [0.2, 0.2, 0.2], view="T")
    assert result["ok"] is True, result.get("error")
    assert result["steps"] == [[pytest.approx(2.2)], [pytest.approx(2.7)]]
    assert result["times"] == pytest.approx([0.0, 0.5])


def test_point_history_outside_reports_distance(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    result = point_history(msh, [5.0, 5.0, 5.0], view="T")
    assert result["ok"] is False
    assert "outside" in result["error"]


@pytest.mark.skipif(not _MPL_AVAILABLE, reason="matplotlib not installed")
def test_plots_are_written(tmp_path):
    msh = _write(tmp_path, _SCALAR_2STEP)
    hist_png = tmp_path / "hist.png"
    hist = field_histogram(msh, view="T", plot_png=hist_png)
    assert hist["ok"] is True, hist.get("error")
    assert hist_png.is_file() and hist_png.stat().st_size > 3000

    history_png = tmp_path / "history.png"
    history = point_history(msh, [0.2, 0.2, 0.2], plot_png=history_png)
    assert history["ok"] is True, history.get("error")
    assert history_png.is_file() and history_png.stat().st_size > 3000


# ----------------------------------------------------------------------
# rendering additions: montage + camera orbit
# ----------------------------------------------------------------------

@pytest.mark.skipif(not _PIL_AVAILABLE, reason="Pillow not installed")
def test_render_montage_grid(tmp_path):
    from PIL import Image

    from radia_mcp.gmsh.render import render_montage

    a = tmp_path / "a.png"
    b = tmp_path / "b.png"
    Image.new("RGB", (60, 40), "red").save(a)
    Image.new("RGB", (50, 45), "blue").save(b)
    out = tmp_path / "grid.png"
    result = render_montage([a, b], out, labels=["case A", "case B"])
    assert result["ok"] is True, result.get("error")
    assert result["grid"] == [1, 2]
    with Image.open(out) as img:
        assert img.size == (120, 45)


def test_export_animation_orbit_mode(tmp_path):
    from radia_mcp.gmsh.render import export_animation

    msh = _write(tmp_path, _SCALAR_2STEP)
    gif = tmp_path / "orbit.gif"
    result = export_animation(msh, gif, orbit_axis="z", orbit_frames=3,
                              width=300, height=250)
    assert result["ok"] is True, result.get("error")
    assert result["num_steps"] == 3
    assert gif.is_file() and gif.stat().st_size > 1000


def test_export_animation_rejects_bad_orbit_axis(tmp_path):
    from radia_mcp.gmsh.render import export_animation

    msh = _write(tmp_path, _SCALAR_2STEP)
    result = export_animation(msh, tmp_path / "x.gif", orbit_axis="q")
    assert result["ok"] is False
    assert "orbit_axis" in result["error"]
