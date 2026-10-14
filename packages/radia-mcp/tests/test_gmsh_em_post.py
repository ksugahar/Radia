"""Analytic goldens + fail-fast contracts for the EM post verbs.

Physics fixtures: a Kuhn 6-tet cube [-0.5, 0.5]^3 carrying NodeData
fields that are UNIFORM or LINEAR, so P1 interpolation is exact and
every reduction below has a closed form:

* uniform B = 2 z^: flux through a tilted patch is |B| A cos(theta);
  the Maxwell stress on a box gives B^2 / (2 mu0) x A per face and a
  net force of exactly zero (a real cancellation of 1e6 N terms, not a
  zero integrand).
* B = x x^ (div B = a, curl B = 0): the stress integral over a box is
  the closed form int div T dV = a^2 c_x V / mu0 -- a NONZERO force
  golden, obtained by putting a source inside the box.
* H = (-y, x, 0): the circulation on a radius-r circle is 2 pi r^2,
  and on an inscribed n-gon it is exactly n r^2 sin(2 pi / n).
* S = x + 2y on a radius-r circle is r cos(theta) + 2 r sin(theta), so
  a_1 = r, b_1 = 2 r and every other bin is zero.
"""

import importlib.util
import math

import pytest
from radia_mcp.gmsh import em_post
from radia_mcp.gmsh.em_post import (MU0, flux_integral, gap_harmonics,
                                    harmonic_series, line_integral,
                                    maxwell_force)

_GMSH_AVAILABLE = importlib.util.find_spec("gmsh") is not None

# Kuhn decomposition of the cube [-0.5, 0.5]^3 (see the particle-trace
# tests): 6 positively oriented tets sharing the diagonal node1->node8.
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

B_UNIFORM = 2.0            # Tesla, along +z


def _nodedata(name, ncomp, rows, time=0.0, step=0):
    lines = ["$NodeData", "1", f'"{name}"', "1", str(time), "3",
             str(step), str(ncomp), str(len(rows))]
    for tag, vals in rows:
        lines.append(str(tag) + " " + " ".join(f"{v:.16g}" for v in vals))
    lines.append("$EndNodeData")
    return "\n".join(lines) + "\n"


def _write(tmp_path, text, name):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


def _uniform_b(tmp_path):
    return _write(tmp_path, _CUBE + _nodedata(
        "B", 3, [(i, [0.0, 0.0, B_UNIFORM]) for i in _P]), "uniform_b.msh")


def _swirl_h(tmp_path):
    """H = (-y, x, 0): circulation 2 pi r^2 on any circle about z."""
    return _write(tmp_path, _CUBE + _nodedata(
        "H", 3, [(i, [-p[1], p[0], 0.0]) for i, p in _P.items()]),
        "swirl_h.msh")


def _linear_bx(tmp_path):
    """B = (x, 0, 0): div B = 1, curl B = 0 -- a source inside the box."""
    return _write(tmp_path, _CUBE + _nodedata(
        "B", 3, [(i, [p[0], 0.0, 0.0]) for i, p in _P.items()]),
        "linear_bx.msh")


def _scalar_s(tmp_path):
    """S = x + 2y: a_1 = r, b_1 = 2 r on the circle of radius r."""
    return _write(tmp_path, _CUBE + _nodedata(
        "S", 1, [(i, [p[0] + 2.0 * p[1]]) for i, p in _P.items()]),
        "scalar_s.msh")


# ----------------------------------------------------------------------
# Fail-fast contracts (no gmsh subprocess involved)
# ----------------------------------------------------------------------

@pytest.fixture
def input_file(tmp_path):
    path = tmp_path / "field.msh"
    path.write_text("placeholder", encoding="ascii")
    return path


@pytest.fixture
def forbid_subprocess(monkeypatch):
    def unexpected(*_args, **_kwargs):
        raise AssertionError("invalid input reached the Gmsh subprocess")

    monkeypatch.setattr(em_post, "_run_em", unexpected)


@pytest.mark.parametrize(
    ("surface", "message"),
    [
        ("rect", "must be a dict naming exactly one of"),
        ({}, "must name exactly ONE of"),
        ({"disk": {}}, "unknown surface kind(s) disk"),
        ({"rect": {"center": [0, 0, 0], "u_vec": [1, 0, 0],
                   "v_vec": [0, 1, 0]},
          "circle": {"center": [0, 0, 0], "normal": [0, 0, 1],
                     "radius": 1.0}}, "must name exactly ONE of"),
        ({"rect": [1, 2, 3]}, "surface.rect must be a dict"),
        ({"rect": {"center": [0, 0, 0]}}, "surface.rect needs keys"),
        ({"rect": {"center": [0, 0, 0], "u_vec": [0, 0, 0],
                   "v_vec": [0, 1, 0]}}, "u_vec must be nonzero"),
        ({"rect": {"center": [0, 0, 0], "u_vec": [1, 0, 0],
                   "v_vec": [2, 0, 0]}}, "must not be parallel"),
        ({"circle": {"center": [0, 0, 0], "normal": [0, 0, 1],
                     "radius": 0.0}}, "radius must be positive"),
        ({"circle": {"center": [0, 0, 0], "normal": [0, 0, 0],
                     "radius": 1.0}}, "normal must be nonzero"),
        ({"circle": {"center": [0, 0], "normal": [0, 0, 1],
                     "radius": 1.0}}, "center needs exactly 3"),
        ({"circle": {"center": [0, 0, math.nan], "normal": [0, 0, 1],
                     "radius": 1.0}}, "center must be a finite number"),
    ],
)
def test_flux_rejects_bad_surface(input_file, forbid_subprocess, surface,
                                  message):
    result = flux_integral(input_file, surface)
    assert result["ok"] is False
    assert message in result["error"]


def test_flux_rejects_bad_grid(input_file, forbid_subprocess):
    good = {"circle": {"center": [0, 0, 0], "normal": [0, 0, 1],
                       "radius": 1.0}}
    r = flux_integral(input_file, good, n_grid=0)
    assert r["ok"] is False and "n_grid must be >= 1" in r["error"]
    r = flux_integral(input_file, good, time_step=-1)
    assert r["ok"] is False and "time_step must be >= 0" in r["error"]


def test_flux_reports_missing_file(tmp_path, forbid_subprocess):
    r = flux_integral(tmp_path / "nope.msh",
                      {"circle": {"center": [0, 0, 0], "normal": [0, 0, 1],
                                  "radius": 1.0}})
    assert r["ok"] is False and "file not found" in r["error"]


@pytest.mark.parametrize(
    ("spec", "message"),
    [
        ({"loop": {}}, "unknown path_spec kind(s) loop"),
        ({}, "must name exactly ONE of"),
        ({"polyline": {"points": [[0, 0, 0]]}}, "needs >= 2 points"),
        ({"polyline": {"points": [[0, 0, 0], [1, 0, 0]], "closed": True}},
         "closed path_spec.polyline needs >= 3"),
        ({"polyline": {"closed": True}}, "path_spec.polyline needs keys"),
        ({"circle": {"center": [0, 0, 0], "normal": [0, 0, 1],
                     "radius": -1.0}}, "radius must be positive"),
    ],
)
def test_line_rejects_bad_path(input_file, forbid_subprocess, spec,
                               message):
    result = line_integral(input_file, spec)
    assert result["ok"] is False
    assert message in result["error"]


def test_line_rejects_bad_controls(input_file, forbid_subprocess):
    good = {"circle": {"center": [0, 0, 0], "normal": [0, 0, 1],
                       "radius": 1.0}}
    r = line_integral(input_file, good, n=1)
    assert r["ok"] is False and "n must be >= 2" in r["error"]
    r = line_integral(input_file, good, expected_ni=0.0)
    assert r["ok"] is False and "expected_ni must be nonzero" in r["error"]


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({}, "box is required"),
        ({"box": [0, 0, 0]}, "box must be a dict"),
        ({"box": {"center": [0, 0, 0]}}, "box needs keys"),
        ({"box": {"center": [0, 0, 0], "half": [0.1, 0.0, 0.1]}},
         "half-extent must be positive"),
        ({"box": {"center": [0, 0, 0], "half": [0.1, 0.1, 0.1]},
          "n_grid": 0}, "n_grid must be >= 1"),
        ({"box": {"center": [0, 0, 0], "half": [0.1, 0.1, 0.1]},
          "mu0": 0.0}, "mu0 must be positive"),
        ({"box": {"center": [0, 0, 0], "half": [0.1, 0.1, 0.1]},
          "torque_about": [0, 0]}, "torque_about needs exactly 3"),
    ],
)
def test_force_rejects_bad_box(input_file, forbid_subprocess, kwargs,
                               message):
    result = maxwell_force(input_file, **kwargs)
    assert result["ok"] is False
    assert message in result["error"]


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"component": "tangent"}, "component must be one of"),
        ({"radius": 0.0}, "radius must be positive"),
        ({"axis": [0.0, 0.0, 0.0]}, "axis must be nonzero"),
        ({"n_samples": 2}, "n_samples must be >= 3"),
        ({"max_harmonic": 0}, "max_harmonic must be >= 1"),
        ({"center": [0.0, 0.0]}, "center needs exactly 3"),
    ],
)
def test_gap_rejects_bad_controls(input_file, forbid_subprocess, kwargs,
                                  message):
    call = {"center": [0.0, 0.0, 0.0], "axis": [0.0, 0.0, 1.0],
            "radius": 0.25}
    call.update(kwargs)
    result = gap_harmonics(input_file, call.pop("center"),
                           call.pop("axis"), call.pop("radius"), **call)
    assert result["ok"] is False
    assert message in result["error"]


# ----------------------------------------------------------------------
# The harmonic reducer is pure -- no gmsh needed
# ----------------------------------------------------------------------

def test_harmonic_series_recovers_planted_orders():
    # cos(2t) + 0.1 cos(10t) on 64 endpoint-excluded samples.
    thetas = [2.0 * math.pi * k / 64 for k in range(64)]
    samples = [math.cos(2.0 * t) + 0.1 * math.cos(10.0 * t)
               for t in thetas]
    out = harmonic_series(samples)
    rows = {r["n"]: r for r in out["harmonics"]}
    assert rows[2]["a_n"] == pytest.approx(1.0, abs=1e-12)
    assert rows[10]["a_n"] == pytest.approx(0.1, abs=1e-12)
    assert rows[2]["b_n"] == pytest.approx(0.0, abs=1e-12)
    # MEASURED leakage 1.666e-16 across every other bin
    leak = max(r["amplitude"] for r in out["harmonics"]
               if r["n"] not in (2, 10))
    assert leak < 1e-12
    # strictly below Nyquist: 64 samples resolve n <= 31
    assert out["n_bins"] == 31
    # no n = 1 content -> THD is undefined, not a roundoff ratio
    assert out["thd"] is None


def test_harmonic_series_phase_and_thd():
    # 2 cos(t - 30 deg) + 0.5 sin(3t): amplitude/phase and THD = 0.25
    thetas = [2.0 * math.pi * k / 128 for k in range(128)]
    samples = [2.0 * math.cos(t - math.radians(30.0))
               + 0.5 * math.sin(3.0 * t) for t in thetas]
    out = harmonic_series(samples)
    rows = {r["n"]: r for r in out["harmonics"]}
    assert rows[1]["amplitude"] == pytest.approx(2.0, rel=1e-12)
    assert rows[1]["phase_deg"] == pytest.approx(30.0, abs=1e-9)
    assert rows[3]["b_n"] == pytest.approx(0.5, abs=1e-12)
    assert out["fundamental"] == pytest.approx(2.0, rel=1e-12)
    assert out["thd"] == pytest.approx(0.25, rel=1e-12)


def test_harmonic_series_max_harmonic_trims_report_only():
    thetas = [2.0 * math.pi * k / 64 for k in range(64)]
    samples = [math.cos(t) + 0.5 * math.cos(5.0 * t) for t in thetas]
    full = harmonic_series(samples)
    trimmed = harmonic_series(samples, max_harmonic=3)
    assert max(r["n"] for r in trimmed["harmonics"]) == 3
    # the 5th harmonic is gone from the report but still in the THD
    assert trimmed["thd"] == pytest.approx(full["thd"], rel=1e-15)
    assert trimmed["thd"] == pytest.approx(0.5, rel=1e-12)


def test_harmonic_series_needs_a_revolution():
    with pytest.raises(ValueError, match="at least 3 samples"):
        harmonic_series([1.0, 2.0])


# ----------------------------------------------------------------------
# Analytic physics (gmsh required)
# ----------------------------------------------------------------------

pytestmark_gmsh = pytest.mark.skipif(not _GMSH_AVAILABLE,
                                     reason="gmsh package not installed")


@pytestmark_gmsh
def test_flux_through_a_tilted_rectangle_is_b_a_cos_theta(tmp_path):
    msh = _uniform_b(tmp_path)
    c30 = math.cos(math.radians(30.0))
    s30 = math.sin(math.radians(30.0))
    result = flux_integral(msh, {"rect": {
        "center": [0.0, 0.0, 0.0],
        "u_vec": [0.4, 0.0, 0.0],
        "v_vec": [0.0, 0.4 * c30, 0.4 * s30]}})
    assert result["ok"] is True, result.get("error")
    # |B| A cos(30 deg) = 2 * 0.16 * cos(30 deg)
    want = B_UNIFORM * 0.16 * c30
    assert want == pytest.approx(0.27712812921102, rel=1e-14)
    # MEASURED rel 2.4e-14 (1024-term summation roundoff; the midpoint
    # rule itself is exact for a uniform field)
    assert result["flux"] == pytest.approx(want, rel=1e-12)
    assert result["area"] == pytest.approx(0.16, rel=1e-12)
    assert result["n_points"] == 1024
    assert result["n_outside"] == 0
    # normal = normalize(u x v): tilted 30 deg about x
    assert result["normal"] == pytest.approx([0.0, -s30, c30], abs=1e-15)


@pytestmark_gmsh
def test_flux_through_a_disc_uses_polar_area_weighting(tmp_path):
    msh = _uniform_b(tmp_path)
    result = flux_integral(msh, {"circle": {
        "center": [0.0, 0.0, 0.0], "normal": [0.0, 0.0, 1.0],
        "radius": 0.2}})
    assert result["ok"] is True, result.get("error")
    want = B_UNIFORM * math.pi * 0.04
    assert want == pytest.approx(0.2513274122871834, rel=1e-14)
    # dA = r dr dtheta makes the midpoint rule exact in r for a constant
    # integrand -- MEASURED rel 1.1e-15
    assert result["flux"] == pytest.approx(want, rel=1e-12)
    assert result["area"] == pytest.approx(math.pi * 0.04, rel=1e-12)


@pytestmark_gmsh
def test_flux_normal_flips_with_the_patch_orientation(tmp_path):
    msh = _uniform_b(tmp_path)
    common = {"center": [0.0, 0.0, 0.0]}
    fwd = flux_integral(msh, {"rect": dict(
        common, u_vec=[0.4, 0.0, 0.0], v_vec=[0.0, 0.4, 0.0])})
    rev = flux_integral(msh, {"rect": dict(
        common, u_vec=[0.0, 0.4, 0.0], v_vec=[0.4, 0.0, 0.0])})
    assert fwd["ok"] and rev["ok"]
    assert fwd["flux"] == pytest.approx(-rev["flux"], rel=1e-14)
    assert fwd["flux"] > 0.0


@pytestmark_gmsh
def test_flux_refuses_a_patch_that_pokes_out_of_the_mesh(tmp_path):
    msh = _uniform_b(tmp_path)
    result = flux_integral(msh, {"rect": {
        "center": [0.0, 0.0, 0.0], "u_vec": [4.0, 0.0, 0.0],
        "v_vec": [0.0, 4.0, 0.0]}}, n_grid=8)
    assert result["ok"] is False
    assert result["n_outside"] > 0
    assert result["n_outside"] < result["n_points"]
    assert "outside the field data" in result["error"]


@pytestmark_gmsh
def test_flux_rejects_a_scalar_view(tmp_path):
    msh = _scalar_s(tmp_path)
    result = flux_integral(msh, {"circle": {
        "center": [0.0, 0.0, 0.0], "normal": [0.0, 0.0, 1.0],
        "radius": 0.2}})
    assert result["ok"] is False
    assert "needs a VECTOR view" in result["error"]


@pytestmark_gmsh
def test_circle_circulation_is_the_analytic_ampere_turns(tmp_path):
    msh = _swirl_h(tmp_path)
    radius = 0.3
    want = 2.0 * math.pi * radius * radius
    assert want == pytest.approx(0.565486677646163, rel=1e-14)
    result = line_integral(msh, {"circle": {
        "center": [0.0, 0.0, 0.0], "normal": [0.0, 0.0, 1.0],
        "radius": radius}}, n=256, expected_ni=want)
    assert result["ok"] is True, result.get("error")
    # parameter-space trapezoid against the analytic tangent --
    # MEASURED rel 3.9e-16 at n = 512
    assert result["integral"] == pytest.approx(want, rel=1e-12)
    assert abs(result["rel_err_vs_expected"]) < 1e-12
    assert result["expected_ni"] == pytest.approx(want, rel=1e-15)
    assert result["n_samples"] == 257
    assert result["n_outside"] == 0


@pytestmark_gmsh
def test_circle_circulation_flips_with_the_normal(tmp_path):
    msh = _swirl_h(tmp_path)
    spec = {"center": [0.0, 0.0, 0.0], "radius": 0.3}
    up = line_integral(msh, {"circle": dict(spec, normal=[0.0, 0.0, 1.0])},
                       n=64)
    down = line_integral(msh,
                         {"circle": dict(spec, normal=[0.0, 0.0, -1.0])},
                         n=64)
    assert up["ok"] and down["ok"]
    # right-hand rule about the normal
    assert up["integral"] == pytest.approx(-down["integral"], rel=1e-12)
    assert up["integral"] > 0.0


@pytestmark_gmsh
def test_polyline_circulation_carries_the_polygon_chord_bias(tmp_path):
    msh = _swirl_h(tmp_path)
    radius, n_seg = 0.3, 256
    points = [[radius * math.cos(2.0 * math.pi * k / n_seg),
               radius * math.sin(2.0 * math.pi * k / n_seg), 0.0]
              for k in range(n_seg)]
    result = line_integral(msh, {"polyline": {"points": points,
                                              "closed": True}},
                           n=2 * n_seg)
    assert result["ok"] is True, result.get("error")
    # exact for the POLYGON: sum of 2 x triangle areas
    # = n r^2 sin(2 pi / n) -- MEASURED rel 4.7e-15
    want_polygon = n_seg * radius * radius * math.sin(2.0 * math.pi / n_seg)
    assert result["integral"] == pytest.approx(want_polygon, rel=1e-12)
    # ... which sits ~1e-4 BELOW the arc answer.  That gap is the
    # polygon-vs-arc GEOMETRY (1 - sinc(2 pi / n) = 1.004e-4 at n = 256),
    # not an integration error: MEASURED 1.0040e-4.
    want_arc = 2.0 * math.pi * radius * radius
    bias = (want_arc - result["integral"]) / want_arc
    assert 5e-5 < bias < 5e-4


@pytestmark_gmsh
def test_open_polyline_matches_its_closed_form(tmp_path):
    msh = _swirl_h(tmp_path)
    # a single straight chord from A to B: int H.dl = |A x B|_z
    a = [0.3, -0.2, 0.0]
    b = [-0.1, 0.35, 0.0]
    result = line_integral(msh, {"polyline": {"points": [a, b]}}, n=16)
    assert result["ok"] is True, result.get("error")
    want = a[0] * b[1] - a[1] * b[0]
    assert result["integral"] == pytest.approx(want, rel=1e-12)


@pytestmark_gmsh
def test_maxwell_force_cancels_a_uniform_field_face_by_face(tmp_path):
    msh = _uniform_b(tmp_path)
    half = 0.4
    result = maxwell_force(msh, box={"center": [0.0, 0.0, 0.0],
                                     "half": [half, half, half]},
                           torque_about=[0.0, 0.0, 0.0])
    assert result["ok"] is True, result.get("error")
    assert result["n_outside"] == 0
    assert result["n_points"] == 6 * 24 * 24
    # every face carries B^2 / (2 mu0) x A -- MEASURED 1.018592e+06 N,
    # so the zero total below is a real cancellation, not a zero
    # integrand
    area = (2.0 * half) ** 2
    face_want = B_UNIFORM ** 2 / (2.0 * MU0) * area
    assert face_want == pytest.approx(1.018592e6, rel=1e-6)
    assert len(result["per_face"]) == 6
    assert {f["face"] for f in result["per_face"]} == {
        "+x", "-x", "+y", "-y", "+z", "-z"}
    for face in result["per_face"]:
        mag = math.hypot(math.hypot(*face["force"][:2]), face["force"][2])
        assert mag == pytest.approx(face_want, rel=1e-9)
        assert face["area"] == pytest.approx(area, rel=1e-12)
    total = math.hypot(math.hypot(*result["force_n"][:2]),
                       result["force_n"][2])
    assert total < 1e-9 * face_want
    assert max(abs(t) for t in result["torque_nm"]) < 1e-9 * face_want * half


@pytestmark_gmsh
def test_maxwell_force_on_an_enclosed_source_is_nonzero_closed_form(
        tmp_path):
    # B = (a x, 0, 0) has div B = a and curl B = 0, so the stress
    # divergence is (div B) B / mu0 and the box integral closes as
    #   F_x = a^2 c_x V / mu0,   F_y = F_z = 0.
    msh = _linear_bx(tmp_path)
    cx, half = 0.1, 0.25
    result = maxwell_force(msh, box={"center": [cx, 0.0, 0.0],
                                     "half": [half, half, half]},
                           torque_about=[cx, 0.0, 0.0])
    assert result["ok"] is True, result.get("error")
    volume = (2.0 * half) ** 3
    want_fx = cx * volume / MU0          # a = 1 for B_x = x
    assert want_fx == pytest.approx(9947.183937828455, rel=1e-14)
    # MEASURED rel 9.1e-16: the +-x faces carry a constant integrand and
    # the transverse faces cancel pairwise, so midpoint is exact here
    assert result["force_n"][0] == pytest.approx(want_fx, rel=1e-9)
    assert abs(result["force_n"][1]) < 1e-9 * want_fx
    assert abs(result["force_n"][2]) < 1e-9 * want_fx
    assert max(abs(t) for t in result["torque_nm"]) < 1e-9 * want_fx * half
    # the pull comes from the +x face being deeper in the field
    faces = {f["face"]: f["force"][0] for f in result["per_face"]}
    assert faces["+x"] > 0.0 > faces["-x"]
    assert faces["+x"] + faces["-x"] == pytest.approx(want_fx, rel=1e-9)


@pytestmark_gmsh
def test_maxwell_force_scales_with_one_over_mu0(tmp_path):
    msh = _linear_bx(tmp_path)
    box = {"center": [0.1, 0.0, 0.0], "half": [0.25, 0.25, 0.25]}
    base = maxwell_force(msh, box=box, n_grid=4)
    doubled = maxwell_force(msh, box=box, n_grid=4, mu0=2.0 * MU0)
    assert base["ok"] and doubled["ok"]
    assert doubled["force_n"][0] == pytest.approx(
        0.5 * base["force_n"][0], rel=1e-12)
    assert doubled["torque_nm"] is None


@pytestmark_gmsh
def test_maxwell_force_refuses_a_box_outside_the_mesh(tmp_path):
    msh = _uniform_b(tmp_path)
    result = maxwell_force(msh, box={"center": [0.0, 0.0, 0.0],
                                     "half": [2.0, 2.0, 2.0]}, n_grid=4)
    assert result["ok"] is False
    assert result["n_outside"] > 0
    assert "outside the field data" in result["error"]


@pytestmark_gmsh
def test_gap_harmonics_recovers_a_planted_linear_scalar(tmp_path):
    # S = x + 2y on r = 0.25 is 0.25 cos(theta) + 0.5 sin(theta)
    msh = _scalar_s(tmp_path)
    radius = 0.25
    result = gap_harmonics(msh, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], radius,
                           n_samples=64)
    assert result["ok"] is True, result.get("error")
    assert result["component"] == "scalar"      # auto on a 1-comp view
    assert result["n_samples"] == 64
    assert result["n_bins"] == 31
    rows = {r["n"]: r for r in result["harmonics"]}
    assert rows[1]["a_n"] == pytest.approx(radius, abs=1e-12)
    assert rows[1]["b_n"] == pytest.approx(2.0 * radius, abs=1e-12)
    assert rows[1]["amplitude"] == pytest.approx(
        math.hypot(radius, 2.0 * radius), rel=1e-12)
    assert rows[1]["phase_deg"] == pytest.approx(
        math.degrees(math.atan2(2.0, 1.0)), abs=1e-9)
    assert abs(rows[0]["a_n"]) < 1e-12
    # MEASURED leakage 8.9e-17 on every bin above the fundamental
    leak = max(r["amplitude"] for r in result["harmonics"] if r["n"] >= 2)
    assert leak < 1e-12
    assert result["thd"] < 1e-12


@pytestmark_gmsh
def test_gap_harmonics_components_split_a_vector_view(tmp_path):
    # H = (-y, x, 0) on a circle about z is purely TANGENTIAL, |H| = r
    msh = _swirl_h(tmp_path)
    radius = 0.3
    tang = gap_harmonics(msh, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], radius,
                         n_samples=32, component="tangential")
    rad = gap_harmonics(msh, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], radius,
                        n_samples=32, component="radial")
    axial = gap_harmonics(msh, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], radius,
                          n_samples=32, component="axial")
    mag = gap_harmonics(msh, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], radius,
                        n_samples=32, component="magnitude")
    assert tang["ok"] and rad["ok"] and axial["ok"] and mag["ok"]
    # tangential and |H| are the constant r; radial and axial vanish
    assert tang["harmonics"][0]["a_n"] == pytest.approx(radius, rel=1e-12)
    assert mag["harmonics"][0]["a_n"] == pytest.approx(radius, rel=1e-12)
    assert abs(rad["harmonics"][0]["a_n"]) < 1e-12
    assert abs(axial["harmonics"][0]["a_n"]) < 1e-12
    # a constant signal has no fundamental -> THD stays undefined
    assert tang["thd"] is None
    assert max(r["amplitude"] for r in tang["harmonics"][1:]) < 1e-12


@pytestmark_gmsh
def test_gap_harmonics_auto_picks_radial_on_a_vector_view(tmp_path):
    msh = _uniform_b(tmp_path)     # B = 2 z^ -> radial component is 0
    result = gap_harmonics(msh, [0.0, 0.0, 0.0], [0.0, 0.0, 1.0], 0.25,
                           n_samples=32)
    assert result["ok"] is True, result.get("error")
    assert result["component"] == "radial"
    assert abs(result["harmonics"][0]["a_n"]) < 1e-12


@pytestmark_gmsh
def test_gap_harmonics_rejects_a_mismatched_component(tmp_path):
    scalar = gap_harmonics(_scalar_s(tmp_path), [0.0, 0.0, 0.0],
                           [0.0, 0.0, 1.0], 0.25, component="radial")
    assert scalar["ok"] is False
    assert "needs a VECTOR view" in scalar["error"]
    vector = gap_harmonics(_uniform_b(tmp_path), [0.0, 0.0, 0.0],
                           [0.0, 0.0, 1.0], 0.25, component="scalar")
    assert vector["ok"] is False
    assert "needs a 1-component view" in vector["error"]


@pytestmark_gmsh
def test_gap_harmonics_refuses_a_circle_outside_the_mesh(tmp_path):
    result = gap_harmonics(_scalar_s(tmp_path), [0.0, 0.0, 0.0],
                           [0.0, 0.0, 1.0], 5.0, n_samples=16)
    assert result["ok"] is False
    assert result["n_outside"] > 0
    assert "outside the field data" in result["error"]


@pytestmark_gmsh
def test_view_selector_names_the_available_views(tmp_path):
    msh = _uniform_b(tmp_path)
    result = flux_integral(msh, {"circle": {
        "center": [0.0, 0.0, 0.0], "normal": [0.0, 0.0, 1.0],
        "radius": 0.2}}, view="Bfield")
    assert result["ok"] is False
    assert "view 'Bfield' not found" in result["error"]
    assert "['B']" in result["error"]
