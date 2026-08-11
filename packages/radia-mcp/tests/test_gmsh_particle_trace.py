"""Analytic validation + fail-fast contracts for gmsh_particle_trace.

Physics fixtures: a Kuhn 6-tet cube [-0.5, 0.5]^3 carrying UNIFORM
NodeData fields (P1 interpolation is exact), so every check has a
closed-form answer:

* uniform B: the orbit is a circle of gyroradius r = m*u_perp/(|q|*B)
  (u = gamma*v).  The Boris polygon inflates the radius by
  (omega*dt)^2/8 ~ 1.2e-3 at 64 steps/gyration -- the tolerances below
  include that known discretization bias.
* uniform E x B: the guiding-center drift is exactly E x B / B^2 (a
  property the Boris pusher preserves independent of dt).
* pure B conserves |v| EXACTLY (the rotation is an isometry), so
  speed_change_rel is a roundoff-level integrator health metric.
"""

import importlib.util
import math

import pytest
from radia_mcp.gmsh import post_process
from radia_mcp.gmsh.post_process import particle_trace, poincare

_GMSH_AVAILABLE = importlib.util.find_spec("gmsh") is not None

C = 299792458.0
QE = 1.602176634e-19
ME = 9.1093837015e-31

# Kuhn decomposition of the cube [-0.5, 0.5]^3 (see the streamlines
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


def _nodedata(name, ncomp, rows, time=0.0, step=0):
    lines = ["$NodeData", "1", f'"{name}"', "1", str(time), "3",
             str(step), str(ncomp), str(len(rows))]
    for tag, vals in rows:
        lines.append(str(tag) + " " + " ".join(f"{v:.16g}" for v in vals))
    lines.append("$EndNodeData")
    return "\n".join(lines) + "\n"


def _electron(ke_ev, r_target):
    """Analytic launch state: gamma, u=gamma*v, v, and the B that makes
    the gyroradius exactly r_target."""
    gamma = 1.0 + ke_ev * QE / (ME * C * C)
    u = C * math.sqrt(gamma * gamma - 1.0)
    v = u / gamma
    b = ME * u / (QE * r_target)
    return gamma, u, v, b


KE_EV = 1.0e4          # 10 keV electron
R_GYRO = 0.2           # meters -- fits the cube with margin
GAMMA0, U0, V0, B0 = _electron(KE_EV, R_GYRO)
T_GYRO = 2.0 * math.pi * GAMMA0 * ME / (QE * B0)


def _uniform_b(bz=B0):
    return _CUBE + _nodedata("B", 3, [(i, [0.0, 0.0, bz]) for i in _P])


def _uniform_b_and_e(bz, ex):
    return (_CUBE
            + _nodedata("B", 3, [(i, [0.0, 0.0, bz]) for i in _P])
            + _nodedata("E", 3, [(i, [ex, 0.0, 0.0]) for i in _P]))


def _write(tmp_path, text, name="case.msh"):
    p = tmp_path / name
    p.write_text(text, encoding="utf-8")
    return p


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

    monkeypatch.setattr(post_process, "_run_post", unexpected)


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"species": "muon"}, "unknown species"),
        ({"charge_e": -1.0}, "BOTH charge_e and mass_amu"),
        ({"mass_amu": 1.0}, "BOTH charge_e and mass_amu"),
        ({"charge_e": 0.0, "mass_amu": 1.0}, "charge_e must be nonzero"),
        ({"charge_e": 1.0, "mass_amu": -1.0}, "mass_amu must be positive"),
        ({"dt_s": 0.0}, "dt_s must be positive"),
        ({"max_time_s": -1.0}, "max_time_s must be positive"),
        ({"steps_per_gyration": 3}, "steps_per_gyration must be >= 4"),
        ({"max_steps": 0}, "max_steps must be positive"),
        ({"color_by": "phase"}, "color_by must be one of"),
        ({"arrows_every": -1}, "arrows_every must be non-negative"),
        ({"time_step": -1}, "time_step must be non-negative"),
        ({"animation_frames": -1}, "animation_frames must be non-negative"),
        ({"animation_mode": "warp"}, "animation_mode must be"),
        ({"comet_window": 0.0}, "comet_window must be in"),
        ({"comet_window": 1.5}, "comet_window must be in"),
    ],
)
def test_rejects_invalid_controls(input_file, forbid_subprocess, kwargs,
                                  message):
    result = particle_trace(input_file, [[0.0, 0.0, 0.0]],
                            [1.0, 0.0, 0.0], 1.0e3, **kwargs)
    assert result["ok"] is False
    assert message in result["error"]


def test_rejects_bad_launch_state(input_file, forbid_subprocess):
    r = particle_trace(input_file, [], [1.0, 0.0, 0.0], 1.0e3)
    assert r["ok"] is False and "at least one point" in r["error"]
    r = particle_trace(input_file, [[0.0, 0.0, 0.0]], [0.0, 0.0, 0.0],
                       1.0e3)
    assert r["ok"] is False and "direction must be nonzero" in r["error"]
    r = particle_trace(input_file, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                       0.0)
    assert r["ok"] is False and "kinetic_energy_ev" in r["error"]
    r = particle_trace(input_file, [[0.0, math.nan, 0.0]],
                       [1.0, 0.0, 0.0], 1.0e3)
    assert r["ok"] is False and "seeds" in r["error"]


# ----------------------------------------------------------------------
# Analytic physics (gmsh required)
# ----------------------------------------------------------------------

pytestmark_gmsh = pytest.mark.skipif(not _GMSH_AVAILABLE,
                                     reason="gmsh package not installed")


@pytestmark_gmsh
def test_uniform_b_gyration_is_the_analytic_circle(tmp_path):
    msh = _write(tmp_path, _uniform_b(), "uniform_b.msh")
    out = tmp_path / "orbit.pos"
    result = particle_trace(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                            KE_EV, species="electron", max_steps=64,
                            return_points=True, out_file=out)
    assert result["ok"] is True, result.get("error")
    assert result["n_tracks"] == 1
    assert result["reasons"] == {"max_steps": 1}
    assert out.is_file()

    track = result["tracks"][0]
    # the reported seed gyroradius IS the analytic construction value
    assert track["gyroradius_seed_m"] == pytest.approx(R_GYRO, rel=1e-12)
    assert track["dt_s"] == pytest.approx(T_GYRO / 64.0, rel=1e-12)
    # pure B: the Boris rotation conserves speed exactly
    assert track["speed_change_rel"] < 1e-12
    assert track["ke_end_ev"] == pytest.approx(KE_EV, rel=1e-12)

    pts = result["polylines"][0]["points"]
    assert len(pts) == 65
    # electron (q < 0), v = +x, B = +z: F = qv x B = +y, center at +y
    xs = [p[0] for p in pts[:64]]
    ys = [p[1] for p in pts[:64]]
    cx, cy = sum(xs) / 64.0, sum(ys) / 64.0
    assert cy > 0.0
    radii = [math.hypot(x - cx, y - cy) for x, y in zip(xs, ys)]
    r_mean = sum(radii) / len(radii)
    # circularity: every sample equidistant from the fitted center
    assert (max(radii) - min(radii)) / r_mean < 2e-3
    # radius = analytic gyroradius (+ the known (w*dt)^2/8 polygon bias)
    assert r_mean == pytest.approx(R_GYRO, rel=5e-3)
    # B || z and no E: the orbit stays exactly planar
    assert max(abs(p[2]) for p in pts) < 1e-12


@pytestmark_gmsh
def test_positron_curls_the_opposite_way(tmp_path):
    msh = _write(tmp_path, _uniform_b(), "uniform_b.msh")
    kwargs = dict(kinetic_energy_ev=KE_EV, max_steps=10,
                  return_points=True)
    r_el = particle_trace(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                          species="electron",
                          out_file=tmp_path / "el.pos", **kwargs)
    r_po = particle_trace(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                          species="positron",
                          out_file=tmp_path / "po.pos", **kwargs)
    assert r_el["ok"] and r_po["ok"]
    y_el = r_el["polylines"][0]["points"][10][1]
    y_po = r_po["polylines"][0]["points"][10][1]
    assert y_el > 0.0 > y_po
    # mirror symmetry: same |y| deflection
    assert y_el == pytest.approx(-y_po, rel=1e-12)


@pytestmark_gmsh
def test_exb_drift_matches_the_analytic_guiding_center(tmp_path):
    ex = 0.01 * V0 * B0            # drift = 1% of the launch speed
    msh = _write(tmp_path, _uniform_b_and_e(B0, ex), "exb.msh")
    result = particle_trace(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                            KE_EV, species="electron", view="B",
                            e_view="E", max_steps=320,
                            return_points=True,
                            out_file=tmp_path / "exb.pos")
    assert result["ok"] is True, result.get("error")
    assert result["tracks"][0]["reason"] == "max_steps"
    pts = result["polylines"][0]["points"]
    dt = result["tracks"][0]["dt_s"]

    def _center(k):
        window = pts[64 * k:64 * (k + 1)]
        return [sum(p[i] for p in window) / len(window) for i in range(3)]

    c0, c4 = _center(0), _center(4)
    drift_y = (c4[1] - c0[1]) / (256.0 * dt)
    v_drift = -ex / B0             # E x B / B^2 = (Ex x^ x Bz z^)/B^2
    assert drift_y == pytest.approx(v_drift, rel=2e-2)
    assert abs(c4[0] - c0[0]) < 0.05 * abs(c4[1] - c0[1])
    # E does work along the cycloid: kinetic energy actually changed
    assert result["tracks"][0]["ke_end_ev"] != pytest.approx(
        KE_EV, rel=1e-6)


@pytestmark_gmsh
def test_orbit_leaving_the_data_terminates_left_data(tmp_path):
    msh = _write(tmp_path, _uniform_b(), "uniform_b.msh")
    # center at (0.4, 0.2): the circle spans x in [0.2, 0.6] and exits
    result = particle_trace(msh, [[0.4, 0.0, 0.0]], [1.0, 0.0, 0.0],
                            KE_EV, max_steps=200,
                            out_file=tmp_path / "exit.pos")
    assert result["ok"] is True, result.get("error")
    track = result["tracks"][0]
    assert track["reason"] == "left_data"
    assert 0 < track["n_steps"] < 64


@pytestmark_gmsh
def test_seed_outside_data_is_skipped_not_fatal(tmp_path):
    msh = _write(tmp_path, _uniform_b(), "uniform_b.msh")
    result = particle_trace(msh, [[2.0, 0.0, 0.0], [0.0, 0.0, 0.0]],
                            [1.0, 0.0, 0.0], KE_EV, max_steps=16,
                            arrows_every=4,
                            out_file=tmp_path / "skip.pos")
    assert result["ok"] is True, result.get("error")
    assert result["skipped_seeds"] == 1
    assert result["n_tracks"] == 1
    assert result["reasons"]["seed_outside_data"] == 1
    assert result["n_arrows"] > 0


@pytestmark_gmsh
def test_zero_b_seed_without_dt_fails_loud(tmp_path):
    # radial field (x, y, z) vanishes at the origin: no gyration period
    rad = _CUBE + _nodedata(
        "B", 3, [(i, list(map(float, p))) for i, p in _P.items()])
    msh = _write(tmp_path, rad, "radial.msh")
    result = particle_trace(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                            KE_EV, max_steps=8,
                            out_file=tmp_path / "b0.pos")
    assert result["ok"] is False
    assert "dt_s" in result["error"]


def _beam_steps(pos_path, anim):
    """Per-step visible-element counts of the beam SL view.

    The threshold sits HALFWAY between the colour-range maximum and the
    sentinel: writing a .pos round-trips the value through ASCII and
    drops the last digits (measured: sentinel 110009.99999999942 comes
    back as 110009.9999999994), so an equality test on the sentinel
    would count every hidden element as visible.
    """
    import gmsh

    hi = anim["color_range"][1]
    threshold = hi + 0.5 * (anim["sentinel"] - hi)

    gmsh.initialize(["-noconfig"])
    try:
        gmsh.option.setNumber("General.Terminal", 0)
        gmsh.open(str(pos_path))
        tag = None
        for t in gmsh.view.getTags():
            i = gmsh.view.getIndex(t)
            if gmsh.option.getString(f"View[{i}].Name").startswith("beam"):
                tag = t
                n_steps = int(gmsh.option.getNumber(f"View[{i}].NbTimeStep"))
        assert tag is not None, "no beam view in the file"
        _dtypes, nels, data = gmsh.view.getListData(tag)
        n_el = int(nels[0])
        block = data[0]
        width = 6 + 2 * n_steps                     # SL: xx yy zz + 2/step
        counts = []
        for s in range(n_steps):
            counts.append(sum(
                1 for e in range(n_el)
                if float(block[e * width + 6 + 2 * s]) < threshold))
        return counts, n_el
    finally:
        gmsh.finalize()


@pytestmark_gmsh
def test_beam_animation_trail_grows_and_hides_with_a_sentinel(tmp_path):
    msh = _write(tmp_path, _uniform_b(), "uniform_b.msh")
    out = tmp_path / "beam.pos"
    result = particle_trace(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                            KE_EV, max_steps=48, animation_frames=8,
                            out_file=out)
    assert result["ok"] is True, result.get("error")
    anim = result["animation"]
    assert anim["n_steps"] == 8 and anim["mode"] == "trail"
    lo, hi = anim["color_range"]
    assert anim["sentinel"] > hi and math.isfinite(anim["sentinel"])
    assert anim["render_hint"]["color"] == {"range": [lo, hi],
                                            "saturate": False}
    counts, n_el = _beam_steps(out, anim)
    assert len(counts) == 8
    assert counts == sorted(counts)          # monotone growth
    assert counts[-1] == n_el                # ends fully drawn
    assert counts[0] < counts[-1]


@pytestmark_gmsh
def test_beam_animation_comet_window_stays_bounded(tmp_path):
    msh = _write(tmp_path, _uniform_b(), "uniform_b.msh")
    out = tmp_path / "comet.pos"
    result = particle_trace(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                            KE_EV, max_steps=64, animation_frames=10,
                            animation_mode="comet", comet_window=0.25,
                            out_file=out)
    assert result["ok"] is True, result.get("error")
    counts, n_el = _beam_steps(out, result["animation"])
    assert max(counts) < n_el                # never the whole track
    assert max(counts) <= 0.4 * n_el
    assert min(counts[1:]) > 0


@pytestmark_gmsh
def test_beam_animation_sentinel_survives_a_constant_colour_quantity(
        tmp_path):
    """A monoenergetic beam coloured by energy has ~1e-13 spread; the
    sentinel must still land far outside, or nothing would be hidden."""
    msh = _write(tmp_path, _uniform_b(), "uniform_b.msh")
    out = tmp_path / "flat.pos"
    result = particle_trace(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                            KE_EV, max_steps=32, color_by="energy",
                            animation_frames=4, out_file=out)
    assert result["ok"] is True, result.get("error")
    anim = result["animation"]
    lo, hi = anim["color_range"]
    assert hi > lo                              # widened, not degenerate
    assert anim["sentinel"] > hi + 0.5 * (hi - lo)
    counts, n_el = _beam_steps(out, anim)
    assert counts[0] < n_el and counts[-1] == n_el


@pytestmark_gmsh
def test_max_time_terminates_at_the_requested_flight_time(tmp_path):
    msh = _write(tmp_path, _uniform_b(), "uniform_b.msh")
    dt = T_GYRO / 64.0
    result = particle_trace(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                            KE_EV, max_steps=500, max_time_s=10.0 * dt,
                            color_by="speed",
                            out_file=tmp_path / "tmax.pos")
    assert result["ok"] is True, result.get("error")
    track = result["tracks"][0]
    assert track["reason"] == "max_time"
    assert track["n_steps"] == 10
    assert track["time_s"] == pytest.approx(10.0 * dt, rel=1e-12)


# ----------------------------------------------------------------------
# Regressions for the audited defects (each FAILS against the old code)
# ----------------------------------------------------------------------

@pytestmark_gmsh
def test_time_step_past_the_last_step_is_a_clean_error_not_a_crash(
        tmp_path):
    """E1: gmsh.view.probe(step=N >= NbTimeStep) FAULTS the child.

    Old behaviour (MEASURED gmsh 4.15.2): "OSError: exception: access
    violation reading 0xFFFFFFFFFFFFFFFF" -- a native out-of-bounds read
    with no Python traceback.  The step must be validated against
    NbTimeStep before the first probe.
    """
    msh = _write(tmp_path, _uniform_b(), "uniform_b.msh")
    result = particle_trace(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                            KE_EV, time_step=3, max_steps=4,
                            out_file=tmp_path / "oob.pos")
    assert result["ok"] is False
    assert "access violation" not in result["error"]
    assert "time_step 3 out of range" in result["error"]
    assert "1 step(s) (0..0)" in result["error"]


@pytestmark_gmsh
def test_scalar_b_view_is_rejected_not_reported_as_seed_outside_data(
        tmp_path):
    """LENSA-02: a scalar view used to yield ok=true / 0 tracks.

    The probe returns fewer than 3 components at every sample, which the
    tracer read as "outside the data" -- a diagnosis pointing at the
    geometry instead of at the view the caller picked (MEASURED:
    n_tracks 0, reasons {"seed_outside_data": 1}).
    """
    scalar = _CUBE + _nodedata("S", 1, [(i, [1.0]) for i in _P])
    msh = _write(tmp_path, scalar, "scalar.msh")
    result = particle_trace(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                            KE_EV, dt_s=1e-11, max_steps=8,
                            out_file=tmp_path / "scalar.pos")
    assert result["ok"] is False
    assert "view 'S' has 1 component(s)" in result["error"]
    assert "3-component vector view" in result["error"]


@pytestmark_gmsh
def test_scalar_e_view_is_rejected_too(tmp_path):
    mixed = (_CUBE
             + _nodedata("B", 3, [(i, [0.0, 0.0, B0]) for i in _P])
             + _nodedata("phi", 1, [(i, [1.0]) for i in _P]))
    msh = _write(tmp_path, mixed, "mixed.msh")
    result = particle_trace(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                            KE_EV, view="B", e_view="phi", max_steps=8,
                            out_file=tmp_path / "mixed.pos")
    assert result["ok"] is False
    assert "view 'phi' has 1 component(s)" in result["error"]
    assert "the E field must be" in result["error"]


@pytestmark_gmsh
def test_beam_sentinel_scales_with_the_data_not_an_absolute_one(tmp_path):
    """LENSA-01: the sentinel floor must be RELATIVE to the values.

    A time-coloured trace carries values in SECONDS.  The old floor
    max(span, |lo|, |hi|, 1.0) took the 1.0 branch and put the sentinel
    at exactly 10.0 -- MEASURED lo=1.67949e-10, hi=6.71778e-06, so
    1.5e6x the data, which poisons every stats verb reading the view.
    The relative floor gives 7.38939e-05, i.e. 11x the data.
    """
    msh = _write(tmp_path, _uniform_b(), "uniform_b.msh")
    out = tmp_path / "long.pos"
    result = particle_trace(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                            KE_EV, max_steps=20000, animation_frames=8,
                            out_file=out)
    assert result["ok"] is True, result.get("error")
    anim = result["animation"]
    lo, hi = anim["color_range"]
    assert lo == pytest.approx(1.67949e-10, rel=1e-4)
    assert hi == pytest.approx(6.71778e-06, rel=1e-4)
    # the old code produced exactly 10.0 here
    assert anim["sentinel"] != pytest.approx(10.0, rel=1e-9)
    assert anim["sentinel"] == pytest.approx(7.38939e-05, rel=1e-4)
    # relative floor => a fixed 11x the range maximum, for any unit
    assert anim["sentinel"] / hi == pytest.approx(11.0, rel=1e-3)
    # ... and hiding still works: it is far outside the colour range
    counts, n_el = _beam_steps(out, anim)
    assert counts == sorted(counts)
    assert counts[0] < counts[-1] == n_el


@pytestmark_gmsh
def test_flat_colour_sentinel_keeps_its_measured_absolute_value(tmp_path):
    """The monoenergetic case pins the OTHER end of the same formula."""
    msh = _write(tmp_path, _uniform_b(), "uniform_b.msh")
    result = particle_trace(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                            KE_EV, max_steps=32, color_by="energy",
                            animation_frames=4,
                            out_file=tmp_path / "flat2.pos")
    assert result["ok"] is True, result.get("error")
    lo, hi = result["animation"]["color_range"]
    # scale = 1e-3 * 1e4 = 10 -> hi = lo + 0.01, sentinel = hi + 100
    assert hi - lo == pytest.approx(0.01, rel=1e-6)
    assert result["animation"]["sentinel"] == pytest.approx(1.0100e4,
                                                            rel=1e-6)


@pytest.fixture
def forbid_gmsh_launch(monkeypatch):
    """Patch the ACTUAL subprocess boundary, not the _run_post wrapper.

    The list-data out_file gate lives inside _run_post (it is keyed on
    the op name), so the coarser forbid_subprocess fixture would mask
    the very thing under test.
    """
    def unexpected(*_args, **_kwargs):
        raise AssertionError("invalid input reached the Gmsh subprocess")

    monkeypatch.setattr(post_process, "run_gmsh_json_subprocess",
                        unexpected)


def test_list_data_out_file_must_be_pos(input_file, forbid_gmsh_launch):
    """LENSA-03: particle_trace(out_file="x.msh") silently wrote v2.2.

    MEASURED: the file began "$MeshFormat / 2.2 0 8" and this lane's own
    field_histogram then rejected it with "MSH v2.2 is not supported
    (ASCII v4.x only)" -- a silent breach of the v4.1-only policy.  The
    rejection happens before gmsh is ever launched.
    """
    result = particle_trace(input_file, [[0.0, 0.0, 0.0]],
                            [1.0, 0.0, 0.0], 1.0e3,
                            out_file=input_file.parent / "tracks.msh")
    assert result["ok"] is False
    assert "must end in '.pos'" in result["error"]
    assert "tracks.msh" in result["error"]
    assert "v2.2" in result["error"]


# ----------------------------------------------------------------------
# poincare: the section of the analytic gyration circle
# ----------------------------------------------------------------------

SPG_POINCARE = 256
DT_POINCARE = T_GYRO / SPG_POINCARE
# The Boris pusher advances the position with the velocity AFTER the
# magnetic rotation, so the sampled polygon runs half a step ahead of
# the analytic circle.  Along the motion that is dt*v/2 -- MEASURED
# 2.454369e-03 m for this launch state.  Every section coordinate below
# carries exactly that offset; it is a known, explainable discretization
# bias, NOT an error.
HALF_STEP_M = DT_POINCARE * V0 / 2.0


def _section(tmp_path, **kwargs):
    msh = _write(tmp_path, _uniform_b(), "uniform_b.msh")
    return poincare(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0], KE_EV,
                    [0.0, R_GYRO, 0.0], [0.0, 1.0, 0.0],
                    steps_per_gyration=SPG_POINCARE,
                    max_steps=int(2.2 * SPG_POINCARE),
                    out_file=tmp_path / "poinc.pos", **kwargs)


@pytestmark_gmsh
def test_poincare_section_of_the_gyration_circle_is_analytic(tmp_path):
    """10 keV electron, r_gyro = 0.2 m, section plane y = +0.2.

    The plane runs through the gyrocentre, so the circle pierces it
    twice per turn at x = +-r.  Tracking 2.2 gyroperiods from phase 0
    (the bottom of the circle) puts crossings at phases 0.25, 0.75,
    1.25, 1.75 -- EXACTLY 4, with the 5th at 2.25 falling outside.
    MEASURED u = +0.19754563 and -0.20245437, v = 0 exactly; both are
    the analytic +-0.2 shifted by the half-step bias (u + dt*v/2 =
    +-0.2 to 6e-8).
    """
    res = _section(tmp_path)
    assert res["ok"] is True, res.get("error")
    # frame convention: n = +y picks the least-aligned world axis (x,
    # lowest index on the x/z tie), so u = +x and v = n x u = -z
    assert res["u_axis"] == [1.0, 0.0, 0.0]
    assert res["v_axis"] == [0.0, 0.0, -1.0]
    assert res["plane_normal"] == [0.0, 1.0, 0.0]

    assert res["n_per_seed"] == [4]
    hits = res["crossings"][0]
    assert len(hits) == 4
    assert [h["direction"] for h in hits] == [
        "positive", "negative", "positive", "negative"]

    # the orbit is planar: the out-of-plane section coordinate is 0
    for h in hits:
        assert abs(h["v"]) < 1e-12

    # |u| = r_gyro, up to the half-step bias (symmetric: the whole
    # circle is shifted by -dt*v/2 in x, so the +x crossing moves IN
    # and the -x crossing moves OUT by the same amount)
    for h in hits:
        assert (R_GYRO - HALF_STEP_M - 1e-4
                <= abs(h["u"]) <= R_GYRO + HALF_STEP_M + 1e-4)
    # the sharper statement: undo the bias and the analytic radius is back
    assert hits[0]["u"] + HALF_STEP_M == pytest.approx(R_GYRO, abs=1e-7)
    assert hits[1]["u"] + HALF_STEP_M == pytest.approx(-R_GYRO, abs=1e-7)
    assert hits[0]["u"] == pytest.approx(0.19754563, abs=1e-7)
    assert hits[1]["u"] == pytest.approx(-0.20245437, abs=1e-7)

    # crossing times: the first pierce is a quarter turn in, minus the
    # same half step (MEASURED 5.33264e-09 s vs T/4 - dt/2 = 5.33237e-09)
    assert hits[0]["t_s"] == pytest.approx(
        T_GYRO / 4.0 - DT_POINCARE / 2.0, rel=1e-3)
    # ... and they are half a gyroperiod apart
    for a, b in zip(hits, hits[1:]):
        assert b["t_s"] - a["t_s"] == pytest.approx(T_GYRO / 2.0, rel=1e-3)
    # a half turn is exactly SPG/2 samples
    assert [h["index"] for h in hits] == [63, 191, 319, 447]


@pytestmark_gmsh
def test_poincare_orbit_closes_on_the_next_turn(tmp_path):
    """Closure: crossings one full gyration apart must coincide.

    MEASURED drift over one turn: |u0 - u2| = 1.49e-08 m and
    |u1 - u3| = 2.49e-08 m (the Boris polygon's per-turn radius drift on
    a 256-step gyration).  The tolerance below is that measurement with
    margin -- a tighter 1e-9 does NOT hold at this step count.
    """
    hits = _section(tmp_path)["crossings"][0]
    assert abs(hits[0]["u"] - hits[2]["u"]) < 1e-7
    assert abs(hits[1]["u"] - hits[3]["u"]) < 1e-7
    assert abs(hits[0]["v"] - hits[2]["v"]) < 1e-12


@pytestmark_gmsh
def test_poincare_crossing_direction_splits_the_section(tmp_path):
    both = _section(tmp_path)["crossings"][0]
    pos = _section(tmp_path, crossing_direction="positive")["crossings"][0]
    neg = _section(tmp_path, crossing_direction="negative")["crossings"][0]
    assert len(pos) == 2 and len(neg) == 2
    assert len(both) == len(pos) + len(neg)
    # +y-going pierces are the +x side, -y-going the -x side
    assert all(h["u"] > 0.0 for h in pos)
    assert all(h["u"] < 0.0 for h in neg)
    assert [h["u"] for h in pos] == [h["u"] for h in both[0::2]]


@pytestmark_gmsh
def test_poincare_plane_missing_the_orbit_yields_no_crossings(tmp_path):
    """The circle spans y in [0, 0.4]; a plane at y = 0.45 is missed."""
    msh = _write(tmp_path, _uniform_b(), "uniform_b.msh")
    res = poincare(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0], KE_EV,
                   [0.0, 0.45, 0.0], [0.0, 1.0, 0.0],
                   steps_per_gyration=64, max_steps=140,
                   out_file=tmp_path / "miss.pos")
    assert res["ok"] is True, res.get("error")
    assert res["n_per_seed"] == [0]
    assert res["n_crossings"] == 0


@pytestmark_gmsh
def test_poincare_writes_a_scatter_png(tmp_path):
    png = tmp_path / "section.png"
    res = _section(tmp_path, png_out=png)
    assert res["ok"] is True, res.get("error")
    assert res.get("plot_error") is None
    assert png.is_file() and png.stat().st_size > 0


@pytestmark_gmsh
def test_poincare_forwards_the_tracer_failure_verbatim(tmp_path):
    """A bad tracer argument must not be masked by the section layer."""
    msh = _write(tmp_path, _uniform_b(), "uniform_b.msh")
    res = poincare(msh, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0], KE_EV,
                   [0.0, 0.2, 0.0], [0.0, 1.0, 0.0], time_step=3,
                   out_file=tmp_path / "bad.pos")
    assert res["ok"] is False
    assert "time_step 3 out of range" in res["error"]


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"crossing_direction": "up"}, "crossing_direction must be one of"),
        ({"plane_normal": [0.0, 0.0, 0.0]}, "plane_normal must be nonzero"),
        ({"plane_normal": [0.0, 1.0]}, "plane_normal needs exactly 3"),
        ({"plane_point": [0.0, math.nan, 0.0]},
         "plane_point must be a finite number"),
    ],
)
def test_poincare_rejects_bad_planes(input_file, forbid_subprocess,
                                     kwargs, message):
    args = {"plane_point": [0.0, 0.0, 0.0], "plane_normal": [0.0, 1.0, 0.0]}
    args.update(kwargs)
    result = poincare(input_file, [[0.0, 0.0, 0.0]], [1.0, 0.0, 0.0],
                      1.0e3, args["plane_point"], args["plane_normal"],
                      **{k: v for k, v in kwargs.items()
                         if k not in ("plane_point", "plane_normal")})
    assert result["ok"] is False
    assert message in result["error"]
