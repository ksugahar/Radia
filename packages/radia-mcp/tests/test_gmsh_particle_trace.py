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
from radia_mcp.gmsh.post_process import particle_trace

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
