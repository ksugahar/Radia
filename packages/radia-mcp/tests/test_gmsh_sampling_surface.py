"""Guards: Gmsh 5 resampling of a planar field-animation record.

Assumes: an analytic affine field J(x, z, t) on a non-uniform 5 x 3 grid with
unequal axis spacing.  Quad4 interpolation reproduces an affine field exactly,
so native Gmsh probes at non-node points must return it to rounding; distinct
signed coefficients per component and step expose axis permutations, sign
flips, id reordering and step/time mix-ups.  The native test needs the qualified
Gmsh 5 runtime (RADIA_GMSH_PROVENANCE) and skips visibly without it.
"""
import os

import numpy as np
import pytest

h5py = pytest.importorskip("h5py")
pytest.importorskip("radia.gmsh_post_export")

from radia_mcp.gmsh import sampling_surface as ss

Y0 = -2.5e-4
XS = np.array([-4.0e-3, -2.5e-3, 0.5e-3, 1.5e-3, 5.0e-3])     # 5 nodes, non-uniform
ZS = np.array([-1.0e-3, 2.0e-3, 6.0e-3])                       # 3 nodes, other spacing
TIMES = np.array([0.0, 0.125, 0.25, 0.5, 0.875, 1.25])


def _affine(t, x, z):
    """Distinct signed coefficients per component; time enters each differently."""
    return np.stack([
        3.0e6 * x - 7.0e5 * z + 11.0 * (1 + t),
        -2.0e5 * x + 9.0e6 * z - 5.0 * t,
        4.0e5 * x * (1 - 2 * t) + 1.3e6 * z + 2.0,
    ], axis=-1)


def _source(path, permute=False):
    pts = np.array([[x, Y0, z] for x in XS for z in ZS])
    ids = np.arange(len(pts))
    J = np.stack([_affine(t, pts[:, 0], pts[:, 2]) for t in TIMES])
    if permute:   # store samples in a shuffled order; ids keep the association
        order = np.random.default_rng(3).permutation(len(ids))
        pts, ids, J = pts[order], ids[order], J[:, order]
    with h5py.File(path, "w") as f:
        f.attrs["schema"] = ss.SCHEMA
        f.attrs["case"] = "analytic_affine"
        f["time_s"] = TIMES
        f["solver_step"] = np.arange(10, 10 + len(TIMES))
        f["magnet/center_m"] = np.tile([0.0, 0.01, 0.0], (len(TIMES), 1))
        f["magnet/rotation"] = np.tile(np.eye(3), (len(TIMES), 1, 1))
        f["magnet/size_m"] = np.full(3, 0.005)
        f["magnet/magnetization_A_per_m"] = np.array([0.0, 1.0e6, 0.0])
        f["conductor/box_min_m"] = np.array([-0.006, -0.0005, -0.006])
        f["conductor/box_max_m"] = np.array([0.006, 0.0, 0.006])
        f["samples/id"] = ids
        f["samples/point_m"] = pts
        f["field/J_A_per_m2"] = J
    return pts, ids


QUERY = np.array([[-3.1e-3, Y0, 0.7e-3], [0.9e-3, Y0, 4.4e-3], [4.2e-3, Y0, -0.3e-3],
                  [-1.0e-3, Y0, 5.5e-3], [1.5e-3, Y0, 2.0e-3]])   # last one is a node
QUERY_IDS = np.array([40, 7, 13, 2, 99])


def test_grid_association_is_exact_and_rejects_scatter(tmp_path):
    pts, ids = _source(tmp_path / "s.h5", permute=True)
    normal, axes, (u, w, c0), shape = ss.grid_association(pts, ids)
    assert (normal, axes, shape, c0) == (1, [0, 2], (5, 3), Y0)
    assert np.array_equal(u, XS) and np.array_equal(w, ZS)
    bad = pts.copy(); bad[0, 0] += 1e-4
    with pytest.raises(ValueError, match="grid|rectilinear"):
        ss.grid_association(bad, ids)
    with pytest.raises(ValueError, match="ids"):
        ss.grid_association(pts, ids + 1)
    tilted = pts.copy(); tilted[:, 1] += 1e-4 * tilted[:, 0]
    with pytest.raises(ValueError, match="plane"):
        ss.grid_association(tilted, ids)


@pytest.mark.parametrize("point,match", [
    ([0.0, Y0 + 1e-6, 0.0], "off the sampling plane"),
    ([6.0e-3, Y0, 0.0], "outside"),
    ([XS[0], Y0, 0.0], "outside or on the grid boundary"),
])
def test_query_points_outside_or_off_plane_are_rejected(tmp_path, point, match):
    pts, ids = _source(tmp_path / "s.h5")
    with pytest.raises(ValueError, match=match):
        ss.check_query_points(pts, ids, np.array([point]))


def test_source_defects_are_rejected_before_gmsh(tmp_path):
    src = tmp_path / "s.h5"
    _source(src)
    with h5py.File(src, "a") as f:
        del f["field/J_A_per_m2"]
    with pytest.raises(ValueError, match="no /field/J_A_per_m2"):
        ss.resample(src, tmp_path / "o.h5", QUERY, QUERY_IDS)
    _source(src)
    with h5py.File(src, "a") as f:
        t = f["time_s"][()]; t[3] = t[2]
        f["time_s"][...] = t
    with pytest.raises(ValueError, match="increasing"):
        ss.resample(src, tmp_path / "o.h5", QUERY, QUERY_IDS)
    _source(src)
    with pytest.raises(ValueError, match="steps"):
        ss.resample(src, tmp_path / "o.h5", QUERY, QUERY_IDS, steps=[3, 1])


def test_surface_msh_is_strict_msh42(tmp_path):
    from radia_mcp.gmsh.msh_inspect import validate_msh
    pts, ids = _source(tmp_path / "s.h5", permute=True)
    with h5py.File(tmp_path / "s.h5", "r") as f:
        J = f["field/J_A_per_m2"][()]
    grid = ss.write_surface_msh(tmp_path / "s.msh", pts, ids, TIMES, J, [0, 2, 5])
    assert grid["n_quad4"] == 8, grid
    report = validate_msh(str(tmp_path / "s.msh"))
    assert report.get("ok", report.get("valid")), report


@pytest.mark.skipif(not os.environ.get("RADIA_GMSH_PROVENANCE"),
                    reason="needs the qualified Gmsh 5 runtime (RADIA_GMSH_PROVENANCE)")
@pytest.mark.parametrize("permute", [False, True])
def test_native_gmsh_reproduces_affine_field_at_non_node_points(tmp_path, permute):
    src = tmp_path / "s.h5"
    _source(src, permute=permute)
    before = ss._sha256(src)
    steps = [0, 2, 5]                       # nonadjacent
    summary = ss.resample(src, tmp_path / "o.h5", QUERY, QUERY_IDS, steps=steps)
    assert ss._sha256(src) == before, "source record must stay unchanged"
    with h5py.File(tmp_path / "o.h5", "r") as g:
        assert g.attrs["schema"] == ss.SCHEMA
        assert np.array_equal(g["time_s"][()], TIMES[steps])
        assert np.array_equal(g["solver_step"][()], 10 + np.array(steps))
        assert np.array_equal(g["samples/id"][()], QUERY_IDS)
        assert g["magnet/rotation"].shape == (3, 3, 3)
        V = g["field/J_A_per_m2"][()]
    expected = np.stack([_affine(t, QUERY[:, 0], QUERY[:, 2]) for t in TIMES[steps]])
    scale = np.abs(expected).max()
    err = np.abs(V - expected).max() / scale
    assert err < 1e-12, f"max rel error {err:.3e}; got {V.tolist()} expected {expected.tolist()}"
    assert summary["gmsh_runtime"]["msh42_roundtrip_probe"] is True


def test_output_never_replaces_the_source_or_an_existing_file(tmp_path):
    """Guards: --out equal to, linked to, or beside an existing file never truncates it."""
    src = tmp_path / "s.h5"
    _source(src)
    before = ss._sha256(src)
    with pytest.raises(ValueError, match="source record"):
        ss.resample(src, src, QUERY, QUERY_IDS)
    with pytest.raises(ValueError, match="source record"):
        ss.resample(src, tmp_path / "." / "s.h5", QUERY, QUERY_IDS)
    try:
        os.link(src, tmp_path / "hard.h5")
    except OSError:
        pass
    else:
        with pytest.raises(ValueError, match="same file"):
            ss.resample(src, tmp_path / "hard.h5", QUERY, QUERY_IDS)
    other = tmp_path / "other.h5"
    other.write_bytes(b"unrelated")
    with pytest.raises(FileExistsError):
        ss.resample(src, other, QUERY, QUERY_IDS)
    assert other.read_bytes() == b"unrelated"
    assert ss._sha256(src) == before, "source record changed"


def test_empty_selection_is_a_clear_value_error(tmp_path):
    src = tmp_path / "s.h5"
    _source(src)
    with pytest.raises(ValueError, match="no solver steps"):
        ss.resample(src, tmp_path / "o.h5", QUERY, QUERY_IDS, steps=[])
    with pytest.raises(ValueError, match="no query points"):
        ss.resample(src, tmp_path / "o.h5", np.empty((0, 3)), np.empty(0, int))


@pytest.mark.parametrize("defect", ["nonfinite", "missing_point", "missing_step"])
def test_malformed_native_values_are_rejected_before_writing(tmp_path, monkeypatch, defect):
    from radia_mcp.gmsh import _gmsh_subprocess
    steps = [0, 2, 5]

    def fake(script, args, **kwargs):
        values = [[[1.0, -2.0, 3.0] for _ in QUERY] for _ in steps]
        if defect == "nonfinite":
            values[1][2][0] = float("nan")
        elif defect == "missing_point":
            values[0].pop()
        else:
            values.pop()
        return {"ok": True, "times": TIMES[steps].tolist(),
                "values": values, "distances": [[0.0] * len(r) for r in values],
                "runtime": {"wheel_sha256": "0" * 64}}

    monkeypatch.setattr(_gmsh_subprocess, "run_gmsh_json_subprocess", fake)
    src = tmp_path / "s.h5"
    _source(src)
    with pytest.raises(RuntimeError, match="shape|non-finite"):
        ss.resample(src, tmp_path / "o.h5", QUERY, QUERY_IDS, steps=steps)
    assert not (tmp_path / "o.h5").exists()
    assert not list(tmp_path.glob(".o.h5.*.tmp")), "temporary output left behind"


def test_publish_does_not_replace_an_output_created_meanwhile(tmp_path, monkeypatch):
    """Guards: publishing is exclusive even if the output appears after the up-front check."""
    from radia_mcp.gmsh import _gmsh_subprocess
    steps = [0, 2]
    monkeypatch.setattr(_gmsh_subprocess, "run_gmsh_json_subprocess", lambda *a, **k: {
        "ok": True, "times": TIMES[steps].tolist(),
        "values": [[[1.0, -2.0, 3.0] for _ in QUERY] for _ in steps],
        "distances": [[0.0] * len(QUERY) for _ in steps], "runtime": {}})
    out = tmp_path / "o.h5"
    real = ss._write_output

    def write_then_race(*args):
        real(*args)
        out.write_bytes(b"created by someone else")

    monkeypatch.setattr(ss, "_write_output", write_then_race)
    src = tmp_path / "s.h5"
    _source(src)
    with pytest.raises(FileExistsError):
        ss.resample(src, out, QUERY, QUERY_IDS, steps=steps)
    assert out.read_bytes() == b"created by someone else"
    assert not list(tmp_path.glob(".o.h5.*.tmp")), "temporary output left behind"
