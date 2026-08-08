"""File-series temporal statistics, and the isosurface "crack" verdict.

``time_series`` is the last ParaView filter the gmsh lane had no answer
for: gmsh's own time steps live inside ONE view, while a transient
solver writes one mesh per step.

The isosurface part records a CORRECTION.  The knowledge used to blame
semi-transparent isosurface "cracks" on hairline T-junctions from
per-element adaptive subdivision.  Measured on a curved p2 field by
counting background pixels enclosed by the surface silhouette, that is
wrong: closed shells give 0 crack pixels at every recursion level 0..3
and at alpha 1.0 and 0.35.  What actually shows through is an OPEN
shell -- one the domain boundary cut -- at 19-29%, independent of the
recursion level.  ``isosurface`` now reports that as a fact.
"""

from __future__ import annotations

import math

import pytest

from radia_mcp.gmsh.post_process import isosurface, time_series

pytest.importorskip("gmsh", reason="gmsh not installed")


def _sphere_field(path, value_at, n=5):
    """Unit box, 5-tet cells, one scalar node view ``f``."""
    nodes, idx, tag = {}, {}, 1
    for i in range(n):
        for j in range(n):
            for k in range(n):
                nodes[tag] = [-1 + 2 * i / (n - 1), -1 + 2 * j / (n - 1),
                              -1 + 2 * k / (n - 1)]
                idx[(i, j, k)] = tag
                tag += 1
    elements, et = {}, 1
    for i in range(n - 1):
        for j in range(n - 1):
            for k in range(n - 1):
                c = [idx[(i, j, k)], idx[(i + 1, j, k)],
                     idx[(i + 1, j + 1, k)], idx[(i, j + 1, k)],
                     idx[(i, j, k + 1)], idx[(i + 1, j, k + 1)],
                     idx[(i + 1, j + 1, k + 1)], idx[(i, j + 1, k + 1)]]
                for a, b, cc, d in ((0, 1, 3, 4), (1, 2, 3, 6), (1, 3, 4, 6),
                                    (1, 4, 5, 6), (3, 4, 6, 7)):
                    elements[et] = [c[a], c[b], c[cc], c[d]]
                    et += 1
    rows = {t: [value_at(p)] for t, p in nodes.items()}
    lo = [min(p[i] for p in nodes.values()) for i in range(3)]
    hi = [max(p[i] for p in nodes.values()) for i in range(3)]
    out = ["$MeshFormat", "4.1 0 8", "$EndMeshFormat",
           "$Entities", "0 0 0 1",
           "1 %g %g %g %g %g %g 0 0" % (*lo, *hi), "$EndEntities",
           "$Nodes", f"1 {len(nodes)} 1 {len(nodes)}",
           f"3 1 0 {len(nodes)}"]
    out += [str(t) for t in nodes]
    out += ["%.15e %.15e %.15e" % tuple(p) for p in nodes.values()]
    out += ["$EndNodes", "$Elements",
            f"1 {len(elements)} 1 {len(elements)}",
            f"3 1 4 {len(elements)}"]
    out += [f"{t} " + " ".join(str(x) for x in ns)
            for t, ns in elements.items()]
    out += ["$EndElements", "$NodeData", "1", '"f"', "1", "0",
            "3", "0", "1", str(len(rows))]
    out += [f"{t} %.9e" % v[0] for t, v in rows.items()]
    out += ["$EndNodeData"]
    path.write_text("\n".join(out) + "\n", encoding="utf-8")
    return path


# --------------------------------------------------------------------
# time series over a FILE series
# --------------------------------------------------------------------

def _series(tmp_path, n_steps=5):
    """A Gaussian bump sweeping along +x: every node's peak occurs at a
    different, ANALYTICALLY KNOWN time.  The bump centres land ON the
    node planes (x = -1, -0.5, 0, 0.5, 1 for n=5), so the sampled peak
    is exactly 1 at every step."""
    paths = []
    for s in range(n_steps):
        x0 = -1.0 + 2.0 * s / (n_steps - 1)
        paths.append(_sphere_field(
            tmp_path / f"step_{s}.msh",
            lambda p, x0=x0: math.exp(-((p[0] - x0) ** 2) / 0.08),
            n=n_steps))
    return paths


def test_time_series_statistics_and_argmax_time(tmp_path):
    paths = _series(tmp_path, 5)
    times = [0.0, 1.0, 2.0, 3.0, 4.0]
    res = time_series(paths, view="f", times=times,
                      out_file=tmp_path / "stats.msh")
    assert res["ok"], res
    assert res["n_steps"] == 5
    assert res["stats_written"][:2] == ["f_min", "f_max"]

    # every node sees the bump pass: max ~ 1 somewhere, min < max
    agg = res["aggregate"]
    assert agg["time"] == times
    assert all(m == pytest.approx(1.0, abs=0.05) for m in agg["max"])

    # argmax_time must MOVE with the bump: read the written view back
    from radia_mcp.gmsh.msh_inspect import read_msh_data

    data = read_msh_data(tmp_path / "stats.msh")
    views = {v["name"]: v for v in data["views"]}
    assert "f_argmax_time" in views and "f_ptp" in views
    nodes = data["nodes"]
    arg = views["f_argmax_time"]["rows"]
    left = [arg[t][0] for t, p in nodes.items() if p[0] < -0.5]
    right = [arg[t][0] for t, p in nodes.items() if p[0] > 0.5]
    assert sum(left) / len(left) < sum(right) / len(right), (
        "argmax_time must increase along the sweep direction")


def test_time_series_point_history_follows_the_bump(tmp_path):
    paths = _series(tmp_path, 5)
    res = time_series(paths, view="f", times=[0, 1, 2, 3, 4],
                      out_file=tmp_path / "s.msh",
                      points=[[-1.0, 0.0, 0.0], [1.0, 0.0, 0.0]])
    assert res["ok"], res
    left, right = res["point_history"]
    assert left["values"].index(max(left["values"])) == 0
    assert right["values"].index(max(right["values"])) == 4


def test_time_series_guards(tmp_path):
    paths = _series(tmp_path, 3)
    with pytest.raises(ValueError, match="at least 2 files"):
        time_series(paths[:1], view="f")
    with pytest.raises(ValueError, match="times has"):
        time_series(paths, view="f", times=[0.0])
    with pytest.raises(ValueError, match="unknown stat"):
        time_series(paths, view="f", stats=("median",))
    missing = time_series(paths, view="nope")
    assert not missing["ok"] and "no view" in missing["error"]


def test_time_series_rejects_a_changed_mesh(tmp_path):
    """A series whose numbering changed is not one time series."""
    a = _sphere_field(tmp_path / "a.msh", lambda p: p[0], n=4)
    b = _sphere_field(tmp_path / "b.msh", lambda p: p[0], n=5)
    res = time_series([a, b], view="f")
    assert not res["ok"]
    assert "tag numbering" in res["error"]


# --------------------------------------------------------------------
# the isosurface "crack" verdict: open shell, not subdivision
# --------------------------------------------------------------------

def test_isosurface_reports_a_closed_shell(tmp_path):
    # f = 1/(0.25 + r^2); level 1.2 closes at r = 0.79, inside the box
    f = _sphere_field(tmp_path / "f.msh",
                      lambda p: 1.0 / (0.25 + p[0] ** 2 + p[1] ** 2
                                       + p[2] ** 2), n=9)
    res = isosurface(f, 1.2, view="f", out_file=tmp_path / "closed.pos")
    assert res["ok"], res
    assert res["n_vertices"] > 0
    assert res["boundary_vertices"] == 0
    assert res["open_surface"] is False
    assert "note" not in res


def test_isosurface_flags_a_shell_cut_by_the_domain(tmp_path):
    # level 0.45 closes at r = 1.44 > the box half-width: cut open
    f = _sphere_field(tmp_path / "f.msh",
                      lambda p: 1.0 / (0.25 + p[0] ** 2 + p[1] ** 2
                                       + p[2] ** 2), n=9)
    res = isosurface(f, 0.45, view="f", out_file=tmp_path / "open.pos")
    assert res["ok"], res
    assert res["open_surface"] is True
    assert res["boundary_vertices"] > 0
    assert "CUT OPEN" in res["note"]


@pytest.mark.parametrize("recur", [0, 1, 2])
def test_closed_shell_stays_closed_at_every_recursion_level(tmp_path, recur):
    """The correction, as a test: adaptive subdivision does not open the
    surface, so it cannot be the source of see-through 'cracks'."""
    f = _sphere_field(tmp_path / "f.msh",
                      lambda p: 1.0 / (0.25 + p[0] ** 2 + p[1] ** 2
                                       + p[2] ** 2), n=9)
    res = isosurface(f, 1.2, view="f", recur_level=recur,
                     out_file=tmp_path / f"c{recur}.pos")
    assert res["ok"], res
    assert res["boundary_vertices"] == 0, (recur, res)
