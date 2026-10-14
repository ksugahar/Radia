"""Cross-mesh field comparison + the beam-animation sentinel report.

Fixtures are Kuhn decompositions of the cube [-0.5, 0.5]^3 at three
refinement levels (6, 48 and 384 tets), carrying analytic NodeData, so
every expected number is a closed form:

* LINEAR  S = x + 2y : P1 interpolation is EXACT on any mesh, so two
  different meshings must agree to roundoff.
* QUADRATIC Q = x^2  : the P1 interpolation error is real and O(h^2),
  so refining the pair halves h and must quarter the difference -- the
  convergence order IS the point of the verb.

The beam-animation half locks the honest reporting of the hide-the-
future sentinel: a magnetic field does no work, so an animation view
whose raw max is 11x the static kinetic energy is a rendering artefact
that field_stats must NAME rather than serve as physics.
"""

import importlib.util
import math

import pytest
from radia_mcp.gmsh.compare import compare_fields
from radia_mcp.gmsh.msh_inspect import field_stats

_GMSH_AVAILABLE = importlib.util.find_spec("gmsh") is not None
_needs_gmsh = pytest.mark.skipif(not _GMSH_AVAILABLE,
                                 reason="gmsh package not installed")

# Kuhn decomposition of one cube: 6 positively oriented tets sharing the
# main diagonal, addressed by the local corner code di + 2*dj + 4*dk.
_LOCAL_TETS = ((0, 1, 3, 7), (0, 1, 7, 5), (0, 2, 7, 3),
               (0, 2, 6, 7), (0, 4, 5, 7), (0, 4, 7, 6))


def _kuhn_cube(n):
    """(msh text, node coordinates) for an n x n x n Kuhn tet cube."""
    m = n + 1
    coords = [(-0.5 + a / n, -0.5 + b / n, -0.5 + c / n)
              for c in range(m) for b in range(m) for a in range(m)]

    def nid(a, b, c):
        return 1 + a + m * b + m * m * c

    tets = []
    for k in range(n):
        for j in range(n):
            for i in range(n):
                corner = [nid(i + (loc & 1), j + ((loc >> 1) & 1),
                              k + ((loc >> 2) & 1)) for loc in range(8)]
                tets.extend(tuple(corner[c] for c in t) for t in _LOCAL_TETS)
    lines = ["$MeshFormat", "4.1 0 8", "$EndMeshFormat", "$Nodes",
             f"1 {len(coords)} 1 {len(coords)}", f"3 1 0 {len(coords)}"]
    lines += [str(i + 1) for i in range(len(coords))]
    lines += [f"{x:.17g} {y:.17g} {z:.17g}" for x, y, z in coords]
    lines += ["$EndNodes", "$Elements", f"1 {len(tets)} 1 {len(tets)}",
              f"3 1 4 {len(tets)}"]
    lines += [f"{e + 1} " + " ".join(str(v) for v in t)
              for e, t in enumerate(tets)]
    lines += ["$EndElements", ""]
    return "\n".join(lines), coords


def _nodedata(name, rows, ncomp=1, step=0, time=0.0):
    lines = ["$NodeData", "1", f'"{name}"', "1", str(time), "3", str(step),
             str(ncomp), str(len(rows))]
    for tag, vals in rows:
        lines.append(str(tag) + " " + " ".join(f"{v:.17g}" for v in vals))
    lines.append("$EndNodeData")
    return "\n".join(lines) + "\n"


def _write_cube(tmp_path, n, views=("S", "Q", "V")):
    mesh, coords = _kuhn_cube(n)
    text = mesh
    fields = {
        "S": lambda x, y, z: [x + 2.0 * y],          # linear: P1-exact
        "Q": lambda x, y, z: [x * x],                # quadratic: O(h^2)
        "V": lambda x, y, z: [x + 2.0 * y, -x, 0.5 * z + x],   # linear vec
    }
    for name in views:
        f = fields[name]
        text += _nodedata(name, [(i + 1, f(*p)) for i, p in enumerate(coords)],
                          ncomp=len(f(0.0, 0.0, 0.0)))
    path = tmp_path / f"cube{n}.msh"
    path.write_text(text, encoding="utf-8")
    return path


# ======================================================================
# compare_fields: analytic goldens
# ======================================================================

@_needs_gmsh
def test_linear_field_is_p1_exact_across_different_meshings(tmp_path):
    """S = x + 2y lives in P1, so 6-tet and 48-tet agree to roundoff."""
    a = _write_cube(tmp_path, 1)
    b = _write_cube(tmp_path, 2)
    r = compare_fields(a, b, view_a="S", view_b="S", n_points=500, seed=0)
    assert r["ok"], r.get("error")
    assert r["ncomp"] == 1
    assert r["view_a"] == "S" and r["view_b"] == "S"
    assert r["n_valid"] == 500 and r["n_skipped"] == 0
    # MEASURED (gmsh 4.15.2, seed 0, 500 pts): linf 4.441e-16, l2 7.5e-17.
    assert r["linf"] < 1e-12
    assert r["l2"] < 1e-12
    assert r["worst_diff"] == r["linf"]
    assert len(r["worst_point"]) == 3


@_needs_gmsh
def test_quadratic_difference_shows_second_order_convergence(tmp_path):
    """Q = x^2 is NOT in P1: the pairwise difference must quarter.

    The difference between two meshings is dominated by the coarser
    mesh's interpolation error, which is O(h^2); halving h therefore
    divides it by ~4.  MEASURED (gmsh 4.15.2, seed 0, 1000 pts):
    linf 0.249323 (6 vs 48 tets) and 0.062491 (48 vs 384 tets), ratio
    3.99 -- the second-order rate, which is exactly what this verb
    exists to expose.  (The exact ceiling is 0.25: at x = 0 the P1
    interpolant of x^2 on the 6-tet cube reads the corner average.)
    """
    c1 = _write_cube(tmp_path, 1)
    c2 = _write_cube(tmp_path, 2)
    c4 = _write_cube(tmp_path, 4)
    coarse = compare_fields(c1, c2, view_a="Q", view_b="Q", n_points=1000)
    fine = compare_fields(c2, c4, view_a="Q", view_b="Q", n_points=1000)
    assert coarse["ok"] and fine["ok"]
    assert coarse["linf"] == pytest.approx(0.2493, abs=2e-3)
    assert fine["linf"] == pytest.approx(0.06249, abs=1e-3)
    ratio = coarse["linf"] / fine["linf"]
    assert 3.5 <= ratio <= 4.5, f"O(h^2) ratio drifted: {ratio}"
    assert coarse["l2"] > fine["l2"] > 0.0
    # relative metrics normalize by the larger field RMS / peak
    assert 0.0 < coarse["linf_rel"] <= 1.0


@_needs_gmsh
def test_same_file_compares_to_exact_zero(tmp_path):
    a = _write_cube(tmp_path, 2)
    r = compare_fields(a, a, view_a="Q", view_b="Q", n_points=200)
    assert r["ok"], r.get("error")
    assert r["linf"] == 0.0 and r["l2"] == 0.0
    assert r["l2_rel"] == 0.0 and r["linf_rel"] == 0.0


@_needs_gmsh
def test_vector_views_use_the_vector_difference_norm(tmp_path):
    """A rotated field has the same |B| everywhere; only ||a-b|| sees it."""
    mesh, coords = _kuhn_cube(2)
    rot = mesh + _nodedata(
        "V", [(i + 1, [-(x + 2.0 * y), x, -(0.5 * z + x)])
              for i, (x, y, z) in enumerate(coords)], ncomp=3)
    flipped = tmp_path / "flipped.msh"
    flipped.write_text(rot, encoding="utf-8")
    straight = _write_cube(tmp_path, 2, views=("V",))

    r = compare_fields(straight, flipped, view_a="V", view_b="V",
                       n_points=200)
    assert r["ok"], r.get("error")
    assert r["ncomp"] == 3
    # b = -a everywhere, so ||a - b|| = 2|a| while |a| - |b| = 0
    assert r["l2"] == pytest.approx(2.0 * r["rms_a"], rel=1e-12)
    assert r["linf"] > 1.0


@_needs_gmsh
def test_grid_sampling_is_a_lattice_and_seed_is_reproducible(tmp_path):
    a = _write_cube(tmp_path, 1)
    b = _write_cube(tmp_path, 2)
    grid = compare_fields(a, b, view_a="Q", view_b="Q", n_points=27,
                          sample="grid")
    assert grid["ok"] and grid["n_valid"] == 27 and grid["n_skipped"] == 0
    # cell centers of a 3x3x3 lattice hit the x = 0 plane exactly, where
    # the 6-tet interpolant of x^2 reads the 0.25 corner average
    assert grid["linf"] == pytest.approx(0.25, abs=1e-12)

    one = compare_fields(a, b, view_a="Q", view_b="Q", n_points=100, seed=7)
    two = compare_fields(a, b, view_a="Q", view_b="Q", n_points=100, seed=7)
    other = compare_fields(a, b, view_a="Q", view_b="Q", n_points=100, seed=8)
    assert one["worst_point"] == two["worst_point"]
    assert one["linf"] == two["linf"]
    assert other["worst_point"] != one["worst_point"]


@_needs_gmsh
def test_bbox_restricts_sampling_and_writes_a_pos_cloud(tmp_path):
    a = _write_cube(tmp_path, 1)
    b = _write_cube(tmp_path, 4)
    out = tmp_path / "diff.pos"
    r = compare_fields(a, b, view_a="Q", view_b="Q", n_points=64,
                       sample="grid",
                       bbox=[[-0.1, -0.1, -0.1], [0.1, 0.1, 0.1]],
                       out_file=out)
    assert r["ok"], r.get("error")
    assert r["bbox"] == [[-0.1, -0.1, -0.1], [0.1, 0.1, 0.1]]
    assert out.is_file()
    text = out.read_text()
    assert text.startswith('View "field difference magnitude" {')
    assert text.count("SP(") == r["n_valid"]


@_needs_gmsh
def test_single_shared_view_name_needs_no_selector(tmp_path):
    a = _write_cube(tmp_path, 1, views=("Q",))
    b = _write_cube(tmp_path, 2, views=("Q",))
    r = compare_fields(a, b, n_points=50)
    assert r["ok"], r.get("error")
    assert r["view_a"] == "Q" and r["view_b"] == "Q"


@_needs_gmsh
def test_ambiguous_or_missing_view_selection_fails_loudly(tmp_path):
    a = _write_cube(tmp_path, 1)
    b = _write_cube(tmp_path, 2)
    ambiguous = compare_fields(a, b, n_points=10)
    assert not ambiguous["ok"]
    assert "share 3 view names" in ambiguous["error"]
    assert "'Q'" in ambiguous["error"] and "'S'" in ambiguous["error"]

    missing = compare_fields(a, b, view_a="B", view_b="Q", n_points=10)
    assert not missing["ok"]
    assert "'B' not found" in missing["error"]
    assert "'Q'" in missing["error"]          # the valid alternatives

    disjoint = _write_cube(tmp_path, 2, views=("V",))
    none_shared = compare_fields(_write_cube(tmp_path, 1, views=("Q",)),
                                 disjoint, n_points=10)
    assert not none_shared["ok"]
    assert "share no view name" in none_shared["error"]


@_needs_gmsh
def test_disjoint_meshes_report_the_empty_overlap(tmp_path):
    a = _write_cube(tmp_path, 1, views=("Q",))
    shifted_text = a.read_text().replace("$Nodes", "$Nodes", 1)
    # move b far away by rewriting its coordinates
    mesh, coords = _kuhn_cube(1)
    moved = mesh.replace("-0.5", "9.5").replace("0.5", "10.5")
    b = tmp_path / "far.msh"
    b.write_text(moved + _nodedata(
        "Q", [(i + 1, [x * x]) for i, (x, y, z) in enumerate(coords)]),
        encoding="utf-8")
    assert shifted_text  # (a is untouched; guards against a silent typo)
    r = compare_fields(a, b, view_a="Q", view_b="Q", n_points=10)
    assert not r["ok"]
    assert "do not overlap" in r["error"]


def test_invalid_arguments_are_rejected_without_gmsh(tmp_path):
    a = _write_cube(tmp_path, 1)
    b = _write_cube(tmp_path, 2)
    assert "not found" in compare_fields("missing.msh", b)["error"]
    bad = compare_fields(a, b, sample="sobol")
    assert "sample must be one of ['random', 'grid']" in bad["error"]
    assert "n_points must be >= 1" in compare_fields(a, b, n_points=0)["error"]
    assert "6 finite numbers" in compare_fields(a, b, bbox=[1, 2, 3])["error"]
    assert "below bbox min" in compare_fields(
        a, b, bbox=[[1, 1, 1], [0, 2, 2]])["error"]


# ======================================================================
# field_stats: the beam-animation sentinel must be NAMED, not served
# ======================================================================

def _beam_msh(tmp_path, values, n_steps=4, name="beam (kinetic energy [eV])",
              static=None):
    """Line-element mesh with a multi-step hide-the-future beam view."""
    n = len(values)
    lo, hi = min(values), max(values)
    span = hi - lo
    scale = max(span, 1e-3 * max(abs(lo), abs(hi)))
    if scale <= 0.0:
        scale = 1.0
    if span < 1e-12 * scale:
        hi = lo + 1e-3 * scale
    sentinel = hi + 10.0 * scale

    nodes = [(i, 0.0, 0.0) for i in range(n + 1)]
    lines = ["$MeshFormat", "4.1 0 8", "$EndMeshFormat", "$Nodes",
             f"1 {n + 1} 1 {n + 1}", f"1 1 0 {n + 1}"]
    lines += [str(i + 1) for i in range(n + 1)]
    lines += [f"{x:.17g} {y:.17g} {z:.17g}" for x, y, z in nodes]
    lines += ["$EndNodes", "$Elements", f"1 {n} 1 {n}", f"1 1 1 {n}"]
    lines += [f"{e + 1} {e + 1} {e + 2}" for e in range(n)]
    lines += ["$EndElements", ""]
    text = "\n".join(lines)

    def _elemdata(view, step, vals):
        rows = ["$ElementData", "1", f'"{view}"', "1", str(float(step)), "3",
                str(step), "1", str(len(vals))]
        rows += [f"{i + 1} {v:.17g}" for i, v in enumerate(vals)]
        rows.append("$EndElementData")
        return "\n".join(rows) + "\n"

    if static is not None:
        text += _elemdata(static, 0, values)
    visible = []
    for k in range(n_steps):
        head = round(n * (k + 1) / n_steps)
        text += _elemdata(name, k,
                          [v if i < head else sentinel
                           for i, v in enumerate(values)])
        visible.extend(values[:head])
    path = tmp_path / "beam.msh"
    path.write_text(text, encoding="utf-8")
    # the pooled VISIBLE entries: the physical statistics are taken over
    # exactly these, with the same all-steps pooling as the raw ones
    return path, sentinel, visible


def test_field_stats_names_the_beam_sentinel_and_reports_physics(tmp_path):
    values = [1.5e7 + 1.0e6 * i for i in range(11)]     # 1.5e7 .. 2.5e7 eV
    path, sentinel, visible = _beam_msh(tmp_path, values,
                                        static="kinetic energy [eV]")
    r = field_stats(path)
    assert r["ok"]
    beam = next(v for v in r["views"] if v["name"].startswith("beam ("))
    static = next(v for v in r["views"]
                  if v["name"] == "kinetic energy [eV]")

    # the raw headline numbers stay exactly as the file has them ...
    assert beam["overall"]["max"] == pytest.approx(sentinel, rel=1e-12)
    assert beam["overall"]["max"] > 4.0 * static["overall"]["max"]
    # ... and the sentinel is NAMED next to them
    assert beam["beam_animation"] is True
    assert beam["sentinel"] == pytest.approx(sentinel, rel=1e-12)
    assert beam["physical_max"] == pytest.approx(2.5e7, rel=1e-12)
    assert beam["physical_min"] == pytest.approx(1.5e7, rel=1e-12)
    # physical stats must match the static per-track view exactly
    assert beam["physical_max"] == static["overall"]["max"]
    assert beam["physical_min"] == static["overall"]["min"]
    assert beam["physical_rms"] == pytest.approx(
        math.sqrt(sum(v * v for v in visible) / len(visible)), rel=1e-12)
    assert beam["physical_samples"] == len(visible)
    assert beam["physical_samples"] + beam["sentinel_entities"] == 44

    warning = beam["warning"]
    assert "SENTINEL" in warning and "NOT physics" in warning
    assert "kinetic energy [eV]" in warning
    assert warning in r["warnings"]
    # a plain static view carries none of this machinery
    assert "sentinel" not in static and "warning" not in static
    assert r["warnings"] == [warning]

    # the hidden count falls to zero on the last, fully drawn frame
    hidden = [s["sentinel_entities"] for s in beam["per_step"]]
    assert hidden == sorted(hidden, reverse=True)
    assert hidden[-1] == 0 and hidden[0] > 0


def test_field_stats_flags_a_time_valued_beam_that_used_to_be_poisoned(
        tmp_path):
    """The regression that motivated the relative sentinel floor.

    Time values are ~1e-6 s; the old ABSOLUTE 1.0 floor put the sentinel
    at 10.0, so field_stats answered "max 10 s" for a 6.7 us flight.
    With the relative floor the sentinel is ~11x the data -- still
    hidden, no longer astronomically off -- and it is reported as a
    sentinel either way.
    """
    values = [1.68e-10 + 6.72e-6 * i / 20 for i in range(21)]
    path, sentinel, _visible = _beam_msh(tmp_path, values,
                                        name="beam (time [s])")
    r = field_stats(path)
    beam = next(v for v in r["views"] if v["name"].startswith("beam ("))
    assert beam["sentinel"] == pytest.approx(sentinel, rel=1e-12)
    assert beam["physical_max"] == pytest.approx(max(values), rel=1e-12)
    assert sentinel < 1e-4                      # not the old 10.0
    assert "time [s]" in beam["warning"]


def test_fully_visible_beam_view_reports_no_sentinel(tmp_path):
    """One frame, nothing hidden: honest 'no sentinel found' note."""
    values = [1.0 + 0.1 * i for i in range(8)]
    path, _sentinel, _visible = _beam_msh(tmp_path, values, n_steps=1)
    r = field_stats(path)
    beam = next(v for v in r["views"] if v["name"].startswith("beam ("))
    assert beam["beam_animation"] is True
    assert beam["sentinel"] is None
    assert "physical_max" not in beam
    assert "no sentinel value was found" in beam["warning"]
    assert beam["overall"]["max"] == pytest.approx(max(values), rel=1e-12)
