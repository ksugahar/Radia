"""Mesh-quality battery: netgen vs cubit on geometries chosen to
exercise known meshing failure modes, judged by the single gmsh minSICN
referee of `cubit_netgen_quality_compare`.

Cases (ALL order 1, per the 2026-08-06 decision: shape quality is a
property of the corner geometry, so first-order elements are the right
scope for simple quality evaluation; the tool's `order` parameter
remains available separately as a curvature-fidelity diagnostic):
  * thin plate (50:1 aspect)      -- sliver risk
  * gapped C-core (lab archetype) -- small feature between pole faces
  * sphere                        -- curvature-limited surface sizing
  * size sweep on the sphere      -- quality must not collapse as h drops

Headless throughout (driving policy).  Requires Cubit + netgen + gmsh
(+ build123d for the geometry authoring); skipped otherwise.  This lane
asserts CORRECTNESS floors and reports the measured numbers -- timing
comparisons use hibino first or an idle-CI mdx fallback, not here.
"""

import json

import pytest

from radia_mcp.cubit import session as cs

pytestmark = pytest.mark.skipif(
    cs.find_cubit_install() is None,
    reason="Coreform Cubit not installed")


def _compare(step, **kw):
    from radia_mcp.cubit.server import cubit_netgen_quality_compare
    out = json.loads(cubit_netgen_quality_compare(str(step), **kw))
    assert out["status"] == "ok", out
    return {r["route"]: r for r in out["rows"]}


def _report(title, by_route):
    print(f"\n[{title}]")
    for route, row in sorted(by_route.items()):
        if not row.get("by_type"):
            print(f"  {route:10s} ERROR {row.get('kind')}: "
                  f"{str(row.get('error'))[:90]}")
            continue
        for bt in row["by_type"]:
            print(f"  {route:10s} {bt['element']:16s} n={bt['n']:6d} "
                  f"min={bt['min']:.3f} mean={bt['mean']:.3f} "
                  f"neg={bt['negative']} below={bt['below_threshold']}")


def _assert_healthy(by_route, routes, min_floor=0.0):
    for route in routes:
        row = by_route[route]
        assert row.get("by_type"), (route, row)
        for bt in row["by_type"]:
            assert bt["negative"] == 0, (route, bt)
            assert bt["min"] > min_floor, (route, bt)


@pytest.fixture(scope="module")
def b3d():
    return pytest.importorskip("build123d")


def test_thin_plate_aspect_50(tmp_path, b3d):
    """A 50 x 50 x 1 plate: tet meshers must not produce inverted or
    sub-threshold slivers at a size comparable to the thickness."""
    from build123d import Box, export_step
    step = tmp_path / "plate.step"
    export_step(Box(50, 50, 1), str(step))
    by_route = _compare(step, netgen_maxh=2.0, cubit_size=2.0)
    _report("thin plate 50x50x1, size 2", by_route)
    _assert_healthy(by_route, ["netgen", "cubit_tet", "cubit_hex"],
                    min_floor=0.1)


def test_gapped_c_core_archetype(tmp_path, b3d):
    """The lab's c_core archetype: an 8 mm air gap between pole faces of
    an 80x60x25 frame -- small-feature meshing on a real magnet shape."""
    from build123d import export_step
    from radia_mcp.build123d.archetypes import c_core
    step = tmp_path / "c_core.step"
    export_step(c_core(width=80, height=60, depth=25, leg=15, gap=8),
                str(step))
    by_route = _compare(step, netgen_maxh=5.0, cubit_size=5.0,
                        schemes=["netgen", "cubit_tet"])
    _report("gapped c_core, size 5", by_route)
    _assert_healthy(by_route, ["netgen", "cubit_tet"], min_floor=0.1)


def test_sphere_curved_surface(tmp_path, b3d):
    """Curved geometry (order 1): both tet meshers must resolve the
    sphere without slivers at a curvature-limited size."""
    from build123d import Sphere, export_step
    step = tmp_path / "sphere.step"
    export_step(Sphere(1.0), str(step))
    by_route = _compare(step, netgen_maxh=0.4, cubit_size=0.4,
                        schemes=["netgen", "cubit_tet"])
    _report("sphere, size 0.4", by_route)
    _assert_healthy(by_route, ["netgen", "cubit_tet"], min_floor=0.1)
    for route in ("netgen", "cubit_tet"):
        assert by_route[route]["by_type"][0]["element"] == "Tetrahedron 4"


def test_sphere_size_sweep_quality_floor(tmp_path, b3d):
    """h-refinement sweep: element counts must grow monotonically and
    the worst minSICN must never collapse below 0.2 for either tet
    mesher."""
    from build123d import Sphere, export_step
    step = tmp_path / "sphere.step"
    export_step(Sphere(1.0), str(step))

    sizes = [0.6, 0.4, 0.25]
    counts = {"netgen": [], "cubit_tet": []}
    mins = {"netgen": [], "cubit_tet": []}
    for h in sizes:
        by_route = _compare(step, netgen_maxh=h, cubit_size=h,
                            schemes=["netgen", "cubit_tet"])
        _report(f"sphere sweep size {h}", by_route)
        for route in counts:
            bt = by_route[route]["by_type"][0]
            counts[route].append(bt["n"])
            mins[route].append(bt["min"])
    for route in counts:
        assert counts[route] == sorted(counts[route]), (route, counts)
        assert min(mins[route]) > 0.2, (route, mins)
