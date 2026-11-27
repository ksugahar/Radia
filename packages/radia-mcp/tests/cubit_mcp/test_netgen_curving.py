"""netgen_curving: GeometricRefit and the periodic coefficient copy outside Netgen.

Both moved from the Radia Netgen fork into
``native/cubit_vol_export/netgen_curving``.  The ``_netgen_curving`` test module
is built from the same sources with the Netgen SDK only (no Cubit):

    cmake -S native/cubit_vol_export -B <build> -DCUBIT_VOL_EXPORT_PLUGIN=OFF
          -DNETGEN_CURVING_PYMODULE=ON -DNetgen_DIR=<venv>/Lib/site-packages/netgen/cmake

Point ``RADIA_NETGEN_CURVING_DIR`` at the directory holding the built module.
``RADIA_REQUIRE_NETGEN_CURVING=1`` turns a missing module into a failure.

While the fork still carries its own copies (``Mesh.GeometricRefit``,
``Mesh.CopyPeriodicCurvedCoefficients``), the parity tests require the moved
code to give the same coefficients bit for bit.

Guards: the refit and periodic copy keep the fork's results and exception
safety after the move out of Netgen.
Assumes: a Netgen runtime with the coefficient accessors used by
netgen_curving/curved_access.hpp (Radia fork r17 or later).
"""
import importlib
import math
import os
import sys

import pytest

# With RADIA_REQUIRE_NETGEN_CURVING=1 nothing here may skip: a missing Netgen,
# OCC or test module fails, and the module must come from RADIA_NETGEN_CURVING_DIR.
REQUIRED = os.environ.get("RADIA_REQUIRE_NETGEN_CURVING") == "1"


def _require(name):
    if REQUIRED:
        return importlib.import_module(name)
    return pytest.importorskip(name)


netgen_meshing = _require("netgen.meshing")


def _load():
    where = os.environ.get("RADIA_NETGEN_CURVING_DIR")
    if REQUIRED and not where:
        raise RuntimeError("RADIA_REQUIRE_NETGEN_CURVING=1 needs RADIA_NETGEN_CURVING_DIR")
    if where:
        sys.path.insert(0, where)
    try:
        module = importlib.import_module("_netgen_curving")
    except ImportError as exc:
        if REQUIRED:
            raise
        pytest.skip(f"_netgen_curving not built: {exc}")
    finally:
        if where:
            sys.path.remove(where)
    if where:
        built = os.path.normcase(os.path.realpath(where))
        loaded = os.path.normcase(os.path.realpath(os.path.dirname(module.__file__)))
        assert loaded == built, (loaded, built)
    return module


nc = _load()
from netgen.csg import CSGeometry, Pnt, Sphere  # noqa: E402
from netgen.meshing import NetgenGeometry, Vec  # noqa: E402

_KEEP = []   # the mesh holds the C++ geometry; Python must keep the subclass alive


def _unit(p):
    r = math.sqrt(p[0] ** 2 + p[1] ** 2 + p[2] ** 2)
    return [p[0] / r, p[1] / r, p[2] / r]


class UnitSphere(NetgenGeometry):
    """Unit sphere in Python; raises on callback number `fail_at`."""

    def __init__(self, fail_at=None):
        super().__init__()
        self.fail_at = fail_at
        self.calls = 0

    def _tick(self):
        self.calls += 1
        if self.calls == self.fail_at:
            raise RuntimeError("injected projection failure")

    def ProjectPointGI(self, surf, p, gi):
        self._tick()
        p[0], p[1], p[2] = _unit([p[0], p[1], p[2]])
        return True

    def GetNormal(self, surf, p, gi):
        self._tick()
        return Vec(*_unit([p[0], p[1], p[2]]))

    def PointBetween(self, p1, p2, s, surf, gi1, gi2, newp, newgi):
        q = _unit([p1[i] + s * (p2[i] - p1[i]) for i in range(3)])
        for i in range(3):
            newp[i] = q[i]


def _set(mesh, geo):
    _KEEP.append(geo)
    mesh.SetGeometry(geo)
    return geo


def _curved_sphere():
    g = CSGeometry()
    g.Add(Sphere(Pnt(0, 0, 0), 1))
    m = g.GenerateMesh(maxh=0.6)
    _set(m, UnitSphere())
    m.Curve(3)
    return m


def _lines(mesh, tmp_path, name):
    f = os.path.join(str(tmp_path), name)
    mesh.Save(f)
    with open(f) as fh:
        return fh.read().splitlines()


# ---- ported from the fork's tests/pytest/test_geometric_refit.py -----------

def test_refit_rejects_second_order_elements_before_changing_coefficients(tmp_path):
    g = CSGeometry()
    g.Add(Sphere(Pnt(0, 0, 0), 1))
    m = g.GenerateMesh(maxh=0.6)
    m.SecondOrder()
    m.Curve(3)
    before = _lines(m, tmp_path, "before.vol")
    with pytest.raises(Exception, match="unsupported surface element type"):
        nc.geometric_refit(m)
    assert _lines(m, tmp_path, "after.vol") == before


def test_refit_exception_leaves_only_original_or_accepted_coefficients(tmp_path):
    original = _lines(_curved_sphere(), tmp_path, "orig.vol")
    full = _curved_sphere()
    geo = _set(full, UnitSphere())
    stats = nc.geometric_refit(full)
    assert stats["edges_accepted"] > 0 and stats["faces_accepted"] > 0, stats
    accepted = _lines(full, tmp_path, "full.vol")
    # Early, middle and late failures: inside an edge trial, then a face trial.
    for k in (45, 47, geo.calls // 3, 2 * geo.calls // 3):
        m = _curved_sphere()
        _set(m, UnitSphere(fail_at=k))
        with pytest.raises(Exception, match="injected"):
            nc.geometric_refit(m)
        got = _lines(m, tmp_path, "k%d.vol" % k)
        assert len(got) == len(original) == len(accepted), (len(got), len(original), len(accepted))
        foreign = [i for i, (a, b, c) in enumerate(zip(got, original, accepted)) if a != b and a != c]
        assert not foreign, "k=%d: %d lines match neither the original nor the accepted value" % (k, len(foreign))


def test_refit_leaves_order_two_unchanged(tmp_path):
    g = CSGeometry()
    g.Add(Sphere(Pnt(0, 0, 0), 1))
    m = g.GenerateMesh(maxh=0.6)
    _set(m, UnitSphere())
    m.Curve(2)
    before = _lines(m, tmp_path, "o2_before.vol")
    assert nc.geometric_refit(m)["edges_tried"] == 0
    assert _lines(m, tmp_path, "o2_after.vol") == before


# ---- parity with the fork's own copy, while it still has one --------------

fork_refit = pytest.mark.skipif(not hasattr(netgen_meshing.Mesh, "GeometricRefit"),
                                reason="Netgen runtime without the fork's GeometricRefit")
fork_periodic = pytest.mark.skipif(
    not hasattr(netgen_meshing.Mesh, "CopyPeriodicCurvedCoefficients"),
    reason="Netgen runtime without the fork's CopyPeriodicCurvedCoefficients")


@fork_refit
@pytest.mark.parametrize("order", [3, 4])
def test_refit_matches_the_fork_bit_for_bit(tmp_path, order):
    def mesh():
        g = CSGeometry()
        g.Add(Sphere(Pnt(0, 0, 0), 1))
        m = g.GenerateMesh(maxh=0.5)
        _set(m, UnitSphere())
        m.Curve(order)
        return m
    a, b = mesh(), mesh()
    sa = a.GeometricRefit()
    sb = nc.geometric_refit(b)
    assert sa == sb
    assert sa["faces_accepted"] > 0
    assert _lines(a, tmp_path, "fork.vol") == _lines(b, tmp_path, "moved.vol")


def _periodic_box(order):
    """Box with a cylindrical hole; the two end faces are identified periodic."""
    occ = _require("netgen.occ")
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(1, 1, 1))
    hole = occ.Cylinder(occ.Pnt(0.5, 0.5, 0), occ.Z, r=0.25, h=1)
    shape = box - hole
    zmin = shape.faces.Min(occ.Z)
    zmax = shape.faces.Max(occ.Z)
    zmin.Identify(zmax, "periodic", occ.IdentificationType.PERIODIC)
    m = occ.OCCGeometry(shape).GenerateMesh(maxh=0.3)
    m.Curve(order)
    return m, zmin, zmax


def _periodic_triangles(mesh):
    """Master/slave triangle pairs (1-based points) of the identified end faces."""
    heights = [p[2] for p in mesh.Points()]

    def z(nr):   # 1-based point number
        return heights[int(nr) - 1]

    pairs = {}   # z = 0 point -> its z = 1 image
    for a, b, _ in mesh.GetIdentifications():
        a, b = a.nr, b.nr
        lo, hi = (a, b) if z(a) < z(b) else (b, a)
        pairs[lo] = hi
    bottom = []
    for el in mesh.Elements2D():
        pts = [int(p.nr) for p in el.vertices]
        if max(z(p) for p in pts) < 1e-9:
            bottom.append(pts)
    master, slave = [], []
    for tri in bottom:
        if all(p in pairs for p in tri):
            master.append(tri)
            slave.append([pairs[p] for p in tri])
    return master, slave


@fork_periodic
@pytest.mark.parametrize("order", [2, 3])
def test_periodic_copy_matches_the_fork_bit_for_bit(tmp_path, order):
    a, _, _ = _periodic_box(order)
    b, _, _ = _periodic_box(order)
    master, slave = _periodic_triangles(a)
    assert master and len(master) == len(slave)
    assert _periodic_triangles(b) == (master, slave)
    ra = a.CopyPeriodicCurvedCoefficients(master, slave)
    rb = nc.copy_periodic_coefficients(b, master, slave)
    assert ra == rb
    assert any(r["entity"] == "edge" for r in rb)
    assert _lines(a, tmp_path, "fork.vol") == _lines(b, tmp_path, "moved.vol")


def test_periodic_copy_rejects_bad_input_without_changing_coefficients(tmp_path):
    m, _, _ = _periodic_box(3)
    master, slave = _periodic_triangles(m)
    before = _lines(m, tmp_path, "before.vol")
    with pytest.raises(Exception, match="triangle count mismatch"):
        nc.copy_periodic_coefficients(m, master, slave[:-1])
    with pytest.raises(Exception, match="source and destination overlap|conflicting"):
        nc.copy_periodic_coefficients(m, master + slave, slave + master)
    assert _lines(m, tmp_path, "after.vol") == before


def test_periodic_copy_makes_slave_coefficients_equal_the_master(tmp_path):
    m, _, _ = _periodic_box(3)
    master, slave = _periodic_triangles(m)
    rows = nc.copy_periodic_coefficients(m, master, slave)
    assert rows
    again = nc.copy_periodic_coefficients(m, master, slave)
    for r in again:   # a second copy finds every slave already equal
        assert r["original"] == r["expected"], r


def _box_with_planted_master(order=3):
    """Curved sphere mesh; returns it with a master triangle whose edge
    coefficients carry known values and a disjoint slave triangle listed so
    that every image edge runs against the master's numbering."""
    g = CSGeometry()
    g.Add(Sphere(Pnt(0, 0, 0), 1))
    m = g.GenerateMesh(maxh=0.6)
    _set(m, UnitSphere())
    m.Curve(order)
    tris = [sorted(int(p.nr) for p in el.vertices) for el in m.Elements2D()]
    master = tris[0]
    slave = next(t for t in tris if not set(t) & set(master))
    slave = slave[::-1]   # master ascending, slave descending: every sign is -1
    planted = {}
    for a in range(3):
        for b in range(a + 1, 3):
            e = nc.edge_number(m, master[a], master[b])
            n = len(nc.edge_coefficients(m, e))
            assert n == order - 1, (e, n)
            vals = [(1.0 + e + 0.25 * k, -2.0 * k - 0.5, 0.125 * (e + k + 1)) for k in range(n)]
            nc.set_edge_coefficients(m, e, vals)
            planted[(a, b)] = vals
    return m, master, slave, planted


def test_periodic_copy_flips_odd_edge_modes_against_the_numbering():
    m, master, slave, planted = _box_with_planted_master(3)
    rows = nc.copy_periodic_coefficients(m, [master], [slave])
    assert any(r["sign"] == -1 for r in rows), rows
    for (a, b), vals in planted.items():
        sign = (1 if slave[a] < slave[b] else -1) * (1 if master[a] < master[b] else -1)
        assert sign == -1
        got = nc.edge_coefficients(m, nc.edge_number(m, slave[a], slave[b]))
        want = [tuple((sign if k % 2 else 1) * x for x in v) for k, v in enumerate(vals)]
        assert got == want, ((a, b), got, want)


@fork_periodic
def test_periodic_copy_sign_matches_the_fork_bit_for_bit(tmp_path):
    a, master, slave, _ = _box_with_planted_master(3)
    b, master_b, slave_b, _ = _box_with_planted_master(3)
    assert (master, slave) == (master_b, slave_b)
    assert a.CopyPeriodicCurvedCoefficients([master], [slave]) == \
        nc.copy_periodic_coefficients(b, [master], [slave])
    assert _lines(a, tmp_path, "fork_sign.vol") == _lines(b, tmp_path, "moved_sign.vol")


def test_coefficient_helpers_reject_bad_edges():
    m = _curved_sphere()
    with pytest.raises(IndexError):
        nc.edge_coefficients(m, -1)
    with pytest.raises(IndexError):
        nc.edge_coefficients(m, 10 ** 7)
    g = CSGeometry()
    g.Add(Sphere(Pnt(0, 0, 0), 1))
    flat = g.GenerateMesh(maxh=0.6)   # never curved
    with pytest.raises(RuntimeError, match="not curved"):
        nc.edge_coefficients(flat, 0)
