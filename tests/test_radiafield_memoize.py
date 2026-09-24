"""RadiaField memoisation and its use by the nonlinear mixed Omega iterations."""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
import radia as rad  # noqa: E402


def _coil():
    rad.UtiDelAll()
    arc = rad.ObjArcCur([-2.5, 0.0, 0.0], [.3, .6], [0.0, 1.2], .5, 40, 'man', 'z', 2.0e5)
    bar = rad.ObjRecCur([-2.2, 0.8, 0.0], [0.3, 0.8, 0.4], [0, 2.0e5, 0])
    return rad.ObjCnt([arc, bar])


def test_memoised_values_equal_direct_evaluation_and_are_reused():
    coil = _coil()
    mesh = ng.Mesh(ng.unit_cube.GenerateMesh(maxh=0.25))
    direct = rad.RadiaField(coil, "h")
    memo = rad.RadiaField(coil, "h")
    memo.SetMemoize(True)
    assert memo.memoize
    integrand = ng.InnerProduct(memo, ng.CF((1.0, 2.0, 3.0)))
    reference = ng.InnerProduct(direct, ng.CF((1.0, 2.0, 3.0)))
    with ng.TaskManager():  # concurrent lookups and inserts from assembly workers
        first = ng.Integrate(integrand, mesh, order=5)
        size = memo.GetCacheStats()["size"]
        second = ng.Integrate(integrand, mesh, order=5)
        expected = ng.Integrate(reference, mesh, order=5)
    stats = memo.GetCacheStats()
    assert size > 0 and stats["size"] == size
    assert stats["hits"] >= size
    # Parallel summation order may differ; the pointwise values below are exact.
    assert first == pytest.approx(expected, rel=1e-13)
    assert second == pytest.approx(expected, rel=1e-13)
    rng = np.random.default_rng(5)
    points = rng.uniform(0.05, 0.95, size=(300, 3))
    located = mesh(points[:, 0], points[:, 1], points[:, 2])
    once = np.asarray(memo(located))
    again = np.asarray(memo(located))
    np.testing.assert_array_equal(once, np.asarray(direct(located)))
    np.testing.assert_array_equal(again, once)
    memo.SetMemoize(False)
    memo.ClearCache()
    assert memo.GetCacheStats()["size"] == 0 and not memo.memoize


@pytest.mark.parametrize("mode", ["prepare", "memoize"])
def test_mirror_symmetric_points_keep_their_own_values(mode):
    """Regression: the former 64-bit FNV key merged (x,-a,-a) with (x,a,a)."""
    from netgen.occ import Box, OCCGeometry, Pnt

    coil = _coil()
    mesh = ng.Mesh(OCCGeometry(Box(Pnt(-1, -1, -1), Pnt(1, 1, 1))).GenerateMesh(maxh=0.5))
    # The colliding pair observed on an order-2 tet rule of the Picard test mesh.
    points = np.array([[0.06909830056250454, -0.9309016994374945, -0.9309016994374945],
                       [0.0690983005625055, 0.9309016994374945, 0.9309016994374955]])
    located = mesh(points[:, 0], points[:, 1], points[:, 2])
    physical = np.asarray(ng.CF((ng.x, ng.y, ng.z))(located)).reshape(-1, 3)
    direct = np.asarray(rad.RadiaField(coil, "h")(located)).reshape(-1, 3)
    assert np.linalg.norm(direct[0] - direct[1]) > 1e-3 * np.linalg.norm(direct[0])
    field = rad.RadiaField(coil, "h")
    if mode == "prepare":
        field.PrepareCache(physical.tolist())
    else:
        field.SetMemoize(True)
        field(located)
    assert field.GetCacheStats()["size"] == 2
    np.testing.assert_array_equal(np.asarray(field(located)).reshape(-1, 3), direct)
