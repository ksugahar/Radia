"""RadiaField.PrepareCache evaluates in parallel blocks with unchanged values."""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
import radia as rad  # noqa: E402


@pytest.mark.parametrize("inside_taskmanager", [True, False])
def test_parallel_cache_values_equal_direct_evaluation(inside_taskmanager):
    rad.UtiDelAll()
    arc = rad.ObjArcCur([0.5, 0.5, 1.3], [.05, .20], [0.0, 1.2], .15, 40, 'man', 'z', 5.44e5)
    bar = rad.ObjRecCur([0.5, 0.5, -0.3], [0.2, 0.6, 0.1], [0, 5.44e5, 0])
    coil = rad.ObjCnt([arc, bar])
    mesh = ng.Mesh(ng.unit_cube.GenerateMesh(maxh=0.3))
    rng = np.random.default_rng(7)
    # More points than one 256-point block, and not a multiple of it.
    points = rng.uniform(0.05, 0.95, size=(1537, 3))
    located = mesh(points[:, 0], points[:, 1], points[:, 2])
    physical = np.asarray(ng.CF((ng.x, ng.y, ng.z))(located)).reshape(-1, 3)
    field = rad.RadiaField(coil, "h")
    if inside_taskmanager:
        with ng.TaskManager():
            field.PrepareCache(physical.tolist())
    else:
        field.PrepareCache(physical.tolist())
    cached = np.asarray(field(located)).reshape(-1, 3)
    stats = field.GetCacheStats()
    assert stats["size"] == len(points)
    assert stats["hits"] >= 0.99 * len(points)
    direct = np.asarray(rad.Fld(coil, "h", physical.tolist()))
    np.testing.assert_allclose(cached, direct, rtol=1e-13, atol=0.0)
