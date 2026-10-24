"""Guards: real sparse storage and operators survive compiler-specific exports.
Assumes: pinned NGSolve/Netgen real-double CPU runtime; caller-owned TaskManager.
"""
import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
from radia.ngsolve_sparse import real_sparse_matrix_type


@pytest.fixture(autouse=True)
def clear_type_cache():
    real_sparse_matrix_type.cache_clear()
    yield
    real_sparse_matrix_type.cache_clear()


@pytest.mark.parametrize("symmetric", [False, True])
def test_type_matches_independently_assembled_storage(symmetric):
    from netgen.geom2d import unit_square
    mesh = ng.Mesh(unit_square.GenerateMesh(maxh=1))
    space = ng.H1(mesh, order=1)
    u, v = space.TnT()
    with ng.TaskManager():
        form = ng.BilinearForm(space, symmetric=symmetric,
                              symmetric_storage=symmetric)
        form += u * v * ng.dx
        form.Assemble()
        actual = real_sparse_matrix_type(symmetric=symmetric)
    assert actual is type(form.mat), (actual, type(form.mat))
    assert real_sparse_matrix_type(symmetric=symmetric) is actual


def test_coo_preserves_nonsymmetric_operator():
    with ng.TaskManager():
        cls = real_sparse_matrix_type()
    matrix = cls.CreateFromCOO([0, 0, 1, 1], [0, 1, 0, 1],
                               [2., 3., 4., 5.], 2, 2)
    source = matrix.CreateColVector()
    source.FV().NumPy()[:] = [7., 11.]
    answer = matrix.CreateRowVector()
    answer.data = matrix * source
    np.testing.assert_array_equal(answer.FV().NumPy(), [47., 83.])


def test_resolution_does_not_require_exported_class_names(monkeypatch):
    with ng.TaskManager():
        types = [real_sparse_matrix_type(symmetric=flag) for flag in (False, True)]
    for cls in types:
        monkeypatch.delattr(ng.la, cls.__name__)
    real_sparse_matrix_type.cache_clear()
    with ng.TaskManager():
        resolved = [real_sparse_matrix_type(symmetric=flag) for flag in (False, True)]
    assert resolved == types


def test_unresolved_storage_raises_without_an_alternate_route(monkeypatch):
    def unavailable(*args, **kwargs):
        raise TypeError("test unavailable storage")
    monkeypatch.setattr(ng, "BilinearForm", unavailable)
    with pytest.raises(RuntimeError, match="cannot resolve NGSolve real-double") as caught:
        real_sparse_matrix_type()
    assert isinstance(caught.value.__cause__, TypeError)
