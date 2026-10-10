"""Native P1 entry ownership, orientation, cache and lifecycle contracts.

Numerical independence from the pre-refactor kernel is checked by
run_hacapk_entry_identity.py's separate baseline/compare workflow.
"""
import gc

import ngsolve as ng
import numpy as np
import pytest

from radia import _radia_pybind as native


def geometry(reverse=False, curved=False):
    tetra = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    faces = np.array([[0, 2, 1], [0, 1, 3], [1, 2, 3], [2, 0, 3]], dtype=np.int64)
    # Disconnected copies exercise all three regular-distance quadrature bins.
    vertices = np.concatenate([tetra+[distance, 0., 0.] for distance in (0., 3., 7., 15.)])
    triangles = np.concatenate([faces+4*k for k in range(4)])
    for k in range(len(triangles)):
        triangles[k] = np.roll(triangles[k], k % 3)
    if reverse:
        triangles = triangles[:, [0, 2, 1]].copy()
    corners = vertices[triangles]
    nodes = np.concatenate([corners, .5*(corners+np.roll(corners, -1, axis=1))], axis=1)
    if curved:
        nodes[:, 3:, 2] += .01
    return vertices, triangles, nodes


@pytest.mark.parametrize('reverse', [False, True])
@pytest.mark.parametrize('curved', [False, True])
def test_owned_entries_and_cache_independence(reverse, curved):
    vertices, triangles, nodes = geometry(reverse, curved)
    ng.SetNumThreads(1)
    with ng.TaskManager():
        expected = native._AssembleSLDL_Galerkin(vertices, triangles, nodes, 7, 6, 1)
        sl, dl, provider = native._CreateP1HACApKGeometry(vertices, triangles, nodes, 7, 6, 0)
        _, _, cached = native._CreateP1HACApKGeometry(vertices, triangles, nodes, 7, 6, 4096)
        # The managers/providers own snapshots, not borrowed NumPy buffers.
        vertices[:] = np.nan; triangles[:] = -1; nodes[:] = np.nan
        del vertices, triangles, nodes
        gc.collect()
        for i in range(16):
            for j in range(16):
                value = provider.Entry(i, j)
                assert value == cached.Entry(i, j)
                m = value[4]
                gamma = (2*m+4)*np.finfo(float).eps/(1-(2*m+4)*np.finfo(float).eps)
                assert abs(value[0]-expected[0][i, j]) <= gamma*value[2]
                assert abs(value[1]-expected[1][i, j]) <= gamma*value[3]
        assert provider.GetStats()['cache_bytes'] == 0
        assert 0 < cached.GetStats()['cache_bytes'] <= 4096
        for handle, matrix in zip((sl, dl), expected):
            with pytest.raises(RuntimeError, match='not been built'):
                handle.MatVec(np.ones(16))
            for _ in range(2):
                assert handle.BuildHMatrix(aca_eps=1e-10, leaf_size=64)
                for transpose in (False, True):
                    result = handle.MatVec(np.arange(16, dtype=float), transpose=transpose)
                    reference = (matrix.T if transpose else matrix) @ np.arange(16)
                    assert np.linalg.norm(result-reference) <= 1e-12*np.linalg.norm(reference)
            assert handle.GetStats()['source_dense_bytes'] == 0
            with pytest.raises(RuntimeError, match='no dense'):
                handle.ReleaseDenseEntries()


def test_bad_geometry_and_entry_indices_fail_loudly():
    vertices, triangles, nodes = geometry()
    with ng.TaskManager():
        _, _, provider = native._CreateP1HACApKGeometry(vertices, triangles, nodes, 7, 6, 0)
        for i, j in ((-1, 0), (0, -1), (16, 0), (0, 16)):
            with pytest.raises(IndexError):
                provider.Entry(i, j)
        invalid = triangles.copy(); invalid[0, 0] = len(vertices)
        with pytest.raises(ValueError, match='outside'):
            native._CreateP1HACApKGeometry(vertices, invalid, nodes)
        invalid = nodes.copy(); invalid[0, 0, 0] += 1
        with pytest.raises(ValueError, match='corners'):
            native._CreateP1HACApKGeometry(vertices, triangles, invalid)
        with pytest.raises(ValueError, match='requires coordinates'):
            native._CreateP1HACApKGeometry(vertices.ravel(), triangles, nodes)
