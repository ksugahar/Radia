"""Exercise rejected inputs and assert no geometry callback executes.

These tests do not observe the internal UpdateTopology call.
"""
import pytest


def _call(volumes, surfaces=(), edges=(), *, coords=None, ids=None):
    import ngsolve
    import numpy as np
    from cubit_mesh_export import cubit_mesh_curver

    points = np.array([
        [0., 0., 0.], [1., 0., 0.], [1., 1., 0.], [0., 1., 0.],
        [0., 0., 1.], [1., 0., 1.], [1., 1., 1.], [0., 1., 1.],
    ])

    def forbidden(*args):
        pytest.fail("Invalid connectivity reached a geometry callback")

    with ngsolve.TaskManager():
        return cubit_mesh_curver.build_curved_mesh(
            points if coords is None else np.asarray(coords, dtype=float),
            np.arange(1, 9, dtype=np.int32) if ids is None else np.asarray(ids, dtype=np.int32), list(volumes),
            list(surfaces), forbidden, forbidden, order=2,
            edge_project_func=forbidden, edge_elements=list(edges),
        )


@pytest.mark.parametrize("kind", ["tet", "hex", "prism", "pyramid", "triangle", "quad", "segment"])
def test_unknown_node_is_rejected_before_geometry(kind):
    volumes = [(0, [1, 2, 3, 4], 1)]
    surfaces, edges = [], []
    volume_types = {"tet": (0, 4), "hex": (1, 8), "prism": (2, 6), "pyramid": (3, 5)}
    if kind in volume_types:
        code, count = volume_types[kind]
        volumes = [(code, [*range(1, count), 99], 1)]
        element_kind = "volume"
    elif kind == "segment":
        surfaces = [dict(id=i, tris=[1, 2, 3], quads=[], uvs={}) for i in (1, 2)]
        edges = [dict(surfnr1=1, surfnr2=2, segments=[[1, 99]])]
        element_kind = "segment"
    else:
        surfaces = [dict(id=7, tris=[1, 2, 99] if kind == "triangle" else [],
                         quads=[1, 2, 3, 99] if kind == "quad" else [], uvs={})]
        element_kind = "triangle" if kind == "triangle" else "quadrilateral"
    with pytest.raises(ValueError, match=rf"{element_kind} element 1 references unknown node ID 99 at local node"):
        _call(volumes, surfaces, edges)


@pytest.mark.parametrize("count", [3, 5])
def test_wrong_volume_connectivity_size_is_not_silently_dropped_or_truncated(count):
    with pytest.raises(ValueError, match=rf"volume element 1 requires 4 node IDs; got {count}"):
        _call([(0, list(range(1, count + 1)), 1)])


def test_unsupported_volume_type_is_rejected_before_geometry():
    with pytest.raises(ValueError, match="unsupported volume element type 99"):
        _call([(99, [1, 2, 3, 4], 1)])


@pytest.mark.parametrize("field,values", [("tris", [1]), ("tris", [1, 2, 3, 4]),
                                         ("quads", [1, 2, 3]), ("quads", [1, 2, 3, 4, 5])])
def test_incomplete_surface_connectivity_is_rejected(field, values):
    surface = dict(id=17, tris=[], quads=[], uvs={})
    surface[field] = values
    with pytest.raises(ValueError, match=rf"surface 17 {field} length must be divisible"):
        _call([], [surface])


@pytest.mark.parametrize("field", ["segments", "dists"])
@pytest.mark.parametrize("count", [0, 1, 3])
def test_edge_pairs_are_neither_dropped_nor_truncated(field, count):
    edge = dict(surfnr1=1, surfnr2=2, segments=[[1, 2]])
    edge[field] = [list(range(1, count + 1))]
    with pytest.raises(ValueError, match="edge 1 .* requires exactly 2"):
        _call([], edges=[edge])


@pytest.mark.parametrize("coords,ids", [([], [1]), ([[0., 0.]], [1]),
    ([[0., 0., 0., 0.]], [1]), ([[0., 0., 0.]], [1, 2]),
    ([[0., 0., 0.], [1., 0., 0.]], [1]), ([[0., 0., 0.]], [[1]])])
def test_node_array_shape_is_checked_before_reading(coords, ids):
    with pytest.raises(ValueError, match="node_coords must have shape"):
        _call([], coords=coords, ids=ids)


def test_duplicate_node_ids_are_rejected():
    with pytest.raises(ValueError, match="duplicate node ID 7"):
        _call([], ids=[1, 2, 3, 4, 5, 6, 7, 7])


@pytest.mark.parametrize("domain", [0, -1])
def test_nonpositive_volume_domain_is_rejected(domain):
    with pytest.raises(ValueError, match="volume element 1 domain must be >= 1"):
        _call([(0, [1, 2, 3, 4], domain)])


@pytest.mark.parametrize("volume", [(), (0,), (0, [1, 2, 3, 4], 1, 2)])
def test_volume_record_arity_is_checked(volume):
    with pytest.raises(ValueError, match="volume element must contain"):
        _call([volume])
