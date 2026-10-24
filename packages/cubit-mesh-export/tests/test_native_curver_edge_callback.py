"""Guards: six-argument edge callback ABI and separate geometric-edge routing.

Assumes: the installed candidate native curver, pinned Netgen, two equatorial
semicircles sharing sphere face descriptors. Source contract, not CAD accuracy.
"""
import numpy as np
import ngsolve as ng
from cubit_mesh_export import cubit_mesh_curver


def test_two_curves_sharing_faces_reach_the_six_argument_python_callback():
    coords = np.array([[0., 0., 0.], [1., 0., 0.], [-1., 0., 0.],
                       [0., 1., 0.], [0., -1., 0.], [0., 0., 1.], [0., 0., -1.]])
    faces = [(x, y, z) for x in (2, 3) for y in (4, 5) for z in (6, 7)]
    elements = []
    for face in faces:
        tet = [1, *face]
        if np.linalg.det(coords[np.array(tet[1:])-1] - coords[0]) < 0:
            tet[2], tet[3] = tet[3], tet[2]
        elements.append((0, tet, 1))
    surfaces = [dict(id=i, tris=[n for f in faces if f[2] == pole for n in f],
                     quads=[], uvs={}) for i, pole in enumerate((6, 7), 1)]
    seen = set()

    def project(surface, x, y, z, u, v, has_hint):
        r = np.linalg.norm([x, y, z])
        return x/r, y/r, z/r, 0., 0.

    def normal(surface, x, y, z):
        r = np.linalg.norm([x, y, z])
        return x/r, y/r, z/r

    def edge_project(edge, face1, face2, x, y, z):
        assert {face1, face2} == {1, 2}, (edge, face1, face2)
        assert edge in (1, 2), edge
        assert y >= 0 if edge == 1 else y <= 0, (edge, x, y, z)
        seen.add(edge)
        r = np.linalg.norm([x, y])
        return x/r, y/r, 0.

    edges = [dict(surfnr1=1, surfnr2=2, segments=[[2, 4], [4, 3]]),
             dict(surfnr1=1, surfnr2=2, segments=[[3, 5], [5, 2]])]
    with ng.TaskManager():
        mesh = cubit_mesh_curver.build_curved_mesh(
            coords, np.arange(1, 8, dtype=np.int32), elements, surfaces,
            project, normal, order=3, edge_project_func=edge_project, edge_elements=edges)
    assert len(list(mesh.Elements3D())) == 8
    assert seen == {1, 2}, seen
