"""Small NGSolve-owned HDiv prerequisites; not MMM or GPU certification."""
from collections import Counter

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")


def connected_mesh(family):
    from netgen.csg import Pnt
    from netgen.meshing import Element2D, Element3D, FaceDescriptor, Mesh, MeshPoint

    points = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0),
              (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1),
              (.5, .5, 2), (-.6, .5, 1.5)]
    hex_faces = [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4),
                 (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]
    pyramid_faces = [(4, 7, 6, 5), (4, 5, 8), (5, 6, 8),
                     (6, 7, 8), (7, 4, 8)]
    tet_faces = [(4, 7, 8), (4, 9, 7), (7, 9, 8), (8, 9, 4)]
    if family == "tet":
        cells = [([4, 7, 8, 9], tet_faces)]
    elif family == "hex":
        cells = [(list(range(8)), hex_faces)]
    elif family == "hex-pair":
        # Two conforming HEX cells share the complete z=1 quadrilateral.
        points.extend([(0, 0, 2), (1, 0, 2), (1, 1, 2), (0, 1, 2)])
        upper = [4, 5, 6, 7, 10, 11, 12, 13]
        cells = [(list(range(8)), hex_faces),
                 (upper, [tuple(upper[i] for i in face) for face in hex_faces])]
    elif family == "pyramid":
        cells = [([4, 5, 6, 7, 8], pyramid_faces)]
    elif family == "tet-pyramid-hex":
        cells = [(list(range(8)), hex_faces), ([4, 5, 6, 7, 8], pyramid_faces),
                 ([4, 7, 8, 9], tet_faces)]
    else:
        raise ValueError(family)
    mesh = Mesh(dim=3)
    mesh.SetMaterial(1, "body")
    fd = mesh.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    mesh.SetBCName(0, "exterior")
    used = sorted({i for vertices, _ in cells for i in vertices})
    ids = {i: mesh.Add(MeshPoint(Pnt(*points[i]))) for i in used}
    counts = Counter(tuple(sorted(f)) for _, faces in cells for f in faces)
    for vertices, faces in cells:
        mesh.Add(Element3D(1, [ids[i] for i in vertices]))
        centre = np.mean([points[i] for i in vertices], axis=0)
        for face in faces:
            if counts[tuple(sorted(face))] != 1:
                continue
            face = list(face)
            xyz = np.array([points[i] for i in face])
            if np.dot(np.cross(xyz[1]-xyz[0], xyz[2]-xyz[0]), xyz.mean(axis=0)-centre) < 0:
                face.reverse()
            mesh.Add(Element2D(fd, [ids[i] for i in face]))
    return ng.Mesh(mesh)


@pytest.mark.parametrize("family", ["tet", "hex", "hex-pair", "tet-pyramid-hex"])
@pytest.mark.parametrize("order", [1, 2])
def test_connected_hdiv_mass_projection_and_flux_balance(family, order):
    """Avoid dual-shape interpolation: solve the independent L2 mass projection."""
    mesh = connected_mesh(family)
    expected = {"tet": {4}, "hex": {8}, "hex-pair": {8}, "tet-pyramid-hex": {4, 5, 8}}[family]
    assert {len(e.vertices) for e in mesh.Elements(ng.VOL)} == expected
    if family == "hex-pair":
        facets = Counter(f.nr for e in mesh.Elements(ng.VOL) for f in e.facets)
        assert len(list(mesh.Elements(ng.VOL))) == 2
        assert sum(count == 2 for count in facets.values()) == 1
        assert len(list(mesh.Elements(ng.BND))) == 10
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=order)
        u, v = fes.TnT()
        a = ng.BilinearForm(fes)
        a += u*v*ng.dx
        a.Assemble()
        assert np.all(np.isfinite(a.mat.COO()[2]))
        constant = ng.CoefficientFunction((.3, -.7, 1.1))
        rhs = ng.LinearForm(fes)
        rhs += constant*v*ng.dx
        rhs.Assemble()
        gf = ng.GridFunction(fes, novisual=True)
        gf.vec.data = a.mat.Inverse(fes.FreeDofs(), inverse="sparsecholesky")*rhs.vec
        assert float(ng.Integrate(ng.InnerProduct(gf-constant, gf-constant), mesh)) < 1e-20
        rng = np.random.default_rng(20261020)
        for _ in range(4):
            gf.vec.FV().NumPy()[:] = rng.standard_normal(fes.ndof)
            volume = float(ng.Integrate(ng.div(gf), mesh, order=10))
            skin = float(ng.Integrate(gf*ng.specialcf.normal(3), mesh, ng.BND, order=10))
            assert abs(volume-skin) <= 1e-10*max(1., abs(volume), abs(skin))
