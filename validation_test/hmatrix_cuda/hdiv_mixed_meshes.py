"""Connected tet, pyramid and hex meshes on a box or a ball, for the HDiv mixed-mesh tests and validation."""
from collections import Counter

import numpy as np

_FACES = {4: [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)],
          5: [(0, 1, 2, 3), (0, 1, 4), (1, 2, 4), (2, 3, 4), (3, 0, 4)],
          8: [(0, 3, 2, 1), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7)]}
_BOX_FACES = ((0, 1, 2, 3), (4, 5, 6, 7), (0, 1, 5, 4), (1, 2, 6, 5), (2, 3, 7, 6), (3, 0, 4, 7))


def build_mesh(points, cells):
    """Conforming NGSolve mesh with material ``body``; faces owned by one cell become boundary ``exterior``."""
    import ngsolve as ng
    from netgen.csg import Pnt
    from netgen.meshing import Element2D, Element3D, FaceDescriptor, Mesh, MeshPoint
    mesh = Mesh(dim=3)
    mesh.SetMaterial(1, "body")
    descriptor = mesh.Add(FaceDescriptor(surfnr=1, domin=1, bc=1))
    mesh.SetBCName(0, "exterior")
    ids = [mesh.Add(MeshPoint(Pnt(*map(float, point)))) for point in points]
    xyz = np.asarray(points, dtype=float)
    owners = Counter(tuple(sorted(cell[i] for i in face)) for cell in cells for face in _FACES[len(cell)])
    for cell in cells:
        mesh.Add(Element3D(1, [ids[i] for i in cell]))
        centre = xyz[cell].mean(axis=0)
        for table in _FACES[len(cell)]:
            face = [cell[i] for i in table]
            if owners[tuple(sorted(face))] != 1:
                continue
            corners = xyz[face]
            if np.dot(np.cross(corners[1] - corners[0], corners[2] - corners[0]), corners.mean(axis=0) - centre) < 0:
                face.reverse()
            mesh.Add(Element2D(descriptor, [ids[i] for i in face]))
    return ng.Mesh(mesh)


def oriented(points, cell):
    """Tet or pyramid vertices ordered so that NGSolve gets a positive map."""
    xyz = np.asarray(points, dtype=float)
    if len(cell) == 4:
        a, b, c, apex = xyz[cell]
        return [cell[1], cell[0], cell[2], cell[3]] if np.dot(np.cross(b - a, c - a), apex - a) < 0 else cell
    a, b, _, d, apex = xyz[cell]
    return [cell[0], cell[3], cell[2], cell[1], cell[4]] if np.dot(np.cross(b - a, d - a), apex - a) < 0 else cell


def pyramid_cube(warp=0.):
    """Unit cube split into six pyramids about its center; ``warp`` shears the top so that some bases are non-planar."""
    corners = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1)]
    points = [(x + warp*x*z*(y - .5), y, z) for x, y, z in corners] + [(.5, .5, .5)]
    return build_mesh(points, [oriented(points, list(face) + [8]) for face in _BOX_FACES])


def tet_cube():
    """Unit cube split into twelve tets about its center."""
    points = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1), (.5, .5, .5)]
    cells = []
    for a, b, c, d in _BOX_FACES:
        cells += [oriented(points, [a, b, c, 8]), oriented(points, [a, c, d, 8])]
    return build_mesh(points, cells)


def transition_cube():
    """Unit cube: a hex on z in [0, 1/2], a pyramid on its top face and ten tets filling the rest."""
    points = [(0, 0, 0), (1, 0, 0), (1, 1, 0), (0, 1, 0), (0, 0, .5), (1, 0, .5), (1, 1, .5), (0, 1, .5),
              (0, 0, 1), (1, 0, 1), (1, 1, 1), (0, 1, 1), (.5, .5, .75)]
    cells = [list(range(8)), oriented(points, [4, 5, 6, 7, 12])]
    for a, b, c, d in ((8, 9, 10, 11), (4, 5, 9, 8), (5, 6, 10, 9), (6, 7, 11, 10), (7, 4, 8, 11)):
        cells += [oriented(points, [a, b, c, 12]), oriented(points, [a, c, d, 12])]
    return build_mesh(points, cells)


def mixed_grid(shape, *, lengths=(1., 1., 1.), origin=(0., 0., 0.), mapping=None):
    """Checkerboard box grid of hexes and center-split cells.

    A split cell is six pyramids, each based on a neighbour's hex face or on the exterior; on an exterior z face the
    base is cut into two tets. All three element types meet, the mesh is conforming, and ``mapping`` moves the points.
    """
    nx, ny, nz = shape
    index = lambda i, j, k: i + (nx + 1)*(j + (ny + 1)*k)
    points = [(origin[0] + lengths[0]*i/nx, origin[1] + lengths[1]*j/ny, origin[2] + lengths[2]*k/nz)
              for k in range(nz + 1) for j in range(ny + 1) for i in range(nx + 1)]
    cells = []
    for k in range(nz):
        for j in range(ny):
            for i in range(nx):
                corners = [index(i, j, k), index(i + 1, j, k), index(i + 1, j + 1, k), index(i, j + 1, k),
                           index(i, j, k + 1), index(i + 1, j, k + 1), index(i + 1, j + 1, k + 1), index(i, j + 1, k + 1)]
                if (i + j + k) % 2:
                    cells.append(corners)
                    continue
                center = len(points)
                points.append(tuple(np.mean([points[c] for c in corners], axis=0)))
                for number, table in enumerate(_BOX_FACES):
                    base = [corners[v] for v in table]
                    if (number == 0 and k == 0) or (number == 1 and k == nz - 1):
                        low = int(np.argmin(base))
                        a, b, c, d = base[low:] + base[:low]
                        cells += [oriented(points, [a, b, c, center]), oriented(points, [a, c, d, center])]
                    else:
                        cells.append(oriented(points, base + [center]))
    if mapping is not None:
        points = [tuple(mapping(*point)) for point in points]
    return build_mesh(points, cells)


def to_ball(radius):
    """Smooth map of the cube [-1, 1]^3 onto a ball; straight cells converge to the sphere as they refine.

    The checkerboard grid of ``mixed_grid`` stays one-signed on it up to n = 6. Its corners flatten as n rises, the
    corner pyramids invert at n = 8, and the mixed route refuses them. The radial max-norm map inverts pyramids from n = 3.
    """
    def mapping(x, y, z):
        a, b, c = x*x, y*y, z*z
        return (radius*x*np.sqrt(1 - b/2 - c/2 + b*c/3), radius*y*np.sqrt(1 - c/2 - a/2 + c*a/3),
                radius*z*np.sqrt(1 - a/2 - b/2 + a*b/3))
    return mapping
