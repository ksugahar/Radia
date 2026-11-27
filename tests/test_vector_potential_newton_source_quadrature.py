"""Reduced-A Newton and linear paths integrate a non-polynomial source alike.

Guards: ``VectorPotentialSolver.solve_nonlinear_newton`` integrates the source
field B_s with the same quadrature as ``solve_linear`` for the same
``bonus_intorder``. In the linear B-H range the two paths must then give the
same per-element iron B. With the default energy quadrature the Newton path
differed by several times the iron field's norm, while the iron mean agreed.
NGSolve lowers the energy order on tetrahedra only, so the hexahedral case
guards against correcting every element type alike (a uniform +2 order gave
0.85 relative L2 there).
Assumes: the B-H table's first segment is linear (mu_r 14874) and the iron
field stays on it; the source is a point dipole outside the iron, so B_s
varies strongly within the coarse iron elements.
"""

import itertools

import numpy as np
import pytest

ng = pytest.importorskip("ngsolve")
occ = pytest.importorskip("netgen.occ")

from radia.vector_potential_solver import MU_0, VectorPotentialSolver

BH = [[0.0, 0.0], [53.5, 1.0], [87.0, 1.5], [130.0, 1.7], [1.0e5, 2.2]]
MU_R_FIRST_SEGMENT = 1.0 / (MU_0 * 53.5)


def _tet_mesh():
    air = occ.Box(occ.Pnt(-1, -1, -1), occ.Pnt(1, 1, 1))
    air.faces.name = "outer"
    iron = occ.Box(occ.Pnt(-0.15, -0.15, -0.4), occ.Pnt(0.15, 0.15, 0.4))
    iron.mat("iron")
    air = air - iron
    air.mat("air")
    return ng.Mesh(occ.OCCGeometry(occ.Glue([air, iron])).GenerateMesh(maxh=0.25))


def _hex_mesh(n=8):
    from netgen.csg import Pnt
    from netgen.meshing import Element2D, Element3D, FaceDescriptor, Mesh, MeshPoint

    m = Mesh(dim=3)
    m.SetMaterial(1, "air")
    m.SetMaterial(2, "iron")
    fd = m.Add(FaceDescriptor(surfnr=1, domin=1, domout=0, bc=1))
    m.SetBCName(0, "outer")
    xs = np.linspace(-1.0, 1.0, n + 1)
    pid = {ijk: m.Add(MeshPoint(Pnt(*(xs[q] for q in ijk))))
           for ijk in itertools.product(range(n + 1), repeat=3)}
    for i, j, k in itertools.product(range(n), repeat=3):
        cx, cy, cz = ((xs[q] + xs[q + 1]) / 2 for q in (i, j, k))
        dom = 2 if abs(cx) < 0.25 and abs(cy) < 0.25 and abs(cz) < 0.5 else 1
        m.Add(Element3D(dom, [pid[i, j, k], pid[i + 1, j, k], pid[i + 1, j + 1, k],
                              pid[i, j + 1, k], pid[i, j, k + 1], pid[i + 1, j, k + 1],
                              pid[i + 1, j + 1, k + 1], pid[i, j + 1, k + 1]]))
    for a, b, axis, side in itertools.product(range(n), range(n), range(3), (0, n)):
        def corner(u, w):
            ijk = [0, 0, 0]
            ijk[axis], ijk[(axis + 1) % 3], ijk[(axis + 2) % 3] = side, u, w
            return pid[tuple(ijk)]
        m.Add(Element2D(fd, [corner(a, b), corner(a + 1, b),
                             corner(a + 1, b + 1), corner(a, b + 1)]))
    return ng.Mesh(m)


def _dipole(position, moment):
    rx, ry, rz = ng.x - position[0], ng.y - position[1], ng.z - position[2]
    r2 = rx * rx + ry * ry + rz * rz
    r = ng.sqrt(r2)
    mr = moment[0] * rx + moment[1] * ry + moment[2] * rz
    c = MU_0 / (4 * np.pi)
    return ng.CF(tuple(c * (3 * mr * q / (r2 * r2 * r) - m / (r2 * r))
                       for q, m in zip((rx, ry, rz), moment)))


@pytest.mark.parametrize("element, order", [("tet", 1), ("tet", 2), ("hex", 1)])
def test_newton_matches_linear_per_iron_element(element, order):
    if element == "tet":
        mesh, position = _tet_mesh(), (0.45, 0.0, 0.0)
    else:
        mesh, position = _hex_mesh(), (0.61, 0.093, 0.071)
    assert {e.type for e in mesh.Elements(ng.VOL)} == {
        ng.ET.TET if element == "tet" else ng.ET.HEX}
    source = _dipole(position, (5.0e3, 0.0, 0.0))
    volume = np.asarray(ng.Integrate(ng.CF(1.0), mesh, element_wise=True))
    iron = np.array([e.mat == "iron" for e in mesh.Elements(ng.VOL)])

    def iron_means(field):
        parts = [np.asarray(ng.Integrate(field[i], mesh, element_wise=True, order=10))
                 for i in range(3)]
        return (np.column_stack(parts) / volume[:, None])[iron]

    with ng.TaskManager():
        linear = VectorPotentialSolver(mesh, mu_r_dict={"iron": MU_R_FIRST_SEGMENT}, order=order)
        linear.set_source_cf(source)
        linear.solve_linear(dirichlet="outer", solver="direct")
        b_linear = iron_means(linear.get_B())
        newton = VectorPotentialSolver(mesh, iron_domains="iron", order=order)
        newton.set_source_cf(source)
        newton.solve_nonlinear_newton(BH, dirichlet="outer", solver="direct",
                                      tol=1e-12, verbose=False)
        b_newton = iron_means(newton.get_B())

    w = volume[iron]
    assert np.abs(b_linear).max() < 0.9, np.abs(b_linear).max()   # stays on the linear segment
    rel = np.sqrt((w * ((b_newton - b_linear) ** 2).sum(1)).sum()
                  / (w * (b_linear ** 2).sum(1)).sum())
    assert rel < 1e-8, (f"{element} order {order}: Newton vs linear iron B relative L2 "
                        f"{rel:.3e}, expected < 1e-8")


def test_newton_rejects_negative_bonus_intorder():
    mesh = ng.Mesh(ng.unit_cube.GenerateMesh(maxh=1.0))
    solver = VectorPotentialSolver(mesh, iron_domains=mesh.GetMaterials()[0], order=1)
    solver.set_source_cf(ng.CF((0.1, 0.0, 0.0)))
    with pytest.raises(ValueError, match="bonus_intorder"):
        solver.solve_nonlinear_newton(BH, bonus_intorder=-1, verbose=False)
