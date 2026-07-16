"""Lock the NGSolve HDiv family convention used by the Mathematica study."""

import pytest

ng = pytest.importorskip("ngsolve")
from netgen.csg import unit_cube
from netgen.geom2d import unit_square
from ngsolve.meshes import MakeStructured2DMesh, MakeStructured3DMesh


def test_local_hdiv_family_dimensions_match_the_symbolic_ledgers():
    meshes = {
        "trig": ng.Mesh(unit_square.GenerateMesh(maxh=2.0)),
        "tet": ng.Mesh(unit_cube.GenerateMesh(maxh=2.0)),
        "quad": MakeStructured2DMesh(quads=True, nx=1, ny=1),
        "hex": MakeStructured3DMesh(hexes=True, nx=1, ny=1, nz=1),
    }
    expected = {
        "trig": {False: [6, 12, 20], True: [8, 15, 24]},
        "tet": {False: [12, 30, 60], True: [15, 36, 70]},
        "quad": {False: [12, 24, 40], True: [12, 24, 40]},
        "hex": {False: [36, 108, 240], True: [36, 108, 240]},
    }

    for name, mesh in meshes.items():
        element = next(iter(mesh.Elements(ng.VOL)))
        for rt in (False, True):
            local_dofs = [
                len(ng.HDiv(mesh, order=p, RT=rt).GetDofNrs(element))
                for p in (1, 2, 3)
            ]
            assert local_dofs == expected[name][rt]

        bare_dofs = [
            len(ng.HDiv(mesh, order=p).GetDofNrs(element))
            for p in (1, 2, 3)
        ]
        assert bare_dofs == expected[name][False]
