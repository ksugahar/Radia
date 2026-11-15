# -*- coding: utf-8 -*-
r"""H1Henrotte on a mesh refined in place (mesh.Refine()).

Guards: after mesh.Refine() NGSolve keeps the coarse-level edges numbered,
but no fine element uses them.  The order-2 space reserves a DOF per edge, so
each such edge is a free DOF with a zero stiffness row unless the space marks
it UNUSED; the direct solve then fails as singular.  Every free DOF must carry
a nonzero stiffness diagonal at both orders, and the number of free DOFs must
equal the count of DOFs referenced by volume elements minus the Dirichlet ones.
Assumes: triangle mesh of an off-axis rectangle, no Dirichlet boundary,
refined once.
"""
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

import numpy as np
import ngsolve as ng
import pytest

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")
from netgen.geom2d import SplineGeometry
from radia.axifem import AxiHenrotteStiffnessBFI, H1Henrotte


def _refined_mesh():
    geo = SplineGeometry()
    geo.AddRectangle((1.0, -0.5), (2.0, 0.5), bc="outer")
    mesh = ng.Mesh(geo.GenerateMesh(maxh=0.3))
    mesh.Refine()
    return mesh


@pytest.mark.parametrize("order", [1, 2])
def test_refined_mesh_has_no_dead_free_dofs(order):
    mesh = _refined_mesh()
    fes = H1Henrotte(mesh, order=order)
    a = ng.BilinearForm(fes)
    a += AxiHenrotteStiffnessBFI(ng.CF(1.0))
    a.Assemble()
    free = np.array([fes.FreeDofs()[i] for i in range(fes.ndof)])
    used = np.zeros(fes.ndof, dtype=bool)
    for el in fes.Elements(ng.VOL):
        used[list(el.dofs)] = True
    diag = np.array([a.mat[i, i] for i in range(fes.ndof)])
    dead = np.where(free & (np.abs(diag) == 0.0))[0]
    assert dead.size == 0, f"{dead.size} free DOFs with a zero stiffness row (first {dead[:5].tolist()})"
    assert int(free.sum()) == int(used.sum()), (
        f"free DOFs {int(free.sum())} vs DOFs used by elements {int(used.sum())}")
