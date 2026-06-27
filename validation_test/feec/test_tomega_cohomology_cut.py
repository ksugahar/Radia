"""Golden: the T-Omega cohomology CUT solve (gmsh-free) reproduces the wire field.

Exercises CohomologyCutSolver.solve() -- not just setup_from_mesh -- the
PRIMARY use of radia.cohomology (the T-Omega total-scalar-potential cut, of
which the Clebsch current-linking demo is one instance).  A straight wire
(current I along z) through an annular air region (b1 = 1):
  * Ampere's law  oint H.dl = NI  (exact, by the unit-circulation cut);
  * the wire field  H_phi = I / (2 pi r)  to FEM accuracy.
"""
import math

import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")


def _solve_tomega_wire(I=100.0, R_out=0.05, r_w=0.005, L=0.10, maxh=0.008,
                       order=2):
    """Solve the straight-wire T-Omega verification problem."""
    import ngsolve as ng
    from netgen.occ import Cylinder, Dir, OCCGeometry, Pnt
    from radia.cohomology import circulation
    from radia.cohomology_cut import CohomologyCutSolver

    ax = Dir(0, 0, 1)
    outer = Cylinder(Pnt(0, 0, -L / 2), ax, r=R_out, h=L)
    inner = Cylinder(Pnt(0, 0, -L / 2), ax, r=r_w, h=L)
    air = outer - inner
    air.mat("air")
    for face in air.faces:
        c = face.center
        rr = math.hypot(c.x, c.y)
        face.name = "wire" if abs(rr - r_w) < 0.3 * r_w else "outer"

    with ng.TaskManager():
        mesh = ng.Mesh(OCCGeometry(air).GenerateMesh(maxh=maxh))
        solver = CohomologyCutSolver()
        b1 = solver.setup_from_mesh(mesh)
        solver.solve([I], dirichlet="outer", order=order)
        H = solver.get_H()

        rs = np.linspace(1.6 * r_w, 0.8 * R_out, 8)
        errs = []
        for r in rs:
            Hy = H(mesh(r, 0.0, 0.0))[1]
            ana = I / (2 * math.pi * r)
            errs.append(abs(Hy - ana) / ana)
        wire_field_err = float(np.mean(errs))

        amp = circulation(H, mesh, 0.0, 0.0, 0.5 * (r_w + R_out))

    return {
        "b1": int(b1),
        "wire_field_err": wire_field_err,
        "ampere_circulation": float(amp),
        "ampere_NI": float(I),
        "ampere_rel_err": float(abs(amp - I) / I),
    }


def test_tomega_wire_cohomology_cut():
    r = _solve_tomega_wire(maxh=0.010)             # slightly coarser for CI
    assert r["b1"] == 1, r                         # one current loop (annulus)
    assert r["ampere_rel_err"] < 1e-3, r           # cut carries NI exactly
    assert r["wire_field_err"] < 0.05, r           # H_phi ~ I/(2 pi r)
