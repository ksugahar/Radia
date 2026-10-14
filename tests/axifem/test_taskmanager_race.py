# -*- coding: utf-8 -*-
r"""Regression: axifem custom assembly is RACE-FREE under NGSolve TaskManager.

The custom FESpace (H1Henrotte -- GetFE placement-news each element on the
per-thread Allocator) and the custom BFI (AxiHenrotteStiffnessBFI::CalcElementMatrix
-- LocalHeap temporaries, const method, only static-const Gauss data) must assemble
a matrix that is INDEPENDENT of the TaskManager thread count.  A data race (e.g. a
future edit that introduces shared mutable state, drops the LocalHeap, or returns a
shared FE) would corrupt the matrix and the checksum would diverge between 1 and 4
threads.

Checksum = (A x).(A x) for a fixed deterministic x; compared at SetNumThreads(1)
vs (4).  Empirically the matrices are BITWISE identical (rel-diff 0.0, NGSolve's
deterministic partitioned assembly); this asserts rel-diff < 1e-12 so it still
tolerates any benign FP-reduction reordering while catching the O(1) corruption a
race produces.

Covers BOTH axifem assembly paths:
  * order-1 V-DOF custom BFI   (AxiHenrotteStiffnessBFI on Q1/P1)
  * order-2 symbolic weak form (custom DiffOp / GetFE via NGSolve's integrator)

Requires the rebuilt axifem extension (Build.ps1 -AxiFemOnly); loads src/radia.
"""
import math
import os
import sys

_SRC = os.path.abspath(os.path.join(os.path.dirname(__file__), "..", "..", "src"))
if _SRC not in sys.path:
    sys.path.insert(0, _SRC)

import numpy as np
import ngsolve
from ngsolve import (BilinearForm, CoefficientFunction as CF, TaskManager, grad,
                     dx, x as r_cf)
from ngsolve.meshes import MakeStructured2DMesh
from radia.axifem import H1Henrotte, AxiHenrotteStiffnessBFI

MU0 = 4e-7 * math.pi
R, Z0, Z1 = 1.0, -0.5, 0.5
NX = NY = 40                      # ~1600 quads / ~3200 tris: parallel work, still fast
THREADS = (1, 4)


def _mesh(quads):
    return MakeStructured2DMesh(quads=quads, nx=NX, ny=NY,
                                mapping=lambda x, y: (R * x, Z0 + (Z1 - Z0) * y))


def _checksum(mat, n):
    """Deterministic (A x).(A x); identical across thread counts iff race-free."""
    x = mat.CreateColVector()
    y = mat.CreateColVector()
    x.FV().NumPy()[:] = 1.0 / (1.0 + np.arange(n))      # fixed vector, no RNG
    y.data = mat * x
    yv = y.FV().NumPy()
    return float(np.dot(yv, yv))


def _assemble_checksum(make_bf, mesh, nthreads):
    saved = ngsolve.ngsglobals.numthreads
    try:
        ngsolve.SetNumThreads(nthreads)
        fes, aK = make_bf(mesh)
        with TaskManager():
            aK.Assemble()
        return _checksum(aK.mat, fes.ndof)
    finally:
        ngsolve.SetNumThreads(saved)


def _bf_vdof(mesh):
    fes = H1Henrotte(mesh, order=1)
    aK = BilinearForm(fes, symmetric=True)
    aK += AxiHenrotteStiffnessBFI(CF(MU0))
    return fes, aK


def _bf_symbolic(mesh):
    fes = H1Henrotte(mesh, order=2)
    u, v = fes.TnT()
    r = r_cf
    aK = BilinearForm(fes, symmetric=True)
    aK += (1.0 / MU0) * (1.0 / r) * (r * grad(u)[0] + u) * (r * grad(v)[0] + v) * dx
    aK += (1.0 / MU0) * r * grad(u)[1] * grad(v)[1] * dx
    return fes, aK


def _reldiff(a, b):
    return abs(a - b) / max(abs(a), abs(b), 1e-300)


def test_vdof_bfi_assembly_race_free():
    """AxiHenrotteStiffnessBFi (order-1 V-DOF) assembles identically at 1 vs 4 threads."""
    mesh = _mesh(quads=True)                       # Q1 axis-aligned V-DOF
    c1 = _assemble_checksum(_bf_vdof, mesh, 1)
    c4 = _assemble_checksum(_bf_vdof, mesh, 4)
    rel = _reldiff(c1, c4)
    print(f"  V-DOF BFI : chk(1t)={c1:.12e}  chk(4t)={c4:.12e}  rel={rel:.2e}")
    assert rel < 1e-12, (
        f"AxiHenrotteStiffnessBFI assembly depends on thread count "
        f"(rel {rel:.2e}) -- data race in the custom BFI/FESpace")


def test_symbolic_assembly_race_free():
    """Order-2 symbolic weak form (custom DiffOp/GetFE) assembles identically 1 vs 4."""
    mesh = _mesh(quads=False)                      # triangles -> P2 symbolic
    c1 = _assemble_checksum(_bf_symbolic, mesh, 1)
    c4 = _assemble_checksum(_bf_symbolic, mesh, 4)
    rel = _reldiff(c1, c4)
    print(f"  symbolic P2: chk(1t)={c1:.12e}  chk(4t)={c4:.12e}  rel={rel:.2e}")
    assert rel < 1e-12, (
        f"H1Henrotte order-2 symbolic assembly depends on thread count "
        f"(rel {rel:.2e}) -- data race in the custom FESpace/DiffOp")


if __name__ == "__main__":
    test_vdof_bfi_assembly_race_free()
    test_symbolic_assembly_race_free()
    print("\nTaskManager race-freedom tests passed.")
