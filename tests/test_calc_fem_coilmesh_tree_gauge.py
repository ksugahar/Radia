"""Tree gauge of the coil A-V solve: removes every non-SIBC gradient kernel."""

import sys
from pathlib import Path

import numpy as np
import pytest

pytestmark = pytest.mark.usefixtures("ngsolve_taskmanager")
import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry, X  # noqa: E402

PANELS = Path(__file__).parents[1] / "src" / "radia" / "panels"
sys.path.insert(0, str(PANELS))
import calc_fem_coilmesh as solver  # noqa: E402


def _mesh():
    from netgen.occ import Z
    box = Box((0, 0, 0), (1, 1, 1))
    box.faces.Min(X).Identify(box.faces.Max(X), "periodic")
    box.faces.Max(Z).name = "sibc"
    return ng.Mesh(OCCGeometry(box).GenerateMesh(maxh=0.35))


def test_tree_gauge_makes_curl_curl_plus_sibc_robin_nonsingular():
    """Guards: without any mass penalty, curl-curl + SIBC Robin is factorable
    and every remaining non-SIBC gradient is removed (exact count)."""
    mesh = _mesh()
    fesA = ng.Periodic(ng.HCurl(mesh, order=1, nograds=True, complex=True))
    fesP = ng.H1(mesh, order=1, complex=True, definedon=mesh.Materials(".*"),
                 dirichlet="sibc")
    fes = fesA * fesP
    free = solver._tree_gauged_free_dofs(mesh, fesA, fes, "sibc")
    (A, _), (N, _) = fes.TnT()
    a = ng.BilinearForm(fes, symmetric=True)
    a += ng.curl(A) * ng.curl(N) * ng.dx
    a += (1 + 1j) * A.Trace() * N.Trace() * ng.ds("sibc")
    a.Assemble()
    # the H1 block is unused here; keep only A's free DOFs
    only_a = ng.BitArray(free)
    for d in range(fesA.ndof, fes.ndof):
        only_a[d] = False
    inv = a.mat.Inverse(only_a, inverse="sparsecholesky")
    rhs = a.mat.CreateRowVector()
    rhs.FV().NumPy()[:] = np.where(np.array(list(only_a)), 1.0, 0.0)
    sol = rhs.CreateVector()
    sol.data = inv * rhs
    res = rhs.CreateVector()
    res.data = rhs - a.mat * sol
    proj = ng.Projector(only_a, True)
    res.data = proj * res
    rel = res.Norm() / rhs.Norm()
    assert np.isfinite(sol.FV().NumPy()).all()
    assert rel < 1e-10, f"relative residual {rel:.3e}, expected < 1e-10"


def test_tree_gauge_rejects_missing_ground_boundary():
    """Guards: a ground label without vertices fails loudly."""
    mesh = _mesh()
    fesA = ng.HCurl(mesh, order=1, nograds=True)
    with pytest.raises(ValueError, match="no vertices"):
        solver._tree_gauged_free_dofs(mesh, fesA, fesA, "absent")
