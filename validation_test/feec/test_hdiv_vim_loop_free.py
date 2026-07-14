"""Production gate for the implicit loop-free HDiv-VIM formulation.

The production solve never constructs a loop/co-loop basis.  This test uses a
small explicit null vector only as a probe: it verifies that a charge-free
HDiv mode is handled by the mass term in the C++ symmetric CG solve, rather
than requiring a separate loop constraint or deflation step.
"""

import numpy as np
import pytest
import scipy.sparse as sp

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")

import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry, Pnt  # noqa: E402

from radia.vim import ChargeGram  # noqa: E402


def test_charge_free_mode_is_solved_without_loop_basis():
    mesh = ng.Mesh(OCCGeometry(Box(Pnt(0, 0, 0), Pnt(1, 1, 1))).GenerateMesh(maxh=0.75))
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        B, gram, mass = ChargeGram(fes, eps=1e-12, ho_far_factor=float("inf"), leafsize=16)

        dense_b = B.toarray()
        _, singular, vtrans = np.linalg.svd(dense_b, full_matrices=True)
        rank = int(np.sum(singular > 1e-10 * singular[0])) if singular.size else 0
        if rank >= fes.ndof:
            pytest.skip("this small mesh has no charge-free HDiv mode")
        loop = np.asarray(vtrans[rank], dtype=float)
        loop /= np.linalg.norm(loop)
        assert np.linalg.norm(B @ loop) < 1e-10

        mcoo = sp.coo_matrix(mass)
        inv_chi = 1.0 / 99.0
        rhs = inv_chi * np.asarray(mass @ loop).ravel()
        result = gram.solve_linear_material_mass_riesz(
            B.indptr.tolist(), B.indices.tolist(), B.data.tolist(), int(fes.ndof),
            mcoo.row.tolist(), mcoo.col.tolist(), mcoo.data.tolist(),
            inv_chi, rhs.tolist(), 1e-10, 200, True)

    assert int(result["iters"]) < 200
    solved = np.asarray(result["m"], dtype=float)
    assert np.linalg.norm(solved - loop) / np.linalg.norm(loop) < 1e-8
