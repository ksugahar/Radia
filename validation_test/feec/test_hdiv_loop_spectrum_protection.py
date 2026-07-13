#!/usr/bin/env python
"""Loop-eigenvalue protection of the HDiv-VIM charge-Gram operator.

Structural fact this golden locks: the loop space K_h = ker(B) of the discrete
charge map is protected from H-matrix compression error BY CONSTRUCTION.  Any
ACA error E enters on the charge side only (the operator is applied as
N v = B^T G_H (B v) with G_H = G + E), so for a loop basis Z with B Z = 0

    Z^T B^T (G + E) B Z = (B Z)^T (G + E) (B Z) = 0        exactly.

Consequently, in the mass generalized eigenproblem of the Hantila operator
A = nu0*M + B^T G_H B, the ENTIRE loop space sits at the single degenerate
eigenvalue nu0 -- at ANY compression tolerance.  A degenerate point is free
for CG, which is why loop deflation gains nothing on this operator (measured:
deflating all loops changes the iteration count by 0; the iteration count is
set by the compact charge-mode spectrum, kappa ~ 5).

This is the quantitative backing for the discussion-section claim of the
2027-01 workshop paper (loops are not merely benign under the SPD mass term:
their eigenvalue is invariant to Gram compression error).

Runtime ~2-4 min (dense generalized eigensolve at 5616 dof, twice) -- the
validation_test lane.  Correctness only; no timing claims (LAB per Benchmark
Policy).  Results saved to loop_spectrum_protection.json.
"""
import json
import platform
import sys
from datetime import datetime
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("ngsolve")

import scipy.linalg as sla                              # noqa: E402
import scipy.sparse as sp                               # noqa: E402
import ngsolve as ng                                    # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh         # noqa: E402

from radia.vim._vim import build_charge_gram            # noqa: E402

HERE = Path(__file__).resolve().parent
NX = 6
X = 1.0 / 200.0                     # mu0 * nu_rev of the goldens' synthetic play (mu_r = 200)
NU0 = X / (1.0 - X)
EPS_TIGHT, EPS_LOOSE = 1e-10, 1e-4


def _densify_N(B, H, n):
    Bc = B.tocsc()
    N = np.zeros((n, n))
    for j in range(n):
        q = np.asarray(Bc[:, j].todense()).ravel()
        N[:, j] = B.T @ np.asarray(H.matvec_sym(q.tolist()), float)
    return 0.5 * (N + N.T)


def test_loop_eigenvalue_protected_from_compression():
    results = dict(
        description="HDiv charge-Gram loop-eigenvalue protection vs ACA compression",
        timestamp=datetime.now().isoformat(), hostname=platform.node(),
        nx=NX, nu0=NU0, eps_cases=[])

    with ng.TaskManager():
        mesh = MakeStructured3DMesh(hexes=True, nx=NX, ny=NX, nz=NX,
                                    mapping=lambda x, y, z: (x - .5, y - .5, z - .5))
        fes = ng.HDiv(mesh, order=1)
        n = fes.ndof

        first = True
        nloop = None
        Z = None
        Mm = None
        for eps in (EPS_TIGHT, EPS_LOOSE):
            B, H, Mm0 = build_charge_gram(fes, eps=eps, leafsize=32, eta=2.0,
                                          far_quad=None, ho_far_factor=None, nonlinear=True)
            B = sp.csr_matrix(B)
            stats = dict(H.stats())
            n_lowrank = int(stats.get("n_lowrank", 0))
            assert n_lowrank > 0, (
                "ACA compression not engaged (n_lowrank=0) -- the eps-invariance "
                "claim would be vacuous at this size; enlarge NX")

            if first:
                Mm = np.asarray(sp.csr_matrix(Mm0).todense())
                Bd = np.asarray(B.todense())
                s = np.linalg.svd(Bd, compute_uv=False)
                rank = int(np.sum(s > s.max() * max(Bd.shape) * np.finfo(float).eps))
                nloop = n - rank
                _, _, vt = np.linalg.svd(Bd, full_matrices=True)
                Z = vt[rank:].T
                assert nloop > 0.3 * n, "unexpectedly small loop space (%d of %d)" % (nloop, n)
                first = False

            N = _densify_N(B, H, n)
            w = sla.eigh(NU0 * Mm + N, Mm, eigvals_only=True)

            at_point = int(np.sum(np.abs(w - NU0) < NU0 * 1e-8))
            below = int(np.sum(w < NU0 * (1.0 - 1e-8)))
            smear = float(np.max(np.abs(
                sla.eigh(Z.T @ N @ Z, Z.T @ Mm @ Z, eigvals_only=True))))
            charge = w[w > NU0 * (1.0 + 1e-8)]
            kappa = float(w[-1] / charge[0])

            # THE protection facts, per eps (incl. the loose one):
            assert below == 0, "eps=%g: %d eigenvalues BELOW the loop point" % (eps, below)
            assert at_point == nloop, (
                "eps=%g: loop cluster not exactly degenerate at nu0 "
                "(%d of %d modes at the point)" % (eps, at_point, nloop))
            assert smear < NU0 * 1e-9, (
                "eps=%g: compression error leaked into ker(B) "
                "(loop-block smear %.3e vs nu0 %.3e)" % (eps, smear, NU0))
            assert 2.0 < kappa < 20.0, (
                "charge-mode spectrum shape changed structurally (kappa_eff %.2f)" % kappa)

            results["eps_cases"].append(dict(
                eps=eps, n_lowrank=n_lowrank,
                compression=float(stats.get("compression", 0.0)),
                at_point=at_point, below_point=below,
                loop_smear_over_nu0=float(smear / NU0),
                first_charge_eig=float(charge[0]), eig_max=float(w[-1]),
                kappa_eff_charge=kappa))

        results.update(ndof=int(n), nloop=int(nloop))

    (HERE / "loop_spectrum_protection.json").write_text(json.dumps(results, indent=1))


if __name__ == "__main__":
    sys.exit(pytest.main([__file__, "-v", "-s"]))
