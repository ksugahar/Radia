"""HEX BDM1 charge-Gram near family: physical spectrum band on elongated distorted cells, exact-anchor
inner for touching pairs, and dense H-matrix leaves for touching hosts.

Background (2026-09-05, ESRF #6 quadrupole): the production chi0-warmstart CG broke at iteration 86
because the HEX charge Gram was indefinite, ``lambda_min(M^-1 N)`` in (-5e-3, -2e-3) against the physical
band [0, 1].  Three near-family defects were measured on that mesh (trilinear-distorted cells, no Q2
curvature):

* touching distorted hosts (centroid ratio 0.56..1.0) were classified "not near" at ``near_grade`` 0.5,
  so their outer clouds were the regular rule against a boundary-singular potential;
* the static-site radial inner of touching pairs underestimated their energies by ~1e-3;
* the graded outer rule shared ``glout_n`` = 4 with the far tensor product, so it could not be raised
  without multiplying the far cost by ``(n/4)^6``.

The same near family overestimates on elongated sector cells: the dense generalized spectrum reached
``lambda_max`` 1.084 (21 modes above the physical bound 1) on the 72-cell sector lattice below with the
legacy rules.  The fixed family (touching classification, exact-anchor radial inner, decoupled near outer
rule, dense H-matrix leaves for touching hosts) keeps it PSD and reduces the error against a fine
reference from -11 %..+18 % to -1 %..+9 % (M-metric); the physical upper bound is still exceeded by
3 % there and is locked as a strict xfail until the near family converges.  The heavy #6 reproduction
lives in ``validation_test/esrf_three_engine/validate_hex_gram_definiteness.py``.
"""

import numpy as np
import pytest

pytest.importorskip("ngsolve")

import ngsolve as ng  # noqa: E402
from ngsolve.meshes import MakeStructured3DMesh  # noqa: E402

from radia.vim import _vim as V  # noqa: E402


def _sector_mesh():
    """Quarter annular sector r in [0.5, 1], theta in [0, pi/2], extruded 0.6: cells with aspect ratio up
    to ~3.5 and tapered (non-affine) trilinear geometry -- the class that breaks the legacy near family."""
    return MakeStructured3DMesh(
        hexes=True, nx=6, ny=4, nz=3,
        mapping=lambda x, y, z: ((0.5 + 0.5 * x) * np.cos(0.5 * np.pi * y),
                                 (0.5 + 0.5 * x) * np.sin(0.5 * np.pi * y), 0.6 * z))


def _generalized_spectrum(mesh, **kwargs):
    sla = pytest.importorskip("scipy.linalg")
    sp = pytest.importorskip("scipy.sparse")
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        B, G, M, _ = V._build_charge_gram_hex(fes, eps=1e-10, **kwargs)
        Bs = sp.csr_matrix(B)
        n = Bs.shape[1]
        N = np.empty((n, n))
        for column in range(n):
            basis = np.zeros(n)
            basis[column] = 1.0
            charge = np.ascontiguousarray(Bs @ basis, dtype=np.float64)
            N[:, column] = Bs.T @ np.asarray(G.matvec_sym(charge))
        N = 0.5 * (N + N.T)
        Md = sp.csr_matrix(M).toarray()
        stats = dict(G.stats())
    return sla.eigvalsh(N, Md), stats


def test_sector_spectrum_is_psd_and_dispatch_is_the_near_family():
    """Elongated distorted cells: the generalized demag spectrum stays PSD, the dispatch reports the
    production near family, and every distorted host pair inside the near band is graded."""
    w, stats = _generalized_spectrum(_sector_mesh())
    assert w[0] > -1e-8, f"demag spectrum lost PSD: min eig {w[0]:.3e}"
    assert stats["hex_near_inner_exact"] == 1.0
    assert stats["hex_cluster_radius_enabled"] == 1.0
    assert stats["hex_glnear_n"] == 8.0 and stats["hex_glout_n"] == 4.0
    assert stats["hex_glin_self_n"] == 12.0 and stats["hex_glin_n"] == 5.0
    assert stats["nonproduction_numerical_override_active"] == 0.0
    # every distorted host pair inside the near band is graded; none took the plain far rule
    assert stats["hex_blk_general_near"] > 0 and stats["hex_blk_general_far"] == 0


@pytest.mark.xfail(strict=True, reason=(
    "open accuracy item (2026-09-05): on this elongated sector lattice the near family is still "
    "+/- 2..9 % off a fine reference in the M-metric (legacy: -11 %..+18 %), so lambda_max reaches "
    "1.03 against the physical bound 1; a converged near family must turn this into a pass"))
def test_sector_spectrum_stays_in_physical_band():
    w, _ = _generalized_spectrum(_sector_mesh())
    assert w[-1] < 1.0 + 1e-3, f"demag spectrum exceeds the physical band: {w[-1]:.6f}"


def test_legacy_site_inner_is_a_flagged_override_and_differs():
    """The static-site radial survives only as a diagnostic A/B path: it is reported as a numerical
    override and its touching-pair entries differ from the exact-anchor family."""
    mesh = _sector_mesh()
    with ng.TaskManager():
        fes = ng.HDiv(mesh, order=1)
        _, G_exact, _, _ = V._build_charge_gram_hex(fes, eps=1e-10, build_hmatrix=False)
        _, G_site, _, _ = V._build_charge_gram_hex(fes, eps=1e-10, build_hmatrix=False, near_inner="site")
        cb = V._charge_basis_hex(fes, cob_quad=2, materialize_mass=False)
    host = np.asarray(cb["host"], int)
    kind = np.asarray(cb["kind"], int)
    stats_exact = dict(G_exact.stats())
    stats_site = dict(G_site.stats())
    assert stats_exact["nonproduction_numerical_override_active"] == 0.0
    assert stats_site["nonproduction_numerical_override_active"] == 1.0
    assert stats_site["hex_near_inner_exact"] == 0.0
    # a cell charge against the charges of a face-sharing neighbour cell: touching, non-self
    cent = np.zeros((mesh.ne, 3))
    for el in mesh.Elements(ng.VOL):
        cent[el.nr] = np.mean([mesh.vertices[v.nr].point for v in el.vertices], axis=0)
    a = int(np.flatnonzero(kind == 0)[0])
    e = int(host[a])
    dist = np.linalg.norm(cent - cent[e], axis=1)
    dist[e] = np.inf
    f = int(np.argmin(dist))
    b = int(np.flatnonzero((kind == 0) & (host == f))[0])
    exact = G_exact.entry(a, b)
    site = G_site.entry(a, b)
    self_scale = np.sqrt(G_exact.entry(a, a) * G_exact.entry(b, b))
    rel = abs(exact - site) / self_scale
    assert 1e-7 < rel < 1e-2, f"exact-anchor and site inners differ by {rel:.3e} of the self scale"
    with pytest.raises(ValueError):
        with ng.TaskManager():
            V._build_charge_gram_hex(fes, eps=1e-10, build_hmatrix=False, near_inner="newton")
