"""Golden test: the production charge-Gram H-matrix is ACCURATE on a NEAR-HEAVY body at the DEFAULT
admissibility (eta=2.0), and the production charge set has NO co-located duplicate charges.

Background (2026-06-27 investigation): an isolated development probe reported a 2-4% ACA
matvec error on a near-heavy closed hex O-ring at eta=2.0, attributed to "eta too lenient".  Re-measured on
the PRODUCTION tet/RT0 path (radia.vim._core.build_demag), the error is ~1e-5 (3-4 orders below the ~1e-4
end-to-end physics error) and END-TO-END eta-INSENSITIVE -- so the default eta=2.0 / leaf=32 was NOT
changed.  The probe's 4% was an ARTIFACT of CO-LOCATED DUPLICATE charges: it emitted all 6 faces of every
hex, so interior shared faces appeared twice at distance 0 (a near-singular block ACA cannot compress).
Production build_demag emits VOLUME charges (cell rho) + BOUNDARY faces only (M.n on ds) -> no duplicates.

These two locks make the near-heavy case no longer golden-invisible (the prior sphere/ellipsoid goldens are
well-separated and miss it) and structurally guard the artifact class:
  (1) the production charge set on a near-heavy annulus has NO duplicate (co-located) charges;
  (2) the production charge-Gram H-matvec at the DEFAULT eta=2.0 reproduces the EXACT analytic Gram
      (entry(), ACA-independent) to well under 1e-3 -- if a future change reintroduced duplicate charges
      (or otherwise corrupted admissibility), this collapses to the percent level and the test fails.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
import radia._radia_pybind as _rp  # noqa: E402
from radia.vim import _core as tet  # noqa: E402
import ngsolve as ng  # noqa: E402
from netgen.csg import CSGeometry, OrthoBrick, Pnt  # noqa: E402


def _annulus(h=0.024, outer=0.06, inner=0.035, zt=0.02):
    """Closed square washer (a flux loop, Betti_1 = 1) -- the canonical near-heavy soft-iron body."""
    o = OrthoBrick(Pnt(-outer, -outer, -zt), Pnt(outer, outer, zt))
    i = OrthoBrick(Pnt(-inner, -inner, -2 * zt), Pnt(inner, inner, 2 * zt))
    g = CSGeometry(); g.Add(o - i)
    with ng.TaskManager():
        return ng.Mesh(g.GenerateMesh(maxh=h))


def test_production_charge_set_has_no_duplicate_charges():
    """The production charge set (cell rho + BOUNDARY-face sigma) must have NO co-located charges.

    Co-located (distance-0) charges are the root cause of the spurious ~4% ACA error: they form
    near-singular blocks the ACA cannot compress.  build_demag must never emit them (it uses ds boundary
    faces, not every face of every cell).  A regression that double-counted interior faces would fail here
    BEFORE it could silently degrade a near-heavy solve."""
    with ng.TaskManager():
        d = tet.build_demag(_annulus())
    cent = np.asarray(d["cent"], float)
    n_charge = int(d["n_charge"])
    assert cent.shape == (n_charge, 3)
    uniq = len({tuple(np.round(c, 9)) for c in cent})
    assert uniq == n_charge, (
        f"production charge set has {n_charge - uniq} co-located duplicate charge(s) "
        f"({uniq} unique of {n_charge}); duplicates are the ACA-error artifact class (mg_aca_tune)")


def test_nearheavy_chargegram_aca_accurate_at_default_eta():
    """On a near-heavy closed annulus, the charge-Gram H-matvec at the DEFAULT eta=2.0 / leaf=32 matches
    the EXACT analytic Gram (entry(), ACA-independent) to well under 1e-3.  near_factor=1e30 so entry() is
    the exact analytic Gram and the H-matvec/entry gap is PURE ACA compression -- this isolates the
    admissibility (eta) accuracy.  Measured ~1e-5; the 1e-3 hard band catches a regression to the percent
    level (e.g. reintroduced duplicate charges) while staying robust to mesh/quadrature noise."""
    with ng.TaskManager():
        d = tet.build_demag(_annulus())
    nch = int(d["n_charge"])
    cv, fv, ne = list(d["cell_verts"]), list(d["face_verts"]), int(d["n_el"])

    # exact analytic Gram via the entry() oracle (build=False, near_factor=1e30 -> ACA-independent)
    Gref = _rp._ChargeGramHMatrix(cell_verts=cv, face_verts=fv, n_el=ne, near_factor=1e30, build=False)
    G = np.empty((nch, nch))
    for i in range(nch):
        for j in range(i, nch):
            v = Gref.entry(i, j); G[i, j] = v; G[j, i] = v

    # production-default H-matrix (eta=2.0, leaf=32); near_factor=1e30 isolates the ACA from the far quad
    H = _rp._ChargeGramHMatrix(cell_verts=cv, face_verts=fv, n_el=ne,
                               eps=1e-9, leaf=32, eta=2.0, near_factor=1e30)
    rng = np.random.default_rng(0)
    X = rng.standard_normal((nch, 8))
    ref = G @ X
    cols = np.column_stack([np.asarray(H.matvec(X[:, k].tolist()), float) for k in range(X.shape[1])])
    rel_max = float((np.linalg.norm(cols - ref, axis=0) / np.maximum(np.linalg.norm(ref, axis=0), 1e-30)).max())
    assert rel_max < 1e-3, (
        f"near-heavy charge-Gram ACA matvec relerr {rel_max:.2e} at default eta=2.0 exceeds 1e-3 "
        f"(measured ~1e-5); a percent-level value indicates duplicate charges or corrupted admissibility")
