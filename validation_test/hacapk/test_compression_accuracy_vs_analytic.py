"""Validation: RT1 ChargeGram ACA accuracy tracks the requested tolerance.

The public ``ChargeGram(HDiv(mesh, order=1))`` result exposes both the
uncompressed quadrature entry oracle and the compressed symmetric H-matrix
apply.  Materializing the former on a small mesh isolates ACA compression
error without retaining the retired RT0 assembly or private constructor.

  1. ACCURACY TRACKS eps -- the compression error decreases with tighter eps (err(eps) <~ 100*eps), so the
     ACA tolerance is a meaningful accuracy knob, not a nominal one.
  2. COMPRESSION IS REAL -- at a loose eps the error is well above machine zero (an UNcompressed / dense
     apply would be exact regardless of eps), proving ACA low-rank compression is actually happening.

Reference numbers (sphere, ~560 charges): err = 6e-4 / 8e-6 / 7e-8 / 1e-9 at eps = 1e-2 / 1e-4 / 1e-6 / 1e-8.
Run explicitly (validation lane): `pytest validation_test/hacapk/`.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.csg")
import ngsolve as ng                                       # noqa: E402
from netgen.csg import CSGeometry, Sphere, Pnt             # noqa: E402
from radia.vim import ChargeGram                           # noqa: E402

_EPS_SWEEP = (1e-2, 1e-4, 1e-6, 1e-8)


def _rt1_fes(maxh=0.35):
    """RT1 tetrahedral sphere with enough charge modes for ACA blocks."""
    g = CSGeometry(); g.Add(Sphere(Pnt(0, 0, 0), 1.0))
    with ng.TaskManager():
        mesh = ng.Mesh(g.GenerateMesh(maxh=maxh))
        return ng.HDiv(mesh, order=1)


def _compression_errors():
    """rel matvec error vs the exact analytic Gram, per ACA eps (Dexact from .entry, built once)."""
    fes = _rt1_fes()
    Dexact = n = X = Yref = None
    errs = {}
    for eps in _EPS_SWEEP:
        with ng.TaskManager():
            _, G, _ = ChargeGram(
                fes, eps=eps, leafsize=32, eta=2.0,
                ho_far_factor=float("inf"),
            )
            if n is None:
                n = G.ndof()
                Dexact = np.array([[G.entry(i, j) for j in range(n)] for i in range(n)])  # EXACT (eps-indep)
                X = np.random.default_rng(0).standard_normal((n, 5))
                Yref = Dexact @ X
            Yh = np.column_stack([np.asarray(G.matvec(X[:, k].tolist()), float) for k in range(X.shape[1])])
        errs[eps] = float(np.linalg.norm(Yh - Yref) / np.linalg.norm(Yref))
    return n, errs


@pytest.fixture(scope="module")
def compression():
    n, errs = _compression_errors()
    assert n > 200, f"need enough charges for far-field compression (got {n})"
    return errs


def test_compression_accuracy_tracks_aca_tolerance(compression):
    """err(eps) <~ 100*eps: the ACA tolerance is a real accuracy knob for the H-matrix vs the exact Gram."""
    for eps, err in compression.items():
        assert err < 100.0 * eps, f"eps={eps:.0e}: compression error {err:.3e} exceeds 100*eps"


def test_compression_accuracy_improves_monotonically(compression):
    """Tighter eps -> smaller error (monotone, allowing a small 3x slack for ACA rank quantisation)."""
    errs = [compression[e] for e in _EPS_SWEEP]           # eps decreasing 1e-2 -> 1e-8
    for a, b, ea, eb in zip(errs, errs[1:], _EPS_SWEEP, _EPS_SWEEP[1:]):
        assert b <= a * 3.0, f"error did not decrease from eps={ea:.0e} ({a:.3e}) to eps={eb:.0e} ({b:.3e})"
    assert errs[-1] < errs[0] * 1e-2, f"tightening eps 1e-2->1e-8 barely helped ({errs[0]:.3e}->{errs[-1]:.3e})"


def test_compression_is_real_not_dense(compression):
    """At a LOOSE eps the error is well above machine zero -- an uncompressed (dense) apply would be exact
    regardless of eps, so a non-trivial eps-dependent error proves ACA low-rank compression is happening."""
    assert compression[1e-2] > 1e-5, \
        f"loose-eps error {compression[1e-2]:.3e} is ~machine-zero -> no compression is occurring (dense?)"
    assert compression[1e-8] < 1e-6, f"tight-eps error {compression[1e-8]:.3e} should be near-exact"
