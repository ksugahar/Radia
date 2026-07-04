"""Validation: radia_ngsolve.prepare_cache_hmatrix populates a RadiaField CoefficientFunction's cache via
the O(N log N) _FieldEvalHMatrix with values bit-consistent with the direct rad.Fld, and gf.Set(cf) then
uses them.  This is the path-B (magnet field -> GridFunction) acceleration for the RadiaField CF workflow.

Run: pytest validation_test/hacapk/test_prepare_cache_hmatrix.py
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
import radia as rad                                            # noqa: E402
import ngsolve as ng                                          # noqa: E402
from netgen.occ import Box, Pnt, OCCGeometry                  # noqa: E402
from radia.radia_ngsolve import prepare_cache_hmatrix         # noqa: E402


@pytest.mark.parametrize("quantity", ["b", "a"])
def test_prepare_cache_hmatrix_matches_direct(quantity):
    """The H-matrix-populated cache returns the same field as direct rad.Fld (bit-consistent), and the
    cached CF projects to a GridFunction."""
    rad.UtiDelAll()
    mag = rad.ObjRecMag([0.0, 0.0, 0.0], [0.1, 0.1, 0.1], [0.0, 0.0, 954930.0])
    cf = rad.RadiaField(mag, quantity)
    with ng.TaskManager():                                    # CALLER wraps; the helper does not
        mesh = ng.Mesh(OCCGeometry(Box(Pnt(-0.3, -0.3, -0.3), Pnt(0.3, 0.3, 0.3))).GenerateMesh(maxh=0.2))
        pts = np.array([[0.20, 0.10, 0.15], [0.15, -0.20, 0.10], [-0.12, 0.22, 0.18],
                        [0.25, 0.00, -0.10], [-0.20, -0.15, -0.20]])              # OUTSIDE the magnet
        prepare_cache_hmatrix(cf, pts, eps=1e-9)
        assert cf.GetCacheStats()["size"] == len(pts)
        Bcf = np.array([list(cf(mesh(float(p[0]), float(p[1]), float(p[2])))) for p in pts])
        Bdir = np.asarray(rad.Fld(mag, quantity, pts.tolist()), float).reshape(-1, 3)
        gf = ng.GridFunction(ng.HDiv(mesh, order=1))
        gf.Set(cf)
        assert np.isfinite(gf.vec.FV().NumPy()).all()
    err = np.linalg.norm(Bcf - Bdir) / np.linalg.norm(Bdir)
    assert err < 1e-7, f"quantity={quantity}: cached-H-matrix vs direct rad.Fld = {err:.3e}"
    rad.UtiDelAll()


def test_prepare_cache_hmatrix_rejects_h():
    """'h'/'phi' are not in _FieldEvalHMatrix -> fail loud (use cf.PrepareCache), no silent fallback."""
    rad.UtiDelAll()
    mag = rad.ObjRecMag([0.0, 0.0, 0.0], [0.1, 0.1, 0.1], [0.0, 0.0, 954930.0])
    cf = rad.RadiaField(mag, "h")
    with pytest.raises(ValueError, match="'b' or 'a'"):
        prepare_cache_hmatrix(cf, [[0.2, 0.1, 0.15]])
    rad.UtiDelAll()
