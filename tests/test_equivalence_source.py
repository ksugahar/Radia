"""Tests for radia.equivalence_source.NearFieldSource (NFS extraction).

Focused on the No-Fallback contract of NearFieldSource.extract_ngsolve: an
unimplemented quadrature_order > 1 (Gauss quadrature) must RAISE rather than
silently sample at order 1 and only warn.
"""
import pytest


def test_extract_ngsolve_quadrature_order_raises():
    """No-Fallback: quadrature_order > 1 (Gauss quadrature) is not implemented,
    so extract_ngsolve RAISES NotImplementedError instead of silently sampling
    at order 1 with a warning (the caller asked for higher accuracy and would
    otherwise get a different, lower-order computation).  The guard fires before
    any mesh work, so a dummy non-None gf_H is enough to exercise it."""
    pytest.importorskip("ngsolve")
    from radia.equivalence_source import NearFieldSource

    with pytest.raises(NotImplementedError, match="quadrature_order"):
        NearFieldSource.extract_ngsolve(
            mesh=None, gf_H=object(), quadrature_order=2)


def test_extract_ngsolve_requires_gf_H():
    """gf_H is mandatory -- a missing magnetic field raises ValueError, not a
    silent zero source."""
    pytest.importorskip("ngsolve")
    from radia.equivalence_source import NearFieldSource

    with pytest.raises(ValueError, match="gf_H"):
        NearFieldSource.extract_ngsolve(mesh=None, gf_H=None)


if __name__ == "__main__":
    pytest.main([__file__, "-v"])
