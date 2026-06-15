"""Per-region LINEAR soft-iron materials for radia.vim.hdiv_demag_solve (mu_r as a dict
{material_name: mu_r}).  N = B^T G B is geometry-only (material-independent), so per-region soft iron
enters ONLY through the (1/chi)-weighted HDiv mass M_invchi = INT (1/chi(x)) u.v dx.  This is the first
per-region productionization increment (docs/hdiv_vim/PRODUCTIONIZATION.md "per-region / mixed").

Locks: (1) a dict with EQUAL mu in every region reproduces the scalar-mu result bit-for-bit (the
weighted mass reduces to (1/chi) M_mass); (2) DIFFERENT mu per region is physical -- the global M_avg
lies between the two single-mu runs and the high-mu region magnetizes more than the low-mu region;
(3) fail-loud on a missing region or mu_r <= 1 (No-Fallbacks)."""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("netgen.occ")
import ngsolve as ng  # noqa: E402
from netgen.occ import Box, OCCGeometry, Pnt, Glue  # noqa: E402

from radia.vim import hdiv_demag_solve  # noqa: E402

H0 = 1000.0
HEXT = ng.CoefficientFunction((0.0, 0.0, H0))


def _two_region_mesh(maxh=0.25):
    """Unit box split at x=0 into materials 'lo' (x<0) and 'hi' (x>0)."""
    lo = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.0, 0.5, 0.5)).mat("lo")
    hi = Box(Pnt(0.0, -0.5, -0.5), Pnt(0.5, 0.5, 0.5)).mat("hi")
    with ng.TaskManager():
        return ng.Mesh(OCCGeometry(Glue([lo, hi])).GenerateMesh(maxh=maxh))


def _region_mean_absMz(mesh, M_el):
    """Mean |M_z| over the 'lo' and 'hi' element sets (M_el is per-element, mesh element order)."""
    mats = [mesh[ng.ElementId(ng.VOL, i)].mat for i in range(mesh.ne)]
    lo = np.array([abs(M_el[i, 2]) for i in range(mesh.ne) if mats[i] == "lo"])
    hi = np.array([abs(M_el[i, 2]) for i in range(mesh.ne) if mats[i] == "hi"])
    return float(lo.mean()), float(hi.mean())


def test_per_region_equal_mu_matches_scalar():
    """A {lo: mu, hi: mu} dict must reproduce the scalar-mu solve exactly (weighted mass == (1/chi)Mm)."""
    mesh = _two_region_mesh()
    with ng.TaskManager():
        rd = hdiv_demag_solve(mesh, mu_r={"lo": 200.0, "hi": 200.0}, H_ext=HEXT)
        rs = hdiv_demag_solve(mesh, mu_r=200.0, H_ext=HEXT)
    assert np.allclose(rd["M"], rs["M"], rtol=1e-8, atol=1e-5)
    assert np.allclose(rd["M_avg"], rs["M_avg"], rtol=1e-8, atol=1e-5)


def test_per_region_different_mu_physics():
    """Different mu per region: the global M_avg_z lies strictly between the all-low and all-high scalar
    runs (monotone response to per-region mu), the mixed run is genuinely distinct from both, and the
    high-mu region develops a larger |M| than the low-mu region.  (The low region's NET M_z can reverse
    sign -- the strong high-mu region's return flux opposes the applied field there -- so we compare
    |M|, not signed M_z, per region.)"""
    mesh = _two_region_mesh()
    with ng.TaskManager():
        r = hdiv_demag_solve(mesh, mu_r={"lo": 50.0, "hi": 500.0}, H_ext=HEXT)
        r_lo = hdiv_demag_solve(mesh, mu_r=50.0, H_ext=HEXT)
        r_hi = hdiv_demag_solve(mesh, mu_r=500.0, H_ext=HEXT)
    az_lo, az_hi = r_lo["M_avg"][2], r_hi["M_avg"][2]
    az = r["M_avg"][2]
    assert az_lo < az < az_hi                               # monotone: bounded by the two scalar runs
    assert abs(az - az_lo) > 0.01 * abs(az_lo)              # genuinely per-region, not == all-low
    assert abs(az - az_hi) > 0.01 * abs(az_hi)              # ... and not == all-high
    m_lo, m_hi = _region_mean_absMz(mesh, r["M"])
    assert m_hi > m_lo                                      # the mu=500 region magnetizes more strongly


def test_per_region_fail_loud():
    """Missing region or mu_r <= 1 must raise (No-Fallbacks), not silently guess."""
    mesh = _two_region_mesh()
    with pytest.raises(ValueError):
        with ng.TaskManager():
            hdiv_demag_solve(mesh, mu_r={"lo": 100.0}, H_ext=HEXT)          # 'hi' missing
    with pytest.raises(ValueError):
        with ng.TaskManager():
            hdiv_demag_solve(mesh, mu_r={"lo": 100.0, "hi": 1.0}, H_ext=HEXT)   # mu_r <= 1
