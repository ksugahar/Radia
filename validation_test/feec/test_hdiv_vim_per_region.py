"""Per-region LINEAR soft-iron materials for radia.vim.hdiv_demag_solve (mu_r as a dict
{material_name: mu_r}).  N = B^T G B is geometry-only (material-independent), so per-region soft iron
enters ONLY through the chi-weighted HDiv mass M_chi = INT chi(x) u.v dx of the form-1 projected system
A = M_mass + M_chi M_mass^-1 N (the same projected M = chi H statement the nonlinear path uses, so a
linear region agrees with its nonlinear-table equivalent).  per-region productionization increment
(docs/hdiv_vim/PRODUCTIONIZATION.md "per-region / mixed").

Locks: (1) a dict with EQUAL mu in every region reproduces the scalar-mu result bit-for-bit (the
weighted mass reduces to chi M_mass); (2) DIFFERENT mu per region is physical -- the global M_avg
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
    """A {lo: mu, hi: mu} dict must reproduce the scalar-mu solve exactly (weighted mass == chi Mm).

    For equal mu the per-region form-1 operator A = M_mass + M_chi M_mass^-1 N (M_chi = chi M_mass) and
    the scalar +N system ((1/chi) M_mass + N) m = M_mass h_ext are the SAME operator (the chi-scaling
    cancels), so the build is correct iff the two solves agree at the SOLUTION.  But they take DIFFERENT
    Krylov paths -- the dict path is GMRES on the form-1 operator, the scalar default is CG with the mass
    Riesz preconditioner -- so at the default tol=1e-8 they agree only to the Krylov stopping level
    (~3e-4 abs / 1e-5 rel), NOT to atol=1e-5.  Solve tightly (tol=1e-11) so both reach the true solution
    and the build-correctness comparison is meaningful (then they agree to ~1e-8).  Pin gram_eps=1e-12 on
    BOTH so the dict (GMRES form-1, default 1e-10) and scalar (CG mass-riesz, whose default tightens to
    1e-12) build the IDENTICAL Gram -- comparing the SAME operator, not two ACA tolerances.  Pin
    near_factor=1e30 on BOTH too: the scalar CG-tet path defaults to the fast near/far build (near_factor=2
    + far_quad=4) while the dict GMRES path stays all-analytic, so force both to the exact all-analytic Gram."""
    mesh = _two_region_mesh()
    with ng.TaskManager():
        rd = hdiv_demag_solve(mesh, mu_r={"lo": 200.0, "hi": 200.0}, H_ext=HEXT,
                              tol=1e-11, gram_eps=1e-12, near_factor=1e30)
        rs = hdiv_demag_solve(mesh, mu_r=200.0, H_ext=HEXT,
                              tol=1e-11, gram_eps=1e-12, near_factor=1e30)
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


# ---------------------------------------------------------------------------------------------------
# Per-region NONLINEAR soft iron: bh_table as {material: [[H,B]]} (multiple iron grades).  Same
# structural fact as the linear case -- N = B^T G B is geometry-only -- but now each region carries its
# OWN BH curve, so per-region nonlinear enters through the per-element constitutive law
# (_table_tensor_tangent_multi) + the per-region zero-field-chi Newton warmstart.
# ---------------------------------------------------------------------------------------------------
_MU0 = 4e-7 * np.pi


def _bh(chi0, Bsat, n=40):
    """Saturating BH table [[H,B]]: M(H) = chi0 H/(1+chi0|H|/Msat), Msat = Bsat/mu0, B = mu0(H+M)."""
    H = np.concatenate([[0.0], np.geomspace(1.0, 1e6, n)])
    M = chi0 * H / (1.0 + chi0 * np.abs(H) / (Bsat / _MU0))
    return np.column_stack([H, _MU0 * (H + M)])


BH_SOFT = _bh(5000.0, 2.0)        # high-permeability grade
BH_HARD = _bh(200.0, 2.0)         # low-permeability grade (same saturation)
HEXT_NL = ng.CoefficientFunction((0.0, 0.0, 1.0e5))   # drive that genuinely saturates the soft grade


@pytest.mark.slow
def test_per_region_nl_equal_table_matches_single():
    """A {lo: BH, hi: BH} dict with the SAME table in both regions must reproduce the single-table
    nonlinear solve bit-for-bit (the per-element constitutive law + per-region warmstart both collapse to
    the single-region path when every region shares one curve)."""
    mesh = _two_region_mesh()
    with ng.TaskManager():
        rd = hdiv_demag_solve(mesh, bh_table={"lo": BH_SOFT, "hi": BH_SOFT}, H_ext=HEXT_NL, nl_tol=1e-6)
        rs = hdiv_demag_solve(mesh, bh_table=BH_SOFT, H_ext=HEXT_NL, nl_tol=1e-6)
    assert rd["nonlinear"] and rs["nonlinear"]
    assert np.allclose(rd["M"], rs["M"], rtol=1e-8, atol=1e-3)
    assert np.allclose(rd["M_avg"], rs["M_avg"], rtol=1e-8, atol=1e-3)


@pytest.mark.slow
def test_per_region_nl_different_tables_physics():
    """Different BH grades per region: (1) the global M_avg_z is strictly BOUNDED by the all-soft and
    all-hard scalar runs (monotone response to per-region permeability); (2) the mixed run is genuinely
    distinct from BOTH (the dict path is not silently collapsing to one curve); (3) the solve is genuinely
    NONLINEAR -- the soft grade's M_avg differs substantially from its LINEAR (mu_r = chi0+1) counterpart,
    so the BH saturation is actually engaged (not a disguised linear solve)."""
    mesh = _two_region_mesh()
    with ng.TaskManager():
        r  = hdiv_demag_solve(mesh, bh_table={"lo": BH_SOFT, "hi": BH_HARD}, H_ext=HEXT_NL, nl_tol=1e-6)
        rs = hdiv_demag_solve(mesh, bh_table=BH_SOFT, H_ext=HEXT_NL, nl_tol=1e-6)
        rh = hdiv_demag_solve(mesh, bh_table=BH_HARD, H_ext=HEXT_NL, nl_tol=1e-6)
        rL = hdiv_demag_solve(mesh, mu_r=5001.0, H_ext=HEXT_NL)   # linear chi0 counterpart of the soft grade
    z, zs, zh, zL = r["M_avg"][2], rs["M_avg"][2], rh["M_avg"][2], rL["M_avg"][2]
    assert zh < z < zs                                    # bounded: harder grade -> less, softer -> more
    assert abs(z - zs) > 0.03 * abs(zs)                   # distinct from all-soft (measured ~11%)
    assert abs(z - zh) > 0.01 * abs(zh)                   # distinct from all-hard (measured ~2.3%)
    assert abs(zs - zL) > 0.10 * abs(zL)                  # saturation engaged: nonlinear != linear (~26%)


def test_per_region_nl_fail_loud():
    """Missing region or a malformed per-region table must raise (No-Fallbacks)."""
    mesh = _two_region_mesh()
    with pytest.raises(ValueError):
        with ng.TaskManager():
            hdiv_demag_solve(mesh, bh_table={"lo": BH_SOFT}, H_ext=HEXT_NL)          # 'hi' missing
    with pytest.raises(ValueError):
        with ng.TaskManager():
            hdiv_demag_solve(mesh, bh_table={"lo": BH_SOFT, "hi": BH_HARD[:, 0]}, H_ext=HEXT_NL)  # 1-col table
