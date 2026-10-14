"""Golden test: the PRODUCTION curved + high-order surface demag Gram via ngsolve.bem (Laplace
single-layer), validated vs the sphere ANALYTIC 1/3.

The HDiv-VIM surface demag Gram IS the Laplace single-layer of sigma = M.n; ngsolve.bem supplies it
high-order + curved + FMM-accelerated.  This locks the combined curved + high-order win on the demag
factor ITSELF (the proper Gram, unlike the crude sub-point method whose quadrature bias masks it):

  - curved + order-2 -> demag EXACT (~1e-4%) at a COARSE mesh;
  - flat is floored (~0.25%) and ORDER-INSENSITIVE (sigma constant per flat face);
  - curved p-convergence: order-2 error is >>10x smaller than order-0.

See validation_test/feec/vim_legacy/hdiv_demag_bem_singlelayer.py.
"""
import os
import sys

import pytest

pytest.importorskip("ngsolve")
pytest.importorskip("ngsolve.bem")
pytest.importorskip("netgen.csg")

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "vim_legacy"))
import hdiv_demag_bem_singlelayer as bem  # noqa: E402


def test_curved_highorder_demag_exact():
    """curved + order-2 single-layer demag is EXACT (within 0.05%) of 1/3 at a coarse mesh."""
    Dz, _ = bem.demag_singlelayer(h=0.6, curve_order=3, charge_order=2)
    assert abs(Dz / (1.0 / 3.0) - 1) < 5e-4, f"curved o2 demag {Dz:.5f} not within 0.05% of 1/3"


def test_flat_is_order_insensitive():
    """On a FLAT mesh sigma=M.n is constant per face, so order does NOT help -- the demag factor is the
    SAME (to 1e-4) at order 0/1/2 and floored ~0.25% by the faceting."""
    d0, _ = bem.demag_singlelayer(0.6, 0, 0)
    d2, _ = bem.demag_singlelayer(0.6, 0, 2)
    assert abs(d0 - d2) < 1e-4, f"flat order should not change demag: o0 {d0:.5f} vs o2 {d2:.5f}"
    assert 1e-3 < abs(d0 / (1.0 / 3.0) - 1) < 0.02, "flat should floor at a small but nonzero faceting error"


def test_curved_p_convergence():
    """The high-order win on CURVED geometry: order-2 error is >10x smaller than order-0 (the n_z that
    varies across a curved face is under-resolved by piecewise-constant charge)."""
    d0, _ = bem.demag_singlelayer(0.6, 3, 0)
    d2, _ = bem.demag_singlelayer(0.6, 3, 2)
    e0 = abs(d0 / (1.0 / 3.0) - 1)
    e2 = abs(d2 / (1.0 / 3.0) - 1)
    assert e2 < 0.1 * e0, f"curved p-convergence weak: o0 err {e0:.2e} -> o2 err {e2:.2e}"


def test_prolate_nonisotropic_shape_exact():
    """The non-isotropic SHAPE test: a 2:1 prolate spheroid has N_z != 1/3 (~0.1736).  curved + order-2
    single-layer nails the ANALYTIC N_z (<0.05%) while the flat faceted ellipsoid floors (order-
    insensitive) -- confirming the single-layer gets the anisotropic demag right, not just isotropy."""
    Nz = bem.prolate_Nz_analytic(2.0)
    assert 0.15 < Nz < 0.20, f"analytic prolate N_z sanity {Nz:.4f}"
    d_flat0, _ = bem.demag_prolate(2.0, 0.6, 0, 0)
    d_flat2, _ = bem.demag_prolate(2.0, 0.6, 0, 2)
    d_curv2, _ = bem.demag_prolate(2.0, 0.6, 3, 2)
    assert abs(d_curv2 / Nz - 1) < 5e-4, f"curved o2 prolate N_z {d_curv2:.5f} not within 0.05% of {Nz:.5f}"
    assert abs(d_flat0 - d_flat2) < 1e-4, "flat order should not change the (faceted) prolate demag"
    assert abs(d_curv2 / Nz - 1) < abs(d_flat0 / Nz - 1), "curved o2 must beat flat on the analytic N_z"


def test_spheroid_full_tensor_and_sum_rule():
    """The full demag TENSOR: for BOTH a prolate (c=2) and an oblate (c=0.5) spheroid, curved + order-2
    single-layer gets the POLAR N_par AND the TRANSVERSE N_perp exact vs analytic (Osborn), and the
    formula-independent SUM RULE N_par + 2 N_perp = 1 holds (N_x=N_y for a body of revolution).  This
    validates the WHOLE anisotropic demag tensor, both axes computed independently from the operator."""
    for c in (2.0, 0.5):
        Npar_an = bem.spheroid_Nz_analytic(c)
        Nperp_an = (1.0 - Npar_an) / 2.0
        Npar, _ = bem.demag_spheroid(c, 0.6, 3, 2, axis=2)     # polar (z) factor
        Nperp, _ = bem.demag_spheroid(c, 0.6, 3, 2, axis=0)    # transverse (x) factor
        assert abs(Npar / Npar_an - 1) < 1e-3, f"c={c}: polar N_par {Npar:.5f} vs analytic {Npar_an:.5f}"
        assert abs(Nperp / Nperp_an - 1) < 1e-3, f"c={c}: transverse N_perp {Nperp:.5f} vs {Nperp_an:.5f}"
        assert abs(Npar + 2.0 * Nperp - 1.0) < 2e-3, f"c={c}: sum rule N_par+2N_perp = {Npar+2*Nperp:.5f} != 1"


def test_triaxial_ellipsoid_full_anisotropy():
    """The canonical fully-anisotropic benchmark: a general triaxial ellipsoid (a!=b!=c) has three
    DISTINCT demag factors.  curved + order-2 single-layer gets each exact vs the analytic Osborn
    integral, and the three sum to 1 -- confirming the operator handles full anisotropy, not just the
    degenerate spheroid symmetry."""
    a, b, c = 1.0, 1.5, 2.0
    tot = 0.0
    for axis in (0, 1, 2):
        Nq, _ = bem.demag_ellipsoid(a, b, c, 0.5, 3, 2, axis=axis)
        Nq_an = bem.ellipsoid_N_analytic(a, b, c, axis)
        assert abs(Nq / Nq_an - 1) < 2e-3, f"axis {axis}: N {Nq:.5f} vs analytic {Nq_an:.5f}"
        tot += Nq
    assert abs(tot - 1.0) < 2e-3, f"triaxial sum rule Na+Nb+Nc = {tot:.5f} != 1"
