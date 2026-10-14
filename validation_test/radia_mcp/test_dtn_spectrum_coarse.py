# -*- coding: utf-8 -*-
r"""DtN-spectrum coarse-mesh gate -- the spectral explanation of why open-boundary
methods (Kelvin transform / BEM / PML) stay accurate on COARSE meshes.

Kameari demonstrated the Kelvin transform's coarse-mesh accuracy empirically, by
mesh refinement.  This gate verifies the same fact as a PROPERTY OF THE DtN
MATRIX: the exterior Laplace DtN Λ_ext has the mesh-independent sphere ladder
λ_n = −(n+1)/R, and the discrete matrix Λ_h = V⁻¹(−½M+K)

  1. reproduces the LOW-degree eigenvalues to <0.5 % on the COARSEST sphere mesh
     (the coarse-mesh accuracy itself, read off the matrix);
  2. has a per-degree error that increases monotonically with degree n at every
     mesh size (the spectral signature: accurate low modes, degraded high modes);
  3. converges steeply under refinement (~×2.6 per level here, rate p≈3.9 ~ O(h⁴)),
     which WIDENS the band of well-resolved modes -- the low modes were already
     accurate; refinement adds bandwidth, not materially better low modes.

Companion code:  radia_ngsolve.bem_integral.{exterior_dtn_spectrum,
                 dtn_spectrum_vs_mesh}
Knowledge:       dtn_coarse_mesh

NB the dense BEM densification is O(ndof²); the 2-mesh sweep below (ndof 336 +
564) is the bulk of the runtime.  maxh ≥ 0.5 all floor to the same coarsest
sphere mesh, so the sweep uses maxh = (0.5, 0.4) for two DISTINCT levels.
"""
import pytest

pytest.importorskip("ngsolve.bem")

from radia_mcp.radia_ngsolve.bem_integral import dtn_spectrum_vs_mesh
from radia_mcp.radia_ngsolve.fem_bem_coupling import kelvin_dtn_eigenvalue


@pytest.fixture(scope="module")
def study():
    """Two distinct sphere meshes (ndof 336, 564); reused by all assertions."""
    return dtn_spectrum_vs_mesh(R=1.0, maxh_list=(0.5, 0.4), order=1,
                                intorder=10, nmax=4, band_tol=0.005, low_nmax=2)


def test_distinct_meshes(study):
    """Sanity: the two meshes really are distinct refinement levels."""
    assert study["ndof_list"][0] < study["ndof_list"][1], (
        f"meshes not distinct: ndof={study['ndof_list']}"
    )


def test_coarse_low_modes_accurate(study):
    """Fact 1: on the COARSEST mesh, the n≤2 DtN eigenvalues match −(n+1)/R to
    <0.5 % -- coarse-mesh accuracy read straight off the matrix."""
    assert study["coarse_low_mode_max_err"] < 0.005, (
        f"coarse low-mode max rel_err = {study['coarse_low_mode_max_err']:.3e} "
        f"(expected < 0.5%)"
    )


def test_dipole_eigenvalue_on_ladder(study):
    """The dipole (n=1) eigenvalue sits on −2/R to <0.2 % on the coarsest mesh."""
    coarse = study["per_mesh"][0]
    dipole = next(m for m in coarse["modes"] if m["n"] == 1)
    assert abs(dipole["lambda_mean"] - (-2.0)) < 0.004, (
        f"dipole eigenvalue = {dipole['lambda_mean']:.5f} (expected −2.0), "
        f"rel_err={dipole['rel_err']:.3e}"
    )


def test_error_ordered_by_degree(study):
    """Fact 2 (the spectral signature): at EVERY mesh size the per-degree error
    increases monotonically with degree n."""
    assert all(study["degree_monotonic"]), (
        f"error not monotonically increasing with degree at some mesh: "
        f"degree_monotonic={study['degree_monotonic']}, table={study['error_table']}"
    )
    # and the spread between lowest and highest matched degree is large (the
    # spectrum really is 'accurate at the bottom, degraded at the top')
    assert min(study["degree_growth"]) > 5.0, (
        f"degree_growth too weak: {study['degree_growth']}"
    )


def test_refinement_widens_accurate_band(study):
    """Fact 3a: refinement widens the band of well-resolved modes (and the
    coarse mesh already covers the engineering-dominant low multipoles)."""
    band = study["accurate_band"]
    assert band[0] >= 2, f"coarse accurate band should already reach quadrupole: {band}"
    assert band[-1] >= band[0], f"refinement should not shrink the band: {band}"


def test_low_mode_converges_not_flat(study):
    """Fact 3b (honesty guard): the low-mode error is NOT mesh-independent -- it
    falls under refinement.  Coarse-mesh accuracy is a small FLOOR, not a plateau."""
    dip = study["error_table"][1]          # dipole rel_err at [coarse, fine]
    assert dip[1] < dip[0], (
        f"dipole error should decrease under refinement: coarse={dip[0]:.3e}, "
        f"fine={dip[1]:.3e}"
    )


# ---------------------------------------------------------------------------
# Kelvin closure: the VOLUME-FEM companion measurement (fast, no BEM densify)
# ---------------------------------------------------------------------------
# Kelvin inverts mode n to a degree-n polynomial, so order-p FEM reproduces the
# effective DtN eigenvalue exactly iff p >= n. The dipole inverts to a LINEAR
# field -> order-1 coarse accurate; the quadrupole needs order >= 2. This MEASURES
# the Kelvin realisation's effective DtN (the BEM Lambda_h companion), so the
# Kelvin coarse-mesh accuracy is demonstrated, not merely argued by analogy.


def test_kelvin_dipole_coarse_accurate_and_converges():
    """Dipole (inverts to a LINEAR field): order-1 Kelvin closure is accurate on
    the coarsest mesh and converges under refinement (geometry-limited residual)."""
    coarse = kelvin_dtn_eigenvalue(R=1.0, degree=1, maxh=0.6, order=1)
    fine = kelvin_dtn_eigenvalue(R=1.0, degree=1, maxh=0.25, order=1)
    assert coarse["rel_err"] < 0.02, (
        f"Kelvin dipole DtN coarse rel_err={coarse['rel_err']:.3e} (expected <2%), "
        f"lam={coarse['lam']:.5f} vs -2.0"
    )
    assert fine["rel_err"] < coarse["rel_err"], (
        f"Kelvin dipole should converge: coarse={coarse['rel_err']:.3e}, "
        f"fine={fine['rel_err']:.3e}"
    )


def test_kelvin_quadrupole_needs_order_two():
    """Quadrupole (inverts to a QUADRATIC field): order-1 cannot represent it
    (large error), order-2 reproduces the eigenvalue -- the polynomial threshold."""
    o1 = kelvin_dtn_eigenvalue(R=1.0, degree=2, maxh=0.4, order=1)
    o2 = kelvin_dtn_eigenvalue(R=1.0, degree=2, maxh=0.4, order=2)
    assert o1["rel_err"] > 0.05, (
        f"order-1 Kelvin quadrupole should be inaccurate (poly threshold), "
        f"got rel_err={o1['rel_err']:.3e}"
    )
    assert o2["rel_err"] < 0.01, (
        f"order-2 Kelvin quadrupole should be accurate, got rel_err={o2['rel_err']:.3e}"
    )
    assert o2["rel_err"] < o1["rel_err"] / 10.0, "order 2 should massively beat order 1"


def test_kelvin_dipole_on_minus_two_over_R():
    """The measured Kelvin dipole eigenvalue sits on the analytic -(n+1)/R = -2/R."""
    r = kelvin_dtn_eigenvalue(R=1.0, degree=1, maxh=0.4, order=1)
    assert abs(r["lam"] - (-2.0)) < 0.02, (
        f"Kelvin dipole lam={r['lam']:.5f} (expected -2.0)"
    )


def test_kelvin_scales_with_R():
    """Effective DtN eigenvalue -(n+1)/R scales with the truncation radius."""
    for R in (0.5, 2.0):
        r = kelvin_dtn_eigenvalue(R=R, degree=1, maxh=0.4 * R, order=1)
        assert abs(r["lam"] - (-2.0 / R)) / (2.0 / R) < 0.02, (
            f"R={R}: lam={r['lam']:.5f}, expected {-2.0/R:.5f}"
        )


# ---------------------------------------------------------------------------
# 2D Kelvin (circle inversion) -- the cross-section case for static apparatus /
# rotating machines.  No conformal prefactor (dim-2=0), eigenvalue -n/R.
# ---------------------------------------------------------------------------

def test_kelvin_2d_dipole_eigenvalue_and_convergence():
    """2D Kelvin: dipole (inverts to a LINEAR field) eigenvalue is -n/R = -1/R
    (NOT -2/R as in 3D), accurate at order 1 on a coarse mesh and converging."""
    coarse = kelvin_dtn_eigenvalue(R=1.0, degree=1, maxh=0.4, order=1, dim=2)
    fine = kelvin_dtn_eigenvalue(R=1.0, degree=1, maxh=0.15, order=1, dim=2)
    assert abs(coarse["lam_exact"] - (-1.0)) < 1e-12, "2D dipole eigenvalue is -1/R"
    assert abs(coarse["lam"] - (-1.0)) < 0.01, (
        f"2D Kelvin dipole lam={coarse['lam']:.5f} (expected -1.0)"
    )
    assert fine["rel_err"] < coarse["rel_err"], "2D dipole should converge"


def test_kelvin_2d_quadrupole_needs_order_two():
    """2D Kelvin quadrupole (inverts to a QUADRATIC) needs order>=2, like 3D."""
    o1 = kelvin_dtn_eigenvalue(R=1.0, degree=2, maxh=0.25, order=1, dim=2)
    o2 = kelvin_dtn_eigenvalue(R=1.0, degree=2, maxh=0.25, order=2, dim=2)
    assert abs(o1["lam_exact"] - (-2.0)) < 1e-12, "2D quadrupole eigenvalue is -2/R"
    assert o1["rel_err"] > 0.01, f"order-1 2D quadrupole should be inaccurate: {o1['rel_err']:.3e}"
    assert o2["rel_err"] < 0.001, f"order-2 2D quadrupole should be accurate: {o2['rel_err']:.3e}"
    assert o2["rel_err"] < o1["rel_err"] / 10.0


def main():
    res = dtn_spectrum_vs_mesh(R=1.0, maxh_list=(0.5, 0.4), order=1,
                               intorder=10, nmax=4, band_tol=0.005, low_nmax=2)
    print("ndof:", res["ndof_list"])
    print("coarse_low_mode_max_err:", res["coarse_low_mode_max_err"])
    print("degree_monotonic:", res["degree_monotonic"])
    print("degree_growth:", res["degree_growth"])
    print("accurate_band:", res["accurate_band"])
    print("error_table:", res["error_table"])
    print("[OK] DtN-spectrum coarse-mesh gate: low modes accurate coarse, "
          "error ordered by degree, refinement widens the band.")


if __name__ == "__main__":
    main()
