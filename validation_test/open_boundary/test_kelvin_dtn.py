# -*- coding: utf-8 -*-
"""Golden tests for radia.open_boundary.kelvin_dtn (the Kelvin-built material-aware /
non-separable DtN -> CLN), ported from the verified research demos
examples/kelvin_transformation/DtN_spectrum/act6_01_kelvin_fem_eddy_dtn.py +
act6_07_cube_eddy_dtn_to_cln.py.

Two tiers:
  * PURE-NUMPY (always run): the radial Kelvin (R/rho')^2-weighted FEM BUILDS the
    eddy DtN and reproduces the closed-form dtn_cln.eddy_dtn (to act6_01's verified
    ~5e-2 over DC->evanescent) with NO DC closure floor; the band-CLN reduces it.
  * NGSOLVE (importorskip): the arbitrary-shape Schur DtN matrix + the generalised
    Steklov ladder of a NON-separable cube is O_h-split (l=2 -> E_g(2)+T_2g(3) = 2+3,
    l=1 dipole degenerate) -- an analytic-free correctness proof -- and band-CLN-reducible.
"""
import numpy as np
import pytest

import radia.open_boundary as ob

BAND = 1j * np.logspace(-4, 2, 40)


# ----------------------------- pure-numpy tier -----------------------------
@pytest.mark.parametrize("n", (1, 2, 3))
def test_kelvin_fem_builds_the_exact_eddy_dtn(n):
    """The radial Kelvin-FEM reproduces the closed-form eddy DtN over DC->evanescent
    (act6_01's verified ~5e-2) and hits the static ladder -(n+1) at DC."""
    G = np.array([ob.kelvin_fem_radial_dtn(n, s) for s in BAND])
    Gex = np.array([ob.eddy_dtn(n, s) for s in BAND])
    nrmse = float(np.sqrt(np.mean(np.abs(G - Gex) ** 2)) / np.sqrt(np.mean(np.abs(Gex) ** 2)))
    assert nrmse < 5e-2, f"n={n}: Kelvin-FEM build != analytic ({nrmse:.1e})"
    dc = abs(ob.kelvin_fem_radial_dtn(n, 1j * 1e-8).real - (-(n + 1)))
    assert dc < 5e-2, f"n={n}: Kelvin-FEM DC != -(n+1) ({dc:.1e})"


def test_kelvin_fem_has_no_dc_closure_floor():
    """The Kelvin compactification has NO DC floor: its DC error is far below the
    (R0/Rfar)^(2n+1) floor a plain truncated (Dirichlet-at-Rfar) FEM would carry."""
    n, Rfar = 1, 3.0
    dc_err = abs(ob.kelvin_fem_radial_dtn(n, 1j * 1e-8) - ob.eddy_dtn(n, 1j * 1e-8)) \
        / abs(ob.eddy_dtn(n, 1j * 1e-8))
    trunc_floor = (1.0 / Rfar) ** (2 * n + 1)            # the floor a truncated FEM cannot beat
    assert dc_err < 1e-3 and dc_err < 0.1 * trunc_floor, \
        f"Kelvin DC err {dc_err:.1e} should be far below the truncated-FEM floor {trunc_floor:.1e}"


def test_band_cln_fit_reduces_eddy_dtn():
    """band_cln_fit reduces a (synthetic) eddy DtN over the band to a few sqrt(s) stages."""
    s_band = 1j * np.logspace(-2, 2, 12)
    G = np.array([ob.eddy_dtn(2, s) for s in s_band])
    _, nrmse = ob.band_cln_fit(s_band, G, 4)
    assert nrmse < 1e-2, f"band-CLN (4 stages) did not reduce the DtN ({nrmse:.1e})"


# ------------------------------ NGSolve tier -------------------------------
def _cube_mesh(maxh, a=1.0, R=3.5):
    pytest.importorskip("netgen")
    from netgen.csg import CSGeometry, Sphere, OrthoBrick, Pnt
    from ngsolve import Mesh
    geo = CSGeometry()
    geo.Add(Sphere(Pnt(0, 0, 0), R).bc("outer")
            - OrthoBrick(Pnt(-a, -a, -a), Pnt(a, a, a)).bc("gamma"))
    return Mesh(geo.GenerateMesh(maxh=maxh))


def test_cube_dtn_is_Oh_split():
    """Kelvin/FEM builds the NON-separable cube eddy DtN; the static Steklov ladder
    splits the sphere's l=2 quintet into E_g(2)+T_2g(3)=2+3 (O_h) while the l=1
    dipole stays a degenerate triplet (T_1u) -- the analytic-free correctness proof."""
    pytest.importorskip("ngsolve")
    from ngsolve import TaskManager
    mesh = _cube_mesh(0.55)
    with TaskManager():
        S, Mg, g_idx = ob.kelvin_dtn_matrix(mesh, 2, 1j * 1e-6)
    w, _ = ob.steklov_spectrum(S, Mg)
    assert np.all(np.diff(w[:9]) >= -1e-6), "Steklov ladder must be real + ordered"
    dip, quad = w[1:4], w[4:9]
    dip_spread = dip.max() - dip.min()
    quad_split = quad.max() - quad.min()
    assert dip_spread < 0.05, "l=1 dipole must stay a degenerate triplet (T_1u)"
    assert quad_split > 4 * dip_spread and quad_split > 0.02, "cube must SPLIT the l=2 quadrupole (O_h)"
    kgap = int(np.argmax(np.diff(quad))) + 1
    assert {kgap, 5 - kgap} == {2, 3}, "O_h must split l=2 into a 2+3 pattern (E_g + T_2g)"


def test_cube_dipole_dtn_band_cln():
    """The built cube dipole DtN interpolates DC->evanescent and a few-stage CLN
    in sqrt(s) reduces it over the band."""
    pytest.importorskip("ngsolve")
    from ngsolve import TaskManager
    mesh = _cube_mesh(0.55)
    with TaskManager():
        S0, Mg, _ = ob.kelvin_dtn_matrix(mesh, 2, 1j * 1e-6)
        w, V = ob.steklov_spectrum(S0, Mg)
        v0 = V[:, 1]                                       # a dipole Steklov mode at DC
        band = 1j * np.logspace(-3, 2, 10)
        Rs = []
        for sv in band:
            Sv, Mgv, _ = ob.kelvin_dtn_matrix(mesh, 2, sv)
            Rs.append(complex(v0 @ Sv @ v0) / complex(v0 @ Mgv @ v0))
    Rs = np.array(Rs)
    assert abs(Rs[0].imag) < 5e-2 and abs(Rs[-1].imag) > 0.1, "Im(DtN) must grow DC->evanescent"
    _, nrmse = ob.band_cln_fit(band, Rs, 6)
    assert nrmse < 2e-2, f"a few-stage CLN must reduce the cube dipole eddy DtN ({nrmse:.1e})"
