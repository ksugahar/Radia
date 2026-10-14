"""Coupled-system physics vs the exact conducting-permeable-sphere alpha(omega).

Protocol (2026-07-28 study): S1 calibrates ONE complex amplitude c at the
smallest a/delta point of the mu_r=1 eddy-only sweep (it absorbs the
reduced-basis drive/extraction amplitude convention); everything downstream is
calibration-frozen.  S2 solves the coupled [[A_M, -K/mu0],[s K^H, Z(s)]]
system at mu_r=100 (the solve_frequency physical scales) and compares the
total dipole m = int M dV + 1/2 int r x J dV against the exact alpha; S3 locks
the vanishing DC current from a static magnetization (the s K^H fix).

Measured anchors (ne=215 sphere, order-2 HCurl, steps=3): |c| = 1.425;
mu_r=1 shape rel err 31% / 65% at a/delta = 1 / 3 (volumetric-basis
saturation ceilings); mu_r=100 static 4.9% (n_M=1 reduction quality),
transition 4.8% / 3.4% / 22% at a/delta_mu = 1 / 3 / 10; DC ratio 0.0100.
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")

import ngsolve as ng  # noqa: E402
from netgen import occ  # noqa: E402

from radia import vim  # noqa: E402
from radia.analytical_formulas import sphere_complex_polarizability  # noqa: E402

MU0 = 4.0e-7 * np.pi
A_R = 0.01
SIGMA = 5.8e7
MU_R = 100.0
MAXH = 0.004


def f_of_ad(ad, mu_r=1.0):
    delta = A_R / ad
    return 1.0 / (np.pi * MU0 * mu_r * SIGMA * delta**2)


@pytest.fixture(scope="module")
def lane():
    sph = occ.Sphere(occ.Pnt(0, 0, 0), A_R)
    sph.mat("cond")
    for face in sph.faces:
        face.name = "skin"
    mesh = ng.Mesh(occ.OCCGeometry(sph).GenerateMesh(maxh=MAXH))

    fes = ng.HCurl(mesh, order=2, nograds=True)
    u, v = fes.TnT()
    stiffness = ng.BilinearForm(fes)
    stiffness += ng.curl(u) * ng.curl(v) * ng.dx + 0.05 * u * v * ng.dx
    mass = ng.BilinearForm(fes)
    mass += u * v * ng.dx
    port = ng.LinearForm(fes)
    port += ng.CoefficientFunction((-ng.y, ng.x, 0.0)) * v * ng.dx
    with ng.TaskManager():
        stiffness.Assemble()
        mass.Assemble()
        port.Assemble()

    def a_ext(points):
        points = np.asarray(points)
        return np.column_stack(
            (-points[:, 1], points[:, 0], np.zeros(points.shape[0])))

    volumetric_gate = vim.EddySIBCApplicability(
        frequency_hz=1.0, sigma=SIGMA, characteristic_thickness_m=2.0 * A_R)
    with ng.TaskManager():
        mixed = vim.NgsolveBDMEddyBubbleVIM(
            mesh, fes, stiffness, mass, port,
            (ng.CoefficientFunction((1.0, 0.0, 0.0)),),
            hdiv_order=1, mu_r=MU_R,
            external_fields=(ng.CoefficientFunction((0.0, 0.0, 1.0)),),
            hdiv_max_modes=3, magnetic_materials="cond",
            steps=3, sigma=SIGMA, conductive_materials="cond",
            response_backend="dense", intorder=2,
            port_vector_potentials=(a_ext,),
            coupling_kernel_epsilon=0.25 * MAXH,
            sibc_applicability=volumetric_gate,
        )

    n_m, n_j = mixed.n_hdiv_mmm_modes, mixed.n_hcurl_vim_modes
    rows = []
    for basis in mixed.eddy_bases:
        if basis.n_modes == 0:
            continue
        r = basis.points
        modes = basis.modes
        cross_z = r[:, 0] * modes[:, :, 1] - r[:, 1] * modes[:, :, 0]
        rows.append(0.5 * np.einsum("mi,i->m", cross_z, basis.weights))
    dip_j = np.concatenate(rows)
    mb = mixed.magnetization_basis
    dip_m = np.einsum("mi,i->m", mb.modes[:, :, 2], mb.weights)
    a_col = np.asarray(mixed.eddy_rhs).reshape(n_j, -1)[:, 0]
    b_m = np.asarray(mixed.magnetic_rhs).reshape(n_m, -1)[:, 0]
    K = np.asarray(mixed.coupling)

    def eddy_Z(s):
        return np.asarray(mixed.eddy_impedance(
            s, surface_impedance=vim.SkinImpedance(s, SIGMA)))

    # S1 calibration at a/delta = 0.3 (mu_r = 1: eddy block only)
    def alpha_eddy(f):
        s = 2j * np.pi * f
        x_j = np.linalg.solve(eddy_Z(s), -s * (0.5 * MU0) * a_col)
        return complex(dip_j @ x_j)

    cal = complex(sphere_complex_polarizability(f_of_ad(0.3), A_R, SIGMA, 1.0)) \
        / alpha_eddy(f_of_ad(0.3))

    red = mixed.hdiv_reduction
    A_M = np.asarray(red.mass) / (MU_R - 1.0) + np.asarray(red.demag)

    def solve_coupled(f):
        s = 2j * np.pi * f
        O = np.block([[A_M, (-1.0 / MU0) * K],
                      [s * K.conj().T, eddy_Z(s)]])
        b = np.concatenate([b_m, -s * (0.5 * MU0) * a_col])
        x = np.linalg.solve(O, b)
        return complex(dip_m @ x[:n_m] + cal * (dip_j @ x[n_m:])), x

    return dict(mixed=mixed, n_m=n_m, cal=cal, alpha_eddy=alpha_eddy,
                solve_coupled=solve_coupled, eddy_Z=eddy_Z, A_M=A_M,
                b_m=b_m, K=K)


def test_calibration_amplitude_is_stable(lane):
    # the basis drive/extraction convention scale; measured |c| = 1.425
    assert 1.2 < abs(lane["cal"]) < 1.7
    assert abs(np.angle(lane["cal"])) < 0.15


def test_mu1_eddy_shape_within_documented_ceilings(lane):
    """Volumetric-basis saturation: quantitative to a/delta ~ 1, ceiling-bound
    beyond (the SIBC surface branch is the designed strong-skin path)."""
    for ad, ceiling in ((1.0, 0.45), (3.0, 0.80)):
        f = f_of_ad(ad)
        num = lane["cal"] * lane["alpha_eddy"](f)
        exact = complex(sphere_complex_polarizability(f, A_R, SIGMA, 1.0))
        assert abs(num - exact) / abs(exact) < ceiling, ad


def test_mu100_static_and_transition(lane):
    static_exact = 4.0 * np.pi * A_R**3 * (MU_R - 1.0) / (MU_R + 2.0)
    m0, x0 = lane["solve_coupled"](1e-6)
    assert abs(m0 - static_exact) / static_exact < 0.10   # n_M=1 reduction band
    assert np.linalg.norm(x0[lane["n_m"]:]) < 1e-3        # no DC circulation

    for ad, band in ((1.0, 0.10), (3.0, 0.10), (10.0, 0.35)):
        f = f_of_ad(ad, MU_R)
        m, _ = lane["solve_coupled"](f)
        exact = complex(sphere_complex_polarizability(f, A_R, SIGMA, MU_R))
        assert abs(m - exact) / abs(exact) < band, ad


def test_static_magnetization_drives_vanishing_current(lane):
    x_m = np.linalg.solve(lane["A_M"], lane["b_m"])
    norms = []
    for f in (1e-6, 1e-4):
        s = 2j * np.pi * f
        x_j = np.linalg.solve(lane["eddy_Z"](s),
                              s * (lane["K"].conj().T @ x_m))
        norms.append(np.linalg.norm(x_j))
    assert norms[0] < 0.02 * norms[1]      # prop. s (volumetric branch)
    assert norms[0] < 1e-3
