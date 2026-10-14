"""Genus-1 (H1 cohomology class) lane: the short-circuited-ring time constant.

A conducting torus in a uniform axial AC field responds through its
topological loop current: (R + sL) I = -s Phi_ext, so the induced dipole is a
single pole m(s) = -s A/(1 + s tau) whose time constant tau = L/R is
CALIBRATION-FREE (the model is invariant under joint (A,R,L) rescaling; tau is
the pure cycle-class observable).

Reference closed forms (thin wire): L = mu0 R_t (ln(8 R_t/r_w) - 7/4),
R = 2 pi R_t/(sigma pi r_w^2).  Two documented mesh compromises (2026-07-28
study): (i) the mesh is held under the 512-tet analytic-interaction cap, which
under-fills the faceted wire by ~25% -- the reference is therefore evaluated
with the MEASURED cross-section A_mesh = V_mesh/(2 pi R_t) (inscribed-polygon
correction); (ii) the loop inductance comes from the SAMPLED mollified
bridge-bridge kernel, whose DEFAULT epsilon heuristic mis-sizes it ~3x at this
coarseness -- kernel_epsilon is set explicitly to ~half the sample spacing
(measured tau/tau_ref: 3.11 default -> 0.93 at 0.5 mm).

Measured anchors: single-pole shape rms 0.33%; tau/tau_ref = 0.93;
flux-freeze plateau |m(3k)|/|m(1k)| = 1.002; DC ratio 0.0100 (prop. s).
"""
import numpy as np
import pytest

pytest.importorskip("ngsolve")
from scipy.optimize import least_squares  # noqa: E402

import ngsolve as ng  # noqa: E402
from netgen import occ  # noqa: E402

from radia import vim  # noqa: E402

MU0 = 4.0e-7 * np.pi
R_T, R_W = 0.02, 0.004
SIGMA = 5.8e7
MU_R = 10.0
KERNEL_EPS = 0.5e-3          # ~half the sample spacing; see module docstring

L_RING = MU0 * R_T * (np.log(8.0 * R_T / R_W) - 1.75)
R_RING = 2.0 * np.pi * R_T / (SIGMA * np.pi * R_W**2)


@pytest.fixture(scope="module")
def lane():
    wp = occ.WorkPlane(occ.Axes(p=(R_T, 0, 0), n=(0, 1, 0), h=(1, 0, 0)))
    torus = wp.Circle(R_W).Face().Revolve(occ.Axis((0, 0, 0), (0, 0, 1)), 360)
    torus.mat("cond")
    for f in torus.faces:
        f.name = "skin"
    mesh = ng.Mesh(occ.OCCGeometry(torus).GenerateMesh(
        maxh=0.010, curvaturesafety=0.8, segmentsperedge=0.3))
    assert 150 < mesh.ne <= 512      # the analytic-interaction cap contract
    v_mesh = float(ng.Integrate(ng.CoefficientFunction(1.0), mesh))

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

    gate = vim.EddySIBCApplicability(
        frequency_hz=1.0, sigma=SIGMA, characteristic_thickness_m=2.0 * R_W)
    with ng.TaskManager():
        mixed = vim.NgsolveBDMEddyBubbleVIM(
            mesh, fes, stiffness, mass, port,
            (ng.CoefficientFunction((1.0, 0.0, 0.0)),),
            hdiv_order=1, mu_r=MU_R,
            external_fields=(ng.CoefficientFunction((0.0, 0.0, 1.0)),),
            hdiv_max_modes=1, magnetic_materials="cond",
            steps=2, sigma=SIGMA, conductive_materials="cond",
            response_backend="dense", intorder=2,
            port_vector_potentials=(a_ext,),
            coupling_kernel_epsilon=1.5e-3,
            kernel_epsilon=KERNEL_EPS,
            sibc_applicability=gate,
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
    a_col = np.asarray(mixed.eddy_rhs).reshape(n_j, -1)[:, 0]
    b_m = np.asarray(mixed.magnetic_rhs).reshape(n_m, -1)[:, 0]
    K = np.asarray(mixed.coupling)

    def eddy_Z(s):
        return np.asarray(mixed.eddy_impedance(
            s, surface_impedance=vim.SkinImpedance(s, SIGMA)))

    def m_of(f):
        s = 2j * np.pi * f
        x_j = np.linalg.solve(eddy_Z(s), -s * (0.5 * MU0) * a_col)
        return complex(dip_j @ x_j)

    # mesh-geometry-consistent ring reference (inscribed-polygon correction)
    a_mesh = v_mesh / (2.0 * np.pi * R_T)
    r_eff = np.sqrt(a_mesh / np.pi)
    tau_ref = (MU0 * R_T * (np.log(8.0 * R_T / r_eff) - 1.75)) \
        / (2.0 * np.pi * R_T / (SIGMA * a_mesh))
    return dict(mixed=mixed, m_of=m_of, eddy_Z=eddy_Z, tau_ref=tau_ref,
                v_mesh=v_mesh, b_m=b_m, K=K, n_m=n_m)


def test_h1_loop_single_pole_and_calibration_free_tau(lane):
    freqs = np.array([10.0, 20.0, 40.0, 70.0, 100.0, 140.0, 200.0, 300.0])
    m = np.array([lane["m_of"](f) for f in freqs])
    w = 2.0 * np.pi * freqs

    def model(p, w):
        A, tau = p
        s = 1j * w
        return -s * A / (1.0 + s * tau)

    def resid(p):
        d = model(p, w) - m
        return np.concatenate([d.real, d.imag]) / np.max(np.abs(m))

    fit = least_squares(resid, x0=[abs(m[0].imag) / w[0], L_RING / R_RING],
                        bounds=([0, 0], [np.inf, np.inf]))
    A_fit, tau_fit = fit.x
    shape_rms = float(np.sqrt(np.mean(np.abs(model(fit.x, w) - m)**2))
                      / np.max(np.abs(m)))
    # the pole EXISTS and is a clean single pole (the H1 class carries it)
    assert shape_rms < 0.05
    # calibration-free topological time constant vs the mesh-corrected ring law
    assert 0.75 < tau_fit / lane["tau_ref"] < 1.15
    # flux-freeze flattening above the pole
    plateau = abs(lane["m_of"](3000.0)) / abs(lane["m_of"](1000.0))
    assert 0.90 < plateau < 1.15


def test_static_magnetization_drives_vanishing_loop_current(lane):
    """The cycle-class DC zero: the s K^H coupling (Faraday on the H1 class)
    makes a static magnetization drive NO persistent loop current."""
    red = lane["mixed"].hdiv_reduction
    A_M = np.asarray(red.mass) / (MU_R - 1.0) + np.asarray(red.demag)
    x_m = np.linalg.solve(A_M, lane["b_m"])
    norms = []
    for f in (1e-6, 1e-4):
        s = 2j * np.pi * f
        x_j = np.linalg.solve(lane["eddy_Z"](s),
                              s * (lane["K"].conj().T @ x_m))
        norms.append(np.linalg.norm(x_j))
    assert norms[0] < 0.02 * norms[1]     # prop. s
