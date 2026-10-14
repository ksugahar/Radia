"""Cross-check the C++ axifem element matrices against the
pure-Python reference implementation in `_reference_python/`.

The pure-Python prototype (axifem_core.py: P1 triangle, axifem_quad.py:
Q1 sigma mass, axifem_quad_q2.py: Q2 quad) is the Henrotte/Meeker
formulation ported directly from FEMM's prob3big.cpp StaticAxisymmetric().
The C++ side under src/ext/axifem/ ships the same formulation as an
NGSolve FESpace (AxiHenrotteFESpace) with assembly bilinear-form
integrators (AxiHenrotteStiffnessBFI / AxiHenrotteSigmaMassBFI).

Q1 stiffness is the one deliberate exception to the historical prototype:
production C++ uses the V-DOF stiffness operator that reproduces a uniform
axial B field at order 1.  The old A=psi closed form is kept as background
reference, but this test compares stiffness against an explicit Python
V-DOF quadrature reference with the same nodal convention.

For a single-element mesh, the global stiffness/mass should equal the
element matrix from the pure-Python prototype (modulo DOF ordering).
This test asserts that equivalence at machine precision for a small
set of triangle / quadrilateral configurations.
"""
from math import pi

import numpy as np
import pytest


pytest.importorskip("ngsolve")

from netgen.occ import OCCGeometry, MoveTo
from ngsolve import (
    BilinearForm, CoefficientFunction, Mesh, TaskManager,
)

_GL4_X = np.array([
    -0.8611363115940526, -0.3399810435848563,
     0.3399810435848563,  0.8611363115940526,
])
_GL4_W = np.array([
    0.3478548451374538, 0.6521451548625461,
    0.6521451548625461, 0.3478548451374538,
])


def _q1_inverse_vandermonde(ra, rb, za, zb):
    sa, sb = ra * ra, rb * rb
    V = np.array([
        [1.0, sa, za, sa * za],
        [1.0, sb, za, sb * za],
        [1.0, sb, zb, sb * zb],
        [1.0, sa, zb, sa * zb],
    ])
    return np.linalg.inv(V)


def _python_q1_quad_vdof_stiffness(ra, rb, za, zb, mu):
    """Reference for production Q1StiffnessVDof.

    K_ij = (2 pi / mu) r_i r_j int grad(psi_i).grad(psi_j) / r dr dz,
    where psi is bilinear in (s=r^2, z).  This is intentionally distinct
    from the older A=psi energy form because it preserves a uniform axial
    B field at order 1.
    """
    inv_V = _q1_inverse_vandermonde(ra, rb, za, zb)
    r_nodes = np.array([ra, rb, rb, ra], dtype=float)
    Ke = np.zeros((4, 4))
    for x, wx in zip(_GL4_X, _GL4_W):
        rq = 0.5 * (ra + rb) + 0.5 * (rb - ra) * x
        if rq <= 1.0e-14:
            continue
        for y, wy in zip(_GL4_X, _GL4_W):
            zq = 0.5 * (za + zb) + 0.5 * (zb - za) * y
            w = wx * wy * 0.25 * (rb - ra) * (zb - za)
            s = rq * rq
            dpr = np.zeros(4)
            dpz = np.zeros(4)
            for i in range(4):
                b, c, d = inv_V[1, i], inv_V[2, i], inv_V[3, i]
                dpr[i] = 2.0 * rq * (b + d * zq)
                dpz[i] = c + d * s
            Ke += (w / rq) * (np.outer(dpr, dpr) + np.outer(dpz, dpz))
    K = (2.0 * pi / mu) * np.outer(r_nodes, r_nodes) * Ke
    return 0.5 * (K + K.T)


def _cpp_q1_quad_matrices(ra, rb, za, zb, mu, sigma):
    """Assemble C++ K, M for a single Q1 quad spanning (ra,rb) x (za,zb)."""
    axifem = pytest.importorskip("radia.axifem")
    box = MoveTo(ra, za).Rectangle(rb - ra, zb - za).Face()
    box.faces.name = "conductor"
    # maxh must dominate both dimensions to coerce Netgen into producing
    # exactly ONE quad element (otherwise it splits along the longer axis).
    maxh = 10 * max(rb - ra, zb - za)
    mesh = Mesh(OCCGeometry(box, dim=2).GenerateMesh(
        maxh=maxh, quad_dominated=True))
    assert mesh.ne == 1, f"expected single-element mesh, got ne={mesh.ne}"

    fes = axifem.H1Henrotte(mesh, order=1)
    assert fes.ndof == 4, f"expected ndof=4 (Q1 quad), got {fes.ndof}"

    a = BilinearForm(fes, symmetric=True)
    a += axifem.AxiHenrotteStiffnessBFI(CoefficientFunction(mu))
    m = BilinearForm(fes, symmetric=True)
    m += axifem.AxiHenrotteSigmaMassBFI(CoefficientFunction(sigma))
    with TaskManager():
        a.Assemble()
        m.Assemble()

    K = np.zeros((4, 4))
    rows, cols, vals = a.mat.COO()
    for r, c, v in zip(rows, cols, vals):
        K[r, c] += v
    M = np.zeros((4, 4))
    rows, cols, vals = m.mat.COO()
    for r, c, v in zip(rows, cols, vals):
        M[r, c] += v
    return K, M


def _python_q1_quad_matrices(ra, rb, za, zb, mu, sigma):
    """Assemble Python ref K, M for the same Q1 quad.

    Stiffness uses the V-DOF reference above, matching production C++.
    Sigma mass reuses axifem_quad.py's closed-form expression.
    """
    from axifem_quad import element_sigma_mass_quad
    K = _python_q1_quad_vdof_stiffness(ra, rb, za, zb, mu)
    M = element_sigma_mass_quad(ra, rb, za, zb, sigma)
    return K, M


@pytest.mark.parametrize(
    "ra,rb,za,zb,label",
    [
        # NOTE: the C++ side asks Netgen to produce one quad with maxh
        # >> element size; tall-aspect rectangles still trip Netgen into
        # splitting along the longer axis even with maxh = 10 * max(side).
        # Keep the parametrised cases close to square so the assertion
        # mesh.ne == 1 holds; the topological coverage we need is interior
        # + axis-touching + a non-zero-z offset.
        (1e-3, 2e-3,  0.0,   1e-3, "interior 1-2 mm"),
        (1e-3, 2e-3,  3e-3,  4e-3, "interior 1-2 mm offset z"),
        (0.0,  1e-3,  0.0,   1e-3, "axis-touching"),
    ],
)
def test_q1_quad_cpp_vs_python(ra, rb, za, zb, label):
    """C++ AxiHenrotteFE_Q1_AxisAligned matches the Python V-DOF reference
    for a single quad, in both stiffness (K) and sigma-mass (M).

    Compared via permutation-invariant spectrum (eigvals) — the local
    DOF ordering may differ between the C++ assembly (NGSolve vertex
    numbering) and the Python prototype (s-midpoint convention)."""
    from scipy.linalg import eigh

    mu0   = 4 * pi * 1e-7
    sigma = 5.8e7

    K_c, M_c = _cpp_q1_quad_matrices(ra, rb, za, zb, mu0, sigma)
    K_p, M_p = _python_q1_quad_matrices(ra, rb, za, zb, mu0, sigma)

    eK_c = np.sort(np.linalg.eigvalsh(K_c))
    eK_p = np.sort(np.linalg.eigvalsh(K_p))
    # K has a constant-mode null space; the smallest eigenvalue is
    # numerically zero (~1e-13).  Use atol scaled to the largest
    # eigenvalue rather than 1e-12 to avoid sign flips on the noise
    # mode tripping the comparison.
    K_atol = max(1e-10, 1e-10 * abs(eK_c).max())
    np.testing.assert_allclose(
        eK_c, eK_p, rtol=1e-10, atol=K_atol,
        err_msg=f"K spectrum mismatch ({label})")

    eM_c = np.sort(np.linalg.eigvalsh(M_c))
    eM_p = np.sort(np.linalg.eigvalsh(M_p))
    M_atol = max(1e-15, 1e-10 * abs(eM_c).max())
    np.testing.assert_allclose(
        eM_c, eM_p, rtol=1e-10, atol=M_atol,
        err_msg=f"M spectrum mismatch ({label})")


@pytest.mark.parametrize(
    "ra,rb,za,zb",
    [
        (1.0e-3, 2.0e-3, -0.5e-3, 0.5e-3),
        (0.0, 1.0e-3, 0.0, 1.0e-3),
    ],
)
def test_q1_native_array_api_matches_independent_python_reference(
        ra, rb, za, zb):
    """The adapter-free array API preserves the validated Q1 formulas."""
    axifem = pytest.importorskip("radia.axifem")
    mu0 = 4 * pi * 1e-7
    sigma = 5.8e7
    actual = axifem.q1_magnetic_element_matrices(
        ra, rb, za, zb, mu0, sigma)
    expected_k, expected_m = _python_q1_quad_matrices(
        ra, rb, za, zb, mu0, sigma)

    np.testing.assert_allclose(
        np.asarray(actual["stiffness"]), expected_k,
        rtol=2e-12, atol=max(1e-12, 2e-12 * np.max(np.abs(expected_k))))
    np.testing.assert_allclose(
        np.asarray(actual["sigma_mass"]), expected_m,
        rtol=2e-12, atol=max(1e-15, 2e-12 * np.max(np.abs(expected_m))))
    assert actual["backend"] == "native-pybind"
    assert actual["node_order"] == "(ra,za),(rb,za),(rb,zb),(ra,zb)"


def test_q1_native_array_api_rejects_invalid_geometry():
    axifem = pytest.importorskip("radia.axifem")
    with pytest.raises(ValueError, match="0 <= ra < rb"):
        axifem.q1_magnetic_element_matrices(
            2.0e-3, 1.0e-3, 0.0, 1.0e-3, 4 * pi * 1e-7, 5.8e7)


def test_python_ref_self_consistency_p1_triangle():
    """The pure-Python P1 triangle prototype is loaded from
    _reference_python/axifem_core.py and its element_matrices()
    function returns symmetric Mr, Mz with finite values for any
    valid triangle.  This test exists to make sure the reference
    implementation itself stays healthy when we move it across
    paths during refactors.
    """
    from axifem_core import element_matrices, MU0

    rn = np.array([1.0e-3, 3.0e-3, 2.0e-3])
    zn = np.array([0.0, 0.0, 1.5e-3])
    M_e, Mr, Mz, *_ = element_matrices(rn, zn, MU0, MU0)
    assert np.all(np.isfinite(M_e))
    assert np.all(np.isfinite(Mr))
    assert np.all(np.isfinite(Mz))
    assert np.allclose(Mr, Mr.T, atol=1e-15)
    assert np.allclose(Mz, Mz.T, atol=1e-15)


def test_p2_triangle_curved_mesh_assembles():
    """C++ order=2 triangle path assembles on a curved OCC mesh.

    This pins the public `radia.axifem` import surface and the production
    P2-triangle/mesh.Curve(2) path.  Q2 curved quads are exercised separately
    in test_q2_curved.py.
    """
    axifem = pytest.importorskip("radia.axifem")
    from netgen.geom2d import SplineGeometry

    geo = SplineGeometry()
    geo.AddCircle((2.0e-3, 0.0), 6.0e-4, leftdomain=1, rightdomain=0,
                  bc="outer")
    geo.SetMaterial(1, "conductor")
    mesh = Mesh(geo.GenerateMesh(maxh=4.0e-4))
    mesh.Curve(2)

    fes = axifem.H1Henrotte(mesh, order=2)
    assert fes.ndof > mesh.nv, "order=2 triangle path did not allocate edge DOFs"

    mu0 = 4 * pi * 1e-7
    sigma = 5.8e7
    a = BilinearForm(fes, symmetric=True, check_unused=False)
    a += axifem.AxiHenrotteStiffnessBFI(CoefficientFunction(mu0))
    m = BilinearForm(fes, symmetric=True, check_unused=False)
    m += axifem.AxiHenrotteSigmaMassBFI(CoefficientFunction(sigma))
    with TaskManager():
        a.Assemble()
        m.Assemble()

    K_vals = np.asarray(a.mat.COO()[2], dtype=float)
    M_vals = np.asarray(m.mat.COO()[2], dtype=float)
    assert K_vals.size > 0
    assert M_vals.size > 0
    assert np.all(np.isfinite(K_vals))
    assert np.all(np.isfinite(M_vals))
