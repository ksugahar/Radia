"""Golden tests for the radia-levitation Mixed-Galerkin polarizability stack.

Locks the numerical invariants of the alpha(s) pipeline (Phase L-3):

  * the closed-form wedge function W(alpha) = (4/pi) cot(alpha/2),
  * the planar SIBC anchor K_SIBC = S sqrt(sigma/mu),
  * the cuboid edge correction c_1 = -16(a+b+c)/(pi mu) (CAD-direct),
  * the perfect-conductor high-frequency limit alpha(s)->V,
  * the Simulink/MATLAB state-space LTI vs the closed-form Y(s),
  * an end-to-end NGSolve cube alpha(s) sweep.

Pure-algebra tests always run.  Tests needing netgen.occ / ngsolve
`importorskip` out when those are absent (e.g. a minimal-dep CI box).

Scope note: the Mixed-Galerkin Y(s) carries a high-frequency Mellin tail
K_SIBC/sqrt(s) + c_1/s which DIVERGES as s->0, so alpha(s)/V is only the
physical polarizability in the mid-to-high band, NOT at DC.  The golden
limit locked here is therefore the PEC limit alpha(s->inf)->V, not the
DC limit.
"""
from __future__ import annotations

import cmath
import math

import numpy as np
import pytest

from radia.levitation.mixed_galerkin import alpha as A
from radia.levitation.mixed_galerkin import rom_fit as RF
from radia.levitation.simulink import export as EX

MU0 = 4.0 * math.pi * 1e-7
SIGMA_CU = 5.8e7


# ---------------------------------------------------------------------------
# Pure-algebra invariants (no external solver)
# ---------------------------------------------------------------------------

def test_wedge_function_special_values():
    """W(pi/2) = 4/pi (right angle), W(pi) = 0 (flat), W(3pi/2) = -4/pi (re-entrant)."""
    assert A.wedge_function(math.pi / 2) == pytest.approx(4.0 / math.pi, rel=1e-12)
    assert A.wedge_function(math.pi * (1 - 1e-9)) == pytest.approx(0.0, abs=1e-5)
    assert A.wedge_function(3 * math.pi / 2) == pytest.approx(-4.0 / math.pi, rel=1e-12)


def test_K_SIBC_total_formula():
    """K_SIBC = S sqrt(sigma/mu)."""
    S = 6.0e-4
    assert A.K_SIBC_total(S, SIGMA_CU, MU0) == pytest.approx(
        S * math.sqrt(SIGMA_CU / MU0), rel=1e-12
    )


def test_c1_polyhedral_cube_closed_form():
    """c_1 = -(1/mu) sum L_e W(alpha_e); a cube has 12 edges of length L at pi/2.

    Closed form: c_1 = -(1/mu) * 12 L * (4/pi) = -48 L / (pi mu)
                     = -16 (a+b+c) / (pi mu)  with a=b=c=L (3L total / axis).
    """
    L = 5e-3
    edges = [(L, math.pi / 2)] * 12
    c1 = A.c1_polyhedral(edges, MU0)
    closed = -16.0 * (3 * L) / (math.pi * MU0)
    assert c1 == pytest.approx(closed, rel=1e-12)


def test_alpha_pec_high_frequency_limit():
    """alpha(s)/V -> 1 as s -> inf (perfect-conductor flux exclusion)."""
    V = 1.25e-7
    lam = np.array([10.0, 40.0, 90.0])
    tau = MU0 * SIGMA_CU / lam
    g_n = np.array([SIGMA_CU * V * 0.30, SIGMA_CU * V * 0.10, SIGMA_CU * V * 0.05])
    K = A.K_SIBC_total(6e-4, SIGMA_CU, MU0)
    c1 = -6.08e4
    a_hi = A.alpha_from_Y(A.Y_mixed(1j * 2 * math.pi * 1e12, lam, tau, g_n, K, c1),
                          V, SIGMA_CU)
    assert (a_hi / V).real == pytest.approx(1.0, abs=2e-3)
    assert abs((a_hi / V).imag) < 1e-2


def test_state_space_lti_matches_closed_form_Y():
    """ss(A,B,C,D) transfer function reproduces Y_mixed(s) to <1% over 1 Hz-10 kHz."""
    V = 1.25e-7
    lam = np.array([10.0, 40.0, 90.0])
    tau = MU0 * SIGMA_CU / lam
    g_n = np.array([SIGMA_CU * V * 0.30, SIGMA_CU * V * 0.10, SIGMA_CU * V * 0.05])
    K = A.K_SIBC_total(6e-4, SIGMA_CU, MU0)
    c1 = -6.08e4

    A_, B_, C_, D_, n_f, n_w, n_i = EX.build_state_space(
        g_n, tau, V, SIGMA_CU, K, c1, n_warburg_rungs=40
    )
    assert n_f == len(g_n)
    assert n_i == 1

    eye = np.eye(A_.shape[0])
    for f in (1.0, 1e2, 1e4):
        s = 1j * 2 * math.pi * f
        y_ss = (C_ @ np.linalg.solve(s * eye - A_, B_) + D_)[0, 0]
        y_cf = A.Y_mixed(s, lam, tau, g_n, K, c1)
        assert abs(y_ss - y_cf) / abs(y_cf) < 1e-2


def test_diffusive_quadrature_residues_positive():
    """Warburg diffusive-quadrature rungs are positive and log-spaced."""
    K = 6.8e6
    xi, r = EX.diffusive_quadrature(K, n_aux=30)
    assert len(xi) == len(r) == 30
    assert np.all(r > 0)
    assert np.all(np.diff(np.log(xi)) > 0)  # strictly log-increasing poles


# ---------------------------------------------------------------------------
# CAD-direct edge extraction (needs netgen.occ)
# ---------------------------------------------------------------------------

def test_cube_cad_topology_c1():
    """CAD-direct: a Box has exactly 12 edges and c_1 == -16(3L)/(pi mu)."""
    occ = pytest.importorskip("netgen.occ")
    from radia.levitation.mixed_galerkin import cad_edges as CE

    L = 5e-3
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(L, L, L))
    c1, L_total, n_edges = CE.cad_topology_c1(box, MU0)

    assert n_edges == 12
    assert L_total == pytest.approx(12 * L, rel=1e-9)
    closed = -16.0 * (3 * L) / (math.pi * MU0)
    assert c1 == pytest.approx(closed, rel=1e-6)


# ---------------------------------------------------------------------------
# End-to-end NGSolve cube alpha(s) (needs ngsolve; fast, ~1 s)
# ---------------------------------------------------------------------------

def test_cube_alpha_sweep_end_to_end():
    """Coarse Cu cube: V exact, PEC high-f limit, c_1 closed form.

    Geometry-exact quantities (V, c_1) are locked tightly; the bulk-Foster
    PEC limit alpha(1 GHz)/V is locked to a band (mesh/eigen-count
    dependent, approaches 1.0 from below).
    """
    pytest.importorskip("ngsolve")
    occ = pytest.importorskip("netgen.occ")
    from ngsolve import Mesh, TaskManager
    from radia.levitation.mixed_galerkin import cad_edges as CE

    L = 5e-3
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(L, L, L))
    for f in box.faces:
        f.name = "outer"
    ng = occ.OCCGeometry(box).GenerateMesh(maxh=L / 4)
    mesh = Mesh(ng)

    with TaskManager():
        lam, tau, g_n, V = A.bulk_foster_via_eigen(mesh, SIGMA_CU, MU0, n_eigen=40)

    # Geometry-exact invariants
    assert V == pytest.approx(L**3, rel=2e-2)          # 125 mm^3
    c1, _, n_edges = CE.cad_topology_c1(box, MU0)
    assert n_edges == 12
    assert c1 == pytest.approx(-16.0 * (3 * L) / (math.pi * MU0), rel=1e-6)

    # Foster spectrum sane
    assert len(lam) == 40
    assert np.all(lam > 0)
    assert np.all(g_n >= 0)

    # PEC high-frequency limit alpha(1 GHz)/V in [0.95, 1.02]
    K = A.K_SIBC_total(6 * L * L, SIGMA_CU, MU0)
    s = 1j * 2 * math.pi * 1e9
    a = A.alpha_from_Y(A.Y_mixed(s, lam, tau, g_n, K, c1), V, SIGMA_CU)
    assert 0.95 < (a / V).real < 1.02


# ---------------------------------------------------------------------------
# Matrix (multi-port) Y(s) + MIMO state-space  (frontier D)
# ---------------------------------------------------------------------------

def _synthetic_matrix_system(P=3):
    """Deterministic rank-1-per-pole matrix system for the algebra tests."""
    V = 1.25e-7
    lam = np.array([10.0, 40.0, 90.0, 160.0, 250.0])
    tau = MU0 * SIGMA_CU / lam
    Bproj = np.array([
        [1.0, 0.2, 0.1],
        [0.1, 1.0, 0.3],
        [0.2, 0.1, 1.0],
        [0.5, 0.4, 0.0],
        [0.0, 0.3, 0.7],
    ])[:, :P]
    G_n = np.array([SIGMA_CU * V * np.outer(Bproj[k], Bproj[k])
                    for k in range(len(lam))])
    K = A.K_SIBC_total(6e-4, SIGMA_CU, MU0)
    Kshape = np.array([[1.0, 0.10, 0.05],
                       [0.10, 0.90, 0.00],
                       [0.05, 0.00, 0.80]])[:P, :P]
    K_mat = K * Kshape
    c1 = -6.08e4
    C1shape = np.array([[1.0, 0.2, 0.0],
                        [0.2, 1.0, 0.1],
                        [0.0, 0.1, 1.0]])[:P, :P]
    C1_mat = c1 * C1shape
    return V, lam, tau, G_n, K_mat, C1_mat


def test_matrix_Y_reduces_to_scalar():
    """P=1 Y_matrix_mixed[0,0] == scalar Y_mixed; alpha likewise (exact)."""
    V = 1.25e-7
    lam = np.array([10.0, 40.0, 90.0])
    tau = MU0 * SIGMA_CU / lam
    g_n = np.array([SIGMA_CU * V * 0.30, SIGMA_CU * V * 0.10, SIGMA_CU * V * 0.05])
    K = A.K_SIBC_total(6e-4, SIGMA_CU, MU0)
    c1 = -6.08e4
    G_n = g_n.reshape(-1, 1, 1)
    for f in (1.0, 1e3, 1e6):
        s = 1j * 2 * math.pi * f
        y_mat = A.Y_matrix_mixed(s, lam, tau, G_n, [[K]], [[c1]])
        y_sca = A.Y_mixed(s, lam, tau, g_n, K, c1)
        assert abs(y_mat[0, 0] - y_sca) <= 1e-12 * abs(y_sca)
        a_mat = A.alpha_matrix_from_Y(y_mat, V, SIGMA_CU)
        a_sca = A.alpha_from_Y(y_sca, V, SIGMA_CU)
        assert abs(a_mat[0, 0] - a_sca) <= 1e-12 * abs(a_sca)


def test_matrix_residue_symmetric_psd():
    """Each Foster residue G_n is symmetric PSD; the SIBC K matrix is PSD."""
    _, _, _, G_n, K_mat, _ = _synthetic_matrix_system()
    for k in range(G_n.shape[0]):
        assert np.allclose(G_n[k], G_n[k].T)
        ev = np.linalg.eigvalsh(G_n[k])
        assert ev.min() > -1e-9 * max(abs(ev).max(), 1e-30)
    assert np.allclose(K_mat, K_mat.T)
    assert np.linalg.eigvalsh(K_mat).min() > 0.0


def test_mimo_state_space_matches_matrix_Y():
    """ss(A,B,C,D) transfer MATRIX reproduces Y_matrix_mixed(s) to <1.5%."""
    V, lam, tau, G_n, K_mat, C1_mat = _synthetic_matrix_system(P=3)
    A_, B_, C_, D_, n_f, n_w, n_i = EX.build_state_space_mimo(
        G_n, tau, V, SIGMA_CU, K_mat, C1_mat, n_warburg_rungs=40)
    assert n_f == G_n.shape[0]    # rank-1 residue -> one state per Foster pole
    assert n_i == 3               # full-rank C1 -> three integrator states
    eye = np.eye(A_.shape[0])
    for f in (1.0, 1e2, 1e4):
        s = 1j * 2 * math.pi * f
        H = C_ @ np.linalg.solve(s * eye - A_, B_) + D_
        Ycf = A.Y_matrix_mixed(s, lam, tau, G_n, K_mat, C1_mat)
        rel = np.linalg.norm(H - Ycf) / np.linalg.norm(Ycf)
        assert rel < 1.5e-2, f"f={f}: MIMO transfer rel err {rel*100:.2f}%"


def test_mimo_reduces_to_scalar_transfer():
    """P=1 MIMO transfer == scalar Y_mixed (different realization, same TF)."""
    V = 1.25e-7
    lam = np.array([10.0, 40.0, 90.0])
    tau = MU0 * SIGMA_CU / lam
    g_n = np.array([SIGMA_CU * V * 0.30, SIGMA_CU * V * 0.10, SIGMA_CU * V * 0.05])
    K = A.K_SIBC_total(6e-4, SIGMA_CU, MU0)
    c1 = -6.08e4
    G_n = g_n.reshape(-1, 1, 1)
    Am, Bm, Cm, Dm, *_ = EX.build_state_space_mimo(
        G_n, tau, V, SIGMA_CU, [[K]], [[c1]], n_warburg_rungs=40)
    eye = np.eye(Am.shape[0])
    for f in (1.0, 1e2, 1e4):
        s = 1j * 2 * math.pi * f
        H = (Cm @ np.linalg.solve(s * eye - Am, Bm) + Dm)[0, 0]
        y = A.Y_mixed(s, lam, tau, g_n, K, c1)
        assert abs(H - y) / abs(y) < 1.5e-2


# ---------------------------------------------------------------------------
# Per-face SIBC envelope + matrix edge term  (frontier C, CAD-direct)
# ---------------------------------------------------------------------------

def test_cad_topology_faces_cube():
    """A Box has 6 faces, each area L^2, unit outward normals; sum == 6 L^2."""
    pytest.importorskip("netgen.occ")
    import netgen.occ as occ
    from radia.levitation.mixed_galerkin import cad_edges as CE

    L = 5e-3
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(L, L, L))
    faces = CE.cad_topology_faces(box)
    assert len(faces) == 6
    for f in faces:
        assert f["area"] == pytest.approx(L * L, rel=1e-9)
        n = np.array(f["normal"])
        assert np.linalg.norm(n) == pytest.approx(1.0, rel=1e-9)
        assert np.max(np.abs(n)) == pytest.approx(1.0, abs=1e-9)  # axis-aligned
    assert sum(f["area"] for f in faces) == pytest.approx(6 * L * L, rel=1e-9)
    cen = np.array([L / 2, L / 2, L / 2])
    for f in faces:  # outward orientation
        assert np.dot(np.array(f["normal"]), np.array(f["center"]) - cen) > 0


def test_edge_moment_matrix_monopole_reduces_to_c1():
    """edge_moment_matrix([1]) reproduces the scalar cad_topology_c1 exactly."""
    pytest.importorskip("netgen.occ")
    import netgen.occ as occ
    from radia.levitation.mixed_galerkin import cad_edges as CE

    L = 5e-3
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(L, L, L))
    C1 = CE.edge_moment_matrix(box, [lambda p: 1.0], MU0)
    c1_scalar, _, _ = CE.cad_topology_c1(box, MU0)
    assert C1.shape == (1, 1)
    assert C1[0, 0] == pytest.approx(c1_scalar, rel=1e-9)
    assert C1[0, 0] == pytest.approx(-16.0 * (3 * L) / (math.pi * MU0), rel=1e-6)


def test_edge_moment_matrix_cube_isotropic():
    """Centered-coordinate edge tensor is symmetric, diagonal-isotropic (cube)."""
    pytest.importorskip("netgen.occ")
    import netgen.occ as occ
    from radia.levitation.mixed_galerkin import cad_edges as CE

    L = 5e-3
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(L, L, L))
    cen = L / 2.0
    drives = [(lambda i: (lambda p: p[i] - cen))(i) for i in range(3)]
    C1 = CE.edge_moment_matrix(box, drives, MU0)
    assert np.allclose(C1, C1.T)
    d = np.diag(C1)
    assert np.allclose(d, d[0], rtol=1e-6)        # three equal diagonal entries
    off = C1 - np.diag(d)
    assert np.max(np.abs(off)) < 1e-6 * abs(d[0])  # vanishing off-diagonals


# ---------------------------------------------------------------------------
# End-to-end NGSolve cube multipole admittance matrix  (frontier C + D)
# ---------------------------------------------------------------------------

def test_cube_admittance_matrix_end_to_end():
    """4-port cube {1, x-c, y-c, z-c}: matrix bulk + surface-moment K + edge.

    Locks the MEANINGFUL matrix-pipeline invariants:
      * the monopole port (0) IS the verified scalar path -- its bulk residue,
        surface K, and edge c_1 equal the scalar values, hence alpha_00(s)
        equals the validated scalar alpha(s);
      * the monopole alpha_00/V reaches the analytic perfect-conductor limit
        (~1) at high frequency -- a NON-trivial check (alpha_00/V is far from
        1, even negative, in the tail-dominated low band);
      * the residues G_n and the SIBC matrix K are symmetric positive-
        semidefinite;
      * the dipole block (ports 1-3) is isotropic for the symmetric cube.

    The scalar-model alpha = V - Y/sigma is NOT a meaningful tensor for the raw
    dipole ports (the V term dominates), so it is exercised only on the
    monopole port; see alpha.alpha_matrix_from_Y / bulk_foster_matrix_via_eigen
    scope notes.
    """
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import netgen.occ as occ
    from ngsolve import Mesh, TaskManager, x, y, z
    from radia.levitation.mixed_galerkin import cad_edges as CE

    L = 5e-3
    box = occ.Box(occ.Pnt(0, 0, 0), occ.Pnt(L, L, L))
    for f in box.faces:
        f.name = "outer"
    mesh = Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=L / 4))
    cen = L / 2.0
    drives = [1.0, x - cen, y - cen, z - cen]          # monopole + 3 dipoles

    with TaskManager():
        lam, tau, G_n, V = A.bulk_foster_matrix_via_eigen(
            mesh, SIGMA_CU, MU0, drives, n_eigen=40)
        K_mat = A.K_SIBC_matrix(mesh, drives, SIGMA_CU, MU0)
        # scalar reference path (monopole only)
        lam_s, tau_s, g_s, V_s = A.bulk_foster_via_eigen(
            mesh, SIGMA_CU, MU0, n_eigen=40)

    assert V == pytest.approx(L**3, rel=2e-2)

    # residues symmetric PSD
    for k in range(0, len(lam), 8):
        assert np.allclose(G_n[k], G_n[k].T)
        ev = np.linalg.eigvalsh(G_n[k])
        assert ev.min() > -1e-9 * max(abs(ev).max(), 1e-30)

    # K symmetric PSD; monopole port == scalar K_SIBC_total
    assert np.allclose(K_mat, K_mat.T)
    ke = np.linalg.eigvalsh(K_mat)
    assert ke.min() > -1e-6 * ke.max()
    S_cad = CE.cad_topology_total_area(box)
    assert K_mat[0, 0] == pytest.approx(A.K_SIBC_total(S_cad, SIGMA_CU, MU0),
                                        rel=2e-2)

    # monopole port reproduces the scalar bulk residue exactly
    assert np.allclose(G_n[:, 0, 0], g_s)

    # dipole block isotropic for the cube (mesh-asymmetry tolerance)
    Yb0 = G_n.sum(axis=0)                           # = Y_bulk(s=0)
    dip = np.diag(Yb0)[1:]
    assert (dip.max() - dip.min()) / dip.mean() < 0.20
    off = Yb0[1:, 1:] - np.diag(np.diag(Yb0[1:, 1:]))
    assert np.max(np.abs(off)) / dip.mean() < 0.10

    # monopole alpha_00(s): scalar identity + analytic PEC limit (non-trivial)
    drive_fns = [lambda p: 1.0] + [(lambda i: (lambda p: p[i] - cen))(i)
                                   for i in range(3)]
    C1_mat = CE.edge_moment_matrix(box, drive_fns, MU0)
    c1_s, _, _ = CE.cad_topology_c1(box, MU0)
    assert C1_mat[0, 0] == pytest.approx(c1_s, rel=1e-6)

    s = 1j * 2 * math.pi * 1e9
    Ymat = A.Y_matrix_mixed(s, lam, tau, G_n, K_mat, C1_mat)
    a00 = A.alpha_matrix_from_Y(Ymat, V, SIGMA_CU)[0, 0] / V
    # equals the scalar path at the same frequency
    K_s = A.K_SIBC_total(S_cad, SIGMA_CU, MU0)
    a_scalar = A.alpha_from_Y(A.Y_mixed(s, lam_s, tau_s, g_s, K_s, c1_s),
                              V_s, SIGMA_CU) / V_s
    assert abs(a00 - a_scalar) < 1e-9 * abs(a_scalar)
    assert 0.95 < a00.real < 1.02     # analytic perfect-conductor limit

    # MIMO LTI transfer matches Y to <1.5%
    Am, Bm, Cm, Dm, n_f, n_w, n_i = EX.build_state_space_mimo(
        G_n, tau, V, SIGMA_CU, K_mat, C1_mat, n_warburg_rungs=30)
    assert n_f == G_n.shape[0]
    eye = np.eye(Am.shape[0])
    for f in (1e3, 1e5):
        sf = 1j * 2 * math.pi * f
        H = Cm @ np.linalg.solve(sf * eye - Am, Bm) + Dm
        Ycf = A.Y_matrix_mixed(sf, lam, tau, G_n, K_mat, C1_mat)
        assert np.linalg.norm(H - Ycf) / np.linalg.norm(Ycf) < 1.5e-2


# ---------------------------------------------------------------------------
# Vector (HCurl) eddy-current Foster bulk  (frontier B)
# ---------------------------------------------------------------------------

def test_vector_bulk_cuboid_eddy_foster():
    """5x2x1 Cu box: HCurl curl-curl GEP eddy modes vs analytic interior-PEC TE.

    Locks the verified vector eddy Foster (de-Rham HCurl partner of the scalar
    H1 bulk): the per-direction leading tau matches the analytic interior-PEC
    TE mode tau = mu sigma / (pi^2 (1/La^2 + 1/Lb^2)); the spectrum is positive;
    the 3x3 residues are symmetric PSD; and a non-cubic box gives three
    DISTINCT leading tau with the shape-split ordering tau_z > tau_y > tau_x.

    Mesh-hungry (h ~ a/20 here): tau_x, tau_z lock tight (<8%); tau_y is
    resolution-sensitive at this h (<15%, ~-0.3% at h=a/28).  See
    mixed_galerkin/vector_bulk.py scope notes.
    """
    pytest.importorskip("ngsolve")
    pytest.importorskip("netgen.occ")
    import netgen.occ as occ
    from ngsolve import Mesh, TaskManager
    from radia.levitation.mixed_galerkin import bulk_foster_vector_via_eigen

    ax, ay, az = 5e-3, 2e-3, 1e-3
    h = 0.25e-3

    def tau_te(L1, L2):
        return MU0 * SIGMA_CU / (math.pi**2 * (1.0 / L1**2 + 1.0 / L2**2))
    tau_ref = {0: tau_te(ay, az), 1: tau_te(ax, az), 2: tau_te(ax, ay)}

    box = occ.Box(occ.Pnt(-ax/2, -ay/2, -az/2), occ.Pnt(ax/2, ay/2, az/2))
    box.mat("conductor").bc("conductor_surface")
    box.maxh = h
    mesh = Mesh(occ.OCCGeometry(box).GenerateMesh(maxh=h))

    with TaskManager():
        lam, tau_n, G_n, V, lead = bulk_foster_vector_via_eigen(
            mesh, SIGMA_CU, MU0, n_per_dir=12, order=2,
            conductor_bnd="conductor_surface")

    assert V == pytest.approx(ax * ay * az, rel=2e-2)
    assert np.all(lam > 0) and np.all(tau_n > 0)        # physical eddy spectrum

    # residues symmetric PSD (rank-1)
    for n in range(0, len(lam), 6):
        assert np.allclose(G_n[n], G_n[n].T)
        ev = np.linalg.eigvalsh(G_n[n])
        assert ev.min() > -1e-9 * max(abs(ev).max(), 1e-30)

    tx, ty, tz = lead[0], lead[1], lead[2]
    # shape split: three distinct leading tau, ordered by cross-section size
    assert tz > ty > tx
    # leading tau vs analytic interior-PEC TE mode
    assert abs(tx - tau_ref[0]) / tau_ref[0] < 0.08    # tight (well resolved)
    assert abs(tz - tau_ref[2]) / tau_ref[2] < 0.08    # tight (well resolved)
    assert abs(ty - tau_ref[1]) / tau_ref[1] < 0.15    # resolution-sensitive at this h


# ---------------------------------------------------------------------------
# Physical-tensor ROM: passive Foster fit of the EXTERIOR-MATCHED alpha(s)
# (rom_fit) -- verified against the analytic sphere Stoll spectrum (pure numpy)
# ---------------------------------------------------------------------------

def _sphere_alpha_causal(omega, a=5e-3, sigma=SIGMA_CU):
    """Causal-convention sphere polarizability alpha(j omega) (LHP poles).

    Landau-Lifshitz: m = (4 pi a^3/mu0) G(x) B0, x=(1+j) a/delta, and the
    FEM-convention alpha = 4 pi a^3 conj(G).  Foster form:
      alpha(s)/(4 pi a^3) = -1/2 + sum_n (3/(n pi)^2)/(1 + s tau_n),
      tau_n = mu0 sigma a^2 / (n pi)^2  (the physical Stoll spectrum).
    """
    delta = math.sqrt(2.0 / (omega * MU0 * sigma))
    x = (1 + 1j) * a / delta
    e = cmath.exp(2j * x)
    cot = 1j * (e + 1.0) / (e - 1.0)
    G = -0.5 * (1.0 - 3.0 / x**2 + 3.0 / x * cot)
    return 4 * math.pi * a**3 * G.conjugate()


def test_rom_fit_sphere_stoll_spectrum():
    """AAA+NNLS passive Foster ROM recovers the analytic sphere Stoll spectrum.

    Locks: the dominant AAA poles == physical tau_n = mu0 sigma a^2/(n pi)^2 to
    ~0%, the passive fit reproduces alpha(s) to <0.5% over 1 Hz..1 GHz, the ROM
    is passive (g_k>=0) and stable, the feedthrough D == the perfect-conductor
    limit -2 pi a^3, and the state-space transfer equals eval(s).
    """
    a = 5e-3
    n = np.arange(1, 4)
    tau_anal = MU0 * SIGMA_CU * a**2 / (n * math.pi) ** 2

    f = np.logspace(0, 9, 400)
    s = 1j * 2 * math.pi * f
    Y = np.array([_sphere_alpha_causal(2 * math.pi * fi, a) for fi in f])

    rom = RF.passive_foster_fit(s, Y, n_filler=20)

    # passive + stable + accurate
    assert np.all(rom.g_n >= 0)
    assert rom.band_fit_relerr < 5e-3
    # feedthrough = perfect-conductor flux exclusion -2 pi a^3
    assert rom.alpha_inf == pytest.approx(-2 * math.pi * a**3, rel=1e-2)

    # AAA-discovered dominant poles == physical Stoll tau_n (top 2 to <1%)
    dt = np.sort(rom.dominant_tau)[::-1]
    assert len(dt) >= 2
    for k in range(2):
        assert abs(dt[k] - tau_anal[k]) / tau_anal[k] < 1e-2

    # state-space transfer == eval (exact realization)
    Amat, Bmat, Cmat, Dmat = rom.state_space()
    Ac = Amat.astype(complex)
    eye = np.eye(rom.n_states)
    sp = 1j * 2 * math.pi * np.array([1e2, 1e3, 1e4, 1e5, 1e6])
    H = np.array([(Cmat @ np.linalg.solve(si * eye - Ac, Bmat))[0, 0] + Dmat[0, 0]
                  for si in sp])
    assert np.allclose(H, rom.eval(sp), rtol=1e-9)
    # stable: all poles in the LHP
    assert np.all(np.diag(Amat) < 0)


def test_rom_fit_diagonal_tensor_mimo():
    """Three scalar ROMs -> a diagonal MIMO LTI (principal-axis tensor).

    Locks: per-axis transfer reproduces each scalar ROM, the off-diagonal is
    identically zero, and the block-diagonal state count is the sum.
    """
    a = 5e-3
    f = np.logspace(0, 9, 300)
    s = 1j * 2 * math.pi * f
    Y = np.array([_sphere_alpha_causal(2 * math.pi * fi, a) for fi in f])
    roms = [RF.passive_foster_fit(s, Y, n_filler=12) for _ in range(3)]

    A_t, B_t, C_t, D_t, ns = RF.diagonal_tensor_state_space(roms)
    assert ns == sum(r.n_states for r in roms)
    assert A_t.shape == (ns, ns) and B_t.shape == (ns, 3) and C_t.shape == (3, ns)

    eye = np.eye(ns)
    sp = 1j * 2 * math.pi * np.array([1e3, 1e4, 1e5])
    Ht = np.array([C_t @ np.linalg.solve(si * eye - A_t.astype(complex), B_t) + D_t
                   for si in sp])
    for p in range(3):
        assert np.allclose(Ht[:, p, p], roms[p].eval(sp), rtol=1e-9)
    assert np.allclose(Ht[:, 0, 1], 0.0)               # no cross coupling
