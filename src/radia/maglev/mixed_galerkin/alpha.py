"""alpha.py -- general-shape polarizability tensor via Mixed Galerkin.

For any conductor shape specified as a `.vol` mesh, compute:

  1. NGSolve scalar diffusion eigenmodes on the conductor volume
     (Foster spectrum: poles tau_n and residues g_n).
  2. Edge dihedral angles -> W(alpha) = (4/pi) cot(alpha/2) per edge,
     summed into c_1 (Wang-Lavers-Zhou 1992 wedge SIBC pinned to the
     cuboid Mellin anchor; Sugahara-Nagamine-Hane 2026 Paper 1 SIII).
  3. Surface area S_total -> K_SIBC = S sqrt(sigma/mu)
     (Mellin leading term).
  4. Schur composition Y_R(s) = Y_bulk_Foster(s) + K_SIBC/sqrt(s)
     + c_1/s.
  5. Polarizability alpha(s) = V - Y(s)/sigma.

Geometry assumption: convex polyhedron with flat faces (cuboid, hex
bar, chamfered bar, L-shape, multi-step bar, ...).  Edge dihedrals
can be arbitrary (sharp); the W(alpha) formula handles the angular
dependence in closed form.  Smooth curved bodies (sphere, ellipsoid)
need a different gamma_k tower (see Paper 1 SIII, smooth-surface
HOIBC).
"""
from __future__ import annotations

import cmath
import math

import numpy as np


# ---------------------------------------------------------------------------
# Bulk Foster spectrum via NGSolve scalar eigenmodes
# ---------------------------------------------------------------------------


def _dirichlet_eigenmodes(mesh, n_eigen: int, dirichlet_label: str):
    """Shared Dirichlet-Laplacian eigensolve for the bulk Foster spectrum.

    Solves (-Laplacian) phi_n = lam_n phi_n on the conductor with phi_n = 0 on
    `dirichlet_label`.  Returns the ascending eigenvalues, the M-normalized
    eigenvectors restricted to the free dofs, the full mass matrix, the
    free-dof boolean mask, the H1 space, and the volume V. The full matrix is
    required because driving functions need not vanish on the boundary. Both the
    scalar
    (`bulk_foster_via_eigen`) and the matrix / multi-port
    (`bulk_foster_matrix_via_eigen`) drivers share this eigensolve so the two
    return bit-identical spectra.
    """
    from ngsolve import H1, BilinearForm, grad, dx, Integrate
    import scipy.sparse as sp
    from scipy.sparse.linalg import eigsh

    fes = H1(mesh, order=2, dirichlet=dirichlet_label)
    u, v = fes.TnT()
    a = BilinearForm(fes, symmetric=True)
    a += grad(u) * grad(v) * dx
    a.Assemble()
    m = BilinearForm(fes, symmetric=True)
    m += u * v * dx
    m.Assemble()

    rows_a, cols_a, vals_a = a.mat.COO()
    rows_m, cols_m, vals_m = m.mat.COO()
    ndof = fes.ndof
    A = sp.csr_matrix(
        (np.asarray(vals_a), (np.asarray(rows_a), np.asarray(cols_a))),
        shape=(ndof, ndof),
    )
    M = sp.csr_matrix(
        (np.asarray(vals_m), (np.asarray(rows_m), np.asarray(cols_m))),
        shape=(ndof, ndof),
    )

    free = np.array([fes.FreeDofs()[i] for i in range(ndof)], dtype=bool)
    A_free = A[free][:, free]
    M_free = M[free][:, free]

    k_eig = min(n_eigen, A_free.shape[0] - 2)
    lam, vecs = eigsh(A_free, k=k_eig, M=M_free, sigma=0, which="LM")
    order = np.argsort(lam)
    lam = lam[order]
    vecs = vecs[:, order]

    # M-normalize each eigenvector
    for k in range(vecs.shape[1]):
        norm = math.sqrt(abs(vecs[:, k] @ M_free.dot(vecs[:, k])))
        if norm > 1e-30:
            vecs[:, k] /= norm

    V = float(Integrate(1.0, mesh))
    return lam, vecs, M, free, fes, V


def _project_drive(fes, M, free, vecs, cf, V):
    """Project a driving CoefficientFunction onto the eigenbasis.

    Returns b[n] = <cf, phi_n>_M / sqrt(V) (length n_eigen). The eigenmode
    is zero on Dirichlet DOFs, but the drive is not, so the product uses the
    free rows and all mass-matrix columns.
    """
    from ngsolve import GridFunction
    gfu = GridFunction(fes)
    gfu.Set(cf)
    full = np.array(gfu.vec.FV().NumPy())
    ip = vecs.T @ (M[free, :] @ full)
    return ip / math.sqrt(V)


def bulk_foster_matrix_via_eigen(mesh, sigma: float, mu: float, drive_cfs,
                                 n_eigen: int = 60,
                                 dirichlet_label: str = "outer"):
    """Matrix (multi-port) bulk Foster spectrum from scalar diffusion eigenmodes.

    Generalizes `bulk_foster_via_eigen` to P driving functions
    {f_0, ..., f_{P-1}} (the "ports"), each projected onto the SAME shared
    Dirichlet-Laplacian eigenbasis.  Returns a per-pole residue MATRIX:

        Y_bulk(s)_{pq} = sum_n G_n[n]_{pq} / (1 + s tau_n)

    with G_n[n] = sigma V b_n b_n^T  (symmetric positive-semidefinite, rank 1
    per mode), b_n[p] = <f_p, phi_n>_M / sqrt(V), tau_n = mu sigma / lam_n.

    This is the Foster-eigenbasis realization of the matrix-form CLN
    (Matsuo 2017/2018c multi-port Kameari: the modal amplitude carries one
    column per port and lambda becomes an N_port x N_port matrix; see
    radia_mcp.mor mor_cln_multiport "multiport_theory").  With drive_cfs the
    set {1, x-x_c, y-y_c, z-z_c, ...} the matrix Y(s)_{pq} is a MULTIPOLE
    expansion of the conductor's scalar eddy response: the monopole port 1
    couples to a uniform external field, the dipole ports x, y, z to field
    gradients.

    Scope note (honest): this is the multipole eddy-ADMITTANCE matrix of the
    SCALAR diffusion model.  It is NOT the physical 3x3 VECTOR (HCurl) eddy
    polarizability tensor -- the alpha = V - Y/sigma relation is monopole-
    specific (for raw dipole drives the V term dominates, see
    `alpha_matrix_from_Y`), and the scalar A_z model has a single field
    component.  The verified physical vector polarizability TENSOR (transverse
    m=1 components + triaxial shape anisotropy) is the full 3D HCurl solve in
    docs/maglev/demos/ellipsoid/ellipsoid_alpha_tensor_3d.py.
    P=1 with drive_cfs=[1.0] reproduces `bulk_foster_via_eigen` exactly.

    Parameters
    ----------
    drive_cfs : sequence
        P driving functions, each acceptable to GridFunction.Set (an NGSolve
        CoefficientFunction or a scalar such as 1.0).

    Returns
    -------
    lam : ndarray (n_eigen,)
    tau_n : ndarray (n_eigen,)
    G_n : ndarray (n_eigen, P, P)
        Symmetric PSD rank-1 residue matrix per Foster pole.
    V : float
    """
    from ngsolve import CoefficientFunction
    lam, vecs, M, free, fes, V = _dirichlet_eigenmodes(
        mesh, n_eigen, dirichlet_label)
    P = len(drive_cfs)
    n_mode = len(lam)
    B = np.zeros((n_mode, P))
    for p, cf in enumerate(drive_cfs):
        B[:, p] = _project_drive(fes, M, free, vecs,
                                 CoefficientFunction(cf), V)
    tau_n = mu * sigma / lam
    G_n = np.zeros((n_mode, P, P))
    for n in range(n_mode):
        G_n[n] = sigma * V * np.outer(B[n], B[n])
    return lam, tau_n, G_n, V


def bulk_foster_via_eigen(mesh, sigma: float, mu: float, n_eigen: int = 60,
                          dirichlet_label: str = "outer"):
    """Compute bulk Foster spectrum from scalar diffusion eigenmodes.

    Solves (-Laplacian) phi_n = lam_n phi_n on conductor V with
    phi_n = 0 on dV.  Projects the constant 1 (driving function for
    the scalar A_z model) onto the eigenbasis and returns the Foster
    representation:

        Y_bulk(s) = sum_n g_n / (1 + s tau_n)

    where g_n = sigma V b_n^2, tau_n = mu sigma / lam_n,
    b_n = <1, phi_n>_M / sqrt(V) (M-normalized eigenvectors).

    This is the single-port (monopole) special case of
    `bulk_foster_matrix_via_eigen` and delegates to it.

    Parameters
    ----------
    mesh : ngsolve.Mesh
        Volumetric mesh of the conductor.
    sigma : float
        Conductivity (S/m).
    mu : float
        Permeability (H/m).
    n_eigen : int
        Number of Dirichlet eigenmodes to compute (default 60).
    dirichlet_label : str
        Boundary label where v=0 is imposed (default "outer").

    Returns
    -------
    lam : ndarray (n_eigen,)
        Eigenvalues, ascending.
    tau_n : ndarray (n_eigen,)
        Foster pole time constants.
    g_n : ndarray (n_eigen,)
        Foster pole residues (units of admittance).
    V : float
        Conductor volume.
    """
    lam, tau_n, G_n, V = bulk_foster_matrix_via_eigen(
        mesh, sigma, mu, [1.0], n_eigen=n_eigen,
        dirichlet_label=dirichlet_label)
    g_n = G_n[:, 0, 0]
    return lam, tau_n, g_n, V


# ---------------------------------------------------------------------------
# Polyhedral wedge function W(alpha) and edge correction c_1
# ---------------------------------------------------------------------------


def wedge_function(alpha: float) -> float:
    """W(alpha) = (4/pi) cot(alpha/2): polyhedral wedge function.

    Closed form from Wang-Lavers-Zhou 1992 (eqs 12-13) pinned to the
    cuboid Mellin anchor W(pi/2) = 4/pi (Sugahara-Nagamine-Hane 2026
    Paper 1, eq wedge_function_closed_form).

    Special values:
        W(pi/2)  = 4/pi      (right-angle cube edge)
        W(pi)    = 0         (flat: no edge contribution)
        W(0+)    -> infty    (knife-edge limit, not physical for solid)
        W(3*pi/2) = -4/pi    (re-entrant 270 deg, e.g. L-shape notch)

    Parameters
    ----------
    alpha : float
        Interior dihedral angle (radians).

    Returns
    -------
    W : float
        Wedge function value.
    """
    return (4.0 / math.pi) * (1.0 / math.tan(alpha / 2.0))


def c1_polyhedral(edges, mu: float) -> float:
    """c_1 = -(1/mu) sum_e L_e W(alpha_e) for a list of (length, dihedral) pairs.

    Parameters
    ----------
    edges : list of (L_e, alpha_e) tuples
        Per-edge length and interior dihedral angle.
    mu : float
        Permeability (H/m).

    Returns
    -------
    c_1 : float
        Mellin asymptote coefficient at order 1/s.
    """
    s = 0.0
    for L_e, alpha_e in edges:
        s += L_e * wedge_function(alpha_e)
    return -s / mu


def K_SIBC_total(area_total: float, sigma: float, mu: float) -> float:
    """K_SIBC = S sqrt(sigma/mu): planar SIBC anchor (Mellin c_0 coefficient)."""
    return area_total * math.sqrt(sigma / mu)


# ---------------------------------------------------------------------------
# Mixed Galerkin admittance Y(s) and polarizability alpha(s)
# ---------------------------------------------------------------------------


def Y_mixed(s, lam, tau, g_n, K_SIBC, c1):
    """Additive composition Y(s) = Y_bulk_Foster(s) + K_SIBC / sqrt(s) + c_1 / s.

    Bulk Foster sum Y_bulk(s) = sum_n g_n / (1 + s tau_n) covers
    The terms are added without a coupling block. For a projected box model
    in which the crossover is selected by a Schur complement, use
    :class:`schur.BoxMixedGalerkin`.
    """
    s_complex = complex(s) if not isinstance(s, complex) else s
    Y_bulk = np.sum(g_n / (1.0 + s_complex * tau))
    Y_planar = K_SIBC / cmath.sqrt(s_complex)
    Y_edge = c1 / s_complex
    return Y_bulk + Y_planar + Y_edge


def alpha_from_Y(Y, V, sigma):
    """alpha(s) = V - Y(s)/sigma  (polarizability from admittance)."""
    return V - Y / sigma


# ---------------------------------------------------------------------------
# Matrix (multi-port) admittance Y(s)_{pq} and polarizability alpha(s)_{pq}
# ---------------------------------------------------------------------------


def surface_moment_matrix(mesh, drive_cfs):
    """Surface moment matrix K_geom[p,q] = integral_{dOmega} f_p f_q dS.

    The matrix generalization of the wetted area S that sets the SIBC
    orthogonal-residual tail (cln_sibc_orthogonal "math": the residual modes
    sum to K_SIBC / sqrt(s) with the area as the geometric factor; for a
    multi-drive projection the area becomes the surface second moment of the
    drives).  For the monopole drive f=1 this is the total boundary area S
    (so K reduces to K_SIBC_total); for centered coordinates {x,y,z} it is
    the surface second-moment tensor.

    Computed by exact boundary integration over the whole surface; for a
    polyhedron this equals the per-face (partition-of-unity) sum
    sum_F integral_F f_p f_q dS because the flat faces are disjoint.

    Returns a symmetric (P, P) ndarray.
    """
    from ngsolve import Integrate, CoefficientFunction
    bnd = mesh.Boundaries(".*")
    cfs = [CoefficientFunction(c) for c in drive_cfs]
    P = len(cfs)
    K = np.zeros((P, P))
    for p in range(P):
        for q in range(p, P):
            val = float(Integrate(cfs[p] * cfs[q], mesh,
                                  definedon=bnd, order=4))
            K[p, q] = val
            K[q, p] = val
    return K


def K_SIBC_matrix(mesh, drive_cfs, sigma, mu):
    """Matrix (multi-port) SIBC tail coefficient K_mat = sqrt(sigma/mu) K_geom.

    Y_matrix_mixed(s) carries a K_mat / sqrt(s) term.  Reduces to the scalar
    K_SIBC_total(S, sigma, mu) for the single monopole drive f=1.
    """
    return math.sqrt(sigma / mu) * surface_moment_matrix(mesh, drive_cfs)


def Y_matrix_mixed(s, lam, tau, G_n, K_mat, C1_mat):
    """Matrix additive admittance Y(s)_{pq} (P x P complex).

        Y(s) = sum_n G_n[n] / (1 + s tau_n) + K_mat / sqrt(s) + C1_mat / s

    This is an additive form without a coupling block. The projected
    multi-port box form is :class:`schur.BoxMixedGalerkin`.

    G_n is the (n_eigen, P, P) residue tensor from
    bulk_foster_matrix_via_eigen; K_mat from K_SIBC_matrix; C1_mat the
    per-edge moment matrix (cad_edges.edge_moment_matrix), all P x P.  The
    P=1 case reproduces Y_mixed(s, lam, tau, g_n, K_SIBC, c1).
    """
    s_complex = complex(s) if not isinstance(s, complex) else s
    P = G_n.shape[1]
    Y = np.zeros((P, P), dtype=complex)
    denom = 1.0 + s_complex * np.asarray(tau)
    for n in range(G_n.shape[0]):
        Y += G_n[n] / denom[n]
    Y += np.asarray(K_mat) / cmath.sqrt(s_complex)
    Y += np.asarray(C1_mat) / s_complex
    return Y


def alpha_matrix_from_Y(Y, V, sigma):
    """alpha(s)_{pq} = V delta_{pq} - Y(s)_{pq} / sigma  (P x P).

    The matrix generalization of `alpha_from_Y`.  It is physically the
    polarizability ONLY for monopole-normalized ports (||f_p||_{L2(V)} =
    sqrt(V), as the constant 1 is): then Y is O(sigma V) and alpha spans
    [0, V].  For raw DIPOLE / coordinate drives the bulk admittance is
    O(sigma V L^2) (Y_ii(0) = sigma * integral (x_i-c)^2 dV = sigma V L^2/12
    for a cube), so the V*I term dominates and alpha ~ V*I at every frequency
    -- in that case work with the admittance matrix Y(s) directly (its eddy
    roll-off is the physical content), not this alpha.  The physical vector
    polarizability tensor is the HCurl solve (see bulk_foster_matrix_via_eigen
    scope note).
    """
    P = Y.shape[0]
    return V * np.eye(P) - np.asarray(Y) / sigma


# ---------------------------------------------------------------------------
# Geometry analysis from mesh
# ---------------------------------------------------------------------------


def measure_total_area(mesh) -> float:
    """Total boundary area S_total of the conductor (mesh-derived, exact)."""
    from ngsolve import Integrate, CoefficientFunction
    return float(Integrate(CoefficientFunction(1), mesh.Boundaries(".*"), order=2))


def measure_total_area_and_edges(mesh):
    """Return (S_total, edge_list=[(length, dihedral), ...]).

    S_total is the boundary area (always mesh-derived, exact).  The edge
    dihedrals come from the OPTIONAL `radia.netgen_mesh_curvature`
    extractor; if the `radia` package is not importable this raises and
    points you at the dependency-free, mesh-independent CAD-direct route
    `mixed_galerkin.cad_edges.cad_topology_edges(shape)` /
    `cad_topology_c1(shape, mu)`.

    There is deliberately NO silent "assume pi/2" fallback (No-Fallbacks
    policy): a guessed dihedral happens to be right only for an all-right-
    angle body (cube) and silently yields a WRONG c_1 for any chamfered /
    L-shape / oblique-edge geometry.  Fail loud instead.
    """
    S_total = measure_total_area(mesh)
    try:
        from radia.netgen_mesh_curvature import mesh_edge_dihedrals
    except Exception as exc:  # radia not installed / module absent
        raise ImportError(
            "measure_total_area_and_edges() needs the optional `radia` "
            "package (radia.netgen_mesh_curvature) for mesh-derived edge "
            "dihedrals; it is not importable.  Use the dependency-free "
            "CAD-direct route instead -- it is mesh-independent and exact "
            "for polyhedra:\n"
            "    from radia.maglev.mixed_galerkin import cad_edges\n"
            "    c1, L_total, n = cad_edges.cad_topology_c1(shape, mu)\n"
            "(pass the OCC primitive / loaded STEP shape, not the mesh)."
        ) from exc
    edge_data = mesh_edge_dihedrals(mesh)
    edges = [(e["length"], e["dihedral"]) for e in edge_data]
    return S_total, edges
