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
    from ngsolve import H1, BilinearForm, grad, dx, GridFunction, Integrate
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

    # Project the constant 1 onto the eigenbasis: <1, phi_n>_M
    one_gfu = GridFunction(fes)
    one_gfu.Set(1.0)
    one_vec_full = np.array(one_gfu.vec.FV().NumPy())
    one_vec_free = one_vec_full[free]
    M_one = M_free.dot(one_vec_free)
    ip = vecs.T @ M_one

    V = float(Integrate(1.0, mesh))
    b_n = ip / math.sqrt(V)
    tau_n = mu * sigma / lam
    g_n = sigma * V * b_n**2

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
    """Y_R(s) = Y_bulk_Foster(s) + K_SIBC / sqrt(s) + c_1 / s.

    Bulk Foster sum Y_bulk(s) = sum_n g_n / (1 + s tau_n) covers
    the low-mid frequency response; the Mellin tail K_SIBC/sqrt(s)
    + c_1/s captures the high-f surface and edge correction beyond
    the truncated bulk sum.
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
            "    from radia.levitation.mixed_galerkin import cad_edges\n"
            "    c1, L_total, n = cad_edges.cad_topology_c1(shape, mu)\n"
            "(pass the OCC primitive / loaded STEP shape, not the mesh)."
        ) from exc
    edge_data = mesh_edge_dihedrals(mesh)
    edges = [(e["length"], e["dihedral"]) for e in edge_data]
    return S_total, edges
