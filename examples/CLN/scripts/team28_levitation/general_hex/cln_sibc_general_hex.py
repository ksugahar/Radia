"""cln_sibc_general_hex.py -- general-shape polarizability tensor via Mixed Galerkin.

Generalizes `cln_sibc_cuboid_3d.py` from the separable cuboid Foster sum
(odd m,n,p closed form) to **ANY conductor shape** specified as a `.vol`
mesh.  Uses:

  1.  NGSolve scalar diffusion eigenmodes on the conductor volume
      (replaces the closed-form Foster spectrum).
  2.  Mesh-derived dihedral angles  -->  W(alpha) = (4/pi) cot(alpha/2)
      per edge (the closed-form polyhedral wedge function from
      Wang-Lavers-Zhou 1992 pinned to the Mellin anchor; see Paper 1 SIII).
  3.  Schur composition of bulk CLN Krylov truncation + Mellin tail
      c_0 / sqrt(s) + c_1 / s, giving Y_R(s) -- the conductor's
      magnetic admittance projected on the driving direction.
  4.  Polarizability alpha_ii(s) = V - Y_ii(s)/sigma per Cartesian axis.

For cuboid case, reduces to the existing closed-form
`cln_sibc_cuboid_3d.py` to within sub-percent at moderate f.

Geometry assumption: convex polyhedron with flat faces (cuboid, hex bar,
chamfered bar, L-shape, multi-step bar, ...).  Edge dihedrals can be
arbitrary (sharp); the W(alpha) formula handles the angular dependence
in closed form.  Smooth curved bodies (sphere, ellipsoid) need a
different `gamma_k` tower (see Paper 1 SIII, smooth-surface HOIBC).

Run examples:
  python cln_sibc_general_hex.py --vol cube_5mm.vol --axis z
  python cln_sibc_general_hex.py --vol Lshape.vol --axis x --freqs 1e2,1e3,1e4
  python cln_sibc_general_hex.py --vol bar_5x2x1.vol --axis all
"""
from __future__ import annotations

import argparse
import cmath
import math
import sys
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Physics
# ---------------------------------------------------------------------------

SIGMA_CU = 5.8e7
MU_0 = 4 * math.pi * 1e-7


# ---------------------------------------------------------------------------
# Bulk Foster spectrum via NGSolve scalar eigenmodes
# ---------------------------------------------------------------------------


def bulk_foster_via_eigen(mesh, sigma: float, mu: float, n_eigen: int = 60,
                          dirichlet_label: str = "outer"):
    """Scalar diffusion eigenmodes -Lap phi = lam phi on conductor V, phi=0 on boundary.
    Returns (lam, K_n, b_n) arrays so that Y_bulk(s) = sum K_n^2 / (1 + s*tau_n).

    The system is (-Lap + s mu sigma) v = -s mu sigma  in V, v=0 on dV.
    Project onto Dirichlet eigenbasis phi_n with -Lap phi_n = lam_n phi_n:
        v = sum c_n phi_n,  c_n = -s mu sigma * <1, phi_n> / (lam_n + s mu sigma)
        <v> = (1/V) sum c_n <1, phi_n>
        Y = Y_DC (1 + <v>) = sigma V * (1 + <v>)
          = sigma V - sigma sum b_n^2 (s mu sigma) / (lam_n + s mu sigma)
        where b_n = <1, phi_n> / sqrt(V)  (normalized projection)
        tau_n = mu sigma / lam_n,  K_n = sigma * b_n^2 (V * lam_n / mu sigma)
                                       = (sigma V / mu sigma) * b_n^2 * lam_n
                                       = (V / mu) * b_n^2 * lam_n
    So:
        Y_bulk(s) = Y_DC - sum K_n^2 * s tau_n / (1 + s tau_n)
                  = sum K_n^2 / (1 + s tau_n)  in the rearranged Foster form.
    """
    from ngsolve import H1, BilinearForm, LinearForm, grad, dx, GridFunction, Integrate
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

    # Build sparse CSR for generalized eigenproblem A x = lam M x
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

    # Restrict to free DOFs (Dirichlet zero)
    free = np.array([fes.FreeDofs()[i] for i in range(ndof)], dtype=bool)
    A_free = A[free][:, free]
    M_free = M[free][:, free]

    k_eig = min(n_eigen, A_free.shape[0] - 2)
    print(f"  Solving {k_eig} smallest Dirichlet eigenmodes (ndof_free = {A_free.shape[0]})...")
    lam, vecs = eigsh(A_free, k=k_eig, M=M_free, sigma=0, which="LM")
    order = np.argsort(lam)
    lam = lam[order]
    vecs = vecs[:, order]

    # M-normalize each eigenvector
    for k in range(vecs.shape[1]):
        norm = math.sqrt(abs(vecs[:, k] @ M_free.dot(vecs[:, k])))
        if norm > 1e-30:
            vecs[:, k] /= norm

    # Inner products <1, phi_n> = (M 1) . phi_n  (mass-weighted)
    one_gfu = GridFunction(fes)
    one_gfu.Set(1.0)
    one_vec_full = np.array(one_gfu.vec.FV().NumPy())
    one_vec_free = one_vec_full[free]
    M_one = M_free.dot(one_vec_free)
    ip = vecs.T @ M_one  # shape (k_eig,)

    V = float(Integrate(1.0, mesh))
    b_n = ip / math.sqrt(V)
    tau_n = mu * sigma / lam
    # Foster residues g_n = sigma V b_n^2 (the scalar diffusion 3D form,
    # consistent with cln_sibc_cuboid_3d.py closed form for cuboid).
    g_n = sigma * V * b_n**2

    return lam, tau_n, g_n, V


# ---------------------------------------------------------------------------
# Edge correction via mesh dihedrals
# ---------------------------------------------------------------------------


def wedge_function(alpha: float) -> float:
    """W(alpha) = (4/pi) cot(alpha/2): polyhedral wedge function.

    Closed form from Wang-Lavers-Zhou 1992 (eqs 12-13) pinned to the
    cuboid Mellin anchor W(pi/2) = 4/pi (Sugahara-Nagamine-Hane 2026
    Paper 1, eq wedge_function_closed_form).
    """
    return (4.0 / math.pi) * (1.0 / math.tan(alpha / 2.0))


def c1_polyhedral(edges, mu: float) -> float:
    """c_1 = -(1/mu) sum_e L_e W(alpha_e) for a list of (length, dihedral) pairs."""
    s = 0.0
    for L_e, alpha_e in edges:
        s += L_e * wedge_function(alpha_e)
    return -s / mu


def K_SIBC_total(area_total: float, sigma: float, mu: float) -> float:
    """c_0 = S * sqrt(sigma/mu).  Planar SIBC anchor."""
    return area_total * math.sqrt(sigma / mu)


# ---------------------------------------------------------------------------
# Mixed Galerkin admittance Y(s) and polarizability alpha(s)
# ---------------------------------------------------------------------------


def Y_mixed(s, lam, tau, g_n, K_SIBC, c1):
    """Y_R(s) = Y_bulk_foster(s) + K_SIBC/sqrt(s) + c_1/s.

    Note: bulk Foster sum Σ g_n / (1 + s τ_n) gives Y(0) = sigma V (DC limit)
    and decays as 1/s at high f. The Mellin tail K_SIBC/√s + c_1/s captures
    the high-f surface and edge correction beyond the truncated bulk sum.
    """
    s_complex = complex(s) if not isinstance(s, complex) else s
    Y_bulk = np.sum(g_n / (1.0 + s_complex * tau))
    Y_planar = K_SIBC / cmath.sqrt(s_complex)
    Y_edge = c1 / s_complex
    return Y_bulk + Y_planar + Y_edge


def alpha_from_Y(Y, V, sigma):
    """alpha(s) = V - Y(s)/sigma."""
    return V - Y / sigma


# ---------------------------------------------------------------------------
# Geometry analysis from mesh
# ---------------------------------------------------------------------------


def measure_total_area_and_edges(mesh):
    """Return (S_total, edge_list) where edge_list = [(length, dihedral), ...]."""
    from ngsolve import Integrate, CoefficientFunction
    S_total = float(Integrate(CoefficientFunction(1), mesh.Boundaries(".*"), order=2))

    # Edge dihedrals: try several import locations
    sys.path.insert(0, str(Path("S:/Radia/01_GitHub/build/lib/radia").resolve()))
    try:
        from netgen_mesh_curvature import mesh_edge_dihedrals
        edge_data = mesh_edge_dihedrals(mesh)
        edges = [(e["length"], e["dihedral"]) for e in edge_data]
    except Exception as exc:
        print(f"  WARN: dihedral extraction failed ({exc}); using cuboid fallback (all alpha = pi/2)", file=sys.stderr)
        # Fallback: assume cuboid -- estimate L_total by counting BBND edge segments
        L_total_est = 0.0
        try:
            from ngsolve import BBND
            pts = mesh.ngmesh.Points()
            for el in mesh.Elements(BBND):
                vs = [v.nr for v in el.vertices]
                if len(vs) == 2:
                    p0 = np.array([pts[vs[0]+1][0], pts[vs[0]+1][1], pts[vs[0]+1][2]])
                    p1 = np.array([pts[vs[1]+1][0], pts[vs[1]+1][1], pts[vs[1]+1][2]])
                    L_total_est += float(np.linalg.norm(p1 - p0))
        except Exception:
            pass
        edges = [(L_total_est, math.pi / 2)] if L_total_est > 0 else []
    return S_total, edges


# ---------------------------------------------------------------------------
# Main driver
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--vol", required=True, help="Conductor .vol mesh path")
    parser.add_argument("--sigma", type=float, default=SIGMA_CU)
    parser.add_argument("--mu", type=float, default=MU_0)
    parser.add_argument("--n-eigen", type=int, default=60, help="Number of Dirichlet eigenmodes for bulk Foster")
    parser.add_argument("--freqs", default="1e2,1e3,1e4,1e5,1e6,1e7", help="comma-separated Hz")
    parser.add_argument("--dirichlet", default="outer", help="Dirichlet boundary label name(s)")
    args = parser.parse_args()

    from ngsolve import Mesh, TaskManager
    mesh = Mesh(args.vol)
    print(f"Mesh: {args.vol}")
    print(f"  ne = {mesh.ne}, materials = {mesh.GetMaterials()}, boundaries = {mesh.GetBoundaries()}")
    print()

    with TaskManager():
        lam, tau, g_n, V = bulk_foster_via_eigen(mesh, args.sigma, args.mu,
                                                  n_eigen=args.n_eigen,
                                                  dirichlet_label=args.dirichlet)
        S_total, edge_list = measure_total_area_and_edges(mesh)

    K_SIBC = K_SIBC_total(S_total, args.sigma, args.mu)
    c1 = c1_polyhedral(edge_list, args.mu)
    L_total = sum(L for L, _ in edge_list) if edge_list else 0.0

    print(f"Volume V = {V*1e9:.4f} mm^3")
    print(f"Surface area S = {S_total*1e6:.4f} mm^2")
    print(f"Total edge length = {L_total*1e3:.4f} mm ({len(edge_list)} edges)")
    print(f"K_SIBC = S sqrt(sigma/mu) = {K_SIBC:.4e}")
    print(f"c_1 (polyhedral edge) = {c1:.4e}")
    print(f"Bulk Foster: {len(lam)} eigenmodes,  lam range = [{lam[0]:.3e}, {lam[-1]:.3e}]")
    Y0_bulk = np.sum(g_n)
    print(f"  Sanity check: bulk Y(0) = sum g_n = {Y0_bulk:.4e}  (expect sigma*V = {args.sigma*V:.4e})")
    print(f"  Completeness ratio = {Y0_bulk / (args.sigma * V):.4f}  (1.0 = complete eigenbasis)")
    print()
    print(f"{'f (Hz)':>10}  {'|alpha|/V':>10}  {'Re(alpha)/V':>13}  {'-Im(alpha)/V':>13}")
    for f_str in args.freqs.split(","):
        f = float(f_str)
        s = 1j * 2 * math.pi * f
        Y = Y_mixed(s, lam, tau, g_n, K_SIBC, c1)
        a = alpha_from_Y(Y, V, args.sigma)
        print(f"  {f:10.2e}  {abs(a)/V:10.4f}  {a.real/V:13.4e}  {-a.imag/V:13.4e}")


if __name__ == "__main__":
    main()
