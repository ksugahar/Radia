"""
BEM inductance via source/sink saddle point EFIE.

Solves the constrained EFIE on the conductor surface:
  [SL  D^T] [J] = [0]
  [D   0  ] [p] = [g]

where:
  SL = LaplaceSL (single layer BEM operator on surface currents)
  D  = divergence matrix (HDivSurface -> SurfaceL2)
  g  = source/sink current injection (+1 at source, -1 at sink)
  J  = surface current (unknowns)
  p  = Lagrange multiplier enforcing current conservation

Then inductance: L = mu_0 * J^T @ SL @ J

This is the standard EFIE approach for partial inductance extraction.
Works for any conductor topology with source/sink ports.

Usage (OCC gapped torus test):
    python inductance_source_sink.py

Usage (Cubit model):
    from inductance_source_sink import compute_inductance_source_sink
    L = compute_inductance_source_sink(mesh)

Requires: NGSolve with ngsolve.bem, scipy
"""

import math
import time
import numpy as np
from scipy.linalg import solve as scipy_solve
from ngsolve import *
from ngsolve.bem import LaplaceSL

MU_0 = 4e-7 * np.pi


def _to_dense(mat):
    """Extract dense NumPy array from NGSolve BaseMatrix via COO.

    NGSolve BEM operators (LaplaceSL etc.) store the full dense matrix
    internally as a SparseMatrix with 100% fill.  The built-in ToDense()
    is O(N) slower than necessary (~144s vs 0.06s at N=5085) because it
    performs N column-by-column MatVecs instead of a direct memory copy.

    This function extracts the COO triplets and converts via scipy,
    which is ~2500x faster.
    """
    from scipy.sparse import coo_matrix
    rows, cols, vals = mat.COO()
    return coo_matrix((vals, (rows, cols)),
                      shape=(mat.height, mat.width)).toarray()


def compute_inductance_source_sink(mesh, source_label="source", sink_label="sink",
                                    fes_order=0):
    """Compute self-inductance via saddle point EFIE with source/sink ports.

    Args:
        mesh: NGSolve Mesh (volume mesh with boundary, or surface-only mesh)
        source_label: Boundary label for current injection face
        sink_label: Boundary label for current extraction face
        fes_order: HDivSurface polynomial order (0=RWG, 1=higher-order)

    Returns:
        dict with L [H], J (current coefficients), diagnostics
    """
    nse = mesh.GetNE(BND)
    nv = mesh.nv

    fes_J = HDivSurface(mesh, order=fes_order)
    fes_L2 = SurfaceL2(mesh, order=0)  # constraint always order=0 (element-wise)
    n_J = fes_J.ndof
    n_f = fes_L2.ndof

    print(f"  DOFs: n_J={n_J} (edges), n_f={n_f} (faces)")
    print(f"  Boundary elements: {nse}, vertices: {nv}")

    # --- Divergence matrix D: n_f x n_J ---
    t0 = time.perf_counter()
    u_J = fes_J.TrialFunction()
    q = fes_L2.TestFunction()
    bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
    bf_D += div(u_J.Trace()) * q * ds
    bf_D.Assemble()
    D = _to_dense(bf_D.mat)
    t_div = time.perf_counter() - t0
    print(f"  Divergence matrix: {t_div:.2f}s")

    # --- LaplaceSL matrix: n_J x n_J ---
    t0 = time.perf_counter()
    jt, jv = fes_J.TnT()
    with TaskManager():
        V_op = LaplaceSL(jt.Trace() * ds, use_fmm=False) * jv.Trace() * ds
    t_asm = time.perf_counter() - t0
    t0 = time.perf_counter()
    SL = _to_dense(V_op.mat)
    t_ext = time.perf_counter() - t0
    print(f"  LaplaceSL assembly: {t_asm:.2f}s, COO extract: {t_ext:.2f}s")

    # --- Source/sink RHS ---
    # g[i] = +1/A_source * A_elem_i  for source elements
    #       = -1/A_sink  * A_elem_i  for sink elements
    #       = 0                       otherwise
    f_src = LinearForm(fes_L2)
    f_src += q * ds(source_label)
    f_src.Assemble()
    g_src = f_src.vec.FV().NumPy().copy()
    A_src = np.sum(g_src)

    f_snk = LinearForm(fes_L2)
    f_snk += q * ds(sink_label)
    f_snk.Assemble()
    g_snk = f_snk.vec.FV().NumPy().copy()
    A_snk = np.sum(g_snk)

    if A_src < 1e-30 or A_snk < 1e-30:
        return {"error": f"Source/sink faces not found (A_src={A_src}, A_snk={A_snk}). "
                f"Check boundary labels: '{source_label}', '{sink_label}'. "
                f"Available: {list(set(mesh.GetBoundaries()))}"}

    # Unit current: total injection = +1, total extraction = -1
    g = g_src / A_src - g_snk / A_snk
    print(f"  Source area: {A_src:.6e} m^2, Sink area: {A_snk:.6e} m^2")
    print(f"  Current balance: sum(g) = {np.sum(g):.2e} (should be ~0)")

    # --- Saddle point system ---
    # D has one redundant row (sum of all rows = 0, charge conservation).
    # Remove last row to make system non-singular.
    D_red = D[:-1, :]
    g_red = g[:-1]
    n_constraint = n_f - 1

    t0 = time.perf_counter()
    K = np.block([
        [SL,              D_red.T],
        [D_red, np.zeros((n_constraint, n_constraint))]
    ])
    rhs = np.zeros(n_J + n_constraint)
    rhs[n_J:] = g_red

    x = scipy_solve(K, rhs)
    J = x[:n_J]
    p = x[n_J:]
    t_solve = time.perf_counter() - t0
    print(f"  Saddle point solve: {t_solve:.2f}s (size {n_J + n_constraint})")

    # --- Inductance: L = mu_0 * J^T @ SL @ J ---
    L = MU_0 * J @ SL @ J

    # Verify: D*J should equal g
    residual = np.max(np.abs(D @ J - g))
    print(f"  Constraint residual: max|D*J - g| = {residual:.2e}")

    return {
        'L': float(L),
        'n_J': n_J,
        'n_f': n_f,
        'A_source': float(A_src),
        'A_sink': float(A_snk),
        't_SL': t_asm + t_ext,
        't_solve': t_solve,
        'residual': float(residual),
        'J': J,
    }


def make_gapped_torus_mesh(R, a, gap_deg=5, maxh=None):
    """Create OCC gapped torus surface mesh with source/sink labels.

    Uses Glue(faces) to create a surface-only mesh (no volume elements),
    which gives correct HDivSurface DOF count for BEM.

    Args:
        R: Major radius [m]
        a: Minor radius [m]
        gap_deg: Gap angle [degrees]
        maxh: Max mesh size (default: a/2)

    Returns:
        NGSolve Mesh with boundary labels "conductor", "source", "sink"
    """
    from netgen.occ import (WorkPlane, Axes, Axis, Pnt, Dir,
                             OCCGeometry, Glue)
    from netgen.meshing import MeshingParameters

    if maxh is None:
        maxh = a / 2

    # Cross-section circle at (R, 0, 0)
    wp = WorkPlane(Axes(p=Pnt(R, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp.Circle(a).Face()

    # Revolve to create gapped torus solid
    torus = circle.Revolve(Axis(Pnt(0, 0, 0), Dir(0, 0, 1)), 360 - gap_deg)

    # Name faces before Glue: identify gap faces by area (smallest = pi*a^2)
    for f in torus.faces:
        f.name = "conductor"
    faces_sorted = sorted(torus.faces, key=lambda f: f.mass)
    faces_sorted[0].name = "source"
    faces_sorted[1].name = "sink"

    # Glue -> surface-only geometry (like inductance_hodge.py)
    torus_surf = Glue(torus.faces)
    geo = OCCGeometry(torus_surf)
    ngmesh = geo.GenerateMesh(
        mp=MeshingParameters(maxh=maxh, curvaturesafety=2.0,
                             segmentsperedge=2))
    mesh = Mesh(ngmesh)
    mesh.Curve(2)

    return mesh


if __name__ == "__main__":
    R, a = 0.05, 0.005
    L_neumann = MU_0 * R * (math.log(8 * R / a) - 2)

    print(f"Torus: R = {R*1e3:.0f} mm, a = {a*1e3:.0f} mm")
    print(f"Neumann formula: L = {L_neumann*1e9:.2f} nH")
    print()

    # Single coarse mesh test first, then convergence
    mesh = make_gapped_torus_mesh(R, a, gap_deg=5, maxh=a)
    nse = mesh.GetNE(BND)
    nv = mesh.nv
    ne_edges = sum(1 for e in mesh.edges)
    euler = nv - ne_edges + nse

    bnd_names = list(set(mesh.GetBoundaries()))
    print(f"Boundary labels: {bnd_names}")
    print(f"Mesh: {nse} faces, {nv} vertices, {ne_edges} edges, Euler={euler}")

    result = compute_inductance_source_sink(mesh)

    if 'error' in result:
        print(f"ERROR: {result['error']}")
    else:
        L = result['L']
        err = (L - L_neumann) / L_neumann * 100
        print(f"\nResult: L = {L*1e9:.4f} nH  (err = {err:+.2f}%)")
        print(f"Neumann reference: {L_neumann*1e9:.2f} nH")
