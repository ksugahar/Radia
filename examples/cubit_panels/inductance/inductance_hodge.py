"""
BEM inductance extraction via Hodge decomposition on genus-1 surface.

Solves the EFIE (not energy method) by decomposing HDivSurface into
solenoidal + harmonic subspaces. The harmonic mode on a genus-1 surface
(torus) represents the circulating current, and its Rayleigh quotient
with the LaplaceSL matrix gives the inductance.

Method (from Lucy Weggler's ngsbem framework):
  1. Build divergence matrix D: HDivSurface -> SurfaceL2
  2. Build vertex-edge incidence C (graph curl)
  3. Hodge decomposition: harmonic = ker([D; C^T M_J])
  4. LaplaceSL on HDivSurface -> project onto harmonic subspace
  5. Generalized eigenvalue: SL_harm * c = lambda * M_harm * c
  6. L = mu_0 * lambda * (2*pi*R) / (2*pi*a)   [3D torus scaling]

CRITICAL: Use Glue(torus.faces) to create a surface-only OCC geometry.
Without Glue, Netgen creates a volume mesh and mesh.edges includes
interior edges, breaking the Euler characteristic (V-E+F != 0).

Usage:
    python inductance_hodge.py

Requires: NGSolve, scipy
"""

import math
import numpy as np
from scipy.linalg import null_space, eigh
from ngsolve import *
from ngsolve.bem import LaplaceSL
from netgen.occ import WorkPlane, Axes, Axis, Pnt, Dir, OCCGeometry, Glue
from netgen.meshing import MeshingParameters

MU_0 = 4e-7 * np.pi


def extract_rect(mat, nrow, ncol):
    M = np.zeros((nrow, ncol))
    ei = mat.CreateRowVector()
    col = mat.CreateColVector()
    for j in range(ncol):
        ei[:] = 0; ei[j] = 1.0
        mat.Mult(ei, col)
        for i in range(nrow):
            M[i, j] = col[i]
    return M


def extract_dense(mat, n):
    M = np.zeros((n, n))
    ei = mat.CreateColVector()
    col = mat.CreateColVector()
    for j in range(n):
        ei[:] = 0; ei[j] = 1.0
        mat.Mult(ei, col)
        for i in range(n):
            M[i, j] = col[i]
    return M


def make_torus_surface_mesh(R, a, curvaturesafety=1.0):
    """Create surface-only OCC torus mesh (genus-1).

    MUST use Glue(torus.faces) to get a surface geometry.
    Without Glue, OCCGeometry generates a volume mesh (dim=3)
    where mesh.edges includes interior edges.
    """
    wp = WorkPlane(Axes(p=Pnt(R, 0, 0), n=Dir(0, 1, 0), h=Dir(0, 0, 1)))
    circle = wp.Circle(a).Face()
    torus = circle.Revolve(Axis(p=Pnt(0, 0, 0), d=Dir(0, 0, 1)), 360)
    torus_surf = Glue(torus.faces)
    geo = OCCGeometry(torus_surf)
    ngmesh = geo.GenerateMesh(
        mp=MeshingParameters(maxh=1.0, curvaturesafety=curvaturesafety,
                             segmentsperedge=2))
    return Mesh(ngmesh)


def compute_inductance_hodge(mesh, R, a):
    """Compute self-inductance via Hodge decomposition + LaplaceSL.

    Args:
        mesh: NGSolve surface mesh (genus-1, Euler=0)
        R: major radius [m]
        a: minor radius [m]

    Returns:
        dict with L_toroidal, L_poloidal, harmonic modes, diagnostics
    """
    nse = mesh.GetNE(BND)
    nv = mesh.nv
    ne = sum(1 for e in mesh.edges)
    euler = nv - ne + nse

    fes_J = HDivSurface(mesh, order=0)
    fes_L2 = SurfaceL2(mesh, order=0)
    n_J, n_f = fes_J.ndof, fes_L2.ndof

    # D: divergence matrix (HDivSurface -> SurfaceL2)
    u_J, q = fes_J.TrialFunction(), fes_L2.TestFunction()
    bf_D = BilinearForm(trialspace=fes_J, testspace=fes_L2)
    bf_D += div(u_J.Trace()) * q * ds
    bf_D.Assemble()
    D = extract_rect(bf_D.mat, n_f, n_J)

    # C: vertex-edge incidence (graph curl)
    C = np.zeros((n_J, nv))
    for i, e in enumerate(mesh.edges):
        vv = list(e.vertices)
        if i < n_J:
            C[i, vv[0].nr] = -1
            C[i, vv[1].nr] = +1

    # M_J: edge mass matrix
    u, v = fes_J.TnT()
    bf_M = BilinearForm(fes_J)
    bf_M += InnerProduct(u.Trace(), v.Trace()) * ds
    bf_M.Assemble()
    M_J = extract_dense(bf_M.mat, n_J)

    # Harmonic subspace: ker([D; C^T M_J])
    A_hodge = np.vstack([D, C.T @ M_J])
    c_h_all = null_space(A_hodge, rcond=1e-10)
    n_harm = c_h_all.shape[1]

    # LaplaceSL (use_fmm=False for reproducibility + dense extraction)
    jt, jv = fes_J.TnT()
    with TaskManager():
        V_op = LaplaceSL(jt.Trace() * ds, use_fmm=False) * jv.Trace() * ds
        SL = V_op.mat.ToDense().NumPy()

    # Project onto harmonic subspace
    SL_harm = c_h_all.T @ SL @ c_h_all
    M_harm = c_h_all.T @ M_J @ c_h_all

    # Generalized eigenvalue: SL c = lambda M c
    eigvals_h, eigvecs_h = eigh(SL_harm, M_harm)

    # 3D torus scaling: L = mu_0 * lambda * R/a
    L_modes = [MU_0 * ev * R / a for ev in eigvals_h]

    return {
        'nse': nse, 'nv': nv, 'ne': ne, 'euler': euler,
        'n_J': n_J, 'n_harm': n_harm,
        'eigvals': eigvals_h,
        'L_modes': L_modes,
    }


if __name__ == "__main__":
    R, a = 0.05, 0.005
    L_neumann = MU_0 * R * (math.log(8 * R / a) - 2)

    print(f"Torus: R = {R*1e3:.0f} mm, a = {a*1e3:.0f} mm, R/a = {R/a:.0f}")
    print(f"Neumann formula: L = {L_neumann*1e9:.2f} nH")
    print()

    header = f"{'cs':>6s} {'nse':>6s} {'n_J':>6s} {'harm':>5s} {'Euler':>6s} {'L_tor(nH)':>12s} {'err':>8s} {'L_pol(nH)':>12s}"
    print(header)
    print("-" * 70)

    for cs in [0.5, 1.0, 1.5, 2.0]:
        mesh = make_torus_surface_mesh(R, a, cs)
        r = compute_inductance_hodge(mesh, R, a)

        if r['n_harm'] >= 2:
            L_pol = r['L_modes'][0]
            L_tor = r['L_modes'][1]
            err = (L_tor - L_neumann) / L_neumann * 100
            print(f"{cs:6.1f} {r['nse']:6d} {r['n_J']:6d} {r['n_harm']:5d} {r['euler']:6d}"
                  f"  {L_tor*1e9:10.4f}  {err:+6.2f}%  {L_pol*1e9:10.4f}")
        elif r['n_harm'] == 1:
            L_0 = r['L_modes'][0]
            print(f"{cs:6.1f} {r['nse']:6d} {r['n_J']:6d} {r['n_harm']:5d} {r['euler']:6d}"
                  f"  {L_0*1e9:10.4f}  (only 1 mode, Euler={r['euler']})")
        else:
            print(f"{cs:6.1f} {r['nse']:6d} {r['n_J']:6d} {r['n_harm']:5d} {r['euler']:6d}  FAILED")

    print(f"\nNeumann reference: {L_neumann*1e9:.2f} nH")
