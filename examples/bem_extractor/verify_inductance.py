"""
verify_inductance.py

Verify BEM inductance calculation against analytical formula
for a circular loop (torus).

Analytical self-inductance (Neumann formula):
    L = mu_0 * R * (ln(8R/a) - 2)

where R = loop radius, a = wire radius.

This script:
  1. Creates a torus surface mesh using Netgen OCC
  2. Assembles the BEM inductance matrix L via ngsolve.bem.LaplaceSL
  3. Computes total inductance via uniform excitation: L = 1/(e^T L^{-1} e)
  4. Compares with analytical formula

Known limitation:
  - Dense surface meshes (curvaturesafety >= 1.0) may have rank-deficient
    BEM matrices due to ngbem surface integration issues on closed surfaces.
  - Use curvaturesafety=0.5 (coarse mesh) for reliable results.
  - For higher accuracy, use p-refinement (order > 0) or quad elements
    from Cubit, NOT finer triangulation.

Usage:
    python verify_inductance.py

Part of Radia project
"""

import sys
import os
import numpy as np
import time

# Add Radia source to path
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../src/radia'))

MU_0 = 4.0 * np.pi * 1e-7


def create_torus_mesh(R, a, curvaturesafety=0.5):
    """Create surface mesh of a torus (circular loop).

    Args:
        R: Major radius (loop radius) [m]
        a: Minor radius (wire radius) [m]
        curvaturesafety: Netgen mesh density control (0.5 = coarse, 1.0 = medium)

    Returns:
        NGSolve Mesh (surface mesh)
    """
    from netgen.occ import WorkPlane, Axes, Axis, Pnt, Dir, OCCGeometry
    from netgen.meshing import MeshingParameters
    from ngsolve import Mesh

    # Cross-section circle at (R, 0, 0), in Y-Z plane
    wp = WorkPlane(Axes(p=Pnt(R, 0, 0),
                        n=Dir(0, 1, 0),
                        h=Dir(0, 0, 1)))
    circle = wp.Circle(a).Face()
    circle.name = "conductor"

    # Revolve around Z axis -> torus
    torus = circle.Revolve(Axis(p=Pnt(0, 0, 0), d=Dir(0, 0, 1)), 360)

    geo = OCCGeometry(torus)
    mp = MeshingParameters(maxh=1.0, curvaturesafety=curvaturesafety,
                           segmentsperedge=2)
    ngmesh = geo.GenerateMesh(mp=mp)

    mesh = Mesh(ngmesh)
    return mesh


def analytical_L(R, a):
    """Neumann formula: L = mu_0 * R * (ln(8R/a) - 2)."""
    return MU_0 * R * (np.log(8.0 * R / a) - 2.0)


def extract_dense(mat, n):
    """Extract dense matrix from NGSolve BaseMatrix."""
    M = np.zeros((n, n))
    ei = mat.CreateColVector()
    col = mat.CreateColVector()
    for i in range(n):
        ei[:] = 0
        ei[i] = 1.0
        mat.Mult(ei, col)
        for j in range(n):
            M[j, i] = col[j]
    return M


def compute_bem_inductance(mesh):
    """Compute self-inductance using BEM (LaplaceSL on HDivSurface).

    L_matrix = mu_0 * LaplaceSL(J, J)

    Total self-inductance = 1 / (e^T @ L^{-1} @ e)
    with uniform current excitation e = ones/n (valid for order=0 RT0 only).

    Args:
        mesh: NGSolve surface mesh

    Returns:
        L_total: Total self-inductance [H]
        n_dof: Number of edge DOFs
        t_assemble: Assembly time [s]
    """
    from ngsolve import HDivSurface, TaskManager, ds
    from ngsolve.bem import LaplaceSL

    t_start = time.perf_counter()

    # HDivSurface function space (RWG/RT0 edge elements on surface)
    fes = HDivSurface(mesh, order=0)
    n_dof = fes.ndof

    # Get actual boundary label from mesh
    bnd_label = list(set(mesh.GetBoundaries()))[0]

    # Assemble LaplaceSL: G(r) = 1/(4*pi*r), MQS kernel
    j_trial = fes.TrialFunction()
    j_test = fes.TestFunction()

    with TaskManager():
        L_op = LaplaceSL(
            j_trial.Trace() * ds(bnd_label)
        ) * j_test.Trace() * ds(bnd_label)

    # Extract dense matrix
    L_dense = extract_dense(L_op.mat, n_dof)
    L_matrix = MU_0 * L_dense

    t_assemble = time.perf_counter() - t_start

    # Check rank
    rank = np.linalg.matrix_rank(L_matrix)
    if rank < n_dof:
        print(f"  WARNING: Rank-deficient matrix ({rank}/{n_dof})")
        print(f"  This occurs on dense surface meshes. Use coarser mesh.")

    # Total inductance with uniform excitation
    e = np.ones(n_dof) / n_dof
    try:
        L_inv_e = np.linalg.solve(L_matrix, e)
        L_total = 1.0 / (e @ L_inv_e)
    except np.linalg.LinAlgError:
        # Fallback: pseudoinverse
        L_inv_e = np.linalg.lstsq(L_matrix, e, rcond=None)[0]
        L_total = 1.0 / (e @ L_inv_e)
        print("  WARNING: L matrix singular, using lstsq")

    return L_total, n_dof, t_assemble


def main():
    print("=" * 60)
    print("BEM Inductance Verification (ngsolve.bem)")
    print("  Circular Loop: L = mu_0*R*(ln(8R/a) - 2)")
    print("=" * 60)

    # Test cases: (R, a, description)
    # Use curvaturesafety=0.5 for reliable BEM assembly
    cases = [
        (0.05, 0.005, "R=50mm, a=5mm (R/a=10)"),
        (0.05, 0.002, "R=50mm, a=2mm (R/a=25)"),
        (0.10, 0.005, "R=100mm, a=5mm (R/a=20)"),
    ]

    results = []
    for R, a, desc in cases:
        print(f"\n--- {desc} ---")
        L_ref = analytical_L(R, a)
        print(f"  Analytical L = {L_ref*1e9:.4f} nH")

        print(f"  Creating mesh (curvaturesafety=0.5)...")
        mesh = create_torus_mesh(R, a, curvaturesafety=0.5)

        n_surf = sum(1 for _ in mesh.Elements())
        print(f"  Surface elements: {n_surf}")

        print(f"  Assembling BEM...")
        L_bem, n_dof, t_asm = compute_bem_inductance(mesh)

        error = (L_bem - L_ref) / abs(L_ref) * 100

        print(f"  BEM L        = {L_bem*1e9:.4f} nH")
        print(f"  Analytical L = {L_ref*1e9:.4f} nH")
        print(f"  Error        = {error:+.2f}%")
        print(f"  DOFs         = {n_dof}")
        print(f"  Time         = {t_asm:.3f} s")

        results.append({
            'R': R, 'a': a, 'desc': desc,
            'L_ref': L_ref, 'L_bem': L_bem,
            'error_pct': error, 'n_dof': n_dof,
            't_asm': t_asm
        })

    # Summary
    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    print(f"{'Case':<30} {'L_ref[nH]':>10} {'L_bem[nH]':>10} {'Error':>8} {'DOFs':>6}")
    print("-" * 70)
    all_pass = True
    for r in results:
        status = "PASS" if abs(r['error_pct']) < 30 else "FAIL"
        if status == "FAIL":
            all_pass = False
        print(f"{r['desc']:<30} {r['L_ref']*1e9:10.4f} {r['L_bem']*1e9:10.4f} "
              f"{r['error_pct']:+7.2f}% {r['n_dof']:6d}")

    print()
    print("NOTE: BEM with uniform excitation (e=ones/n) and coarse mesh")
    print("gives ~15-25% error. For higher accuracy:")
    print("  - Use higher-order elements (order > 0 with proper excitation)")
    print("  - Use quad elements from Cubit (better geometry on curves)")
    print("  - Do NOT simply refine mesh (rank-deficiency issue on closed surfaces)")

    if all_pass:
        print("\nAll cases PASS (error < 30%)")
    else:
        print("\nSome cases FAIL")

    return results


if __name__ == "__main__":
    main()
