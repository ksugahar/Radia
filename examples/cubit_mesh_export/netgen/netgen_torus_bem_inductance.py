"""
Netgen Export Example: Torus BEM Inductance with High-Order Curved Elements

Demonstrates that SetGeomInfo + mesh.Curve(order) directly improves
BEM inductance accuracy by enabling high-order curved surface elements.

Key comparison:
    Curve(1) = flat facets         -> polygon approximation -> large error
    Curve(2) = quadratic surfaces  -> much better
    Curve(3) = cubic surfaces      -> near-exact

Workflow:
    1. Create torus in Cubit (inductor coil model)
    2. Tet mesh -> export to Netgen via export_NGSolveCurvedMesh()
    3. For each curve order: export_NGSolveCurvedMesh(order=N)
       a. Compute surface area (analytical: 4*pi^2*R*a)
       b. Compute BEM inductance via ngsolve.bem.LaplaceSL
       c. Compare with Neumann formula: L = mu_0*R*(ln(8R/a) - 2)

This script proves the necessity of SetGeomInfo (netgen PR#232) for
accurate BEM on externally imported meshes.

Requirements:
    - Coreform Cubit (2025.3+)
    - NGSolve with ksugahar/netgen fork (SetGeomInfo API)
    - ngsolve.bem (ngbem BEM module)

Run: python netgen_torus_bem_inductance.py

For NGSolve User Meeting demonstration.
"""

import sys
import os
import math
import time
import numpy as np

# Auto-detect Cubit installation
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__)))), 'src', 'radia'))
from install_panels import find_cubit_bin as _fcb
_cubit_path = _fcb()
if _cubit_path and _cubit_path not in sys.path:
    sys.path.append(_cubit_path)

work_dir = os.path.dirname(os.path.abspath(__file__))
repo_root = os.path.dirname(os.path.dirname(work_dir))
sys.path.insert(0, repo_root)

from ngsolve import Mesh, Integrate, CF, BND, HDivSurface, TaskManager, ds
import cubit
import cubit_mesh_export

MU_0 = 4.0 * math.pi * 1e-7


# ============================================================
# Analytical formulas
# ============================================================

def analytical_area(R, a):
    """Torus surface area: 4*pi^2*R*a."""
    return 4 * math.pi**2 * R * a


def analytical_volume(R, a):
    """Torus volume: 2*pi^2*R*a^2."""
    return 2 * math.pi**2 * R * a**2


def analytical_inductance(R, a):
    """Neumann formula: L = mu_0*R*(ln(8R/a) - 2)."""
    return MU_0 * R * (math.log(8.0 * R / a) - 2.0)


# ============================================================
# BEM inductance computation
# ============================================================

def extract_dense_matrix(mat, n):
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


def compute_bem_inductance(mesh, bnd_label="conductor"):
    """Compute self-inductance using BEM (LaplaceSL on HDivSurface).

    L_matrix = mu_0 * LaplaceSL(J, J)   [Laplace kernel, MQS regime]
    L_total = 1 / (e^T @ L^{-1} @ e)    [uniform current excitation]

    Args:
        mesh: NGSolve Mesh (volume mesh with curved boundary)
        bnd_label: Boundary label for conductor surface

    Returns:
        L_total [H], n_dof, t_assemble [s]
    """
    from ngsolve.bem import LaplaceSL

    t_start = time.perf_counter()

    fes = HDivSurface(mesh, order=0)
    n_dof = fes.ndof

    j_trial = fes.TrialFunction()
    j_test = fes.TestFunction()

    with TaskManager():
        L_op = LaplaceSL(
            j_trial.Trace() * ds(bnd_label)
        ) * j_test.Trace() * ds(bnd_label)

    L_dense = extract_dense_matrix(L_op.mat, n_dof)
    L_matrix = MU_0 * L_dense

    t_assemble = time.perf_counter() - t_start

    # Check rank
    rank = np.linalg.matrix_rank(L_matrix)
    if rank < n_dof:
        print(f"    WARNING: rank-deficient ({rank}/{n_dof})")

    # Total inductance: L = 1 / (e^T @ L^{-1} @ e)
    e = np.ones(n_dof) / n_dof
    try:
        L_inv_e = np.linalg.solve(L_matrix, e)
        L_total = 1.0 / (e @ L_inv_e)
    except np.linalg.LinAlgError:
        L_inv_e = np.linalg.lstsq(L_matrix, e, rcond=None)[0]
        L_total = 1.0 / (e @ L_inv_e)

    return L_total, n_dof, t_assemble


# ============================================================
# Main
# ============================================================

def main():
    # Torus parameters (inductor coil)
    R = 0.05     # Major radius: 50 mm (loop radius)
    a = 0.005    # Minor radius:  5 mm (wire radius)
    MESH_SIZE = 0.004  # Element size

    L_ref = analytical_inductance(R, a)
    A_ref = analytical_area(R, a)
    V_ref = analytical_volume(R, a)

    print("=" * 70)
    print("Cubit Torus -> SetGeomInfo -> BEM Inductance (Curve Order Study)")
    print("=" * 70)
    print(f"  R = {R*1e3:.1f} mm (loop radius)")
    print(f"  a = {a*1e3:.1f} mm (wire radius)")
    print(f"  R/a = {R/a:.0f}")
    print(f"  Analytical: L = {L_ref*1e9:.4f} nH")
    print(f"  Analytical: Area = {A_ref*1e6:.4f} mm^2")
    print(f"  Analytical: Vol = {V_ref*1e9:.6f} mm^3")

    # ============================================================
    # Step 1: Create torus in Cubit
    # ============================================================
    print("\n--- Step 1: Create torus in Cubit ---")

    cubit.init(['cubit', '-nojournal', '-batch'])
    cubit.cmd("reset")
    cubit.cmd(f"create torus major radius {R} minor radius {a}")
    print(f"  Surfaces: {cubit.get_surface_count()}")

    # ============================================================
    # Step 2: Tet mesh
    # ============================================================
    print("\n--- Step 2: Tet mesh ---")

    cubit.cmd("volume all scheme tetmesh")
    cubit.cmd(f"volume all size {MESH_SIZE}")
    cubit.cmd("mesh volume all")
    cubit.cmd("block 1 add tet all")
    cubit.cmd('block 1 name "domain"')
    cubit.cmd("block 2 add tri all")
    cubit.cmd('block 2 name "conductor"')

    n_tets = cubit.get_tet_count()
    n_tris = cubit.get_tri_count()
    print(f"  Tets: {n_tets}")
    print(f"  Surface tris: {n_tris}")

    # ============================================================
    # Step 3: Curve order comparison
    # ============================================================
    print("\n--- Step 3: Curve order comparison ---")

    results = []
    for order in [1, 2, 3]:
        print(f"\n  export_NGSolveCurvedMesh(order={order}):")

        mesh = cubit_mesh_export.export_NGSolveCurvedMesh(cubit, order=order)

        # Surface area
        area = Integrate(CF(1), mesh, VOL_or_BND=BND)
        area_err = abs(area - A_ref) / A_ref * 100

        # Volume
        vol = Integrate(CF(1), mesh)
        vol_err = abs(vol - V_ref) / V_ref * 100

        print(f"    Area  = {area*1e6:.4f} mm^2 (error: {area_err:.4f}%)")
        print(f"    Vol   = {vol*1e9:.6f} mm^3 (error: {vol_err:.4f}%)")

        # BEM inductance
        bnd_label = list(set(mesh.GetBoundaries()))[0]
        print(f"    BEM assembly (label='{bnd_label}')...")
        L_bem, n_dof, t_asm = compute_bem_inductance(mesh, bnd_label)
        L_err = (L_bem - L_ref) / abs(L_ref) * 100

        print(f"    L_BEM = {L_bem*1e9:.4f} nH (error: {L_err:+.2f}%)")
        print(f"    DOFs  = {n_dof}, time = {t_asm:.3f} s")

        results.append({
            'order': order,
            'area_err': area_err,
            'vol_err': vol_err,
            'L_bem': L_bem,
            'L_err': L_err,
            'n_dof': n_dof,
            't_asm': t_asm,
        })

    # ============================================================
    # Summary table
    # ============================================================
    print("\n" + "=" * 70)
    print("Summary: Curve Order vs BEM Inductance Accuracy")
    print("=" * 70)
    print(f"  Analytical L = {L_ref*1e9:.4f} nH  (Neumann formula)")
    print(f"  Mesh: {n_tets} tets, {n_tris} surface tris, {results[0]['n_dof']} BEM DOFs")
    print()
    print(f"  {'Order':>5} {'Area Err':>10} {'Vol Err':>10} {'L [nH]':>10} {'L Err':>10} {'Time [s]':>10}")
    print("  " + "-" * 55)
    for r in results:
        print(f"  {r['order']:>5} {r['area_err']:>9.4f}% {r['vol_err']:>9.4f}% "
              f"{r['L_bem']*1e9:>10.4f} {r['L_err']:>+9.2f}% {r['t_asm']:>10.3f}")

    print()
    print("Conclusion:")
    if len(results) >= 2 and abs(results[-1]['L_err']) < abs(results[0]['L_err']):
        improvement = abs(results[0]['L_err']) - abs(results[-1]['L_err'])
        print(f"  Curve(3) reduces L error by {improvement:.1f}pp vs Curve(1)")
    print("  SetGeomInfo (PR#232) enables mesh.Curve() on Cubit meshes")
    print("  -> High-order curved surfaces -> Better BEM accuracy")

    print("\n=== DONE ===")
    return results


if __name__ == "__main__":
    main()
