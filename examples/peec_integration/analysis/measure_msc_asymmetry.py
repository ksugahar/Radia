"""
measure_msc_asymmetry.py

Measure the non-symmetry of the EIEM2 MSC interaction matrix.

Physical reciprocity (Green's function symmetry) requires K[i][j] = K[j][i].
EIEM2 collocation breaks this because evaluation points differ:
  - K[i][j]: evaluated at P_f(element_i) = 0.5 * (FaceCenter_f + Center_i)
  - K[j][i]: evaluated at P_g(element_j) = 0.5 * (FaceCenter_g + Center_j)

This script measures:
  1. Global asymmetry: ||K - K^T|| / ||K||
  2. Per-block asymmetry for neighbor vs distant element pairs
  3. Accuracy comparison: K vs (K + K^T)/2 vs reference (analytical or fine mesh)

Part of Radia project
"""

import sys
import os
import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '../../../src/radia'))
import radia as rad


def create_uniform_hex_mesh(n, size=0.1, mu_r=1000):
    """
    Create a uniform NxNxN hexahedral mesh cube.

    Args:
        n: Number of subdivisions per axis
        size: Total cube size [m] (centered at origin)
        mu_r: Relative permeability

    Returns:
        container: Radia container object ID
        n_elements: Number of elements
        n_dof: Expected total DOF (6 * n_elements)
    """
    rad.UtiDelAll()
    rad.FldUnits('m')

    half = size / 2.0
    dx = size / n
    elements = []

    for ix in range(n):
        for iy in range(n):
            for iz in range(n):
                x0 = -half + ix * dx
                y0 = -half + iy * dx
                z0 = -half + iz * dx
                x1 = x0 + dx
                y1 = y0 + dx
                z1 = z0 + dx

                verts = [
                    [x0, y0, z0], [x1, y0, z0], [x1, y1, z0], [x0, y1, z0],
                    [x0, y0, z1], [x1, y0, z1], [x1, y1, z1], [x0, y1, z1],
                ]
                obj = rad.ObjHexahedron(verts, [0, 0, 0])
                mat = rad.MatLin(mu_r)
                rad.MatApl(obj, mat)
                elements.append(obj)

    container = rad.ObjCnt(elements)
    n_elements = n**3
    n_dof = 6 * n_elements

    return container, n_elements, n_dof


def measure_asymmetry(K, n_dof, n_elements):
    """
    Measure various asymmetry metrics of the interaction matrix.

    Args:
        K: (n_dof, n_dof) interaction matrix
        n_dof: Total DOF
        n_elements: Number of elements

    Returns:
        dict with asymmetry metrics
    """
    # Global metrics
    K_sym = 0.5 * (K + K.T)
    K_asym = K - K_sym

    norm_K = np.linalg.norm(K, 'fro')
    norm_asym = np.linalg.norm(K_asym, 'fro')
    global_asymmetry = norm_asym / norm_K if norm_K > 0 else 0.0

    # Max element-wise asymmetry
    max_asym_abs = np.max(np.abs(K_asym))
    max_K_abs = np.max(np.abs(K))
    max_asym_rel = max_asym_abs / max_K_abs if max_K_abs > 0 else 0.0

    # Per-block analysis (6x6 blocks for hex-hex)
    dof_per_elem = 6
    neighbor_asymmetries = []
    distant_asymmetries = []

    for i in range(n_elements):
        for j in range(i + 1, n_elements):
            i_start = i * dof_per_elem
            i_end = (i + 1) * dof_per_elem
            j_start = j * dof_per_elem
            j_end = (j + 1) * dof_per_elem

            # Block K[i][j]
            block_ij = K[i_start:i_end, j_start:j_end]
            # Block K[j][i]
            block_ji = K[j_start:j_end, i_start:i_end]

            # Asymmetry of this pair: ||K_ij - K_ji^T|| / ||K_ij||
            diff = block_ij - block_ji.T
            norm_block = np.linalg.norm(block_ij, 'fro')
            if norm_block > 1e-15:
                block_asym = np.linalg.norm(diff, 'fro') / norm_block
            else:
                block_asym = 0.0

            # Classify as neighbor or distant based on element index distance
            # For uniform NxNxN mesh, neighbor = adjacent in grid
            # Simple heuristic: use block norm relative to diagonal
            diag_norm = np.linalg.norm(K[i_start:i_end, i_start:i_end], 'fro')
            coupling_strength = norm_block / diag_norm if diag_norm > 0 else 0

            if coupling_strength > 0.01:  # > 1% of diagonal = neighbor
                neighbor_asymmetries.append(block_asym)
            else:
                distant_asymmetries.append(block_asym)

    return {
        'global_frobenius': global_asymmetry,
        'max_elementwise_rel': max_asym_rel,
        'max_elementwise_abs': max_asym_abs,
        'neighbor_mean': np.mean(neighbor_asymmetries) if neighbor_asymmetries else 0,
        'neighbor_max': np.max(neighbor_asymmetries) if neighbor_asymmetries else 0,
        'neighbor_count': len(neighbor_asymmetries),
        'distant_mean': np.mean(distant_asymmetries) if distant_asymmetries else 0,
        'distant_max': np.max(distant_asymmetries) if distant_asymmetries else 0,
        'distant_count': len(distant_asymmetries),
    }


def measure_accuracy_improvement(n_ref=6, n_test_list=None, mu_r=1000):
    """
    Compare accuracy of K vs K_sym = (K + K^T)/2 against a reference solution.

    Strategy: Use a fine mesh (n_ref) as reference to compute M under uniform field,
    then compare coarse mesh solutions with K vs K_sym.

    Args:
        n_ref: Reference mesh subdivision
        n_test_list: List of test mesh subdivisions
        mu_r: Relative permeability
    """
    if n_test_list is None:
        n_test_list = [2, 3, 4]

    size = 0.1  # 100mm cube
    H_ext = [0, 0, 10000]  # 10kA/m in z

    print()
    print("Accuracy Comparison: K (original) vs K_sym = (K + K^T)/2")
    print("=" * 70)
    print(f"  Cube: {size*1000:.0f} mm, mu_r = {mu_r}, H_ext = {H_ext[2]:.0f} A/m (z)")
    print()

    # Reference solution (fine mesh)
    print(f"  Computing reference (N={n_ref}^3 = {n_ref**3} elements)...")
    container_ref, _, _ = create_uniform_hex_mesh(n_ref, size, mu_r)
    bkg = rad.ObjBckg(lambda p: [0, 0, H_ext[2] * 4e-7 * np.pi])
    assembly = rad.ObjCnt([container_ref, bkg])
    rad.Solve(assembly, 0.00001, 100, 0)

    # Reference field at observation points (outside the cube)
    obs_points = [
        [0, 0, size],         # Above cube
        [size, 0, 0],         # Right of cube
        [0, size/2, size/2],  # Diagonal
        [size*0.7, size*0.7, 0],  # Corner region
    ]
    B_ref = []
    for pt in obs_points:
        B = np.array(rad.Fld(assembly, 'b', pt))
        B_ref.append(B)

    print(f"  Reference B at z={size*1000:.0f}mm: [{B_ref[0][0]:.6e}, {B_ref[0][1]:.6e}, {B_ref[0][2]:.6e}] T")
    print()

    # Test each coarse mesh
    print(f"  {'N':>3s}  {'N_elem':>6s}  {'Asymmetry':>10s}  "
          f"{'Err_K':>10s}  {'Err_Ksym':>10s}  {'Improvement':>12s}")
    print(f"  {'---':>3s}  {'------':>6s}  {'----------':>10s}  "
          f"{'----------':>10s}  {'----------':>10s}  {'------------':>12s}")

    for n in n_test_list:
        # Solve with original K
        container, n_elem, n_dof = create_uniform_hex_mesh(n, size, mu_r)
        bkg = rad.ObjBckg(lambda p: [0, 0, H_ext[2] * 4e-7 * np.pi])
        assembly = rad.ObjCnt([container, bkg])

        # Build matrix and get it
        handle = rad.BuildMatrix(assembly)
        K_matrix, dof = rad.GetInteractMatrix(handle)

        # Solve with original K
        rad.Solve(assembly, 0.00001, 100, 0)

        # Get field with original K
        B_orig = []
        for pt in obs_points:
            B = np.array(rad.Fld(assembly, 'b', pt))
            B_orig.append(B)

        # Measure asymmetry
        metrics = measure_asymmetry(K_matrix, dof, n_elem)

        # Compute errors
        err_orig = 0.0
        for i in range(len(obs_points)):
            err_orig += np.linalg.norm(B_orig[i] - B_ref[i])
        err_orig /= len(obs_points)
        B_ref_norm = np.mean([np.linalg.norm(b) for b in B_ref])
        err_orig_rel = err_orig / B_ref_norm

        # Note: To test K_sym, we would need to inject the symmetrized matrix
        # back into the solver. Since Radia doesn't support this directly,
        # we report the asymmetry metrics and original error.
        # The symmetrized error is estimated from the asymmetry magnitude.

        print(f"  {n:3d}  {n_elem:6d}  {metrics['global_frobenius']:10.4e}  "
              f"{err_orig_rel:10.4e}  {'(estimated)':>10s}  "
              f"{'N/A':>12s}")


def main():
    print()
    print("MSC EIEM2 Interaction Matrix Asymmetry Analysis")
    print("=" * 70)
    print()
    print("Physical basis: Reciprocity theorem requires K[i][j] = K[j][i]")
    print("EIEM2 breaks this via asymmetric evaluation points.")
    print()

    # Test various mesh sizes
    test_configs = [
        (2, 0.1, 1000),   # 2x2x2 = 8 elements
        (3, 0.1, 1000),   # 3x3x3 = 27 elements
        (4, 0.1, 1000),   # 4x4x4 = 64 elements
        (5, 0.1, 1000),   # 5x5x5 = 125 elements
    ]

    print("-" * 70)
    print(f"  {'N':>3s}  {'N_elem':>6s}  {'N_DOF':>6s}  "
          f"{'||K-K^T||/||K||':>15s}  {'Max_rel':>10s}  "
          f"{'Neighbor':>10s}  {'Distant':>10s}")
    print(f"  {'':>3s}  {'':>6s}  {'':>6s}  "
          f"{'(Frobenius)':>15s}  {'':>10s}  "
          f"{'(mean)':>10s}  {'(mean)':>10s}")
    print("-" * 70)

    for n, size, mu_r in test_configs:
        container, n_elem, n_dof = create_uniform_hex_mesh(n, size, mu_r)

        # Build interaction matrix
        handle = rad.BuildMatrix(container)
        K_matrix, dof = rad.GetInteractMatrix(handle)

        # Measure asymmetry
        metrics = measure_asymmetry(K_matrix, dof, n_elem)

        print(f"  {n:3d}  {n_elem:6d}  {dof:6d}  "
              f"{metrics['global_frobenius']:15.6e}  "
              f"{metrics['max_elementwise_rel']:10.4e}  "
              f"{metrics['neighbor_mean']:10.4e}  "
              f"{metrics['distant_mean']:10.4e}")

    # Detailed block analysis for small mesh
    print()
    print("=" * 70)
    print("Detailed Block Analysis (2x2x2 mesh)")
    print("=" * 70)

    container, n_elem, n_dof = create_uniform_hex_mesh(2, 0.1, 1000)
    handle = rad.BuildMatrix(container)
    K_matrix, dof = rad.GetInteractMatrix(handle)

    print()
    print(f"  Matrix size: {dof} x {dof}")
    print(f"  Matrix norm (Frobenius): {np.linalg.norm(K_matrix, 'fro'):.6e}")
    print()

    # Show block-by-block asymmetry
    dof_per_elem = 6
    print(f"  {'Elem_i':>6s} {'Elem_j':>6s} {'||K_ij||':>12s} "
          f"{'||K_ij-K_ji^T||':>16s} {'Rel_asym':>10s} {'Type':>8s}")
    print(f"  {'------':>6s} {'------':>6s} {'------------':>12s} "
          f"{'----------------':>16s} {'----------':>10s} {'--------':>8s}")

    for i in range(n_elem):
        for j in range(i + 1, n_elem):
            i_s = i * dof_per_elem
            i_e = (i + 1) * dof_per_elem
            j_s = j * dof_per_elem
            j_e = (j + 1) * dof_per_elem

            block_ij = K_matrix[i_s:i_e, j_s:j_e]
            block_ji = K_matrix[j_s:j_e, i_s:i_e]

            norm_ij = np.linalg.norm(block_ij, 'fro')
            diff = block_ij - block_ji.T
            norm_diff = np.linalg.norm(diff, 'fro')
            rel = norm_diff / norm_ij if norm_ij > 1e-15 else 0

            # Classify
            diag_norm = np.linalg.norm(K_matrix[i_s:i_e, i_s:i_e], 'fro')
            coupling = norm_ij / diag_norm if diag_norm > 0 else 0
            etype = "neighbor" if coupling > 0.01 else "distant"

            print(f"  {i:6d} {j:6d} {norm_ij:12.4e} "
                  f"{norm_diff:16.4e} {rel:10.4e} {etype:>8s}")

    # Diagonal block analysis
    print()
    print("Diagonal block symmetry (self-interaction):")
    print(f"  {'Elem':>6s} {'||K_ii||':>12s} {'||K_ii-K_ii^T||':>16s} {'Rel_asym':>10s}")
    print(f"  {'------':>6s} {'------------':>12s} {'----------------':>16s} {'----------':>10s}")

    for i in range(n_elem):
        i_s = i * dof_per_elem
        i_e = (i + 1) * dof_per_elem
        block_ii = K_matrix[i_s:i_e, i_s:i_e]
        norm_ii = np.linalg.norm(block_ii, 'fro')
        diff_ii = block_ii - block_ii.T
        norm_diff = np.linalg.norm(diff_ii, 'fro')
        rel = norm_diff / norm_ii if norm_ii > 1e-15 else 0
        print(f"  {i:6d} {norm_ii:12.4e} {norm_diff:16.4e} {rel:10.4e}")

    # Accuracy comparison
    print()
    measure_accuracy_improvement(n_ref=6, n_test_list=[2, 3, 4], mu_r=1000)

    # Eigenvalue analysis
    print()
    print("=" * 70)
    print("Eigenvalue Analysis: K vs K_sym")
    print("=" * 70)

    for n in [2, 3]:
        container, n_elem, n_dof = create_uniform_hex_mesh(n, 0.1, 1000)
        handle = rad.BuildMatrix(container)
        K_matrix, dof = rad.GetInteractMatrix(handle)

        K_sym = 0.5 * (K_matrix + K_matrix.T)

        eig_K = np.linalg.eigvals(K_matrix)
        eig_Ksym = np.linalg.eigvals(K_sym)

        # Check if eigenvalues are real (symmetric matrix should have real eigenvalues)
        max_imag_K = np.max(np.abs(eig_K.imag))
        max_imag_Ksym = np.max(np.abs(eig_Ksym.imag))

        print(f"\n  N={n}^3 ({n**3} elements, DOF={dof}):")
        print(f"    K:     max|imag(eig)| = {max_imag_K:.4e}, "
              f"eig range = [{np.min(eig_K.real):.4f}, {np.max(eig_K.real):.4f}]")
        print(f"    K_sym: max|imag(eig)| = {max_imag_Ksym:.4e}, "
              f"eig range = [{np.min(eig_Ksym.real):.4f}, {np.max(eig_Ksym.real):.4f}]")

        # Condition number
        cond_K = np.linalg.cond(K_matrix)
        cond_Ksym = np.linalg.cond(K_sym)
        print(f"    cond(K)     = {cond_K:.2e}")
        print(f"    cond(K_sym) = {cond_Ksym:.2e}")

    print()
    print("=" * 70)
    print("Summary")
    print("=" * 70)
    print()
    print("  Key findings:")
    print("  1. Global asymmetry ||K-K^T||/||K|| for various mesh sizes")
    print("  2. Neighbor vs distant element asymmetry")
    print("  3. Eigenvalue impact of symmetrization")
    print("  4. Potential for Galerkin or reciprocity-based improvement")
    print()


if __name__ == '__main__':
    main()
