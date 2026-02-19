"""
Verify U = V for symmetric Lanczos transformation

This script verifies that for symmetric matrices K and N,
the Lanczos transformation produces U = V (or U proportional to V).

Mathematical background:
- For symmetric A: A = Q * T * Q^T where T is symmetric tridiagonal
- In the generalized case with K, N both symmetric:
  K' = U^H * K * V (tridiagonal)
  N' = U^H * N * V (diagonal)

The question is: Does U = V for symmetric K, N?

Answer: In general, U != V for the generalized Lanczos with two matrices.
However, for the special case where N = I (identity), we get U = V = Q.

For PRIMA (Cauer Ladder Network) transformation where:
  K = L (inductance, tridiagonal)
  N = R (resistance, diagonal)

The algorithm uses a different recurrence, and U, V are constructed differently.
Let's verify numerically what the relationship is.
"""

import sys
import numpy as np

# Add PRIMA library path
sys.path.insert(0, r'W:\30_CauerLadderNetwork\2021_01_22_CauerI_to_CauerII\Python')

try:
    from prima import lanczos
except ImportError:
    print("ERROR: 'prima' module not found.")
    print("This example requires the PRIMA library from:")
    print("  W:\\30_CauerLadderNetwork\\2021_01_22_CauerI_to_CauerII\\Python")
    print("Skipping this example.")
    sys.exit(0)


def test_uv_equality():
    """Test if U = V for symmetric Lanczos."""
    print("=" * 60)
    print("Testing U vs V in Lanczos Transformation")
    print("=" * 60)

    # Test 1: Simple symmetric case - L tridiagonal, R diagonal
    print("\n--- Test 1: PRIMA I form (L tridiagonal, R diagonal) ---")

    n = 5
    L_diag = np.array([1.0, 0.5, 0.3, 0.2, 0.1])  # Inductance
    R_diag = np.array([0.1, 0.1, 0.1, 0.1, 0.1])  # Resistance

    # Build tridiagonal L matrix
    L = np.diag(L_diag)
    for i in range(n - 1):
        L[i, i+1] = -L_diag[i+1]
        L[i+1, i] = -L_diag[i+1]

    R = np.diag(R_diag)

    print(f"L matrix (tridiagonal):\n{L}")
    print(f"R matrix (diagonal):\n{R}")

    result = lanczos(L, R, n_iter=n)

    print(f"\nU:\n{result.U}")
    print(f"\nV:\n{result.V}")

    # Check if U = V
    diff_uv = np.linalg.norm(result.U - result.V, 'fro')
    norm_u = np.linalg.norm(result.U, 'fro')
    rel_diff = diff_uv / norm_u if norm_u > 1e-15 else diff_uv

    print(f"\n||U - V||_F = {diff_uv:.6e}")
    print(f"||U||_F = {norm_u:.6e}")
    print(f"Relative diff = {rel_diff:.6e}")

    if rel_diff < 1e-10:
        print("=> U = V (identical within numerical precision)")
    else:
        print("=> U != V (different matrices)")

        # Check if U and V span the same subspace
        # by checking if U^T V is close to identity (up to sign/permutation)
        print("\n  Checking if U and V span the same subspace...")
        UtV = result.U.T @ result.V
        print(f"  U^T V =\n{UtV}")

        # Check if columns are proportional
        print("\n  Checking column-by-column proportionality...")
        for j in range(min(result.U.shape[1], 3)):
            u_col = result.U[:, j]
            v_col = result.V[:, j]
            ratio = u_col / (v_col + 1e-30)  # Avoid division by zero
            print(f"  Column {j}: U/V ratios = {ratio[:3]}...")

    # Test 2: Dense symmetric L - PRIMA I型
    # lanczos(K, N): K->LL(三重対角), N->RR(対角)
    # PRIMA I型 (s=0展開): R対角, L三重対角
    # よって: lanczos(K=L, N=R) -> LL=L'(三重対角), RR=R'(対角)
    print("\n\n--- Test 2: PRIMA I型 (K=L dense, N=R diagonal) ---")
    print("lanczos(K, N): K->LL(tridiag), N->RR(diag)")
    print("For PRIMA I: lanczos(K=L, N=R) -> LL=L'(tridiag), RR=R'(diag)")

    # Dense L with mutual inductance
    L_dense = np.array([
        [1.0, 0.3, 0.1, 0.05, 0.02],
        [0.3, 1.0, 0.3, 0.1, 0.05],
        [0.1, 0.3, 1.0, 0.3, 0.1],
        [0.05, 0.1, 0.3, 1.0, 0.3],
        [0.02, 0.05, 0.1, 0.3, 1.0]
    ])

    print(f"L_dense (symmetric):\n{L_dense}")
    print(f"R matrix (diagonal):\n{R}")

    # PRIMA I型: lanczos(K=L, N=R) -> LL=L'(三重対角), RR=R'(対角)
    result2 = lanczos(L_dense, R, n_iter=n)

    diff_uv2 = np.linalg.norm(result2.U - result2.V, 'fro')
    norm_u2 = np.linalg.norm(result2.U, 'fro')
    rel_diff2 = diff_uv2 / norm_u2 if norm_u2 > 1e-15 else diff_uv2

    print(f"\n||U - V||_F = {diff_uv2:.6e}")
    print(f"Relative diff = {rel_diff2:.6e}")

    if rel_diff2 < 1e-10:
        print("=> U = V (identical)")
    else:
        print("=> U != V (different)")

    # Test 3: Standard symmetric eigenvalue (K symmetric, N = I)
    print("\n\n--- Test 3: Standard symmetric case (N = I) ---")

    I = np.eye(n)
    result3 = lanczos(L_dense, I, n_iter=n)

    diff_uv3 = np.linalg.norm(result3.U - result3.V, 'fro')
    norm_u3 = np.linalg.norm(result3.U, 'fro')
    rel_diff3 = diff_uv3 / norm_u3 if norm_u3 > 1e-15 else diff_uv3

    print(f"||U - V||_F = {diff_uv3:.6e}")
    print(f"Relative diff = {rel_diff3:.6e}")

    if rel_diff3 < 1e-10:
        print("=> U = V (identical)")
    else:
        print("=> U != V (different)")

    # Verify the transformation preserves structure
    print("\n\n--- Verification of PRIMA I Transformation ---")
    print("lanczos(K=R, N=L) function:")
    print("  R -> R_diag (diagonalized)")
    print("  L -> L_tridiag (tridiagonalized)")

    # PRIMA I型: lanczos(K=R, N=L) -> R_diag (対角), L_tridiag (三重対角)
    print("\n--- Test 4: PRIMA I型 (K=R, N=L_dense) ---")
    result4 = lanczos(R, L_dense, n_iter=n)

    print(f"R_diag (should be DIAGONAL):\n{result4.R_diag}")
    print(f"L_tridiag (should be TRIDIAGONAL):\n{result4.L_tridiag}")

    # Check DIAGONAL structure for R_diag
    is_diag_R = np.allclose(result4.R_diag, np.diag(np.diag(result4.R_diag)), atol=1e-10)
    print(f"\nR_diag is diagonal: {is_diag_R}")

    # Check TRIDIAGONAL structure for L_tridiag
    is_tridiag_L = True
    for i in range(result4.L_tridiag.shape[0]):
        for j in range(result4.L_tridiag.shape[1]):
            if abs(i - j) > 1 and abs(result4.L_tridiag[i, j]) > 1e-10:
                is_tridiag_L = False
                break

    print(f"L_tridiag is tridiagonal: {is_tridiag_L}")

    if is_diag_R and is_tridiag_L:
        print("\n=> PRIMA I型 structure confirmed!")

    print("\n" + "=" * 60)
    print("Conclusion:")
    print("=" * 60)
    print("""
PRIMA I型モデル縮約:
   lanczos(K=R, N=L) -> R_diag (対角), L_tridiag (三重対角)

PEEC model order reduction:
   lanczos(K=R_dc, N=L_ext) -> R_diag, L_tridiag (PRIMA I型)
""")


if __name__ == '__main__':
    test_uv_equality()
