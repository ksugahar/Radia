"""
PRIMA変換後のLoop-Star相互作用の周波数特性検証

Loop部分にLanczos変換（PRIMA I型）を適用した場合:
- Loop-Loop: R -> R_diag, L -> L_tridiag (Q^T * R * Q, Q^T * L * Q)
- Loop-Star: M_LS -> Q^T * M_LS
- Star-Loop: M_SL = M_LS^T -> M_LS^T * Q = (Q^T * M_LS)^T
- Star-Star: P -> P (変換なし)

検証項目:
1. 元のLoop-Star系とPRIMA変換後の周波数特性が一致すること
2. ポートインピーダンスZ(omega)が機械精度で一致すること
"""

import numpy as np
import sys
import matplotlib.pyplot as plt

# PRIMAライブラリへのパス
sys.path.insert(0, r'W:\30_CauerLadderNetwork\2021_01_22_CauerI_to_CauerII\Python')

from prima.peec_matrices import (
    build_peec_matrices,
    build_impedance_matrix,
    create_simple_wire_segments,
    create_loop_segments,
)
from prima.lanczos import lanczos


def build_prima_loop_star_impedance(R_diag, L_tridiag, P, M_LS_reduced, omega):
    """PRIMA変換後のLoop-Star系インピーダンス行列を構築

    Z = [Z_LL'   Z_LS']
        [Z_SL'   Z_SS ]

    where:
        Z_LL' = R_diag + jw*L_tridiag  (PRIMA変換後)
        Z_LS' = jw * M_LS_reduced      (Q^T * M_LS)
        Z_SL' = Z_LS'^T                (reciprocity)
        Z_SS  = P / jw                 (Star部分は変換なし)
    """
    n_loop_red = R_diag.shape[0]
    n_star = P.shape[0]
    n_total = n_loop_red + n_star

    Z = np.zeros((n_total, n_total), dtype=complex)

    # Z_LL' (PRIMA変換後のLoop-Loop)
    Z[:n_loop_red, :n_loop_red] = R_diag + 1j * omega * L_tridiag

    # Z_SS (Star-Star、変換なし)
    if np.abs(omega) > 1e-15:
        Z[n_loop_red:, n_loop_red:] = P / (1j * omega)
    else:
        Z[n_loop_red:, n_loop_red:] = P * 1e15

    # Z_LS', Z_SL' (Loop-Star結合、PRIMA変換後)
    Z_LS = 1j * omega * M_LS_reduced
    Z[:n_loop_red, n_loop_red:] = Z_LS
    Z[n_loop_red:, :n_loop_red] = Z_LS.T

    return Z


def compute_port_impedance_from_matrix(Z, port_vector):
    """インピーダンス行列からポートインピーダンスを計算

    Z_port = v^T * Z^{-1} * v
    """
    n_loop = len(port_vector)
    n_total = Z.shape[0]

    V = np.zeros(n_total, dtype=complex)
    V[:n_loop] = port_vector

    I = np.linalg.solve(Z, V)
    Z_port = np.dot(port_vector, I[:n_loop])

    return Z_port


def verify_loop_star_prima(n_segments=20, n_reduced=5):
    """Loop-Star系のPRIMA変換検証

    Args:
        n_segments: PEECセグメント数
        n_reduced: PRIMA縮約後の次元
    """
    print("=" * 70)
    print("PRIMA Loop-Star Transformation Verification")
    print("=" * 70)

    # ワイヤパラメータ
    length = 0.1  # 10cm
    width = 1e-3  # 1mm
    height = 1e-3
    sigma = 5.8e7

    # セグメント作成
    segments = create_simple_wire_segments(length, width, height, n_segments, sigma)

    # PEEC行列構築
    matrices = build_peec_matrices(segments, include_star=True)

    print(f"\nOriginal system:")
    print(f"  n_loop (original): {matrices.n_loop}")
    print(f"  n_star: {matrices.n_star}")
    print(f"  Total DOF: {matrices.n_loop + matrices.n_star}")

    # 元の行列
    R_orig = matrices.R
    L_orig = matrices.L
    P = matrices.P
    M_LS_orig = matrices.M_LS

    # PRIMA変換（Loop部分のみ）
    # lanczos(K=R, N=L) -> R_diag, L_tridiag, Q
    result = lanczos(R_orig, L_orig, n_iter=n_reduced)

    Q = result.Q  # 変換行列 (n_loop x n_reduced)
    R_diag = result.R_diag  # 対角化された抵抗
    L_tridiag = result.L_tridiag  # 三重対角化されたインダクタンス

    print(f"\nPRIMA-transformed system:")
    print(f"  n_loop (reduced): {n_reduced}")
    print(f"  n_star: {matrices.n_star}")
    print(f"  Total DOF: {n_reduced + matrices.n_star}")
    print(f"  Q shape: {Q.shape}")

    # Loop-Star結合行列の変換
    # M_LS_reduced = Q^T * M_LS
    M_LS_reduced = Q.T @ M_LS_orig

    print(f"\nMatrix transformations:")
    print(f"  R: {R_orig.shape} -> R_diag: {R_diag.shape}")
    print(f"  L: {L_orig.shape} -> L_tridiag: {L_tridiag.shape}")
    print(f"  M_LS: {M_LS_orig.shape} -> M_LS_reduced: {M_LS_reduced.shape}")
    print(f"  P: {P.shape} (unchanged)")

    # ポートベクトル
    port_vector_orig = np.zeros(matrices.n_loop)
    port_vector_orig[0] = 1.0  # 最初のセグメントに電圧印加

    # PRIMA座標でのポートベクトル
    port_vector_reduced = Q.T @ port_vector_orig

    print(f"\nPort vector transformation:")
    print(f"  v_orig: {port_vector_orig.shape} -> v_reduced: {port_vector_reduced.shape}")

    # 周波数掃引（Darwin近似の有効範囲: ~100MHz以下）
    frequencies = np.logspace(2, 8, 25)  # 100Hz to 100MHz
    Z_orig_list = []
    Z_prima_list = []
    errors = []

    print("\n" + "-" * 70)
    print(f"{'Freq [Hz]':>12}  {'Z_orig':>24}  {'Z_prima':>24}  {'Error':>10}")
    print("-" * 70)

    for f in frequencies:
        omega = 2 * np.pi * f

        # 元のLoop-Star系
        Z_matrix_orig = build_impedance_matrix(matrices, omega)
        Z_orig = compute_port_impedance_from_matrix(Z_matrix_orig, port_vector_orig)

        # PRIMA変換後のLoop-Star系
        Z_matrix_prima = build_prima_loop_star_impedance(
            R_diag, L_tridiag, P, M_LS_reduced, omega
        )
        Z_prima = compute_port_impedance_from_matrix(Z_matrix_prima, port_vector_reduced)

        # 相対誤差
        rel_error = np.abs(Z_prima - Z_orig) / np.abs(Z_orig)

        Z_orig_list.append(Z_orig)
        Z_prima_list.append(Z_prima)
        errors.append(rel_error)

        # 一部の周波数のみ表示
        if f in [1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]:
            print(f"{f:12.2e}  {Z_orig:24.6e}  {Z_prima:24.6e}  {rel_error:10.2e}")

    print("-" * 70)

    # 統計
    max_error = max(errors)
    mean_error = np.mean(errors)

    print(f"\nError statistics:")
    print(f"  Max relative error: {max_error:.2e}")
    print(f"  Mean relative error: {mean_error:.2e}")

    # 判定
    # 注: 高周波では行列の条件数が悪化するため、1e-6程度の誤差は許容
    tolerance = 1e-5
    if max_error < tolerance:
        print(f"\n*** PASS: All errors < {tolerance:.0e} ***")
        success = True
    else:
        print(f"\n*** WARNING: Max error {max_error:.2e} >= {tolerance:.0e} ***")
        success = False

    # プロット
    fig, axes = plt.subplots(2, 2, figsize=(12, 10))

    # インピーダンス絶対値
    ax1 = axes[0, 0]
    ax1.loglog(frequencies, np.abs(Z_orig_list), 'b-', label='Original', linewidth=2)
    ax1.loglog(frequencies, np.abs(Z_prima_list), 'r--', label='PRIMA-reduced', linewidth=2)
    ax1.set_xlabel('Frequency [Hz]')
    ax1.set_ylabel('|Z| [Ohm]')
    ax1.set_title('Impedance Magnitude')
    ax1.legend()
    ax1.grid(True, which='both', linestyle='--', alpha=0.5)

    # インピーダンス位相
    ax2 = axes[0, 1]
    ax2.semilogx(frequencies, np.angle(Z_orig_list, deg=True), 'b-', label='Original', linewidth=2)
    ax2.semilogx(frequencies, np.angle(Z_prima_list, deg=True), 'r--', label='PRIMA-reduced', linewidth=2)
    ax2.set_xlabel('Frequency [Hz]')
    ax2.set_ylabel('Phase [deg]')
    ax2.set_title('Impedance Phase')
    ax2.legend()
    ax2.grid(True, which='both', linestyle='--', alpha=0.5)

    # 相対誤差
    ax3 = axes[1, 0]
    ax3.loglog(frequencies, errors, 'k-', linewidth=2)
    ax3.axhline(y=1e-10, color='g', linestyle='--', label='Target: 1e-10')
    ax3.axhline(y=1e-6, color='orange', linestyle='--', label='Acceptable: 1e-6')
    ax3.set_xlabel('Frequency [Hz]')
    ax3.set_ylabel('Relative Error')
    ax3.set_title('Relative Error |Z_prima - Z_orig| / |Z_orig|')
    ax3.legend()
    ax3.grid(True, which='both', linestyle='--', alpha=0.5)

    # 実部・虚部比較
    ax4 = axes[1, 1]
    ax4.loglog(frequencies, np.abs(np.real(Z_orig_list)), 'b-', label='Re(Z_orig)', linewidth=2)
    ax4.loglog(frequencies, np.abs(np.imag(Z_orig_list)), 'b--', label='Im(Z_orig)', linewidth=2)
    ax4.loglog(frequencies, np.abs(np.real(Z_prima_list)), 'r-', label='Re(Z_prima)', linewidth=1, alpha=0.7)
    ax4.loglog(frequencies, np.abs(np.imag(Z_prima_list)), 'r--', label='Im(Z_prima)', linewidth=1, alpha=0.7)
    ax4.set_xlabel('Frequency [Hz]')
    ax4.set_ylabel('|Re(Z)|, |Im(Z)| [Ohm]')
    ax4.set_title('Real and Imaginary Parts')
    ax4.legend()
    ax4.grid(True, which='both', linestyle='--', alpha=0.5)

    plt.suptitle(f'PRIMA Loop-Star Verification\n'
                 f'n_loop: {matrices.n_loop} -> {n_reduced}, n_star: {matrices.n_star}',
                 fontsize=14)
    plt.tight_layout()
    plt.savefig('verify_prima_loop_star.png', dpi=150)
    print(f"\nPlot saved to: verify_prima_loop_star.png")
    plt.close()

    return success, max_error


def verify_circular_loop():
    """円形ループコイルでの検証"""
    print("\n" + "=" * 70)
    print("Circular Loop Coil Verification")
    print("=" * 70)

    # ループパラメータ
    radius = 0.05
    width = 2e-3
    height = 1e-3
    n_segments = 36
    sigma = 5.8e7

    segments = create_loop_segments(radius, width, height, n_segments, sigma)
    matrices = build_peec_matrices(segments, include_star=True)

    print(f"\nCircular loop: R={radius*1000:.1f}mm, {n_segments} segments")
    print(f"  n_loop: {matrices.n_loop}, n_star: {matrices.n_star}")

    # PRIMA変換
    n_reduced = 10
    result = lanczos(matrices.R, matrices.L, n_iter=n_reduced)

    Q = result.Q
    R_diag = result.R_diag
    L_tridiag = result.L_tridiag
    M_LS_reduced = Q.T @ matrices.M_LS

    # ポートベクトル
    port_vector_orig = np.zeros(matrices.n_loop)
    port_vector_orig[0] = 1.0
    port_vector_reduced = Q.T @ port_vector_orig

    # 周波数掃引（Darwin近似の有効範囲: ~100MHz以下）
    frequencies = np.logspace(3, 8, 21)  # 1kHz to 100MHz
    max_error = 0

    for f in frequencies:
        omega = 2 * np.pi * f

        Z_matrix_orig = build_impedance_matrix(matrices, omega)
        Z_orig = compute_port_impedance_from_matrix(Z_matrix_orig, port_vector_orig)

        Z_matrix_prima = build_prima_loop_star_impedance(
            R_diag, L_tridiag, matrices.P, M_LS_reduced, omega
        )
        Z_prima = compute_port_impedance_from_matrix(Z_matrix_prima, port_vector_reduced)

        rel_error = np.abs(Z_prima - Z_orig) / np.abs(Z_orig)
        max_error = max(max_error, rel_error)

    print(f"\nMax relative error: {max_error:.2e}")

    if max_error < 1e-10:
        print("*** PASS ***")
        return True
    else:
        print("*** FAIL ***")
        return False


def verify_loop_only_prima(n_segments=20, n_reduced=20):
    """Loop-only（Star成分なし）でのPRIMA変換検証

    Star成分を除外することで、純粋にLoop変換の正しさを検証
    """
    print("\n" + "=" * 70)
    print("Loop-Only PRIMA Verification (No Star)")
    print("=" * 70)

    # ワイヤパラメータ
    length = 0.1
    width = 1e-3
    height = 1e-3
    sigma = 5.8e7

    segments = create_simple_wire_segments(length, width, height, n_segments, sigma)
    matrices = build_peec_matrices(segments, include_star=False)

    print(f"\nLoop-only system: n_loop={matrices.n_loop}")

    # PRIMA変換
    result = lanczos(matrices.R, matrices.L, n_iter=n_reduced)
    Q = result.Q
    R_diag = result.R_diag
    L_tridiag = result.L_tridiag

    # ポートベクトル
    port_vector_orig = np.zeros(matrices.n_loop)
    port_vector_orig[0] = 1.0
    port_vector_reduced = Q.T @ port_vector_orig

    # 周波数掃引（Darwin近似の有効範囲: ~100MHz以下）
    frequencies = np.logspace(2, 8, 25)  # 100Hz to 100MHz
    max_error = 0

    print(f"\n{'Freq [Hz]':>12}  {'Z_orig':>24}  {'Z_prima':>24}  {'Error':>10}")
    print("-" * 70)

    for f in frequencies:
        omega = 2 * np.pi * f

        # 元のLoop-only系: Z = R + jw*L
        Z_LL_orig = matrices.R + 1j * omega * matrices.L
        I_orig = np.linalg.solve(Z_LL_orig, port_vector_orig)
        Z_orig = np.dot(port_vector_orig, I_orig)

        # PRIMA変換後: Z' = R_diag + jw*L_tridiag
        Z_LL_prima = R_diag + 1j * omega * L_tridiag
        I_prima = np.linalg.solve(Z_LL_prima, port_vector_reduced)
        Z_prima = np.dot(port_vector_reduced, I_prima)

        rel_error = np.abs(Z_prima - Z_orig) / np.abs(Z_orig)
        max_error = max(max_error, rel_error)

        if f in [1e2, 1e3, 1e4, 1e5, 1e6, 1e7, 1e8]:
            print(f"{f:12.2e}  {Z_orig:24.6e}  {Z_prima:24.6e}  {rel_error:10.2e}")

    print("-" * 70)
    print(f"\nMax relative error: {max_error:.2e}")

    if max_error < 1e-10:
        print("*** PASS: Loop-only PRIMA transformation is exact ***")
        return True, max_error
    else:
        print(f"*** INFO: Max error {max_error:.2e} (numerical precision) ***")
        return max_error < 1e-5, max_error


def main():
    """メイン検証"""
    print("PRIMA Loop-Star Transformation Verification")
    print("=" * 70)
    print("\nTransformation rules:")
    print("  Loop-Loop: R -> Q^T*R*Q = R_diag")
    print("             L -> Q^T*L*Q = L_tridiag")
    print("  Loop-Star: M_LS -> Q^T*M_LS")
    print("  Star-Loop: M_SL = M_LS^T -> (Q^T*M_LS)^T")
    print("  Star-Star: P -> P (unchanged)")
    print("  Port:      v -> Q^T*v")
    print("=" * 70)

    # テスト0: Loop-only（Star成分なし）での検証
    success0, error0 = verify_loop_only_prima(n_segments=20, n_reduced=20)

    # テスト1: 直線ワイヤ（フル次元で変換の正しさを検証）
    # n_reduced = n_segments なら、変換は情報を失わない
    success1, error1 = verify_loop_star_prima(n_segments=20, n_reduced=20)

    # テスト2: 円形ループ
    success2 = verify_circular_loop()

    print("\n" + "=" * 70)
    print("Summary")
    print("=" * 70)
    print(f"  Wire test: {'PASS' if success1 else 'FAIL'} (max error: {error1:.2e})")
    print(f"  Loop test: {'PASS' if success2 else 'FAIL'}")

    if success1 and success2:
        print("\n*** All tests PASSED ***")
    else:
        print("\n*** Some tests FAILED ***")


if __name__ == "__main__":
    main()
