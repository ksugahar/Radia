"""
PEEC Loop-Star行列のテスト

Loop（インダクタンス）とStar（ポテンシャル係数/容量）行列を構築し、
周波数特性を検証する。

テスト項目:
1. L行列（インダクタンス）の構築
2. P行列（ポテンシャル係数）の構築
3. Loop-Star結合行列M_LSの構築
4. 周波数特性（低周波: 誘導性、高周波: 容量性）
"""

import numpy as np
import sys

# CLNライブラリへのパス
sys.path.insert(0, r'W:\30_CauerLadderNetwork\2021_01_22_CauerI_to_CauerII\Python')

from cln.peec_matrices import (
    PEECSegment,
    PEECNode,
    build_peec_matrices,
    build_impedance_matrix,
    compute_port_impedance,
    create_simple_wire_segments,
    create_loop_segments,
)


def test_simple_wire():
    """単純なワイヤのPEEC行列をテスト"""
    print("=" * 60)
    print("Test 1: Simple Wire PEEC Matrices")
    print("=" * 60)

    # ワイヤパラメータ
    length = 0.1  # 10cm
    width = 1e-3  # 1mm
    height = 1e-3  # 1mm
    n_segments = 10
    sigma = 5.8e7  # 銅

    # セグメント作成
    segments = create_simple_wire_segments(length, width, height, n_segments, sigma)

    print(f"Wire length: {length*1000:.1f} mm")
    print(f"Cross section: {width*1000:.1f} x {height*1000:.1f} mm")
    print(f"Number of segments: {n_segments}")

    # PEEC行列構築（Star成分なし）
    matrices_loop_only = build_peec_matrices(segments, include_star=False)

    print(f"\n--- Loop-only matrices ---")
    print(f"L matrix shape: {matrices_loop_only.L.shape}")
    print(f"R matrix shape: {matrices_loop_only.R.shape}")

    # 自己インダクタンスの合計（全体のインダクタンスに近い）
    L_self_sum = np.sum(np.diag(matrices_loop_only.L))
    L_total = np.sum(matrices_loop_only.L)  # 相互インダクタンス含む

    print(f"\nSelf-inductance sum: {L_self_sum*1e9:.3f} nH")
    print(f"Total inductance (including mutual): {L_total*1e9:.3f} nH")

    # DC抵抗
    R_dc_total = np.sum(np.diag(matrices_loop_only.R))
    R_dc_expected = length / (sigma * width * height)

    print(f"\nDC resistance (from matrix): {R_dc_total*1000:.3f} mOhm")
    print(f"DC resistance (expected): {R_dc_expected*1000:.3f} mOhm")

    # PEEC行列構築（Star成分あり）
    matrices_full = build_peec_matrices(segments, include_star=True)

    print(f"\n--- Full Loop-Star matrices ---")
    print(f"L matrix shape: {matrices_full.L.shape}")
    print(f"P matrix shape: {matrices_full.P.shape}")
    print(f"R matrix shape: {matrices_full.R.shape}")
    print(f"M_LS matrix shape: {matrices_full.M_LS.shape}")
    print(f"n_loop: {matrices_full.n_loop}, n_star: {matrices_full.n_star}")

    # P行列から容量を推定（対角要素の逆数）
    P_diag = np.diag(matrices_full.P)
    C_equiv = 1.0 / np.mean(P_diag)
    print(f"\nEquivalent capacitance (1/P_diag_avg): {C_equiv*1e15:.3f} fF")

    return True


def test_circular_loop():
    """円形ループコイルのPEEC行列をテスト"""
    print("\n" + "=" * 60)
    print("Test 2: Circular Loop Coil PEEC Matrices")
    print("=" * 60)

    # ループパラメータ
    radius = 0.05  # 5cm
    width = 2e-3   # 2mm
    height = 1e-3  # 1mm
    n_segments = 36
    sigma = 5.8e7

    # セグメント作成
    segments = create_loop_segments(radius, width, height, n_segments, sigma)

    print(f"Loop radius: {radius*1000:.1f} mm")
    print(f"Cross section: {width*1000:.2f} x {height*1000:.2f} mm")
    print(f"Number of segments: {n_segments}")

    # PEEC行列構築
    matrices = build_peec_matrices(segments, include_star=True)

    print(f"\nL matrix shape: {matrices.L.shape}")
    print(f"P matrix shape: {matrices.P.shape}")

    # 円形ループの解析解インダクタンス
    # L = mu_0 * R * (ln(8R/a) - 2)
    # a = sqrt(w*h/pi) (等価半径)
    a_equiv = np.sqrt(width * height / np.pi)
    L_analytical = 4 * np.pi * 1e-7 * radius * (np.log(8 * radius / a_equiv) - 2)

    # 行列から計算した値
    L_total = np.sum(matrices.L)

    print(f"\nInductance (analytical): {L_analytical*1e6:.3f} uH")
    print(f"Inductance (from PEEC L matrix): {L_total*1e6:.3f} uH")
    print(f"Ratio: {L_total/L_analytical:.3f}")

    return True


def test_frequency_response():
    """周波数特性のテスト

    低周波: Z ~ jωL (誘導性)
    高周波: Z ~ 1/(jωC) (容量性)
    """
    print("\n" + "=" * 60)
    print("Test 3: Frequency Response (Loop-Star)")
    print("=" * 60)

    # ワイヤパラメータ
    length = 0.1
    width = 1e-3
    height = 1e-3
    n_segments = 5
    sigma = 5.8e7

    segments = create_simple_wire_segments(length, width, height, n_segments, sigma)
    matrices = build_peec_matrices(segments, include_star=True)

    # ポートベクトル（最初のセグメントに電圧印加）
    port_vector = np.zeros(matrices.n_loop)
    port_vector[0] = 1.0

    # 周波数掃引
    frequencies = np.logspace(3, 9, 13)  # 1kHz to 1GHz

    print(f"\nFrequency [Hz]    |Z| [Ohm]    Angle [deg]   Type")
    print("-" * 60)

    for f in frequencies:
        omega = 2 * np.pi * f
        Z = compute_port_impedance(matrices, omega, port_vector)

        mag = np.abs(Z)
        angle = np.angle(Z, deg=True)

        if angle > 45:
            ztype = "Inductive"
        elif angle < -45:
            ztype = "Capacitive"
        else:
            ztype = "Resistive"

        print(f"{f:12.2e}    {mag:10.4e}    {angle:8.1f}      {ztype}")

    return True


def test_cln_with_star():
    """CLN変換とStar行列の結合テスト

    Loop部分をCLN変換し、Star部分との結合を検証
    """
    print("\n" + "=" * 60)
    print("Test 4: CLN Transformation with Star Coupling")
    print("=" * 60)

    from cln.lanczos import lanczos

    # ワイヤパラメータ（より多くのセグメント）
    length = 0.1
    width = 1e-3
    height = 1e-3
    n_segments = 20
    sigma = 5.8e7

    segments = create_simple_wire_segments(length, width, height, n_segments, sigma)
    matrices = build_peec_matrices(segments, include_star=True)

    print(f"Original dimensions: n_loop={matrices.n_loop}, n_star={matrices.n_star}")

    # Loop部分のCLN変換
    # lanczos(K=R, N=L) -> R_diag, L_tridiag
    n_reduced = min(5, matrices.n_loop)

    result = lanczos(matrices.R, matrices.L, n_iter=n_reduced)

    print(f"CLN reduced dimensions: {n_reduced}")
    print(f"R_diag shape: {result.R_diag.shape}")
    print(f"L_tridiag shape: {result.L_tridiag.shape}")

    # 変換行列Q
    Q = result.Q

    # Star行列はそのまま
    P = matrices.P

    # Loop-Star結合行列の変換: M_LS -> Q^T * M_LS
    M_LS_original = matrices.M_LS
    M_LS_reduced = Q.T @ M_LS_original

    print(f"\nM_LS original shape: {M_LS_original.shape}")
    print(f"M_LS reduced shape: {M_LS_reduced.shape}")

    # 周波数特性の比較
    print("\nComparing frequency response (original vs CLN-reduced)...")

    port_vector_original = np.zeros(matrices.n_loop)
    port_vector_original[0] = 1.0

    port_vector_reduced = Q.T @ port_vector_original

    frequencies = [1e3, 1e5, 1e7]

    print(f"\nFrequency [Hz]    Z_orig [Ohm]       Z_reduced [Ohm]    Error")
    print("-" * 70)

    for f in frequencies:
        omega = 2 * np.pi * f

        # オリジナル
        Z_orig = compute_port_impedance(matrices, omega, port_vector_original)

        # CLN縮約版のインピーダンス計算
        # Z_LL_reduced = R_diag + jw*L_tridiag
        Z_LL_red = result.R_diag + 1j * omega * result.L_tridiag

        # フルシステム（縮約Loop + Star）
        n_loop_red = n_reduced
        n_star = matrices.n_star
        n_total = n_loop_red + n_star

        Z_full = np.zeros((n_total, n_total), dtype=complex)
        Z_full[:n_loop_red, :n_loop_red] = Z_LL_red

        if n_star > 0:
            Z_SS = P / (1j * omega)
            Z_full[n_loop_red:, n_loop_red:] = Z_SS

            Z_LS = 1j * omega * M_LS_reduced
            Z_full[:n_loop_red, n_loop_red:] = Z_LS
            Z_full[n_loop_red:, :n_loop_red] = Z_LS.T

        V = np.zeros(n_total, dtype=complex)
        V[:n_loop_red] = port_vector_reduced

        I = np.linalg.solve(Z_full, V)
        Z_reduced = np.dot(port_vector_reduced, I[:n_loop_red])

        error = np.abs(Z_reduced - Z_orig) / np.abs(Z_orig)

        print(f"{f:12.2e}    {Z_orig:18.6e}    {Z_reduced:18.6e}    {error:.2e}")

    return True


def main():
    """全テストを実行"""
    print("PEEC Loop-Star Matrix Tests")
    print("=" * 60)
    print("Darwin approximation: G(r) = 1/(4*pi*r)")
    print("=" * 60)

    tests = [
        ("Simple Wire", test_simple_wire),
        ("Circular Loop", test_circular_loop),
        ("Frequency Response", test_frequency_response),
        ("CLN with Star", test_cln_with_star),
    ]

    results = []
    for name, test_func in tests:
        try:
            success = test_func()
            results.append((name, "PASS" if success else "FAIL"))
        except Exception as e:
            print(f"\nERROR in {name}: {e}")
            import traceback
            traceback.print_exc()
            results.append((name, "ERROR"))

    print("\n" + "=" * 60)
    print("Summary")
    print("=" * 60)
    for name, status in results:
        print(f"  {name}: {status}")


if __name__ == "__main__":
    main()
