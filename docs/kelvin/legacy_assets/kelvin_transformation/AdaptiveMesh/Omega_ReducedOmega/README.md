# Omega-Reduced Omega Formulation with Kelvin Transformation

3D磁場解析のための適応メッシュ細分化手法の研究。Omega法（スカラーポテンシャル）とKelvin変換を組み合わせ、無限遠境界条件を有限領域で扱う。

## 概要

本研究では、以下の3つのメッシュ細分化手法を比較検討した：

1. **Uniform Refinement (Refine_all_elements)**: 全要素を一様に細分化
2. **Adaptive Refinement (Refine_with_zz_estimator)**: ZZ誤差推定器に基づくDörflerマーキングによる適応細分化
3. **Metric-based Remeshing (metric_based)**: ZZ誤差推定器＋Dörflerマーキング＋メトリックベース再メッシュのハイブリッド手法

## 問題設定

### 物理パラメータ
- 透磁率: μ₀ = 4π × 10⁻⁷ H/m（真空）
- 比透磁率: μᵣ = 100（磁性体）
- 外部磁場: H_s = (0, 0, 1) A/m（z方向一様場）

### 形状
- **Sphere_3D**: 半径0.3mの磁性球（1/8モデル、対称性利用）
- **Cylinder_3D**: 半径0.3m、半高さ0.5mの磁性円筒（1/8モデル、対称性利用）

### Kelvin変換
- Kelvin半径: R_K = 1.5m
- 内部領域（r < R_K）: 通常の空気領域
- 外部領域（r > R_K）: Kelvin変換により r' = R_K²/r にマッピング
- 変換後の透磁率: μ_kelvin = μ₀ × (R_K/r)⁶

## 定式化

### Omega法（スカラーポテンシャル）
磁場 H = -∇Ω + H_s として、スカラーポテンシャル Ω を求める。

弱形式:
∫ μ ∇Ω · ∇v dx = ∫ μ H_s · ∇v dx

### 境界条件
- 対称面: ∂Ω/∂n = 0（Neumann）
- Kelvin境界: 周期境界条件（内側・外側をペアリング）

## ZZ誤差推定器

### 実装の要点

1. **L2射影によるフラックス回復**（補間ではない）
   ```python
   fes_flux = HDiv(mesh, order=order-1)
   # L2射影: (σ, τ) = (B, τ) for all τ in H(div)
   ```

2. **B_bounded の定義**（air_outerでμ₀を使用）
   ```python
   B_bounded_dict = {
       "magnetic": (mu_r * mu0) * H_pert,
       "air_inner": mu0 * H_pert,
       "air_outer": mu0 * H_pert  # μ_kelvinではなくμ₀
   }
   ```
   air_outerでμ_kelvinを使用すると、r→R_Kでμ_kelvin→∞となり特異性が生じるため、誤差推定にはμ₀を使用。

3. **誤差計算**
   ```python
   element_errors = Integrate(|B - B_recovered|², mesh, element_wise=True)
   total_error = sqrt(sum(element_errors))
   ```

## エネルギー計算

摂動エネルギー W = (1/2) ∫ μ |H_pert|² dV を各領域で計算：

```python
# 磁性体領域
energy_magnetic = 0.5 * (mu_r * mu0) * |grad(Ω) - H_s|²

# 空気領域（内部）
energy_air_inner = 0.5 * mu0 * |grad(Ω_reduced)|²

# 空気領域（外部、Kelvin変換後）
energy_air_outer = 0.5 * mu_kelvin * |grad(Ω_kelvin)|²
```

## 結果

### Sphere_3D（磁性球）

| Order | Method | Final DOF | Final Error | Iterations |
|-------|--------|-----------|-------------|------------|
| p=2 | Uniform | ~27,000 | ~2.0e-06 | 3 |
| p=2 | Adaptive (ZZ) | ~137,000 | ~2.2e-06 | 7 |
| p=2 | Metric-based | ~117,000 | ~2.6e-06 | 30 |
| p=3 | Uniform | ~88,000 | ~1.3e-06 | 3 |
| p=3 | Adaptive (ZZ) | ~160,000 | ~1.1e-06 | 6 |
| p=3 | Metric-based | ~126,000 | ~7.8e-07 | 11 |
| p=4 | Uniform | ~27,000 | ~1.3e-06 | 2 |
| p=4 | Adaptive (ZZ) | ~155,000 | ~6.9e-07 | 5 |
| p=4 | Metric-based | ~133,000 | ~6.6e-07 | 10 |

### Cylinder_3D（磁性円筒）

| Order | Method | Final DOF | Final Error | Iterations |
|-------|--------|-----------|-------------|------------|
| p=2 | Uniform | 27,122 | 1.97e-06 | 3 |
| p=2 | Adaptive (ZZ) | 136,727 | 2.22e-06 | 7 |
| p=2 | Metric-based | 117,116 | 2.63e-06 | 30 |
| p=3 | Uniform | 87,698 | 1.35e-06 | 3 |
| p=3 | Adaptive (ZZ) | 159,866 | 1.11e-06 | 6 |
| p=3 | Metric-based | 126,036 | 7.76e-07 | 11 |
| p=4 | Uniform | 27,122 | 1.26e-06 | 2 |
| p=4 | Adaptive (ZZ) | 155,034 | 6.88e-07 | 5 |
| p=4 | Metric-based | 132,964 | 6.59e-07 | 10 |

### エネルギー収束値

両形状とも、全手法で以下のエネルギー値に収束：
- Total Energy: ~1.659e-05 J（1/8モデル、8倍で全体）
- magnetic: ~1.582e-05 J
- air_inner: ~7.51e-07 J
- air_outer: ~2.48e-08 J

## 理論的収束率

3D問題における有限要素法の理論的収束率：
- Error ∝ N^(-p/3)（Nは自由度数、pは多項式次数）

| Order | 理論的傾き |
|-------|-----------|
| p=2 | O(N^(-2/3)) |
| p=3 | O(N^(-1)) |
| p=4 | O(N^(-4/3)) |

## ファイル構造

Per-order runner scripts under `order=*/.../*.py` have been promoted to the
result-bearing archive notebook `docs/kelvin/kelvin_adaptive_mesh_archive.ipynb`
with full source text and SHA-256 hashes in
`docs/kelvin/kelvin_adaptive_mesh_archive_results.json`.  The aggregate
`compare_convergence.py` scripts and `CubeMesh.py` were then promoted into the
final remaining-examples archive:
`docs/kelvin/kelvin_remaining_examples_archive.ipynb` and
`docs/kelvin/kelvin_remaining_examples_archive_results.json`.

```
Omega_ReducedOmega/
├── README.md                    # 本ファイル
├── debug.md                     # デバッグ記録・実装詳細
├── Sphere_3D/
│   └── order=*/                 # archived historical runner outputs
└── Cylinder_3D/
    └── order=*/                 # archived historical runner outputs
```

## 出力ファイル

各シミュレーションは以下のファイルを出力：
- `*.mat`: MATLAB形式のデータファイル（収束履歴、エネルギー値等）
- `*.vtu`: VTK形式のメッシュ・解データ（ParaViewで可視化可能）
- `*.png`: 収束プロット画像

## 実行方法

For retired Python source, use
`docs/kelvin/kelvin_adaptive_mesh_archive.ipynb` and
`docs/kelvin/kelvin_adaptive_mesh_archive_results.json` for per-order runners,
and `docs/kelvin/kelvin_remaining_examples_archive_results.json` for the final
aggregate scripts.

## 依存ライブラリ

- NGSolve/Netgen
- NumPy
- SciPy
- Matplotlib

## 参考文献

1. Zienkiewicz, O.C. and Zhu, J.Z., "A simple error estimator and adaptive procedure for practical engineering analysis", Int. J. Numer. Methods Eng., 1987
2. Dörfler, W., "A convergent adaptive algorithm for Poisson's equation", SIAM J. Numer. Anal., 1996
3. Kelvin transformation for unbounded domain problems
