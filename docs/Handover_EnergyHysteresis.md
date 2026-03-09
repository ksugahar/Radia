# Handover: Table-Based Energy Hysteresis Model

## 概要

B-input Play モデルの形状関数をテーブルデータとして C++ Energy-Based Hysteresis モデルに取り込む作業。インフラは完成、変換ロジックが未解決。

## 完了した作業

### 1. C++ テーブル補間エンジン (解析関数を完全削除)

`A_s * tan(pi/2 * r / J_s)` を完全に削除し、テーブル補間に置換。

**OperatorTable 構造体** (`rad_material_def.h`):
```cpp
struct OperatorTable {
    std::vector<double> r;     // |J| grid [0, r_max]
    std::vector<double> f;     // f_k(r) = U_k'(r)
    std::vector<double> U;     // 台形積分で前計算
    std::vector<double> df;    // 中央差分で前計算
    int n;
    double r_max;
};
```

**InterpolateTable** (`rad_material_impl.cpp`): 二分探索 + 線形補間。Uk/dUk/d2Uk 全てがテーブル参照に変更。

### 2. API チェーン

```
Python                          C interface                    C++ core
rad.MatEnergyHysteresis(    ->  RadMatEnergyHysteresis(    ->  SetEnergyHysteresisMaterial(
  K, chi, f_k_tables, eps)       K, chi, r_flat, f_flat,       K, chi, r_tables,
                                  table_sizes, eps)              f_tables, eps)
```

- `f_k_tables`: Python list of K tuples `(r_array, f_array)`
- pybind11 で flatten して C interface に渡す
- C interface で `std::vector<std::vector<double>>` に再構築

### 3. Forward/Inverse 命名規則

B-input Play の自然な方向に合わせて統一済み:
- **Forward(B)** = B -> H (Schur complement, O(K))
- **Inverse(H)** = H -> B (各 J_k 独立, O(K))
- `MatMvsH(mat, 'm', H)` は Inverse 方向 (H -> M)

### 4. テストインフラ

実データ (Jiles-Atherton 40ループ, BH.mat) からテスト fixture を作成:
- `examples/hysteresis/binput_play_fixture.npz` (K=80 full + K=20 subset)
- `examples/hysteresis/verify_cpp_hysteresis.py` (5 tests, ALL PASS)
- `examples/hantila_solver/test_binput_cpp.py` (5 tests, ALL PASS)

### 5. 変更したファイル一覧

| ファイル | 変更内容 |
|---------|---------|
| `src/core/rad_material_def.h` | OperatorTable, m_tables, コンストラクタ変更 |
| `src/core/rad_material_impl.cpp` | InterpolateTable, PrecomputeOperatorTable, Uk/dUk/d2Uk |
| `src/core/rad_application.h` | SetEnergyHysteresisMaterial 宣言 |
| `src/core/rad_c_interface.cpp` | flat array -> vector 再構築 |
| `src/lib/radentry.h` | RadMatEnergyHysteresis 宣言 |
| `src/lib/radentry.cpp` | ラッパー関数 |
| `src/radia/radia_pybind.cpp` | MatEnergyHysteresis pybind |
| `src/radia/hysteresis_io.py` | convert_play_to_energy (As/Js 削除) |

## 未解決: B-input Play -> Energy 変換

### 問題

`convert_play_to_energy()` が B-input Play の形状関数 f_k をそのまま Energy モデルに渡しているが、二つのモデルでは f_k の物理的意味が異なる:

```
B-input Play model:
  H = sum_k f_k(|p_k|) * p_k/|p_k|
  → f_k は H (磁場) への寄与

Energy-based model:
  M = sum_k J_k    (J_k は各オペレータの磁化)
  各 J_k は U_k'(|J_k|) = f_k(|J_k|) で決まるエネルギーランドスケープで最適化
  → f_k は U_k' (エネルギー微分) であり、M への間接的な寄与
```

### 症状

K=20 operators, r_max=2.0 T の場合:
- 各 J_k が最大 2.0 T → 合計 J ≤ 40 T (現実: J_s ~ 2 T)
- H=100 A/m で mu_r ~ 67,000 (現実: 100-1000)
- BEM ソルバで NaN 発生

### 考えられるアプローチ

1. **スケーリング**: f_k を K で割る? r_max を J_s/K に制限?
   - 単純だが物理的根拠が薄い

2. **逆変換**: B-H 関係 H(B) = sum f_k(|p_k|) から H(M) 関係を導出
   - B = mu_0*(H + M) の関係を使って変換
   - 非線形逆問題になる

3. **B-input Play を C++ に直接実装** (Energy モデルを迂回)
   - 最も確実だが、Energy モデルの利点 (Jacobian の解析的計算) を失う

4. **Egger 論文の再精読**
   - Energy-based formulation が H-input Play 用に設計されている可能性
   - B-input 用の Energy formulation が別途必要かもしれない

### 参考文献

- Egger, Engertsberger, Schafelner: "Efficient evaluation of forward and inverse energy-based magnetic hysteresis operators", MAGCON-25-07-0171
- Potter, Schmulian: B-input Play model (MATLAB コード: `VectorPlayModel/`)
- Sugahara, Hane: B-input Play model (近畿大学)

## データファイル

| ファイル | 内容 | 場所 |
|---------|------|------|
| `BH.mat` | 40 Jiles-Atherton ループ (1x40 struct) | `W:\...\VectorPlayModel\` |
| `HdataLUT.mat` | H lookup table (61x30) | 同上 |
| `Hfunc.mat` | Shape functions (109x54) | 同上 |
| `B.mat` | Reference trajectory (3770 pts) | 同上 |
| `B_input_case_*.mat` | B-input test cases | `W:\...\2024_03_08_H-input_B-input\` |
| `binput_play_fixture.npz` | テスト fixture (K=80/20) | `examples/hysteresis/` |

### BH.mat の読み方

```python
from scipy.io import loadmat
data = loadmat(path)  # squeeze_me=False
BH = data['BH']       # shape (1, 40), structured array
for i in range(40):
    B = np.array(BH[0, i]['B'][0, :], dtype=float)  # (n_pts,)
    H = np.array(BH[0, i]['H'][0, :], dtype=float)  # (n_pts,)
```

### .hys ファイルの問題

`load_hys()` は現在 JilesAtherton.hys 等の読み込みに失敗する。原因: 各ループ間にヘッダ行 (2行) が挿入されており、`np.loadtxt()` がパースできない。修正が必要。

## pip vs local .pyd

Build.ps1 は `src/radia/` にコピーするが、Python は pip インストール版を優先的にロードする:
```
cp build-msvc/_radia_pybind.cp312-win_amd64.pyd \
   "C:/Program Files/Python312/Lib/site-packages/radia/_radia_pybind.pyd"
```

## 設計方針

- **Python 前処理** (build_shape_functions, convert_play_to_energy) → C++ はテーブルを受け取るだけ
- C++ は B-input Play モデルの詳細を知る必要なし
- `hysteresis_io.py` が全ての I/O と変換を担当
