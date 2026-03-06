# MMM (Magnetic Moment Method) pybind11 リファクタリング

## ステータス: 実装完了 (2026-01-14)

pybind11ベースのmmm_coreモジュールが実装完了しました。

## 1. 実装完了した機能

### 1.1 mmm_coreモジュール

| クラス | 機能 | ステータス |
|--------|------|-----------|
| `MMMBuilder` | 要素追加、相互作用行列構築 | ✅ 完了 |
| `MMMSolver` | LU/BiCGSTABソルバー | ✅ 完了 |
| `MMMHACApKSolver` | H行列加速ソルバー | ✅ 完了 |
| `MMMFieldComputer` | B/H磁場計算 | ✅ 完了 |
| `compute_chi_from_bh` | B-H曲線からχ計算 | ✅ 完了 |
| `check_convergence` | 収束判定 | ✅ 完了 |

### 1.2 APIサマリー

```python
import mmm_core
import numpy as np

# === 要素追加 ===
builder = mmm_core.MMMBuilder()

# 四面体メッシュから追加
tet_vertices = np.array([[0,0,0], [1,0,0], [0.5,0.866,0], [0.5,0.289,0.816]])
tet_elements = np.array([[0, 1, 2, 3]])
builder.add_tetrahedra_from_mesh(tet_vertices, tet_elements)

# 六面体メッシュから追加
hex_vertices = np.array([...])  # 8頂点
hex_elements = np.array([[0,1,2,3,4,5,6,7]])
builder.add_hexahedra_from_mesh(hex_vertices, hex_elements)

# === 相互作用行列構築 ===
N, dof_offset = builder.build()
# N: (total_dof, total_dof) 相互作用行列
# dof_offset: (n_elem+1,) DOFオフセット

# === 線形ソルブ ===
solver = mmm_core.MMMSolver()
solver.set_matrix(N, dof_offset)

# LU直接法
mu_r = 1000
inv_chi = np.full(n_elem, 1.0 / (mu_r - 1))
H_ext = np.array([0, 0, 1e5, ...])  # total_dof
M = solver.solve_lu(inv_chi, H_ext, chi_per_element=True)

# BiCGSTAB反復法
M, iterations = solver.solve_bicgstab(inv_chi, H_ext, tol=1e-8, max_iter=1000, chi_per_element=True)

# === H行列加速ソルブ (HACApK) ===
hacapk_solver = mmm_core.MMMHACApKSolver()
hacapk_solver.set_from_builder(builder, dof_offset)
success = hacapk_solver.build_hmatrix(inv_chi, eps=1e-4, leaf_size=10, eta=2.0, print_level=0)
y = hacapk_solver.matvec(x)  # H行列-ベクトル積 (内部で置換処理済み)
stats = hacapk_solver.get_stats()  # 圧縮統計

# === 磁場計算 ===
field_computer = mmm_core.MMMFieldComputer()
field_computer.set_elements_from_builder(builder)
obs_points = np.array([[0.1, 0.0, 0.0], [0.0, 0.1, 0.0]])
B = field_computer.compute_b_field(M, obs_points)
H = field_computer.compute_h_field(M, obs_points)

# === 非線形ヘルパー ===
chi = mmm_core.compute_chi_from_bh(M_flat, H_flat, bh_curve)
max_change = mmm_core.check_convergence(B_new, B_old, B_sat)
```

## 2. ファイル構成

```
src/
├── core/
│   ├── rad_mmm_matrices.h       # MMMBuilder, MMMSolver
│   ├── rad_mmm_matrices.cpp     # 実装
│   ├── rad_mmm_hacapk.h         # MMMHACApKSolver
│   └── rad_mmm_hacapk.cpp       # H行列ソルバー実装
├── lib/
│   └── rad_mmm_matrices_api.cpp # pybind11バインディング
└── ext/HACApK/
    ├── cHACApK_base.h           # HACApK基本構造体
    ├── cHACApK_base.c           # HACApK ACA+実装
    ├── cHACApK_cpp.h            # C++ opaqueポインタAPI
    └── cHACApK_cpp_impl.c       # Cラッパー実装
```

## 3. HACApK H行列ソルバー

### 3.1 特徴

- **ACA+圧縮**: 低ランク近似によるメモリ削減 (O(N log N))
- **高速matvec**: H行列-ベクトル積の高速計算
- **置換処理**: HACApK内部のlod配列による置換を自動処理

### 3.2 パラメータ

| パラメータ | デフォルト | 説明 |
|-----------|-----------|------|
| `eps` | 1e-4 | ACA許容誤差 |
| `leaf_size` | 10 | リーフブロック最小サイズ |
| `eta` | 2.0 | 許容性パラメータ |

### 3.3 統計情報

```python
stats = hacapk_solver.get_stats()
# stats = {
#     'n_dof': int,          # 総自由度
#     'n_leaves': int,       # リーフブロック数
#     'n_lowrank': int,      # 低ランクブロック数
#     'n_dense': int,        # 密ブロック数
#     'max_rank': int,       # 最大ランク
#     'memory_mb': float,    # H行列メモリ使用量 [MB]
#     'dense_memory_mb': float,  # 密行列換算メモリ [MB]
#     'compression': float,  # 圧縮率
#     'build_time': float,   # 構築時間 [s]
# }
```

## 4. 削除されたレガシーコード

以下のファイルは削除されました：

- `src/ext/HACApK/cHACApK_radia.c` - 旧Radia専用ラッパー
- `src/ext/HACApK/cHACApK_radia.h` - 旧Radia専用ヘッダー
- `test_hacapk_quick.py` - 開発用テスト
- `test_hacapk_simple.py` - 開発用テスト
- `test_import.py` - 開発用テスト

## 5. 将来拡張 (Phase 2)

### 5.1 PEEC連携
- H行列をPEEC導体問題にも適用
- CLN縮約との統合

### 5.2 Netgen OCC連携
- 要素可視化機能
- メッシュビューワー統合

---

**実装完了日**: 2026-01-14
**対象バージョン**: Radia 1.6
