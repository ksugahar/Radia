# Python APIリファレンス

## インポート

```python
import ngsolve  # Must be imported first (loads shared libraries)
from sparsesolv_ngsolve import (
    # Preconditioners (IC/SGS)
    ICPreconditioner,      # Incomplete Cholesky preconditioner
    SGSPreconditioner,     # Symmetric Gauss-Seidel preconditioner

    # Compact AMG/AMS Preconditioners
    CompactAMGPreconditioner,            # Classical AMG (for H1)
    CompactAMSPreconditioner,            # Real HCurl (magnetostatics)
    CompactAMSPreconditionerImpl,        # Real AMS type (with Update())
    ComplexCompactAMSPreconditioner,     # Complex Re/Im fused AMS
    ComplexCompactAMSPreconditionerImpl, # Complex AMS type (with Update())
    has_compact_ams,                     # Check Compact AMG/AMS availability

    # Iterative Solvers
    SparseSolvSolver,      # Unified iterative solver (ICCG/SGSMRTR/CG/COCR)
    SparseSolvResult,      # Solve result
    COCRSolver,            # COCR (complex symmetric systems, native C++)
    GMRESSolver,           # GMRES (non-symmetric systems, left-preconditioned)
)
```

ファクトリ関数は `mat.IsComplex()` により行列の型（実数/複素数）を自動判定する。

---

## ICPreconditioner

不完全コレスキー（IC）前処理。

### コンストラクタ

```python
pre = ICPreconditioner(mat, freedofs=None, shift=1.05)
```

| パラメータ | 型 | デフォルト | 説明 |
|-----------|------|---------|-------------|
| `mat` | `SparseMatrix` | - | SPD行列（実数/複素数） |
| `freedofs` | `BitArray` or `None` | `None` | 自由度。`None` の場合は全自由度を対象とする |
| `shift` | `float` | `1.05` | シフトパラメータ |

### プロパティ

| プロパティ | 型 | 説明 |
|----------|------|-------------|
| `shift` | `float` | シフトパラメータ（読み書き可） |
| `use_abmc` | `bool` | ABMCオーダリングの有効化（読み書き可） |
| `abmc_block_size` | `int` | ABMCブロックサイズ（読み書き可） |

### メソッド

| メソッド | 説明 |
|--------|-------------|
| `Update()` | 前処理を再計算する（行列変更後に呼び出すこと） |

### 使用例

```python
pre = ICPreconditioner(a.mat, freedofs=fes.FreeDofs(), shift=1.05)
inv = CGSolver(a.mat, pre, tol=1e-10)
gfu.vec.data = inv * f.vec
```

---

## SGSPreconditioner

対称ガウス・ザイデル（SGS）前処理。

### コンストラクタ

```python
pre = SGSPreconditioner(mat, freedofs=None)
```

| パラメータ | 型 | デフォルト | 説明 |
|-----------|------|---------|-------------|
| `mat` | `SparseMatrix` | - | SPD行列 |
| `freedofs` | `BitArray` or `None` | `None` | 自由度 |

### メソッド

| メソッド | 説明 |
|--------|-------------|
| `Update()` | 前処理を再計算する |

### 使用例

```python
pre = SGSPreconditioner(a.mat, freedofs=fes.FreeDofs())
inv = CGSolver(a.mat, pre, tol=1e-10)
gfu.vec.data = inv * f.vec
```

---

## SparseSolvSolver

統合反復ソルバー。ICCG、SGSMRTR、CG、COCRメソッドに対応。
`BaseMatrix` として使用可能であり、`gfu.vec.data = solver * f.vec` の形式で利用できる。

### コンストラクタ

```python
solver = SparseSolvSolver(mat, method="ICCG", freedofs=None,
                           tol=1e-10, maxiter=1000, shift=1.05,
                           save_best_result=True,
                           save_residual_history=False,
                           printrates=False, conjugate=False,
                           use_abmc=False, abmc_block_size=4,
                           abmc_num_colors=4, abmc_reorder_spmv=False,
                           abmc_use_rcm=False)
```

### パラメータ

| パラメータ | 型 | デフォルト | 説明 |
|-----------|------|---------|-------------|
| `mat` | `SparseMatrix` | - | SPD行列 |
| `method` | `str` | `"ICCG"` | `"ICCG"`, `"SGSMRTR"`, `"CG"`, `"COCR"` |
| `freedofs` | `BitArray` | `None` | 自由度 |
| `tol` | `float` | `1e-10` | 収束判定閾値 |
| `maxiter` | `int` | `1000` | 最大反復回数 |
| `shift` | `float` | `1.05` | IC分解のシフト量（ICCG用） |
| `save_best_result` | `bool` | `True` | 最良解の追跡 |
| `save_residual_history` | `bool` | `False` | 残差履歴の記録 |
| `printrates` | `bool` | `False` | 収束情報の出力 |
| `conjugate` | `bool` | `False` | 共役内積（エルミート系用） |

### ABMC関連パラメータ

| パラメータ | 型 | デフォルト | 説明 |
|-----------|------|---------|-------------|
| `use_abmc` | `bool` | `False` | ABMCオーダリングの有効化 |
| `abmc_block_size` | `int` | `4` | ブロックあたりの行数 |
| `abmc_num_colors` | `int` | `4` | 目標カラー数 |
| `abmc_reorder_spmv` | `bool` | `False` | ABMC空間でのSpMV実行 |
| `abmc_use_rcm` | `bool` | `False` | 事前にRCMバンド幅縮小を適用 |

### プロパティ（構築後に設定可能）

| プロパティ | 型 | 説明 |
|----------|------|-------------|
| `auto_shift` | `bool` | IC分解の自動シフト調整 |
| `diagonal_scaling` | `bool` | 対角スケーリング |
| `divergence_check` | `bool` | 停滞検出時の早期終了 |
| `divergence_threshold` | `float` | 発散検出の閾値 |
| `divergence_count` | `int` | 発散と判定するまでの連続不良反復回数 |
| `last_result` | `SparseSolvResult` | 直近の求解結果 |

### メソッド

| メソッド | 引数 | 戻り値 | 説明 |
|--------|-----------|---------|-------------|
| `Solve(rhs, sol)` | `BaseVector`, `BaseVector` | `SparseSolvResult` | 初期推定値付きで求解 |

### 使用例

```python
# Use as BaseMatrix (zero initial guess)
solver = SparseSolvSolver(a.mat, method="ICCG",
                           freedofs=fes.FreeDofs(), tol=1e-10)
gfu.vec.data = solver * f.vec

# Use the Solve method to obtain detailed results
result = solver.Solve(f.vec, gfu.vec)
print(f"Converged: {result.converged}, Iters: {result.iterations}")
```

---

## SparseSolvResult

求解結果を格納する構造体。

| フィールド | 型 | 説明 |
|-------|------|-------------|
| `converged` | `bool` | ソルバーが収束したかどうか |
| `iterations` | `int` | 反復回数 |
| `final_residual` | `float` | 最終相対残差 |
| `residual_history` | `list[float]` | 各反復における残差（`save_residual_history=True` の場合） |

```python
result = solver.Solve(f.vec, gfu.vec)
if result.converged:
    print(f"Converged in {result.iterations} iterations")
    print(f"Final residual: {result.final_residual:.2e}")
```

---

## CompactAMSPreconditioner

実数HCurlシステム（静磁界のcurl-curl + 質量項）用のCompact AMS前処理。

ヘッダオンリーC++実装。外部ライブラリ不要。
非線形ソルバー（Newton反復）に対応: `Update()` は幾何情報を保持しつつ、行列依存部分のみを再構築する。

### コンストラクタ

```python
pre = CompactAMSPreconditioner(
    mat, grad_mat, freedofs=None,
    coord_x=[], coord_y=[], coord_z=[],
    cycle_type=1, print_level=0,
    subspace_solver=0, num_smooth=1)
```

| パラメータ | 型 | デフォルト | 説明 |
|-----------|------|---------|-------------|
| `mat` | `SparseMatrix` (real) | - | 実数SPD行列（curl-curl + 質量項） |
| `grad_mat` | `SparseMatrix` (real) | - | 離散勾配 G (HCurl -> H1) |
| `freedofs` | `BitArray` or `None` | `None` | 自由度 |
| `coord_x` | `list[float]` | `[]` | 頂点のx座標 |
| `coord_y` | `list[float]` | `[]` | 頂点のy座標 |
| `coord_z` | `list[float]` | `[]` | 頂点のz座標 |
| `cycle_type` | `int` | `1` | AMSサイクル型（1="01210", 7="0201020"） |
| `print_level` | `int` | `0` | 出力の詳細レベル |
| `subspace_solver` | `int` | `0` | 0=CompactAMG, 1=SparseCholesky |
| `num_smooth` | `int` | `1` | l1-Jacobiスムージングのステップ数 |

### メソッド

| メソッド | 説明 |
|--------|-------------|
| `Update()` | 前処理を再構築する（幾何情報を保持し、行列依存部分のみ再計算） |
| `Update(new_mat)` | 新しい行列で前処理を再構築する |

### 使用例

```python
import sparsesolv_ngsolve as ssn
from ngsolve import *
from ngsolve.krylovspace import CGSolver

# Magnetostatic problem
fes = HCurl(mesh, order=1, nograds=True, dirichlet="outer")
G_mat, h1_fes = fes.CreateGradient()
coord_x = [mesh.ngmesh.Points()[i+1][0] for i in range(mesh.nv)]
coord_y = [mesh.ngmesh.Points()[i+1][1] for i in range(mesh.nv)]
coord_z = [mesh.ngmesh.Points()[i+1][2] for i in range(mesh.nv)]

a = BilinearForm(fes)
a += nu * curl(u) * curl(v) * dx + 1e-6 * u * v * dx
a.Assemble()

pre = ssn.CompactAMSPreconditioner(
    a.mat, G_mat, freedofs=fes.FreeDofs(),
    coord_x=coord_x, coord_y=coord_y, coord_z=coord_z)

inv = CGSolver(a.mat, pre, tol=1e-10, maxiter=500)
gfu.vec.data = inv * f.vec
```

### Newton反復での使用方法（非線形問題）

```python
pre = ssn.CompactAMSPreconditioner(a.mat, G_mat,
    freedofs=fes.FreeDofs(),
    coord_x=coord_x, coord_y=coord_y, coord_z=coord_z)

for k in range(max_newton):
    a.Assemble()          # Reassemble matrix with B-H curve
    pre.Update(a.mat)     # Rebuild preconditioner (geometry retained)
    inv = CGSolver(a.mat, pre, tol=1e-8, maxiter=500)
    delta = inv * rhs
    gfu.vec.data += delta
```

---

## ComplexCompactAMSPreconditioner

複素渦電流システム用のCompact AMS前処理（Re/Im融合処理）。

ヘッダオンリーC++実装。外部ライブラリ不要。
Re/Im融合SpMVにより行列データを一度だけロードし、実部と虚部を同時に処理する。
対称前処理（l1-Jacobiスムーザ）-- **COCRSolverの使用を推奨**。

### コンストラクタ

```python
pre = ComplexCompactAMSPreconditioner(
    a_real_mat, grad_mat, freedofs=None,
    coord_x=[], coord_y=[], coord_z=[],
    ndof_complex=0, cycle_type=1, print_level=0,
    correction_weight=1.0, subspace_solver=0, num_smooth=1)
```

| パラメータ | 型 | デフォルト | 説明 |
|-----------|------|---------|-------------|
| `a_real_mat` | `SparseMatrix` (real) | - | 実数SPD補助行列 (K + eps*M + \|omega\|*sigma*M) |
| `grad_mat` | `SparseMatrix` (real) | - | 離散勾配 G (HCurl -> H1) |
| `freedofs` | `BitArray` or `None` | `None` | 自由度 |
| `coord_x` | `list[float]` | `[]` | 頂点のx座標 |
| `coord_y` | `list[float]` | `[]` | 頂点のy座標 |
| `coord_z` | `list[float]` | `[]` | 頂点のz座標 |
| `ndof_complex` | `int` | `0` | 複素自由度数（0=行列から自動検出） |
| `cycle_type` | `int` | `1` | AMSサイクル型（1="01210", 7="0201020"） |
| `print_level` | `int` | `0` | 出力の詳細レベル |
| `correction_weight` | `float` | `1.0` | 補正重み |
| `subspace_solver` | `int` | `0` | 0=CompactAMG, 1=SparseCholesky |
| `num_smooth` | `int` | `1` | l1-Jacobiスムージングのステップ数 |

### メソッド

| メソッド | 説明 |
|--------|-------------|
| `Update()` | 前処理を再構築する（幾何情報を保持し、行列依存部分のみ再計算） |
| `Update(new_a_real)` | 新しい実数補助行列で前処理を再構築する |

### 使用例

```python
import sparsesolv_ngsolve as ssn

pre = ssn.ComplexCompactAMSPreconditioner(
    a_real_mat=a_real.mat, grad_mat=G_mat,
    freedofs=fes_real.FreeDofs(),
    coord_x=cx, coord_y=cy, coord_z=cz,
    ndof_complex=fes.ndof, cycle_type=1, print_level=0)

with TaskManager():
    inv = ssn.COCRSolver(a.mat, pre, freedofs=fes.FreeDofs(),
                          maxiter=500, tol=1e-10)
    gfu.vec.data = inv * f.vec
```

### Newton反復での使用方法（非線形渦電流問題）

```python
pre = ssn.ComplexCompactAMSPreconditioner(a_real.mat, G_mat, ...)

for k in range(max_newton):
    a_real.Assemble()
    pre.Update(a_real.mat)
    inv = ssn.COCRSolver(a_complex.mat, pre, tol=1e-10)
    delta = inv * rhs
    gfu.vec.data += delta
```

詳細は [compact_ams_cocr.md](compact_ams_cocr.md) を参照。

---

## has_compact_ams

Compact AMG/AMSサポートが利用可能かどうかを確認する関数。

```python
ssn.has_compact_ams()  # -> True if Compact AMG/AMS support is available
```

---

## COCRSolver

COCR（Conjugate Orthogonal Conjugate Residual）ソルバー。ネイティブC++実装。
複素対称系（A^T=A、エルミートではない）に最適な短漸化式Krylovソルバー。

非共役内積（x^T y）を使用。||A*r~||_2 を最小化し、COCGよりも滑らかな収束を実現する。
1反復あたりのコスト: 1回のMatVec + 1回の前処理適用（CGと同等）。

**参考文献**: Sogabe & Zhang (2007), J. Comput. Appl. Math., 199(2), 297-303.

### 使用法1: COCRSolver（外部前処理との組み合わせ）

NGSolveの `CGSolver` と同じインタフェース。AMS前処理等と組み合わせて使用する。

```python
import sparsesolv_ngsolve

inv = sparsesolv_ngsolve.COCRSolver(mat, pre, freedofs=fes.FreeDofs(),
                                    maxiter=500, tol=1e-8, printrates=False)
gfu.vec.data = inv * f.vec
print(f"COCR converged in {inv.iterations} iterations")
```

| パラメータ | 型 | デフォルト | 説明 |
|-----------|------|---------|-------------|
| `mat` | `BaseMatrix` | - | 複素対称行列 |
| `pre` | `BaseMatrix` | - | 前処理行列 |
| `freedofs` | `BitArray` or `None` | `None` | 自由度。`None` の場合は全自由度を対象とする |
| `maxiter` | `int` | `500` | 最大反復回数 |
| `tol` | `float` | `1e-8` | 収束判定閾値（相対残差） |
| `printrates` | `bool` | `False` | 収束情報の出力 |

### 使用法2: SparseSolvSolver(method="COCR")

SparseSolvSolver統合インタフェース経由。内部IC前処理を使用する。

```python
solver = sparsesolv_ngsolve.SparseSolvSolver(mat, method="COCR",
    freedofs=fes.FreeDofs(), tol=1e-10, maxiter=1000)
gfu.vec.data = solver * f.vec
result = solver.last_result
```

### プロパティ

| プロパティ | 型 | 説明 |
|----------|------|-------------|
| `iterations` | `int` | 実行された反復回数（COCRSolverのみ） |

### 収束判定基準

```
sqrt(|rt^T * r|) / sqrt(|rt0^T * r0|) < tol
```
ここで `rt = M^{-1} * r`（前処理付き残差）。
NGSolveの `CGSolver(conjugate=False)` と同等の収束判定基準。

### COCGについて

COCG（Conjugate Orthogonal CG）は `CGSolver(conjugate=False)` と数学的に等価である。
別個のクラスは提供していない。

```python
from ngsolve.krylovspace import CGSolver
inv = CGSolver(a.mat, pre, conjugate=False, maxiter=500, tol=1e-8)
```

---

## GMRESSolver

左前処理付きGMRES（Generalized Minimal Residual）ソルバー。ネイティブC++実装。
非対称行列に対応。前処理が非対称な場合に使用する。

**注意**: 対称前処理（IC、SGS、Compact AMS）の場合はCGまたはCOCRを推奨する。
GMRESは前処理が非対称な場合にのみ使用すること。

### コンストラクタ

```python
inv = sparsesolv_ngsolve.GMRESSolver(mat, pre, freedofs=None,
                                      maxiter=500, tol=1e-8,
                                      restart=30, printrates=False)
gfu.vec.data = inv * f.vec
print(f"GMRES converged in {inv.iterations} iterations")
```

| パラメータ | 型 | デフォルト | 説明 |
|-----------|------|---------|-------------|
| `mat` | `BaseMatrix` | - | 係数行列（対称/非対称） |
| `pre` | `BaseMatrix` | - | 左前処理行列 |
| `freedofs` | `BitArray` or `None` | `None` | 自由度。`None` の場合は全自由度を対象とする |
| `maxiter` | `int` | `500` | 最大反復回数 |
| `tol` | `float` | `1e-8` | 収束判定閾値（相対残差） |
| `restart` | `int` | `30` | リスタート周期 |
| `printrates` | `bool` | `False` | 収束情報の出力 |

### プロパティ

| プロパティ | 型 | 説明 |
|----------|------|-------------|
| `iterations` | `int` | 実行された反復回数 |

### 使用例

```python
import sparsesolv_ngsolve as ssn
from ngsolve.krylovspace import CGSolver

# SGS preconditioner (non-symmetric) + GMRES
pre = ssn.SGSPreconditioner(a.mat, freedofs=fes.FreeDofs())
inv = ssn.GMRESSolver(a.mat, pre, freedofs=fes.FreeDofs(),
                       maxiter=500, tol=1e-10, restart=30)
gfu.vec.data = inv * f.vec
```

---

## CompactAMGPreconditioner

古典的AMG（代数的マルチグリッド）前処理。H1 Poissonシステム用に設計されている。
PMIS粗視化 + 古典的補間 + l1-Jacobiスムージングを使用する。

CompactAMSの内部で部分空間ソルバーとしても使用される。
H1有限要素によるSPD問題に対して単独で使用可能。

### コンストラクタ

```python
pre = CompactAMGPreconditioner(mat, freedofs=None,
                                theta=0.25, max_levels=25,
                                min_coarse=50, num_smooth=1,
                                print_level=0)
```

| パラメータ | 型 | デフォルト | 説明 |
|-----------|------|---------|-------------|
| `mat` | `SparseMatrix` (real) | - | SPD行列 |
| `freedofs` | `BitArray` or `None` | `None` | 自由度 |
| `theta` | `float` | `0.25` | 強結合閾値 |
| `max_levels` | `int` | `25` | 最大レベル数 |
| `min_coarse` | `int` | `50` | 最粗レベルの最小自由度数 |
| `num_smooth` | `int` | `1` | スムージングステップ数 |
| `print_level` | `int` | `0` | 出力の詳細レベル |

### 使用例

```python
import sparsesolv_ngsolve as ssn
from ngsolve.krylovspace import CGSolver

pre = ssn.CompactAMGPreconditioner(a.mat, freedofs=fes.FreeDofs())
inv = CGSolver(a.mat, pre, tol=1e-10, maxiter=500)
gfu.vec.data = inv * f.vec
```
