# Compact AMS + COCR ソルバー

ヘッダーオンリーのC++前処理+Krylovソルバーで、HCurl問題に対応する。
外部依存なし。すべての並列処理はNGSolve TaskManagerで行う。

**実数**（静磁場）と**複素数**（渦電流）の両方の系に対応する。
`Update()` メソッドにより非線形問題のNewton反復に使用可能。

## 問題設定

渦電流有限要素離散化により以下の連立方程式が得られる：

```
A x = b,   A = K + jw * sigma * M
```

ここで K はcurl-curl剛性行列、M は質量行列（導体領域）、
sigma は電気伝導率である。系行列 A は**複素対称**（A^T = A、エルミートではない）
であり、COCRの短漸化式Krylov法が適用可能となる。

## アーキテクチャ

```
                        COCR Solver
                  (complex symmetric Krylov)
                           |
                  ComplexCompactAMS
              (fused Re/Im preconditioner)
                           |
           +---------+---------+---------+
           |         |         |         |
      FineSmooth  GradCorr  NodalCorr  CompactAMG
      (l1-Jacobi)  (G^T AG)  (Pi^T A Pi)  (V-cycle)
```

3つの層が連携して動作する：

1. **COCR**（外側Krylov法）- A^T = A を利用しO(n)メモリで動作、リスタート不要
2. **ComplexCompactAMS**（前処理）- 全レベルでRe/Im融合処理
3. **CompactAMG**（粗グリッドソルバー）- DualMult融合付き古典的AMG

## クイックスタート

### 実数系（静磁場）

```python
import sparsesolv_ngsolve as ssn
from ngsolve import *
from ngsolve.krylovspace import CGSolver

# Real HCurl system (static magnetic field)
fes = HCurl(mesh, order=1, nograds=True, dirichlet="outer")
u, v = fes.TnT()

a = BilinearForm(fes)
a += curl(u) * curl(v) * dx + 1e-6 * u * v * dx
a.Assemble()

f = LinearForm(fes)
# ... assemble right-hand side ...
f.Assemble()

gfu = GridFunction(fes)

# Gradient and coordinates
G_mat, h1_fes = fes.CreateGradient()
coord_x = [mesh.ngmesh.Points()[i+1][0] for i in range(mesh.nv)]
coord_y = [mesh.ngmesh.Points()[i+1][1] for i in range(mesh.nv)]
coord_z = [mesh.ngmesh.Points()[i+1][2] for i in range(mesh.nv)]

pre = ssn.CompactAMSPreconditioner(
    mat=a.mat, grad_mat=G_mat,
    freedofs=fes.FreeDofs(),
    coord_x=coord_x, coord_y=coord_y, coord_z=coord_z)

with TaskManager():
    inv = CGSolver(a.mat, pre, tol=1e-10, maxiter=500)
    gfu.vec.data = inv * f.vec
```

### 複素数系（渦電流）

```python
import sparsesolv_ngsolve as ssn
from ngsolve import *

omega = 2 * pi * 30000          # 30 kHz
sigma_cu = 5.96e7               # copper conductivity [S/m]
nu = 1.0                        # reluctivity 1/mu_0 (air/copper region)

# Complex HCurl space
fes = HCurl(mesh, order=1, nograds=True, dirichlet="outer", complex=True)
u, v = fes.TnT()

# Complex system matrix: K + jw*sigma*M
a = BilinearForm(fes)
a += nu * curl(u) * curl(v) * dx
a += 1j * omega * sigma_cu * u * v * dx("conductor")
a += 1e-6 * u * v * dx                       # regularization
a.Assemble()

f = LinearForm(fes)
# ... assemble right-hand side (e.g., source current) ...
f.Assemble()

gfu = GridFunction(fes)

# Real auxiliary space (non-complex)
fes_real = HCurl(mesh, order=1, nograds=True, dirichlet="outer")
u_r, v_r = fes_real.TnT()

# Real SPD auxiliary matrix: nu*curl-curl + |omega|*sigma*M + eps*M
a_real = BilinearForm(fes_real)
a_real += nu * curl(u_r) * curl(v_r) * dx
a_real += omega * sigma_cu * u_r * v_r * dx("conductor")
a_real += 1e-6 * u_r * v_r * dx
a_real.Assemble()

# Gradient and coordinates
G_mat, h1_fes = fes_real.CreateGradient()
coord_x = [mesh.ngmesh.Points()[i+1][0] for i in range(mesh.nv)]
coord_y = [mesh.ngmesh.Points()[i+1][1] for i in range(mesh.nv)]
coord_z = [mesh.ngmesh.Points()[i+1][2] for i in range(mesh.nv)]

# Preconditioner (uses real auxiliary matrix, applied to complex system)
pre = ssn.ComplexCompactAMSPreconditioner(
    a_real_mat=a_real.mat,
    grad_mat=G_mat,
    freedofs=fes_real.FreeDofs(),
    coord_x=coord_x, coord_y=coord_y, coord_z=coord_z,
    ndof_complex=fes.ndof,
    cycle_type=1,        # 1=01210 (default), 7=0201020
    print_level=0)

# Solve with COCR (complex symmetric Krylov)
with TaskManager():
    inv = ssn.COCRSolver(
        a.mat, pre,
        freedofs=fes.FreeDofs(),
        maxiter=500, tol=1e-10)
    gfu.vec.data = inv * f.vec

print(f"Converged in {inv.iterations} iterations")
```

### なぜ a_real が必要か

AMS前処理は、実際の系が複素数であっても**実数SPD（対称正定値）** 行列上で動作する。
その理由は以下の通りである：

1. **AMG粗視化には実数行列が必要。** PMIS粗視化と古典的補間は、実数行列要素に基づく
   結合強度（strength-of-connection）を使用する。複素数要素には強度閾値のための
   自然な順序が存在しない。

2. **実数補助行列はスペクトル特性を反映する。** 複素系
   `A = K + jw*sigma*M` に対して、実数補助行列 `A_real = K + |omega|*sigma*M + eps*M`
   は同じスパースパターンを持ち、類似のスペクトル特性を有する。`|omega|*sigma`
   の項（`j` なし）により、導体領域からの質量行列の寄与が前処理に反映される。

3. **Re/Im融合適用。** 前処理は同一の実数AMS V-cycleを残差の実部と虚部の両方に
   同時に適用する（DualMult）。これによりメモリ帯域コストが半減する。

複素系行列 `a.mat` は外側のCOCRソルバーによる行列ベクトル積にのみ使用される。
前処理が複素行列を直接参照することはない。

### よくある間違い

**複素前処理に誤ったfreedofsを渡す：**
```python
# WRONG: passing complex freedofs to preconditioner
pre = ssn.ComplexCompactAMSPreconditioner(
    a_real_mat=a_real.mat, ...,
    freedofs=fes.FreeDofs())        # fes is complex -> wrong size

# CORRECT: use real freedofs
pre = ssn.ComplexCompactAMSPreconditioner(
    a_real_mat=a_real.mat, ...,
    freedofs=fes_real.FreeDofs())    # fes_real is non-complex
```

**TaskManagerの欠落：**
```python
# WRONG: no TaskManager -> single-threaded, 5-10x slower
inv = ssn.COCRSolver(a.mat, pre, ...)
gfu.vec.data = inv * f.vec

# CORRECT: wrap solve in TaskManager
with TaskManager():
    inv = ssn.COCRSolver(a.mat, pre, ...)
    gfu.vec.data = inv * f.vec
```

**nograds=Trueの欠落：**
```python
# WRONG: without nograds, HCurl includes gradient DOFs that
# make the system much larger and harder to precondition
fes = HCurl(mesh, order=1, dirichlet="outer")

# CORRECT: nograds=True removes gradient DOFs (kernel of curl)
fes = HCurl(mesh, order=1, nograds=True, dirichlet="outer")
```

**ComplexCompactAMSPreconditionerに誤った行列型を渡す：**
```python
# WRONG: passing complex matrix
pre = ssn.ComplexCompactAMSPreconditioner(
    a_real_mat=a.mat, ...)          # a is complex BilinearForm

# CORRECT: pass real auxiliary matrix
pre = ssn.ComplexCompactAMSPreconditioner(
    a_real_mat=a_real.mat, ...)     # a_real is real BilinearForm
```

**order > 1 は未対応：**
```python
# WRONG: order=2 not supported by CompactAMS
fes = HCurl(mesh, order=2, nograds=True, dirichlet="outer")

# CORRECT: use order=1
fes = HCurl(mesh, order=1, nograds=True, dirichlet="outer")
```

## パラメータ

**ComplexCompactAMSPreconditioner**:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `a_real_mat` | (required) | 実数SPD補助行列 (K + eps*M + \|omega\|*sigma*M) |
| `grad_mat` | (required) | 離散勾配 G (HCurl -> H1) |
| `freedofs` | None | HCurl空間の自由度マスク |
| `coord_x/y/z` | (required) | 頂点座標（長さ = ndof_h1） |
| `ndof_complex` | 0 | 複素自由度数（0 = 行列から自動導出） |
| `cycle_type` | 1 | AMSサイクル: 1 = "01210", 7 = "0201020" |
| `print_level` | 0 | 出力レベル（0 = 無出力） |
| `correction_weight` | 1.0 | 部分空間補正の重み |
| `subspace_solver` | 0 | 0 = CompactAMG, 1 = SparseCholesky（診断用） |
| `num_smooth` | 1 | l1-Jacobiの平滑化回数 |

**COCRSolver**:

| Parameter | Default | Description |
|-----------|---------|-------------|
| `mat` | (required) | 系行列（実数または複素数のBaseMatrix） |
| `pre` | (required) | 前処理（BaseMatrix） |
| `freedofs` | None | 自由度マスク |
| `maxiter` | 500 | 最大反復回数 |
| `tol` | 1e-8 | 相対収束判定閾値 |
| `printrates` | False | 収束情報の出力 |

## ソースファイル

| File | Lines | Description |
|------|-------|-------------|
| `compact_amg.hpp` | ~900 | 古典的AMG: PMIS粗視化、l1-Jacobi、V-cycle、DualMult |
| `compact_ams.hpp` | ~750 | AMS前処理: Hiptmair-Xu補助空間法 |
| `complex_compact_ams.hpp` | ~340 | 複素系向けRe/Im融合ComplexCompactAMS |
| `sparsesolv_solvers.hpp` | ~300 | 融合ベクトル演算付きCOCRSolverNGS |

すべてのファイルは `include/sparsesolv/preconditioners/`（AMG/AMS）および
`include/sparsesolv/ngsolve/`（COCRソルバー、Pythonエクスポート）に配置されている。

## アルゴリズムの詳細

### CompactAMG (compact_amg.hpp)

以下のコンポーネントを持つ古典的代数的マルチグリッド法：

- **粗視化**: PMIS (Parallel Modified Independent Set)、De Sterck et al. (2006) に基づく
- **補間**: 古典的直接補間 (Ruge-Stueben 1987)
- **平滑化**: l1-Jacobi（打ち切りl1ノルム: `diag + 0.5*|off_diag|`、`4/3*diag` で上限）
- **最粗レベル**: NGSolve SparseCholesky直接法ソルバー
- **結合強度閾値**: 0.25（`amg_theta` で設定可能）

**DualMult**: すべてのAMGレベルで2つの右辺を同時に処理する。
SpMVはメモリ帯域律速であるため、行列の各行を1回のロードで2つのベクトルに適用すれば
帯域コストが半減する：

```cpp
// DualResidual: res1 = b1 - A*x1, res2 = b2 - A*x2  (one matrix pass)
ParallelFor(n, [&](size_t i) {
    auto cols = A.GetRowIndices(i);
    auto vals = A.GetRowValues(i);
    double d1 = 0, d2 = 0;
    for (int j = 0; j < cols.Size(); j++) {
        int c = cols[j]; double v = vals[j];
        d1 += v * x1[c]; d2 += v * x2[c];
    }
    res1[i] = b1[i] - d1;
    res2[i] = b2[i] - d2;
});
```

DualMultには以下の融合版が含まれる: V-cycle、l1-Jacobi平滑化、残差計算、SpMV、MultAdd。
各AMGレベルは二重作業ベクトル（`residual2`、`tmp2`、`correction2`）を保持する。

### CompactAMS (compact_ams.hpp)

Hiptmair-Xu (2007) のHCurl系向け補助空間前処理を実装する。

**AMS V-cycle** (cycle_type=1, パターン "01210"):

```
1. FineSmooth     : l1-Jacobi on A_bc (fine HCurl grid)
2. GradCorrect    : G^T * A_bc * G  solved by CompactAMG (H1 gradient space)
3. NodalCorrect   : Pi_d^T * A_bc * Pi_d  solved by CompactAMG x3 (H1 nodal, d=x,y,z)
4. GradCorrect    : repeat
5. FineSmooth     : repeat
```

**Pi行列の構成**（Kolev-Vassilevskiの公式）：

```
Gd[e] = (G * coord_d)[e]           // edge vector in d-direction
Pi_d[e,v] = |G[e,v]| * 0.5 * Gd[e] // nodal interpolation operator
```

**節点補正**は加法的に行われる: Pix、Piy、Piz補正は同一の残差から計算され、
加算される。これにより乗法的アプローチと比較して、サイクル当たり2回の
細レベルSpMVが節約される。

**セットアップ**では `ParallelFor(4, ...)` により4つのCompactAMGインスタンス
（G、Pix、Piy、Piz）を並列に構築する。

### ComplexCompactAMS (complex_compact_ams.hpp)

複素渦電流系向けにCompactAMSをラップする。主要な最適化は
**Re/Im融合処理**である：細レベルのSpMV演算で行列データを1回ロードし、
実部と虚部を同時に処理する。

**融合演算**:

| Operation | Matrix | Conventional | Fused |
|-----------|--------|-------------|-------|
| FineSmooth | A_bc (HCurl) | 2 passes | 1 pass |
| Residual | A_bc (HCurl) | 2 passes | 1 pass |
| Restrict | G^T, Pi_d^T | 2 passes each | 1 pass each |
| Prolongate | G, Pi_d | 2 passes each | 1 pass each |
| AMG V-cycle | coarse levels | 2 calls | DualMult (fused) |

**FusedFineSmooth** は競合状態を回避するため2フェーズで実行する必要がある：

```
Phase 1: res[i] = b[i] - (A*x)[i]   for ALL i   (reads OLD x)
Phase 2: x[i] += res[i] / l1[i]     for ALL i   (no data dependency)
```

1フェーズ方式（残差計算と更新を同時に行う）ではカオス的Jacobiとなる：
ParallelForが別スレッドの更新前後で x[j] を読む可能性があり、発散を引き起こす。

**AMG DualMult** はすべての粗レベルソルブ（勾配 + 3節点）に使用される。
部分空間ソルバーがCompactAMGでない場合（例：デバッグ用SparseCholesky）は、
逐次的なRe/Im処理にフォールバックする。

### COCR ソルバー (sparsesolv_solvers.hpp)

複素対称系のための共役直交共役残差法（Conjugate Orthogonal Conjugate Residual）。

COCRは複素対称性（A^T = A）を短漸化式（3項）で利用し、
反復回数によらず6本の作業ベクトルのみで動作する。

**融合ベクトル演算**によりメモリトラフィックを削減する：

```cpp
// 3 updates -> 1 ParallelFor
ParallelFor(n, [=](size_t i) {
    sol_d[i] += alpha * p_d[i];    // solution update
    r_d[i]   -= alpha * q_d[i];    // residual update
    rt_d[i]  -= alpha * qt_d[i];   // preconditioned residual update
});

// 2 updates -> 1 ParallelFor
ParallelFor(n, [=](size_t i) {
    p_d[i] = rt_d[i] + beta * p_d[i];   // direction update
    q_d[i] = t_d[i]  + beta * q_d[i];   // A*direction update
});
```

`GetVectorData<SCAL>()` による生ポインタアクセスにより、NGSolveの
BaseVector仮想ディスパッチのオーバーヘッドを回避する。

## 性能最適化の履歴

Hirumaメッシュ mesh1_3.5T（197,395自由度、168反復）での計測結果：

| Optimization | ms/iter | Cumulative Speedup |
|-------------|---------|-------------------|
| Baseline (sequential Re/Im AMS, NGSolve vector ops) | 49.2 | 1.00x |
| + Fused Re/Im fine-level SpMV | 41.1 | 1.20x |
| + Fused COCR vector updates (ParallelFor) | 36.6 | 1.34x |
| + Fused AMG DualMult (all coarse levels) | 32.7 | **1.50x** |

## Newton反復用 Update()

`CompactAMSPreconditioner`（実数）と `ComplexCompactAMSPreconditioner`（複素数）の
両方が `Update()` をサポートする。各ステップで系行列が変化する非線形ソルバーに対応する。

### Update() で保持・再構築される要素

| Category | Components | Rebuilt? |
|----------|-----------|----------|
| **幾何情報**（初回のみ） | Pi行列、G転置、作業ベクトル | NO |
| **行列依存**（Update毎） | A_bc、Galerkin射影 (A_G, A_Pi)、AMG階層、l1ノルム | YES |

この分離により、メッシュと座標にのみ依存し行列値に依存しない
幾何データの冗長な再計算を回避する。

### API

```python
# Update with current matrix (matrix was modified in-place)
pre.Update()

# Update with a new matrix object
pre.Update(new_mat)
```

### 実数系（静磁場 + CG）

```python
import sparsesolv_ngsolve as ssn
from ngsolve import *
from ngsolve.krylovspace import CGSolver

fes = HCurl(mesh, order=1, nograds=True, dirichlet="outer")
G_mat, h1_fes = fes.CreateGradient()
coord_x = [mesh.ngmesh.Points()[i+1][0] for i in range(mesh.nv)]
coord_y = [mesh.ngmesh.Points()[i+1][1] for i in range(mesh.nv)]
coord_z = [mesh.ngmesh.Points()[i+1][2] for i in range(mesh.nv)]

a = BilinearForm(fes)
a += nu_cf * curl(u) * curl(v) * dx + 1e-6 * u * v * dx
a.Assemble()

pre = ssn.CompactAMSPreconditioner(a.mat, G_mat,
    freedofs=fes.FreeDofs(),
    coord_x=coord_x, coord_y=coord_y, coord_z=coord_z)

# Newton loop
for k in range(max_newton):
    a.Assemble()          # Reassemble with updated nu(B)
    pre.Update(a.mat)     # Rebuild preconditioner (geometry preserved)
    inv = CGSolver(a.mat, pre, tol=1e-8, maxiter=500)
    delta = inv * rhs
    gfu.vec.data += delta
```

### 複素数系（渦電流 + COCR）

```python
pre = ssn.ComplexCompactAMSPreconditioner(a_real.mat, G_mat,
    freedofs=fes_real.FreeDofs(),
    coord_x=coord_x, coord_y=coord_y, coord_z=coord_z,
    ndof_complex=fes.ndof)

# Newton loop
for k in range(max_newton):
    a_real.Assemble()
    pre.Update(a_real.mat)
    inv = ssn.COCRSolver(a_complex.mat, pre, tol=1e-10)
    delta = inv * rhs
    gfu.vec.data += delta
```

## 主要な設計判断

### なぜCOCRか

渦電流系 A = K + jw*sigma*M は A^T = A（複素対称）を満たす。
COCRはこの対称性を短漸化式（3項）で利用し、
反復回数によらず6本の作業ベクトルのみで動作する。

### なぜRe/Im融合か

現代のハードウェアではSpMVはメモリ帯域律速である。行列 A_bc はReとImの
両方の部分で同一であるため、各行を1回ロードして両方の内積を計算すれば
帯域コストが半減する。適用対象は以下の通りである：
- 細レベル平滑化（l1-Jacobi）
- 残差計算（b - A*x）
- 制限（P^T * residual）
- 延長（x += P * correction）
- すべての粗レベルでのAMG V-cycle（DualMult経由）

### なぜGauss-SeidelではなくI1-Jacobiか

l1-Jacobiは完全に並列化可能（行間のデータ依存性なし）であり、
TaskManager並列化に最適である。Gauss-Seidelは逐次的またはマルチカラー
順序付けが必要であり、複雑性が増す。l1ノルム打ち切り
（`diag + 0.5*|off_diag|`、`4/3*diag` で上限）により、
過減衰なしにロバストな平滑化を実現する。

### なぜ加法的節点補正か

Pix、Piy、Piz補正はすべて同一の細レベル残差を使用し、1回だけ計算する。
これにより乗法的アプローチ（各Pi補正の間で残差を再計算する）と比較して、
AMSサイクル当たり2回の細レベルSpMV評価が節約される。

## ベンチマーク結果

**テスト問題**: Hiruma渦電流モデル (SA-26-001)。
銅導体 (sigma = 5.96e7 S/m) と強磁性鉄心 (mu_r = 1000)。
周波数 30 kHz、HCurl次数 1、収束判定閾値 1e-10。

**環境**: Windows Server 2022, Intel Xeon (8 cores), MSVC 2022, MKL 2024.2。

### メッシュサイズに対するスケーリング

| Mesh | Elements | HCurl DOFs | H1 DOFs | Iters | Setup [s] | Solve [s] | Total [s] | ms/iter | Memory [MB] |
|------|----------|-----------|---------|-------|-----------|-----------|-----------|---------|-------------|
| mesh1_2.5T | 130,460 | 155,527 | 23,731 | 144 | 0.7 | 3.7 | 4.5 | 25.8 | 368 |
| mesh1_3.5T | 166,198 | 197,395 | 29,901 | 168 | 0.9 | 6.3 | 7.3 | 37.7 | 460 |
| mesh1_4.5T | 203,403 | 241,129 | 36,473 | 210 | 1.3 | 9.2 | 10.5 | 43.8 | 543 |
| mesh1_5.5T | 280,278 | 331,595 | 49,643 | 249 | 1.9 | 14.3 | 16.2 | 57.3 | 725 |
| mesh1_20.5T | 1,227,241 | 1,441,102 | 211,337 | 499 | 24.9 | 197.7 | 222.6 | 396.2 | 2,933 |

すべてのケースで収束を確認（`true ||b-Ax||/||b|| < 2e-10`）。

### IC前処理との比較 (ABMC-ICCG)

IC（不完全コレスキー分解）前処理にABMC（代数的ブロックマルチカラー）
並列順序付けを適用。30 kHz、tol=1e-10、maxiter=20,000。

| Mesh | DOFs | Method | Iters | Total [s] | Status |
|------|------|--------|-------|-----------|--------|
| mesh1_3.5T | 197k | Compact AMS + COCR | 168 | 7.2 | converged |
| mesh1_3.5T | 197k | ABMC-ICCG | 17,178 | 438.4 | not converged (res=2.8e-10) |

IC前処理はHCurl離散化に固有のcurl-curl零空間を処理できない。
AMSは勾配および節点補助空間補正によりこれを解決する（Hiptmair-Xu 2007）。

### AMS vs ICCGのスケーリング（鉄心渦電流問題）

材料コントラストを含む鉄心渦電流問題における Compact AMS + COCR と ICCG の比較
（mu_r=1000 鉄心、sigma=1e6 S/m、エアギャップ付き、30 kHz）。

このベンチマークはAMSの主要な利点である**メッシュサイズに依存しない反復回数**を示す。

| Mesh | HCurl DOFs | Method | Iters | Solve [s] | Notes |
|------|-----------|--------|-------|-----------|-------|
| Small | 2,728 | ICCG | 97 | 0.03 | |
| Small | 2,728 | AMS + COCR | 46 | 0.11 | 1.5x fewer iters |
| Medium | 6,382 | ICCG | 147 | 0.09 | |
| Medium | 6,382 | AMS + COCR | 52 | 0.13 | 2.8x fewer iters |
| Large | 19,357 | ICCG | 234 | 0.59 | |
| Large | 19,357 | AMS + COCR | 52 | 0.23 | **4.5x fewer iters, 2.6x faster** |

**主な観察結果**:

1. **AMS反復回数は安定**（46 → 52 → 52）、メッシュサイズに依存しない
2. **ICCG反復回数は増加**（97 → 147 → 234）、メッシュ細分化に伴いO(h^{-1})スケーリング
3. **AMSの計算時間面での優位性は問題サイズとともに拡大** -- 19K自由度でAMSは既に2.6倍高速
4. より大規模な問題（>100K自由度）では、ICCGの反復増加によりAMSが不可欠となる（上記IC前処理との比較を参照：197K自由度で17,178反復）

設定: `maxh` = 0.08/0.06/0.04、`order=1`、`nograds=True`、`tol=1e-8`、`maxiter=10000`。

### EMD前処理との比較 (Hiruma SA-26-001)

EMD（辺ベース磁場分解法）の論文結果
(Hiruma, SA-26-001, 3,670,328自由度, 30 kHz)：

| Method | Iters | Time [s] | ms/iter |
|--------|-------|----------|---------|
| IC only | 15,838 | 5964.8 | 376.7 |
| EMD (IC + AMG V-cycle) | 4,069 | 1716.8 | 422.0 |
| EMD (IC + AMG W-cycle) | 2,935 | 1552.9 | 529.1 |
| EMD (IC + GenEO-DDM, 24 domains) | 1,004 | 550.8 | 548.6 |

本手法の 1.44M自由度 (mesh1_20.5T) での結果: 499反復、222.6秒（8 CPUコア）。
自由度数が異なるため（1.44M vs 3.67M）、直接比較はできない。

## 参考文献

- R. Hiptmair, J. Xu. "Nodal Auxiliary Space Preconditioning in H(curl) and H(div) Spaces."
  SIAM J. Numer. Anal. 45(6), 2007.
- T. Kolev, P. Vassilevski. "Parallel Auxiliary Space AMG Solver for H(div) Problems."
  J. Comput. Math. 27(5), 2009.
- J. Ruge, K. Stueben. "Algebraic Multigrid." In Multigrid Methods, 1987.
- H. De Sterck, U. Yang, J. Heys. "Reducing Complexity in Parallel Algebraic Multigrid
  Preconditioners." SIAM J. Matrix Anal. Appl. 27(4), 2006.
- T. Sogabe, S.-L. Zhang. "A COCR method for solving complex symmetric linear systems."
  J. Comput. Appl. Math. 199(2), 2007.

## ベンチマークの再現方法

```bash
cd examples/hiruma

# Compact AMS + COCR (single mesh)
python bench_compact_ams.py mesh1_3.5T

# Compact AMS + COCR (all standard meshes)
python bench_compact_ams.py --all

# Comparison with ABMC-ICCG
python bench_ams_vs_abmc.py mesh1_3.5T

# Results saved to results_compact_ams.json, results_ams_vs_abmc.json
```
