# ELF と Radia のメッシュ規約の比較

**作成日**: 2026-01-30
**更新日**: 2026-01-30
**目的**: ELF_MAGIC と Radia の六面体（ヘキサヘドロン）要素における規約の違いを文書化

## 概要

Radia と ELF_MAGIC は同じ MSC（Magnetic Surface Charge）法を使用していますが、いくつかの規約の違いがあります。これらの違いを理解することで、両ソルバー間で結果を正しく比較できます。

### 主要な違いの要約

| 項目 | ELF (Fortran) | Radia (C++) |
|------|---------------|-------------|
| 行列格納順序 | Column-major | Row-major |
| 対角成分符号 | 負 | 負（v1.4.4以降同一）|
| 面の順序 | Face 0-5 (kkh array) | Face 0-5 (kkh array) |

### 変換結果

fullモデル（52要素）での検証:
- **最大絶対誤差**: 3.59e-4
- **最大相対誤差**: 0.0053%
- **結論**: 数値積分精度の範囲内で一致

## 行列格納順序（Row-Major vs Column-Major）

### 結論: 転置が必要

ELF（Fortran）とRadia（C++）は行列の格納順序が異なります。

| 言語 | 格納順序 | A[i,j]のメモリ位置 |
|------|----------|-------------------|
| Fortran (ELF) | Column-major | A[j*M + i] |
| C++ (Radia) | Row-major | A[i*N + j] |

### ELF行列ファイルの読み込み

ELFの`.mat`ファイルをPython/NumPyで読み込む際、転置が必要です:

```python
# 正しい読み込み方法
elf_matrix_raw = read_elf_matrix(mat_file)
elf_matrix = elf_matrix_raw.T  # 転置して使用
```

### 転置なしの場合の問題

転置せずに使用すると:
- 対角成分は正しい
- 非対角成分が転置される
- 結果として行列全体が転置された状態になる

## 頂点順序（Vertex Ordering）

### 結論: 同一

両ソルバーとも**標準 CHEXA 頂点順序**（Netgen/Nastran 互換）を使用しています。

```
      v7 -------- v6
      /|         /|
     / |        / |
   v4 -------- v5 |
    |  v3 -----|-- v2
    | /        | /
    |/         |/
   v0 -------- v1

頂点位置:
  v0: (x-, y-, z-)  下面左前
  v1: (x+, y-, z-)  下面右前
  v2: (x+, y+, z-)  下面右後
  v3: (x-, y+, z-)  下面左後
  v4: (x-, y-, z+)  上面左前
  v5: (x+, y-, z+)  上面右前
  v6: (x+, y+, z+)  上面右後
  v7: (x-, y+, z+)  上面左後
```

**注意**: ELF の MEG ファイルでは節点番号が異なる順序で定義されていますが、要素の連結性（connectivity）を通じて標準 CHEXA 順序に変換されます。

## 面順序（Face Ordering）

### 結論: 異なる

| Radia面 | Radia方向 | ELF面 | ELF方向 |
|---------|----------|-------|--------|
| 0 | z- (底面) | 4 | z- |
| 1 | z+ (上面) | 5 | z+ |
| 2 | y- (前面) | 0 | y- |
| 3 | y+ (後面) | 2 | y+ |
| 4 | x- (左面) | 1 | x- |
| 5 | x+ (右面) | 3 | x+ |

### 置換行列

Radia から ELF への変換:

```
置換: Radia面 -> ELF面
  0 -> 4
  1 -> 5
  2 -> 0
  3 -> 2
  4 -> 1
  5 -> 3

置換ベクトル: [4, 5, 0, 2, 1, 3]

置換行列 P:
[[0 0 1 0 0 0]
 [0 0 0 0 1 0]
 [0 0 0 1 0 0]
 [0 0 0 0 0 1]
 [1 0 0 0 0 0]
 [0 1 0 0 0 0]]
```

## 符号規約（Sign Convention）

### 結論: 同一（v1.4.4以降）

| 項目 | Radia | ELF |
|------|-------|-----|
| 対角成分 | 負 (-1.5684) | 負 (-1.5684) |
| 行列定義 | N = 相互作用行列 | A = N - diag(1/chi) |

### 物理的解釈

- **Radia**: 負の対角成分（ELF互換、自己減磁場は磁化を減少させる）
- **ELF**: 負の対角成分（エネルギー定式化）

**注意**: v1.4.3以前のRadiaは正の対角成分を使用していました。

## 変換式

ELF 行列 A と Radia 行列 N の関係:

```
A_ELF_raw = P @ N_Radia^T @ P^T

ここで:
  A_ELF_raw = ELFファイルから直接読み込んだ行列（転置なし）
  N_Radia^T = Radiaの相互作用行列の転置
  P = 置換行列（上記参照）
  @ = 行列乗算
```

**注意**: ELF行列を転置して使用する場合は `A_ELF_T = P @ N_Radia @ P^T` となります。

### 検証結果

**単一立方体要素** (100mm角、mu_r=1000):
- 最大絶対差: 1.7e-05
- 結論: 数値精度内で一致

**単一直方体要素** (100x90x120mm、mu_r=1000):
- 最大絶対差: 2.8e-05
- 結論: 数値精度内で一致

**fullモデル** (52要素、直方体):
- 最大絶対差: 3.6e-04 (0.005%)
- 結論: 数値精度内で一致

**重要**: 置換は要素形状（立方体/直方体）に依存せず、面の命名規約のみに依存します。

## コード例

### Python による行列変換

```python
import numpy as np

def create_permutation_matrix():
    """Radia -> ELF の置換行列を作成"""
    perm = [4, 5, 0, 2, 1, 3]  # Radia面 -> ELF面
    P = np.zeros((6, 6))
    for radia_i, elf_i in enumerate(perm):
        P[elf_i, radia_i] = 1.0
    return P

def radia_to_elf_matrix(radia_N):
    """Radia の N 行列を ELF の A 行列形式に変換（raw形式）"""
    P = create_permutation_matrix()
    # 変換: A_raw = P @ N^T @ P^T (v1.4.4以降、符号は同一)
    A_elf_raw = P @ radia_N.T @ P.T
    return A_elf_raw

def elf_raw_to_radia_matrix(elf_A_raw):
    """ELF の raw A 行列を Radia の N 行列形式に変換"""
    P = create_permutation_matrix()
    # 逆変換: N^T = P^T @ A_raw @ P
    # よって N = (P^T @ A_raw @ P)^T = P^T @ A_raw^T @ P
    N_radia = P.T @ elf_A_raw.T @ P
    return N_radia

def read_elf_matrix_and_convert(mat_file):
    """ELFファイルを読み込み、Radia形式に変換"""
    from elf_nastran_reader import read_elf_matrix
    elf_raw = read_elf_matrix(mat_file)
    # 6x6ブロックごとに変換
    n = elf_raw.shape[0]
    n_elem = n // 6
    radia_N = np.zeros_like(elf_raw)
    for i in range(n_elem):
        for j in range(n_elem):
            i0, i1 = i*6, (i+1)*6
            j0, j1 = j*6, (j+1)*6
            radia_N[i0:i1, j0:j1] = elf_raw_to_radia_matrix(
                elf_raw[i0:i1, j0:j1])
    return radia_N
```

## 推奨事項

### Radia ユーザー向け

1. Radia は Netgen 標準に従っています
2. ELF との結果比較時は上記変換を適用してください
3. 面順序の違いは DOF（自由度）の並び替えに影響します

### 将来の互換性

- Radia は Netgen/NGSolve との統合を優先
- ELF 形式への変換は Python ユーティリティで提供
- 新規開発では Radia（Netgen）規約を推奨

## 参考ファイル

- **解析スクリプト**: `examples/electromagnet/mu=1000/analyze_elf_face_ordering.py`
- **単一要素比較**: `examples/electromagnet/mu=1000/compare_single_element.py`
- **ELF MEG ファイル例**: `S:\ELF_MAGIC\...\single\ELF_MAGIC.meg`

## 更新履歴

| 日付 | 変更内容 |
|------|----------|
| 2026-01-30 | 初版作成。面順序と符号規約の違いを特定。 |
| 2026-01-30 | 対角成分符号をELF互換に変更（正→負）。変換式から負号を削除。 |
