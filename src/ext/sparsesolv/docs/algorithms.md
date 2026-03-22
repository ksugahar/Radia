# アルゴリズムリファレンス

## 1. BDDC (Balancing Domain Decomposition by Constraints)

### 1.1 概要

BDDCは大規模有限要素問題に対する**領域分割前処理**である [Dohrmann 2003]。
メッシュの各要素を「部分領域」として扱い、部分領域間の結合を**粗空間**を通じて処理する。
代数的定式化については [Mandel, Dohrmann, Tezaur 2005] を参照。

**特徴**:
- 反復回数がメッシュサイズ *h* にほぼ依存しない（スケーラブル）
- 要素行列から直接構築可能（組み立て済みの全体行列が不要）

**注意**: SparseSolvの独自BDDC実装は削除された。以下の記述はアルゴリズムの参考資料として残す。
代わりにNGSolve組込みBDDC (`a.mat.Inverse(fes.FreeDofs(), inverse="bddc")`) を使用すること。

### 1.2 自由度の分類

有限要素空間の自由度 (DOF) は以下の2種類に分類される:

| カテゴリ | 名称 | 役割 | NGSolve CouplingType |
|----------|------|------|---------------------|
| **Wirebasket** | 粗空間自由度 | 頂点および辺の自由度。全体的に結合 | `WIREBASKET_DOF` |
| **Interface** | 局所自由度 | 面および内部の自由度。各部分領域内で消去 | その他全て |

高次有限要素（次数2以上）では、wirebasket自由度は少数の低次自由度であり、
interface自由度が大多数を占める。BDDCの効率はこの階層構造に依存する。

### 1.3 要素レベルSchur補行列

各要素 *e* の剛性行列をwirebasket (w) とinterface (i) のブロックに分割する:

```
K^(e) = | K_ww  K_wi |
        | K_iw  K_ii |
```

**Schur補行列**（wirebasket上の等価行列）:

```
S^(e) = K_ww - K_wi * K_ii^{-1} * K_iw
```

**調和拡張**（wirebasketの値からinterfaceの値を復元）:

```
he^(e) = -K_ii^{-1} * K_iw
```

疑似コード:
```
K_ii_inv = LU_factorize(K_ii)   // PA=LU（部分ピボット付き）
he = -K_ii_inv * K_iw            // 調和拡張
S = K_ww + K_wi * he             // Schur補行列
```

### 1.4 重み付け

共有自由度（複数の要素に属する自由度）からの寄与を正しく分配するために、
**重み付け（平均化）**を適用する。

要素 *e* におけるinterface自由度 *k* の重み:

```
w_k^(e) = |K_ii^(e)(k,k)|   （K_ii対角成分の絶対値）
```

重みの正規化: 全要素の寄与を合計し、逆数をとる:

```
W_k = 1 / Sigma_e w_k^(e)
```

これらの正規化された重みを用いて、調和拡張、内部ソルブ、
調和拡張の転置の行列成分をスケーリングする。


### 1.5 粗空間ソルバー

全要素のSchur補行列から組み立てられた**wirebasket行列** S_global を
直接法で解く。

直接ソルバー（SparseCholesky、MKL PARDISOなど）を使用する。
wirebasket自由度の数は全自由度に比べて少ないため、
直接解法のコストは許容範囲内である。
NGSolve組込みBDDCでは、SparseCholeskyまたはPARDISOが粗空間ソルバーとして利用可能。

### 1.6 適用アルゴリズム（5ステップ）

BDDC前処理の適用は以下の5ステップで構成される。
入力 *r* は残差ベクトル、出力 *y* は前処理済みベクトルである。

```
ステップ1: y = r                              （コピー）
ステップ2: y += H^T * r                       （調和拡張の転置）
ステップ3: y_wb = S^{-1} * y_wb              （wirebasket成分に対する粗空間ソルブ）
ステップ4: y += K_inner * r                   （内部ソルブ: interface成分）
ステップ5: y = y + H * y                      （調和拡張）
```

ここで:
- H: 調和拡張行列（全要素のhe^(e)を重み付きで組み立て）
- H^T: その転置
- K_inner: 全要素のK_ii^{-1}を重み付きで組み立て

NGSolveでの使用例:
```python
# NGSolve組込みBDDC
inv = a.mat.Inverse(fes.FreeDofs(), inverse="bddc")
gfu.vec.data = inv * f.vec
```

---

## 2. IC分解 (不完全コレスキー分解)

### 2.1 概要

対称正定値行列 *A* の近似分解を計算する [Meijerink, van der Vorst 1977]:

```
A ~ L * D * L^T
```

- *L*: 下三角行列（Aと同じ疎性パターン）
- *D*: 対角行列

前処理の適用: y = (LDL^T)^{-1} * x を前進代入、対角スケーリング、後退代入の3段階で解く。

**実装ファイル**: `preconditioners/ic_preconditioner.hpp`

### 2.2 シフトパラメータ

安定化のためにIC分解の対角に**シフト**を加える:

```
d_i = alpha * a_ii - Sigma_{k<i} l_ik^2 * d_k^{-1}
```

ここでalphaはシフトパラメータ（デフォルト: 1.05）。

- alpha = 1.0: 標準IC(0)。不定値行列では破綻する可能性あり
- alpha = 1.05: デフォルト。やや安定性が高い
- alpha = 1.1〜1.2: 困難な問題向け

### 2.3 自動シフト

半正定値行列（curl-curl問題）では、対角成分が過小となりIC分解が破綻する場合がある。
自動シフトはこれを自動検出し、シフトを増加させて分解をやり直す:

```cpp
// ic_preconditioner.hpp: compute_ic_factorization()
if (abs_s < config_.min_diagonal_threshold && abs_orig > 0.0) {
    shift += increment;     // シフトを増加
    increment *= 2;         // 指数バックオフ
    restart = true;          // 分解をやり直す
}
```

### 2.4 対角スケーリング

行列の条件数を改善するために対角スケーリングを適用可能:

```
scaling[i] = 1 / sqrt(|A[i,i]|)
A_scaled[i,j] = scaling[i] * A[i,j] * scaling[j]
```

実装: `ICPreconditioner::compute_scaling_factors()`, `apply_scaling_to_L()`

### 2.5 三角ソルブの並列化

三角ソルブの適用（前進代入と後退代入）は本質的に逐次処理であるが、
データ依存性の解析により部分的に並列化できる。

2つのアプローチを提供する:

| アプローチ | ファイル | 並列粒度 | 特徴 |
|----------|------|---------|------|
| レベルスケジューリング | `level_schedule.hpp` | 行レベル | 単純、オーバーヘッドが少ない |
| ABMC順序付け | `abmc_ordering.hpp` | ブロックレベル | より高い並列性 |

詳細は第4節（ABMC順序付け）を参照。

---

## 3. SGS-MRTR

### 3.1 概要

SGS-MRTRは**対称ガウス-ザイデル (SGS)** 前処理を**MRTR反復法**に組み込んだソルバーである。
分割公式を用いて、前処理の前進部分 *L* と後退部分 *L^T* を個別に適用する。

**特徴**:
- IC分解が不要（分解コストなし）
- DAD変換（対角スケーリング）を内蔵
- 複素対称行列をサポート（精度は限定的）

**実装ファイル**: `solvers/sgs_mrtr_solver.hpp`

**参考文献**:
- T. Tsuburaya, Y. Mifune, T. Iwashita, E. Takahashi,
  "Numerical Experiments on Preconditioned Iterative Methods Based on the MRTR Method",
  *IEEJ Technical Meeting on Static Apparatus*, SA-12-64, 2012.
- T. Tsuburaya, Y. Okamoto, K. Fujiwara, S. Sato,
  "Improvement of the Preconditioned MRTR Method With Eisenstat's Technique
  in Real Symmetric Sparse Matrices",
  *IEEE Trans. Magnetics*, Vol. 49, No. 5, pp. 1641-1644, 2013.
  [DOI: 10.1109/TMAG.2013.2240283](https://doi.org/10.1109/TMAG.2013.2240283)

### 3.2 DAD変換

行列 *A* の対角スケーリング:

```
D = diag(|A|)^{-1/2}
A' = D * A * D
b' = D * b
x  = D * x'
```

スケーリング後、A'の全対角成分は1に近くなり、条件数が改善される。

### 3.3 分割公式

SGSの前進部分Lと後退部分L^Tを個別に適用する。
MRTR反復の各ステップにおいて:

```
rd = L^{-1} * r          （前進ソルブ）
u  = L^{-T} * rd         （後退ソルブ）
ARd = u + L^{-1}(rd - u) （M^{-1}Ardの近似）
```

### 3.4 MRTR反復

MRTRは二項再帰を用いて解を更新する最小残差型反復法である:

```
p_{k+1} = u_k + (eta_k zeta_{k-1} / zeta_k) p_k
x_{k+1} = x_k + zeta_k p_k
```

ここでzeta_kとeta_kは最適化パラメータである。

### 3.5 CSR対角位置の前提条件

SGS-MRTRの前進・後退代入は、CSR行列の各行において**対角要素が行インデックス位置に
存在すること**を前提とする。NGSolveの`SparseMatrix`はこの条件を満たすが、
任意のCSR行列では保証されない。

具体的には、行 *i* の列インデックス `col_idx[row_ptr[i]:row_ptr[i+1]]` は
ソートされており、対角成分 `A[i,i]` が存在しなければならない。

### 3.6 複素数サポートに関する注意

SGS-MRTRのzeta_k計算で除算が発生する場合、複素数の比較に`std::real()`を使用する:

```cpp
// sgs_mrtr_solver.hpp
if (std::abs(denom) < constants::DENOMINATOR_BREAKDOWN) {
    denom = (std::real(denom) >= 0) ? ... : ...;
}
```

複素対称行列 (A^T = A) の渦電流問題では、DAD変換が最適でない場合があり、
収束精度が約5%に限定されることがある。
そのような場合はICCG+ABMC、COCR、またはCompact AMS+COCRを推奨する。

---

## 4. ABMC順序付け (Algebraic Block Multi-Color Ordering)

本節ではSparseSolvで使用されるAlgebraic Block Multi-Color (ABMC) 順序付け法の
数学的背景とC++実装の詳細を記述する。ABMCは不完全コレスキー (IC) 前処理における
前進代入・後退代入の並列化手法である。

**実装ファイル:**
- `core/abmc_ordering.hpp` -- ABMCスケジュール構築（BFS集約 + グラフ彩色）
- `core/level_schedule.hpp` -- レベルスケジューリング（比較ベースライン）
- `core/parallel.hpp` -- 並列プリミティブ (`parallel_for`, `SpinBarrier`, `get_num_threads`)
- `core/solver_config.hpp` -- ABMCパラメータ設定
- `preconditioners/ic_preconditioner.hpp` -- 三角ソルブ実行、適用パス選択

**参考文献**:
- T. Iwashita, H. Nakashima, Y. Takahashi,
  "Algebraic Block Multi-Color Ordering Method for Parallel
  Multi-Threaded Sparse Triangular Solver in ICCG Method",
  *Proc. IEEE IPDPS*, 2012.

### 4.1 理論的背景

#### 4.1.1 並列三角ソルブの課題

ICCGなどの反復ソルバーでは、前処理の適用に三角系の求解が必要となる:
- 前進代入: $L y = x$ を解く
- 後退代入: $L^T z = y$ を解く

$y_i$ の計算は $L_{ij} \neq 0$ かつ $j < i$ であるすべての $y_j$ に依存する。このデータ依存性が並列化のボトルネックとなる。

#### 4.1.2 レベルスケジューリング

依存関係を解決する基本的な手法。各行のレベル（依存の深さ）を計算し、同一レベルの行を並列処理する:

```
level[i] = max(level[j] for j in L[i,:] where j < i) + 1
```

**問題点**: FEM行列はバンド幅が広く依存チェーンが深い。レベル数が数百に達し、各レベルの行数が少ないため並列効率が低い。

#### 4.1.3 マルチカラーリング

グラフ彩色により依存関係を断ち切るアプローチ。隣接するノードが同色にならないよう彩色し、同色のノードを並列処理する。ただし、古典的な点単位のマルチカラーリングは、色割り当てに基づくメモリアクセスが遠隔の位置にジャンプするため、キャッシュ局所性が悪化する。

#### 4.1.4 ABMC: ブロック + マルチカラー

ABMC [Iwashita, Nakashima, Takahashi 2012] は並列性とキャッシュ局所性を両立させるため、2レベルの階層を導入する:

1. **ブロッキング**: 近傍の行をブロックにまとめる -- キャッシュ局所性を維持
2. **カラーリング**: ブロック間の依存関係を彩色する -- 少数の同期点で並列性を実現

幾何情報は不要で、CSR行列のパターンのみで動作する。

### 4.2 実装パイプライン

ABMCスケジュールは`ICPreconditioner::setup()`中に構築される。`ABMCSchedule::build()`は4つのステージを順次実行する (`abmc_ordering.hpp`)。

#### ステージ1: 代数的ブロッキング（BFS集約）

**関数**: `make_algebraic_blocks()`

行列グラフ上のBFS（幅優先探索）を用いて近傍の行をブロックにまとめる。

```
入力: CSR行列 (row_ptr, col_idx), ブロックサイズB
出力: raw_block_list[b] = {元の行インデックスの集合}
      block_assign[i] = 行iのブロックID

1. block_assignを-1（未割当）で初期化
2. 最初の未割当行をシードとして新しいブロックを開始
3. シードの隣接行をBFSキューに追加（-2: キュー内を示す）
4. ブロックサイズBに達するまでキューから取り出してブロックに追加
5. キュー内の残りのノードを-1にリセット（次のブロック用）
6. すべての行が割り当てられるまで繰り返す
```

**番兵値の使い方:**
| 値 | 意味 |
|-------|---------|
| `-1` | 未割当 |
| `-2` | BFSキュー内（重複挿入を防止） |
| `>= 0` | ブロックIDに割り当て済み |

**計算量**: O(nnz) -- 各非零要素は最大2回走査される

ブロックサイズ `B`（デフォルト: 4、`SolverConfig::abmc_block_size`）は、ブロック内の行が逐次処理されるため、キャッシュ効率と並列粒度のバランスをとる。

#### ステージ2: ブロック隣接グラフの構築

**関数**: `build_block_graph()`

ブロック $A$ の任意の行がブロック $B$ の任意の行と非零結合を持つ場合、2つのブロックは隣接する。

```
入力: raw_block_list, block_assign, CSR行列
出力: ブロック隣接グラフ (blk_row_ptr, blk_col_idx) -- CSR形式

for each block bi:
    for each row row in bi:
        for each nonzero column j in row:
            bj = block_assign[j]
            if bj != bi and last_seen[bj] != bi:
                add bj to neighbors[bi]
                last_seen[bj] = bi  <-- 重複排除用番兵
```

**計算量**: O(nnz)。`last_seen`番兵により各ブロック対は1回のみ追加される。

#### ステージ3: 貪欲彩色

**関数**: `color_graph()`

下三角依存関係のあるブロックが異なる色を持つように、ブロック隣接グラフを彩色する。

**2段階アルゴリズム:**

```
フェーズ1: 最小色数の決定
    num_colors = target_colors  (SolverConfig::abmc_num_colors, デフォルト: 4)
    for each block i:
        lower_count = |{j : j < i and (i,j) are adjacent}|
        num_colors = max(num_colors, lower_count + 1)

フェーズ2: 禁止色集合による貪欲彩色
    forbidden[c] = i は「色cがブロックiの下三角隣接で使用済み」を意味
    for each block i (順番に):
        下三角隣接ブロック j < i の色を禁止としてマーク
        最小の未使用色を選択
        すべての色が衝突する場合: ++num_colors（新しい色を追加）
```

**重要**: `target_colors` は**下限値**であり、保証ではない。実際の色数は:

```
num_colors >= max(target_colors, 1 + max_i |{j < i : (i,j) in E}|)
```

フェーズ2でさらに増加する場合がある（フェーズ1は下三角次数のみに基づくが、彩色順序により追加の色が必要になることがある）。

**計算量**: O(次数の合計) -- 各ブロックの隣接を1回走査

#### ステージ4: 行レベル順序付けの構築

**関数**: `build_row_ordering()`

ブロック彩色の結果からフラットな行レベルの配列を構築する。

```
出力:
  color_offsets[c] -- 色cの開始ブロック位置（CSR）
  color_blocks[k] -- グローバルブロックID
  block_offsets[b] -- ブロックbの開始行位置（CSR）
  ordering[old] = new -- 旧→新の順列
  reverse_ordering[new] = old -- 新→旧の逆順列

走査順: 色0 -> 色1 -> ... -> 色C-1
  各色内: ブロック順
    各ブロック内: 元のインデックスでソート (std::sort)
```

**連続番号付け**: 走査順に新しい行インデックスが連続的に割り当てられる。ブロック $b$ の行は `block_offsets[b]` から `block_offsets[b+1]-1` まで連続するため、別途行インデックス配列は不要。三角ソルブでは `for (i = block_offsets[blk]; i < block_offsets[blk+1]; ++i)` と直接ループできる。

ブロック内ソートの目的は、同一ブロック内の行が元の行列で近接したインデックスを持つようにし、前処理適用時のTLBミスを削減してキャッシュライン利用を向上させることである。

### 4.3 並列三角ソルブのロジック

#### 4.3.1 データ構造

ABMCスケジューリング後の三角ソルブに関わるデータ:

| データ | 空間 | 説明 |
|------|-------|------|
| `L_`, `Lt_`, `inv_diag_` | ABMC空間 | IC分解結果（下三角/上三角因子、逆対角） |
| `abmc_x_perm_` | ABMC空間 | 入力ベクトル順列バッファ (サイズn) |
| `abmc_y_perm_` | ABMC空間 | 出力ベクトル順列バッファ (サイズn) |
| `work_temp_` | ABMC空間 | 中間結果ベクトル (サイズn、ABMC/レベルスケジューリング共用) |
| `work_temp2_` | ABMC空間 | 第2中間結果（対角スケーリング用のみ、サイズn） |
| `composite_perm_` | 元空間 -> ABMC | 合成順列 |
| `composite_scaling_` | 元空間 | 対角スケーリング係数（元空間） |

**ワークベクトルの共有**: `work_temp_`はABMCパスとレベルスケジューリングパスで共有される。`setup()`時に1回割り当てられ、`apply()`時のヒープ確保を排除する。

#### 4.3.2 前進代入 ($L y = x$)

`forward_substitution_abmc()` -- 3レベルループ:

```cpp
for (color c = 0; c < num_colors; ++c) {            // 逐次: 色間依存
    parallel_for(num blocks in color c, [&](bidx) {  // 並列: 同色ブロックは独立
        blk = color_blocks[blk_begin + bidx];
        row_begin = block_offsets[blk];
        row_end = block_offsets[blk + 1];
        for (i = row_begin; i < row_end; ++i) {      // 逐次: ブロック内依存
            s = x[i];                                 // ABMC空間の行インデックスを直接使用
            for (k : off-diagonal entries in row i of L_)
                s -= L_.values[k] * y[L_.col_idx[k]];
            y[i] = s / L_.values[diag];
        }
    });
    // 暗黙の同期バリア (parallel_for完了 = 次の色へ進む)
}
```

`block_offsets`が連続したインデックス範囲を定義するため、ブロック内の行は`row_begin`から`row_end`まで直接ループできる。間接参照配列は不要。

#### 4.3.3 後退代入 ($L^T z = D^{-1} y$)

`backward_substitution_abmc()` -- 前進代入の逆順:
- 色を逆順に走査: C-1 -> 0
- ブロック内の行も逆順 (`row_end -> row_begin`)
- $L^T$（上三角）を使用

```cpp
for (c = nc; c-- > 0;) {                              // 色を逆順に
    parallel_for(num blocks in color c, [&](bidx) {
        for (i = row_end; i-- > row_begin;) {          // 行を逆順に
            s = 0;
            for (k : off-diagonal entries in row i of Lt_)
                s -= Lt_.values[k] * y[Lt_.col_idx[k]];
            y[i] = s * inv_diag_[i] + x[i];           // D^{-1}を適用 + ソースを加算
        }
    });
}
```

#### 4.3.4 正当性の根拠

ABMC彩色は同色ブロック間に下三角依存関係がないことを保証する。したがって:
- 色cのブロックの前進代入は、色 < c の結果のみに依存
- 同色cのブロックは互いに独立 -- 並列実行可能
- 元の行列からのブロック内依存関係は保存される -- 逐次実行で正しい

### 4.4 適用パスの選択

`ICPreconditioner`は3つの適用パスを提供し、入出力ベクトルの「空間」が異なる。

#### 4.4.1 標準パス: `apply_abmc()`（推奨）

CGが元空間で動作する場合に使用。前処理適用時に順列を行う。

```
x (元空間) -> composite_perm_で順列 -> abmc_x_perm_ (ABMC空間)
  -> forward_substitution_abmc -> work_temp_
  -> backward_substitution_abmc -> abmc_y_perm_ (ABMC空間)
y (元空間) <- composite_perm_で逆順列 <- abmc_y_perm_
```

対角スケーリングが有効な場合、順列と同時にスケーリングを適用:
```cpp
abmc_x_perm_[perm[i]] = x[i] * composite_scaling_[i];
y[i] = abmc_y_perm_[perm[i]] * composite_scaling_[i];
```

`composite_scaling_[i] = scaling_[composite_perm_[i]]` -- 事前計算済み。

#### 4.4.2 ABMC全空間パス: `apply_in_reordered_space()`

CGがABMC空間で動作する場合に使用 (`abmc_reorder_spmv=True`)。入出力が既にABMC空間なので順列が不要。

#### 4.4.3 RCM+ABMC分離パス: `apply_rcm_abmc()`

CGがRCM空間で動作する場合に使用。RCMからABMCへの順列のみ必要:
```
x (RCM空間) -> composite_perm_rcm_で順列 -> ABMC空間
  -> 三角ソルブ
y (RCM空間) <- composite_perm_rcm_で逆順列 <- ABMC空間
```

#### 4.4.4 パス選択基準

| パス | SpMV行列 | 前処理 | 設定 |
|------|-------------|----------------|---------------|
| 標準 (4.4.1) | 元の行列 | 合成順列 + ABMC三角ソルブ | `use_abmc=True` |
| RCM+ABMC (4.4.3) | RCM行列 | RCM-to-ABMC順列 | `use_abmc=True, abmc_use_rcm=True` |
| ABMC全空間 (4.4.2) | ABMC行列 | 順列なし | `abmc_reorder_spmv=True` |

NGSolve統合 (`sparsesolv_ngsolve`モジュール) では、SpMVがNGSolveの`BaseMatrix::Mult`で処理され元の行列空間で動作するため、**標準パス (4.4.1)** を使用する。

### 4.5 RCMとの組み合わせ

#### 4.5.1 合成順列

RCM+ABMCを使用する場合、2段階の順列が必要となる:

```
元の行i -> RCM空間: rcm_ordering_[i]
RCM空間 -> ABMC空間: abmc_ord[rcm_ordering_[i]]
```

毎回2段階の間接参照を行う代わりに、セットアップ時に合成順列を事前計算する:

```cpp
// build_composite_permutations()
composite_perm_[i] = abmc_ord[rcm_ordering_[i]];  // 元 -> ABMC（1ステップ）
composite_perm_rcm_ = abmc_ord;                    // RCM -> ABMC（1ステップ）
```

#### 4.5.2 行列コピーの流れ

RCM+ABMCパスには2つの行列コピーが伴う:

```
元の行列A
  | reorder_matrix_with_perm(A, rcm_ordering_)
RCM行列 rcm_csr_
  | reorder_matrix(rcm_view)  [ABMC順序で再順序付け]
ABMC行列 reordered_csr_
  | extract_lower_triangular()
L_ (IC分解用)
```

ABMCのみのパスではコピーは1回のみ (`A -> reordered_csr_ -> L_`)。

#### 4.5.3 メモリ管理

分解完了後、不要な行列コピーを解放する:
```cpp
if (!config_.abmc_reorder_spmv)
    reordered_csr_.clear();  // SpMVは元の行列を使用
if (!config_.abmc_use_rcm)
    rcm_csr_.clear();
```

### 4.6 レベルスケジューリングとの比較

#### 4.6.1 並列化アプローチの相違

SparseSolvは2つの並列三角ソルブ手法を提供する。ABMCを使用しない場合、**永続並列領域付きレベルスケジューリング**が採用される。

| 特性 | レベルスケジューリング（永続） | ABMC |
|---------------|-------------------------------|------|
| 並列粒度 | 行レベル | ブロックレベル |
| 同期回数 | レベル数（数百） | 色数（4〜数十） |
| 同期メカニズム | SpinBarrier（単一parallel_for内） | parallel_for呼び出し（色ごと） |
| FEMでの並列性 | 低い（深い依存チェーン） | 高い |
| セットアップコスト | O(nnz) -- レベル計算のみ | O(nnz) -- BFS + 彩色 + 行列再順序付け |
| 追加メモリ | レベル配列のみ | ワークベクトル2*n + 行列コピー |
| 適用時のメモリ確保 | なし（事前確保済み） | なし（事前確保済み） |

**追加メモリの内訳（ABMCパス）**:
- `abmc_x_perm_` (n), `abmc_y_perm_` (n) -- ABMC固有
- `work_temp_` (n), `work_temp2_` (n) -- 両パスで共有
- `composite_perm_` (n), `composite_scaling_` (n) -- 順列/スケーリング
- `reordered_csr_` -- 再順序付け行列（`abmc_reorder_spmv=false`の場合、分解後に解放）

#### 4.6.2 永続並列領域

素朴なレベルスケジューリング実装では各レベルで`parallel_for`を呼び出す。レベル数が数百に達すると、スレッドプール起動のオーバーヘッドが問題となる。

`forward_substitution_persistent()` はすべてのレベルを単一の`parallel_for(nthreads, ...)`内で処理し、レベル間同期にSpinBarrierを使用する:

```cpp
parallel_for(nthreads, [&](index_t thread_id) {
    for (int lev = 0; lev < num_levels; ++lev) {
        // 各スレッドが割り当てられた行範囲を処理
        const index_t my_start = level_size * thread_id / nthreads;
        const index_t my_end = level_size * (thread_id + 1) / nthreads;
        for (idx = my_start; idx < my_end; ++idx) { ... }
        barrier.wait();  // レベル間同期
    }
});
```

このパスはスレッド数が1より大きい場合に自動選択される (`get_num_threads() > 1`)。

#### 4.6.3 SpinBarrierの実装

`SpinBarrier` (`parallel.hpp`) はセンス反転アルゴリズムを使用する:

```cpp
class SpinBarrier {
    alignas(64) std::atomic<int> count_;   // 到着済みスレッド数
    alignas(64) std::atomic<int> sense_;   // 現在の世代番号
    int num_threads_;
public:
    void wait() {
        int my_sense = sense_.load(acquire);
        if (count_.fetch_add(1, acq_rel) == num_threads_ - 1) {
            // 最後のスレッド: カウンタをリセット + センスを反転
            count_.store(0, relaxed);
            sense_.fetch_add(1, release);
        } else {
            // スピン待ち（4096回反復後にyield）
            while (sense_.load(acquire) == my_sense) { ... }
        }
    }
};
```

**キャッシュライン分離**: `count_`と`sense_`は`alignas(64)`で別々のキャッシュラインに配置し、偽共有を防止する。

**スピン-to-yieldストラテジー**: 4096回のスピン反復後に`std::this_thread::yield()`を呼び出す。短い同期にはスピンし、長い同期にはCPUを譲る。

### 4.7 性能特性

#### 4.7.1 ベンチマーク条件

- 8スレッド (NGSolve TaskManager)
- 3D HCurl curl-curl問題 (order=2, `nograds=True`)
- 単位立方体メッシュ、シフト付きICCG (`shift=1.05, auto_shift=True, diagonal_scaling=True`)
- ベースライン: 永続並列領域付きレベルスケジューリング
- 3回実行の最小値

#### 4.7.2 スケーリング結果

| 自由度数 | レベルスケジューリング | ABMC（最良） | 高速化率 | 反復回数 |
|------|-------------|-------------|---------|------------|
| 11K | 0.044s | 0.051s | **0.86x** | 31 |
| 27K | 0.158s | 0.154s | **1.02x** | 39 |
| 82K | 0.763s | 0.650s | **1.17x** | 59 |
| 186K | 2.505s | 1.962s | **1.28x** | 79 |

#### 4.7.3 NB02トロイダルコイルでの比較

[02_performance_comparison.ipynb](02_performance_comparison.ipynb)のトロイダルコイル問題 (148K自由度) において:

| ソルバー | 反復回数 | 実行時間 | 高速化率 |
|--------|-----------|-----------|---------|
| ICCG (レベルスケジューリング) | 463 | 10.0s | 1.0x |
| ICCG + ABMC (8色) | 415 | 6.0s | **1.7x** |

単位立方体 (186K自由度) での高速化率は1.28xであるのに対し、トロイダルコイル (148K自由度) では1.7xとなる。この差異の要因は:
- **反復回数**: 463 vs 79。反復回数が多いほど三角ソルブに費やされる時間割合が増え、ABMCの効果が顕著になる
- **行列構造**: トロイダルコイルはバンド幅が広く、レベルスケジューリングのレベル数が増加し、同期オーバーヘッドが増大するため、ABMCの優位性が大きくなる

#### 4.7.4 損益分岐点

ABMCには以下のオーバーヘッドがある:
1. **セットアップ**: BFS集約、ブロックグラフ構築、彩色、行列再順序付け
2. **適用ごと**: ベクトル順列（scatter/gatherの2回の`parallel_for`呼び出し）

これらのオーバーヘッドが三角ソルブの高速化を上回る領域では、ABMCは逆効果となる。

**損益分岐点: 約25K〜30K自由度**（8スレッド、単位立方体、HCurl次数2）

- 自由度 < 25K: レベルスケジューリングが高速
- 自由度 > 30K: ABMCが有利。問題サイズに比例して効果が増大

#### 4.7.5 有効性に影響する要因

| 要因 | ABMCに有利 | ABMCに不利 |
|--------|-------------------|---------------------|
| 問題サイズ | 大規模 (>30K自由度) | 小規模 (<25K自由度) |
| 反復回数 | 多い (>100) | 少ない (<30) |
| 行列バンド幅 | 広い（深い依存チェーン） | 狭い |
| ABMC色数 | 少ない (4〜8) | 多い (>20) |
| スレッド数 | 多い (>=4) | 少ない (1〜2) |

### 4.8 設定パラメータ

`SolverConfig`で制御されるABMC関連パラメータ:

| パラメータ | 型 | デフォルト | 説明 |
|-----------|------|---------|------|
| `use_abmc` | bool | false | ABMCを有効化 |
| `abmc_block_size` | int | 4 | BFS集約のブロックサイズ（行/ブロック） |
| `abmc_num_colors` | int | 4 | 目標色数（下限値、自動拡張あり） |
| `abmc_reorder_spmv` | bool | false | ABMC空間でCGを実行（SpMVも再順序付け） |
| `abmc_use_rcm` | bool | false | ABMC前にRCMバンド幅削減を適用 |

**推奨設定**: `use_abmc=True`のみ（他はデフォルト）で十分。

- `abmc_block_size`: 4〜16の範囲で性能差は数%以内
- `abmc_num_colors`: ほぼ常に自動拡張されるため、指定値の影響は小さい
- `abmc_reorder_spmv`: FEMメッシュ順序のキャッシュ局所性を保持するため通常false
- `abmc_use_rcm`: 追加の行列コピーとセットアップコストのため通常false

### 4.9 ABMCScheduleデータ構造

`ABMCSchedule`構造体のメンバー (`abmc_ordering.hpp`):

```cpp
struct ABMCSchedule {
    // 色 -> ブロック（CSR形式）
    std::vector<index_t> color_offsets;    // サイズ: num_colors + 1
    std::vector<index_t> color_blocks;     // 色ごとのブロックIDリスト

    // ブロック -> 行（CSR形式、連続番号付けにより実質オフセットのみ）
    std::vector<index_t> block_offsets;    // サイズ: num_blocks + 1

    // 行の順列
    std::vector<index_t> ordering;         // ordering[old_row] = new_row
    std::vector<index_t> reverse_ordering; // reverse_ordering[new_row] = old_row
};
```

**三角ソルブでの走査パターン**:
```
for c = 0 to num_colors - 1:
    blk_begin = color_offsets[c]
    blk_end = color_offsets[c + 1]
    parallel_for(blk_begin .. blk_end):
        blk = color_blocks[bidx]
        for i = block_offsets[blk] to block_offsets[blk+1] - 1:
            // ABMC空間の行iを処理
```

**メモリ効率**: `block_offsets`が連続した行範囲を定義するため、行-to-ブロックマッピング配列は不要。ブロック $b$ の行は `[block_offsets[b], block_offsets[b+1])` で完全に決定される。

### 4.10 CGカーネル融合 (v2.3.0)

#### 4.10.1 動機

CG反復はメモリバンド幅律速であり（算術強度 ~0.33 FLOP/byte）、カーネル間でのベクトル再読み込みがボトルネックとなる。融合前のCG反復は7カーネル起動で、合計メモリトラフィックはnnz + 12nリードと4nライト（IC適用を除く）。

#### 4.10.2 融合カーネル

**フェーズ1: SpMV + 内積融合** -- `cg_solver.hpp`

```cpp
Scalar pAp = parallel_reduce_sum<Scalar>(n, [&](index_t i) -> Scalar {
    Scalar s = Scalar(0);
    for (index_t k = A_rowptr[i]; k < A_rowptr[i + 1]; ++k)
        s += A_vals[k] * p[A_colidx[k]];
    Ap[i] = s;
    return p[i] * s;  // dot(p, Ap)の部分和
});
```

SpMV (`Ap = A*p`) と内積 (`pAp = dot(p, Ap)`) を1回のパスで計算。p[]とAp[]の再読み込み (2nリード) を排除する。

**フェーズ2: AXPY + ノルム融合**

```cpp
double norm_r_sq = parallel_reduce_sum<double>(n, [&](index_t i) -> double {
    x[i] += alpha * p[i];
    r[i] -= alpha * Ap[i];
    return std::norm(r[i]);  // |r[i]|^2
});
```

ベクトル更新 (`x += alpha*p`, `r -= alpha*Ap`) と残差ノルム計算を1回のパスで実行。r[]の再読み込み (nリード) を排除する。

**フェーズ3: 前処理適用 + dot(r, z)融合** -- `preconditioner.hpp`, `ic_preconditioner.hpp`

```cpp
// Preconditioner::apply_fused_dot() -- 基底クラスの仮想メソッド
Scalar apply_fused_dot(const Scalar* r_for_dot, const Scalar* x,
                       Scalar* y, index_t size, bool conjugate) const {
    apply(x, y, size);
    return parallel_reduce_sum<Scalar>(size, [&](index_t i) -> Scalar {
        return conjugate ? std::conj(r_for_dot[i]) * y[i]
                         : r_for_dot[i] * y[i];
    });
}
```

各IC前処理パスは`apply_fused_dot()`をオーバーライドし、後退代入の出力フェーズで`dot(r, z)`を同時に計算する。ABMC空間パス（パス1）では、永続並列領域内でスレッドローカルの部分和を蓄積し、最後にシリアルリダクションする。

**効果**: 反復あたり7カーネルが4に削減され、メモリトラフィックが約25%削減。

#### 4.10.3 ABMC並列auto_shift

v2.3.0以前、auto_shiftが有効な場合はABMC並列IC分解が利用できず、逐次パスにフォールバックしていた。v2.3.0ではアトミックフラグ (`std::atomic<bool>`) によるリスタートメカニズムを導入:

1. ABMC色並列IC分解を実行
2. いずれかのスレッドが対角の破綻を検出した場合、アトミックフラグをセット
3. 色ループの終了時にフラグを確認 -- シフトを増加（指数バックオフ）-- 元の値からやり直し

シフト増加は指数バックオフ (`increment *= 2`) を用い、少ない回数のリスタートで適切なシフト値に到達する。
ABMC順序付けはパターンベースのため、リスタート時に再計算は不要。

**効果**: 比留間渦電流問題 (HCurl p=1) において、並列スケーリングが1.5xから2.6x（8コア）に改善。

#### 4.10.4 ABMC空間CG (`abmc_reorder_spmv=True`)

CGをABMC再順序付け空間で完全に実行するモード。SpMV、IC適用、すべてのベクトル演算がABMC空間で行われる。

利点:
- IC適用時の入出力ベクトル順列（scatter/gather）を排除
- SpMVのキャッシュ局所性がABMCブロック構造と一致
- パラメータ`abmc_reorder_spmv=True`で有効化

**効果**: 比留間問題で2.7xから2.89x（8コア）。

#### 4.10.5 永続並列領域 + dot(r,z)融合 (`apply_fused_dot`)

`apply_in_reordered_space()`（ABMC空間CGパス）は元々各色ごとに`forward_substitution_abmc()` / `backward_substitution_abmc()`を呼び出し、合計2*nc回の`parallel_for`ディスパッチが発生していた。これを永続並列領域（1回のディスパッチ + 2*nc回のバリア）に変換し、`apply_abmc()`と同様に、さらに後退代入の最終出力で`dot(r, z)`を計算する。

`apply_in_reordered_space_fused_dot()`の構造:

```cpp
parallel_for(nthreads, [&](index_t tid) {
    // 前進代入（色は逐次、ブロックは並列）
    for (index_t c = 0; c < nc; ++c) {
        // ... 色c内の自分のブロック ...
        barrier.wait();
    }
    // 後退代入 + 内積蓄積
    Scalar local_dot = 0;
    for (index_t c = nc; c-- > 0;) {
        // ... 自分のブロックの後退ソルブ ...
        // 最終書き込み時: local_dot += r[i] * y[i];
        barrier.wait();
    }
    partial_dots[tid] = local_dot;
});
// 部分和のシリアルリダクション
```

全3パスでサポート:
- **パス1** (`abmc_reorder_spmv=True`): `apply_in_reordered_space_fused_dot()`
- **パス2** (RCM+ABMC): `apply_rcm_abmc_fused_dot()`
- **パス3** (標準ABMC): `apply_abmc_fused_dot()` -- フェーズ4（出力順列）でのdot融合

**効果**: 比留間問題で2.89xから3.14x（8コア、理論最大4xの79%）。

### 4.11 既知の制限事項

#### 4.11.1 色数の予測困難性

`target_colors`は貪欲彩色の下限値であり、ブロック間依存が密な場合は自動拡張される。複雑なメッシュ（例: ヘリカルコイル）では、実際の色数が指定値を大幅に超える場合がある。

**潜在的改善策**: 実際の色数をログ出力し、ユーザーが確認できるようにする。

#### 4.11.2 ブロックサイズの自動チューニングなし

現在、`block_size`（デフォルト: 4）はユーザー指定の固定値である。行列構造に基づく自動チューニングは実装されていない。

ベンチマーク結果では、ブロックサイズ4〜16の差は小さい（数%以内）。

#### 4.11.3 RCM事前再順序付けの効果が限定的

RCMはバンド幅削減によりABMC色数を減少させる可能性があるが、追加の行列コピーコストとセットアップ時間を伴う。現在のベンチマークでは、`abmc_use_rcm=False`と`True`の差は誤差範囲内であるため、デフォルトは無効。

### 4.12 改善履歴

本節では実装レビューに基づいて行われた改善を記録する。

#### 4.12.1 `block_rows`配列の除去

**問題**: `build_row_ordering()`で構築された`block_rows[ridx]`配列は常に`ridx`に等しかった（新しいインデックスは連続的に割り当てられるため、恒等写像であった）。`i = block_rows[ridx]`の代わりに単に`i = ridx`と書ける。

**解決**: `ABMCSchedule`から`block_rows`メンバーを除去。三角ソルブのループを`for (i = row_begin; i < row_end; ++i)`に簡素化。O(n)のメモリ削減。

#### 4.12.2 ワークベクトルの事前確保と共有

**問題**: レベルスケジューリングパスの`apply_level_schedule()`は毎回`std::vector<Scalar> temp(size)`をヒープ確保していた。ABMCパスは`setup()`時に事前確保していたが、異なる変数名 (`abmc_temp_`) を使用。

**解決**: `work_temp_`（および`work_temp2_`）を両パスの共有メンバー変数として統一。`setup()`時に1回確保し、`apply()`時のヒープ確保を完全に排除。

#### 4.12.3 `build_row_ordering`の未使用パラメータ除去

**問題**: `color_graph()`に未使用の出力パラメータがあり、`build_row_ordering()`にも未使用のパラメータがあった。

**解決**: `color_graph()`から不要な出力パラメータを除去。`build_row_ordering()`のシグネチャを簡素化。

### 4.13 謝辞

ABMC順序付けの実装は、鶴谷祐紀氏（福岡大学）が[JP-MARs/SparseSolv](https://github.com/JP-MARs/SparseSolv)で提供したコードに基づく。このコードは岩下武史教授（京都大学）のグループによるABMC法 [5] およびSGS-MRTR前処理の研究を実装したものである。本リポジトリはヘッダオンリー再構成、NGSolve統合、auto_shift IC、Compact AMS前処理、COCRソルバーなどの独自拡張を加えている。

---

## 5. CG法 (共役勾配法)

### 5.1 概要

対称正定値 (SPD) 行列に対する**前処理付き共役勾配 (PCG)** 法 [Hestenes, Stiefel 1952]。

**実装ファイル**: `solvers/cg_solver.hpp`

### 5.2 アルゴリズム

```
r_0 = b - A*x_0
z_0 = M^{-1} * r_0
p_0 = z_0

for k = 0, 1, 2, ...:
    alpha_k = (r_k, z_k) / (p_k, A*p_k)
    x_{k+1} = x_k + alpha_k * p_k
    r_{k+1} = r_k - alpha_k * A*p_k
    z_{k+1} = M^{-1} * r_{k+1}
    beta_k = (r_{k+1}, z_{k+1}) / (r_k, z_k)
    p_{k+1} = z_{k+1} + beta_k * p_k
```

### 5.3 複素対称行列のサポート

FEM行列は複素対称 (A^T = A) であり、**エルミートではない** (A^H != A)。
例: 渦電流方程式 curl-curl + i*sigma*mass。

この場合、内積は**非共役**でなければならない:

```
(a, b) = Sigma a_i * b_i       （非共役: 複素対称）
(a, b) = Sigma conj(a_i) * b_i  （共役: エルミート）
```

実装での切り替え:

```cpp
// iterative_solver.hpp: dot_product()
if (config_.conjugate)
    sum += std::conj(a[i]) * b[i];  // エルミート
else
    sum += a[i] * b[i];  // 複素対称（デフォルト）
```

これはNGSolveの`CGSolver(conjugate=False)`に対応する。

### 5.4 収束判定

相対残差ノルムにより収束を判定する:

```
||r_k|| / ||r_0|| < tol
```

追加機能:
- **最良結果の追跡**: 最良の解を保持し、収束しない場合にそれを返す
- **発散検出**: 残差が最良値の一定倍を超えた場合に早期終了
- **残差履歴**: 各反復での残差を記録

---

## 6. DenseMatrix LU逆行列

### 6.1 概要

小規模密行列のLU分解と逆行列計算。COO/CSR行列ビルダーの内部で使用される。
部分ピボット付きPA = LU分解を実行する。

**実装ファイル**: `core/dense_matrix.hpp`

### 6.2 置換行列の構築

部分ピボット付きLU分解は、行の交換を記録する置換ベクトル`piv[]`を生成する。
置換行列Pは以下のように構築される:

```
P[k, piv[k]] = 1   （正しい: PA = LU）
```

**注意**: `P[piv[k], k] = 1`はP^T（転置）であり、誤りとなる。
この区別は対角優位行列 (H1) ではほとんど影響しないが、
非自明なピボットを持つ行列 (HCurl) では決定的に重要となる。

```cpp
// dense_matrix.hpp: invert()
DenseMatrix inv(n, n);
for (index_t k = 0; k < n; ++k) {
    inv(k, piv[k]) = Scalar(1);  // P[k, piv[k]] = 1
}
```

### 6.3 逆行列計算の3段階

```
1. LU分解: PA = LU（部分ピボット）
2. 前進代入: L * Y = P（Pは置換行列）
3. 後退代入: U * X = Y -> X = A^{-1}
```

---

## 7. コード起源と貢献者

ngsolve-sparsesolv は複数の研究者のコードを統合・拡張している。

### 7.1 起源マッピング

| コンポーネント | 起源 | 原著者 | 元のファイル (JP-MARs/SparseSolv) |
|---|---|---|---|
| CG反復ソルバー | 佐藤 | 佐藤 (tsato) | `MatSolvers_ICCG.cpp` |
| IC分解および逐次前進/後退代入 | 佐藤 | 佐藤 (tsato) | `MatSolvers_ICCG.cpp` |
| ABMC順序付け | 比留間 | 比留間慎吾 | `MatSolvers_ABMCICCG.cpp` |
| ABMC並列IC分解および並列前進/後退代入 | 比留間 | 比留間慎吾 | `MatSolvers_ABMCICCG.cpp` |
| IC-MRTR反復 | 円谷 -> 佐藤 | 円谷友紀（理論）、佐藤（実装） | `MatSolvers_ICMRTR.cpp` |
| SGS-MRTR反復（Eisenstatテクニック） | 円谷 -> 佐藤 | 円谷友紀（理論）、佐藤（実装） | `MatSolvers_SGSMRTR.cpp` |

### 7.2 元のコード

- **円谷のコード**: IC-MRTRおよびSGS-MRTR反復公式のリファレンス実装（C言語、1始まりインデックス）。
  シフトパラメータ付きIC分解のauto_shiftループ（gamma = 1.05、+0.05刻み）はこのコードに由来する

- **JP-MARs/SparseSolv**: `https://github.com/JP-MARs/SparseSolv`
  - 佐藤 (tsato) による疎行列ソルバーフレームワーク（C++、Eigenベース）
  - 比留間によりABMC並列化を追加
  - ngsolve-sparsesolvのfork元

### 7.3 ngsolve-sparsesolvでの拡張

fork後に追加された機能:

| 機能 | 説明 |
|---|---|
| IC自動シフト | 半正定値行列に対する自動シフト調整 |
| 対角スケーリング | DAD変換による条件数改善 |
| 零対角処理 | curl-curl行列の零対角自由度のサポート |
| Localized IC (Block IC) | Fukuhara 2009に基づくパーティションレベルIC |
| 永続並列領域 | SpinBarrier同期によるレベルスケジューリングの高速化 |
| RCM順序付け | ABMC前のReverse Cuthill-McKeeバンド幅削減 |
| 複素対称サポート | COCG / 複素対称CG |
| NGSolve統合 | pybind11モジュール、BaseMatrix互換性 |

---

## 8. Compact AMS + COCR（複素渦電流問題向け）

### 8.1 概要

Compact AMSは、Hiptmair-Xu (2007) による**Auxiliary Space Preconditioning**法の
ヘッダオンリーC++実装である。外部ライブラリは不要で、NGSolve TaskManagerを用いて並列化する。

**特徴**:
- HCurl有限要素によるcurl-curl + mass問題に特化（実数・複素数の両方をサポート）
- 離散勾配行列Gと頂点座標のみが必要（要素行列は不要）
- 補助空間ソルバーとしてCompactAMG（古典的AMG）を使用
- `Update()`は非線形ソルバー（Newton反復）をサポート: 幾何情報を保持し、行列依存部分のみを再構築
- 実数の磁気静解析: `CompactAMSPreconditioner` + CG
- 複素渦電流: `ComplexCompactAMSPreconditioner` + COCR

**実装ファイル**: `preconditioners/compact_amg.hpp`, `compact_ams.hpp`, `complex_compact_ams.hpp`

詳細は[compact_ams_cocr.md](compact_ams_cocr.md)を参照。

### 8.2 AMS理論

HCurl空間のHelmholtz分解:

```
H(curl) = grad(H^1) + (H(curl) intersection ker(div))
```

curl-curl行列Kの核はgrad(H^1)勾配空間に含まれる。
AMSはこの空間分解を利用して、2つの補助空間で前処理を行う:

1. **G補正**（勾配補正）: grad(H^1)成分を`G * A_G^{-1} * G^T`で処理
2. **Pi補正**（Nedelec補間）: curl成分を`Pi * A_Pi^{-1} * Pi^T`で処理

ここで:
- G: 離散勾配行列 (HCurl -> H1)、`fes.CreateGradient()`で取得
- Pi: Nedelec補間行列、頂点座標 (x, y, z) から自動構築
- A_G, A_Pi: 補助空間上の粗格子行列（CompactAMGによる近似逆行列）

### 8.3 複素渦電流への適用

複素対称系 A = K + jw*sigma*M に対しては、Re/Im分離を用いる。
AMS前処理は実数SPD補助行列 `A_real = K + eps*M + |w|*sigma*M_cond` から構築され、
同じAMSインスタンスが複素ベクトルの実部と虚部に適用される。

**Re/Im融合処理**: SpMVはメモリバンド幅律速であるため、行列データを1回ロードして
Re/Imを同時に処理することでバンド幅コストを半減させる。この最適化は
細レベルのスムーザ、残差計算、制限/延長、AMG Vサイクルの全レベル (DualMult) に適用される。

### 8.4 AMG粗格子構築の計算量

`BuildClassicalInterp()`（古典的補間行列構築）は各行について強結合列の
非零パターンを走査するため、最悪計算量は**O(nnz^2/n)**となる。
密な結合パターンを持つ行列（3D高次要素）では、粗視化のセットアップが支配的になりうる。

現時点ではこのコストは許容範囲内であるが、大規模問題ではPMIS粗視化の品質
（粗格子サイズ）とのトレードオフに注意が必要である。

### 8.5 推奨選択

| 条件 | 推奨 | 理由 |
|-----------|---------------|--------|
| H1低次 (p=1-2) | **ABMC ICCG** | セットアップが最速、メモリ最小 |
| H1/HCurl高次 (p>=3) | **BDDC**（NGSolve組込み） | h非依存の反復回数、高次で卓越 |
| HCurl磁気静解析（実数、p=1） | **Compact AMS+CG** | AMS補助空間がcurl-curl零空間を処理 |
| HCurl非線形（Newton反復） | **Compact AMS+CG** | `Update()`で前処理を再構築 |
| HCurl渦電流（大規模、p=1） | **Compact AMS+COCR** | Re/Im融合 + AMS補助空間 |
| HCurl渦電流（中規模） | **BDDC+CG**（NGSolve組込み） | メモリ効率が良い、対称前処理 |
| メモリ制約が厳しい場合 | **ABMC ICCG** | CG/COCR: 約5〜6本のワークベクトル |

---

## 参考文献

1. C. R. Dohrmann,
   "A Preconditioner for Substructuring Based on Constrained Energy Minimization",
   *SIAM J. Sci. Comput.*, Vol. 25, No. 1, pp. 246-258, 2003.
   [DOI: 10.1137/S1064827502412887](https://doi.org/10.1137/S1064827502412887)

2. J. Mandel, C. R. Dohrmann, R. Tezaur,
   "An Algebraic Theory for Primal and Dual Substructuring Methods
   by Constraints",
   *Appl. Numer. Math.*, Vol. 54, No. 2, pp. 167-193, 2005.
   [DOI: 10.1016/j.apnum.2004.09.022](https://doi.org/10.1016/j.apnum.2004.09.022)

3. J. A. Meijerink, H. A. van der Vorst,
   "An Iterative Solution Method for Linear Systems of Which the
   Coefficient Matrix is a Symmetric M-Matrix",
   *Math. Comp.*, Vol. 31, No. 137, pp. 148-162, 1977.
   [DOI: 10.1090/S0025-5718-1977-0438681-4](https://doi.org/10.1090/S0025-5718-1977-0438681-4)

4. M. R. Hestenes, E. Stiefel,
   "Methods of Conjugate Gradients for Solving Linear Systems",
   *J. Research of the National Bureau of Standards*,
   Vol. 49, No. 6, pp. 409-436, 1952.
   [DOI: 10.6028/jres.049.044](https://doi.org/10.6028/jres.049.044)

5. T. Iwashita, H. Nakashima, Y. Takahashi,
   "Algebraic Block Multi-Color Ordering Method for Parallel
   Multi-Threaded Sparse Triangular Solver in ICCG Method",
   *Proc. IEEE IPDPS*, 2012.
   [DOI: 10.1109/IPDPS.2012.51](https://doi.org/10.1109/IPDPS.2012.51)

6. T. Tsuburaya, Y. Mifune, T. Iwashita, E. Takahashi,
   "Numerical Experiments on Preconditioned Iterative Methods Based on the MRTR Method",
   *IEEJ Technical Meeting on Static Apparatus*, SA-12-64, 2012.

7. T. Tsuburaya, Y. Okamoto, K. Fujiwara, S. Sato,
   "Improvement of the Preconditioned MRTR Method With Eisenstat's Technique
   in Real Symmetric Sparse Matrices",
   *IEEE Trans. Magnetics*, Vol. 49, No. 5, pp. 1641-1644, 2013.
   [DOI: 10.1109/TMAG.2013.2240283](https://doi.org/10.1109/TMAG.2013.2240283)

8. R. Hiptmair, J. Xu,
   "Nodal Auxiliary Space Preconditioning in H(curl) and H(div) Spaces",
   *SIAM J. Numer. Anal.*, Vol. 45, No. 6, pp. 2483-2509, 2007.
   [DOI: 10.1137/060660588](https://doi.org/10.1137/060660588)

9. T. V. Kolev, P. S. Vassilevski,
   "Parallel Auxiliary Space AMG for H(curl) Problems",
   *J. Comput. Math.*, Vol. 27, No. 5, pp. 604-623, 2009.

10. T. Sogabe, S.-L. Zhang,
    "A COCR method for solving complex symmetric linear systems",
    *J. Comput. Appl. Math.*, Vol. 199, No. 2, pp. 297-303, 2007.

11. E. Cuthill, J. McKee,
    "Reducing the Bandwidth of Sparse Symmetric Matrices",
    *Proc. 24th Nat. Conf. ACM*, pp. 157-172, 1969.
    [DOI: 10.1145/800195.805928](https://doi.org/10.1145/800195.805928)
