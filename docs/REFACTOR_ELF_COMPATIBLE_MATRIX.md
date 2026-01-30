# Radia行列定式化のELF互換化リファクタリング

## 1. 概要

### 1.1 目的

RadiaのMSC (Magnetic Surface Charge) 行列定式化をELF/MAGICと互換にする。
これにより:
- ELFの検証済み結果との直接比較が可能
- TrfMlt（対称変換）の正確な実装が可能
- 既存ELFユーザーの移行が容易

### 1.2 変更範囲

| コンポーネント | 変更内容 |
|---------------|---------|
| 相互作用行列構築 | 符号規約の変更 |
| 対角項の計算 | 1/χ → -1/χ |
| ソルバー | RHS符号の調整 |
| テスト | ELF検証データとの比較 |

## 2. 現状と目標

### 2.1 現在のRadia定式化

```
行列方程式: A · M = H_ext
行列定義:   A = -N + diag(1/χ)

対角要素:   A_ii = 1/χ - N_ii  (正の値)
非対角要素: A_ij = -N_ij
```

**コード箇所**: `src/core/rad_hacapk.h` line 84
```cpp
// The H-matrix approximates A = -N + diag(1/chi)
```

### 2.2 目標: ELF互換定式化

```
行列方程式: A · M = H_ext
行列定義:   A = -N - diag(1/χ)

対角要素:   A_ii = -1/χ - N_ii  (負の値)
非対角要素: A_ij = -N_ij
```

### 2.3 数学的関係

変換は単純:
```
A_ELF = A_Radia - 2·diag(1/χ)
```

または:
```
A_ELF[i,i] = A_Radia[i,i] - 2/χ_i
A_ELF[i,j] = A_Radia[i,j]  (i ≠ j)
```

## 3. 影響を受けるファイル

### 3.1 コアファイル

| ファイル | 変更内容 |
|---------|---------|
| `rad_interaction.cpp` | 行列構築時の対角項符号 |
| `rad_interaction.h` | コメント更新 |
| `rad_hacapk.cpp` | UpdateDiagonal()の符号 |
| `rad_hacapk.h` | コメント更新 |
| `rad_relaxation_methods.cpp` | ソルバーの調整 |
| `rad_mmm_hacapk.cpp` | MMM用H-matrix |

### 3.2 テストファイル

| ファイル | 変更内容 |
|---------|---------|
| `tests/test_msc_matrix.py` | 期待値の符号変更 |
| `examples/electromagnet/mu=1000/verify_matrix.py` | ELF比較テスト |

## 4. 実装計画

### Phase 1: 対角項の符号変更

#### 4.1.1 rad_interaction.cpp

**現在のコード** (SetupInteractMatrix内):
```cpp
// 対角ブロックの計算
for(int i = 0; i < AmOfMainElem; i++)
{
    // Self-interaction term
    // 現在: N_ii が計算される
    // 変更不要（Nの符号はそのまま）
}
```

**変更箇所** (UpdateDiagonalWithChi):
```cpp
// 現在:
FlatInteractMatrix[diagIdx] = N_ii + 1.0 / chi;

// 変更後:
FlatInteractMatrix[diagIdx] = N_ii - 1.0 / chi;  // ELF互換
```

#### 4.1.2 rad_hacapk.cpp

**UpdateDiagonal関数の変更**:
```cpp
void RadHACApKManager::UpdateDiagonal(const std::vector<double>& inv_chi) {
    // 現在:
    // diag[i] = N_ii + inv_chi[i]

    // 変更後:
    // diag[i] = N_ii - inv_chi[i]  // ELF互換: -1/χ
}
```

### Phase 2: ソルバーの調整

#### 4.2.1 rad_relaxation_methods.cpp

RHSの符号は変更不要（A·M = H_ext のまま）。
ただし、残差計算と収束判定の確認が必要。

**確認項目**:
- `r = H_ext - A·M` の計算
- Jacobi前処理の対角抽出
- 収束判定 `|r| / |H_ext|`

### Phase 3: テスト更新

#### 4.3.1 ELF検証テスト

```python
# verify_matrix.py の更新

def compare_matrices(radia_matrix, elf_matrix):
    # 現在: 符号反転と順列を試行
    # 変更後: 順列のみで一致するはず

    # 対角要素の符号確認
    assert np.all(np.diag(radia_matrix) < 0), "対角は負であるべき"

    # ELFとの直接比較（順列適用後）
    P = find_permutation(radia_matrix, elf_matrix)
    radia_permuted = P @ radia_matrix @ P.T

    np.testing.assert_allclose(radia_permuted, elf_matrix, rtol=1e-4)
```

## 5. 面順序の統一

### 5.1 調査結果 (2025-01-30)

**ELFとRadiaの面定義は同一である。**

両者とも以下の面順序を使用（1-indexed頂点）:
```
面0/1: 頂点 1-4-3-2 (z-)
面1/2: 頂点 5-6-7-8 (z+)
面2/3: 頂点 1-2-6-5 (y-)
面3/4: 頂点 3-4-8-7 (y+)
面4/5: 頂点 1-5-8-4 (x-)
面5/6: 頂点 2-3-7-6 (x+)
```

### 5.2 コサイン類似度

Phase 1の符号変更後の比較結果:
- **コサイン類似度: 0.918** (以前は0.35)
- 最適な並び替え: `(2, 4, 3, 5, 0, 1)`

### 5.3 残差の原因

行列が完全に一致しない理由:
1. **Nastranファイルの頂点順序**: ELFのNastranファイルで定義された頂点順序が
   標準CHEXAと異なる場合がある（メッシュ生成ツール依存）
2. **評価点の違い**: 面上の積分評価点の選択が異なる可能性
3. **数値精度**: 立体角計算の実装の違い

### 5.4 結論

面順序の変更は**不要**。ELFとRadiaは同じ面定義を使用している。
並び替え `(2, 4, 3, 5, 0, 1)` はNastranファイル固有の頂点順序に起因する。

**推奨対応**:
- 面順序の変更は行わない
- 頂点順序が異なるメッシュに対しては、読み込み時に正規化を行う
- Phase 3のELF検証テストで物理結果（磁場）を比較する

## 6. 検証計画

### 6.1 単体テスト

| テスト | 内容 | 期待結果 |
|--------|------|---------|
| 単一立方体 | 10cm立方体, μ_r=1000 | 対角要素が負 |
| 対称性 | 行列 A = A^T | 対称行列 |
| ELF比較 | 同一メッシュ | 順列後一致 |

### 6.2 統合テスト

| テスト | 内容 | 期待結果 |
|--------|------|---------|
| 線形解析 | μ_r=1000電磁石 | ELFと同じ磁場 |
| 非線形解析 | BHカーブ | 収束、ELFと一致 |
| TrfMlt | 1/4モデル対称 | 全体モデルと一致 |

## 7. 互換性への影響

### 7.1 破壊的変更

- `GetInteractMatrix()` の戻り値が符号変更
- 行列要素を直接使用するユーザーコードは修正が必要

### 7.2 非破壊的項目

- `Solve()` API は変更なし
- `Fld()` 計算結果は変更なし（同じ物理結果）
- Pythonユーザー向けAPIは変更なし

### 7.3 バージョニング

- マイナーバージョンアップ: v1.x.0 → v1.(x+1).0
- 変更履歴に明記

## 8. 実装順序

```
Week 1: Phase 1 - 対角項符号変更
  ├── rad_interaction.cpp の修正
  ├── rad_hacapk.cpp の修正
  └── 単体テスト作成

Week 2: Phase 2 - 面順序統一
  ├── 面順序マッピングの調査（ELF側）
  ├── rad_polyhedron.cpp の修正
  └── 順列テスト

Week 3: Phase 3 - 統合テスト
  ├── ELF検証テストの実行
  ├── TrfMlt対称テスト
  └── ドキュメント更新

Week 4: リリース準備
  ├── CHANGELOG更新
  ├── バージョンアップ
  └── PyPIリリース
```

## 9. ロールバック計画

変更が問題を引き起こした場合:

1. `#define RADIA_ELF_COMPAT` フラグで切り替え可能に
2. デフォルトは新方式（ELF互換）
3. レガシーモードも維持

```cpp
#ifdef RADIA_ELF_COMPAT
    // ELF互換: 負の対角
    diag = N_ii - inv_chi;
#else
    // レガシー: 正の対角
    diag = N_ii + inv_chi;
#endif
```

## 10. まとめ

### 変更の本質

```
現在:  A = -N + diag(+1/χ)  → 対角が正
変更後: A = -N + diag(-1/χ)  → 対角が負（ELF互換）
```

### 期待される効果

1. **ELF検証**: 同一メッシュで行列が一致
2. **TrfMlt修正**: 対称変換が正しく動作
3. **保守性向上**: 単一の標準定式化

---

*作成日: 2025-01-30*
*ステータス: 計画中*
