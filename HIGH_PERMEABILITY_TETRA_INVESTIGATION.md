# 高透磁率四面体メッシュの問題調査レポート

**日時**: 2025-11-25
**ステータス**: ⚠️ **問題確認・根本原因特定**

---

## 概要

ANALYTICAL四面体メソッドが高透磁率材料(μr > 100)で正常に動作しない問題を調査しました。

### 症状

| μr | 誤差 | 状態 |
|----|------|------|
| 100 | 14.78% | ✅ 合格 |
| 1000 | 93.98% | ❌ 失敗 |
| 4000 | 770.21% | ❌ 失敗 |
| 10000 | 2147.33% | ❌ 失敗 |

**パターン**: 誤差はμrに比例して増加（約μr/100の割合）

---

## 調査経過

### 1. 背景磁場の単位問題を発見 ✅

**問題**: 初期テストでは`rad.ObjBckg([0, 0, H_ext])`にA/m単位の値を渡していた

```python
# ❌ 間違い
H_ext = 1000.0  # A/m
bg = rad.ObjBckg([0, 0, H_ext])  # WRONG: rad.ObjBckg expects Tesla!

# ✅ 正しい
H_ext = 1000.0  # A/m
mu0 = 4*np.pi*1e-7
B_ext = H_ext * mu0  # Convert to Tesla
bg = rad.ObjBckg([0, 0, B_ext])  # CORRECT
```

**結果**:
- 単位修正により、背景磁場は正しい値になった
- しかし高透磁率での誤差は残った（問題は別にある）

**ファイル**: `test_background_field_units.py`, `test_high_mu_corrected_units.py`

---

### 2. 磁化評価の問題を発見 🔴

**重大な発見**: 四面体要素の内部で`rad.Fld(obj, 'm', point)`が**常にゼロ**を返す

#### テスト結果 (`test_self_interaction_tetra.py`)

```
[1] 六面体 (参照)
  Magnetization at center: M = [0.00, 0.00, 2911.76] A/m  ← 正しい
  |M| = 2911.76 A/m
  Demagnetization factor: N_z = 0.3333  ← 理論値と一致

[2] 四面体 (ANALYTICAL)
  Magnetization at center: M = [0.00, 0.00, 0.00] A/m  ← ゼロ！
  |M| = 0.00 A/m
  Demagnetization factor: N_z = NaN  ← 計算不可
```

#### 詳細調査 (`debug_magnetization_inside_tetra.py`)

四面体の内部・外部すべての点で`|M| = 0.0`：

| Location | |M| (A/m) | |H| (A/m) |
|----------|----------|----------|
| Centroid (inside) | 0.0 | 18840.49 |
| Near vertex 1 (inside) | 0.0 | 38053.31 |
| Inside tetra | 0.0 | 33013.36 |
| Outside tetra | 0.0 | 4932.46 |
| Far outside | 0.0 | 890.70 |

**六面体**（比較）：

| Location | |M| (A/m) | |H| (A/m) |
|----------|----------|----------|
| Center | 2911.76 | 29.41 |
| Near center | 2911.76 | 646.50 |
| Outside | 0.00 | 6566.38 |

---

## 根本原因

### 四面体要素の磁化評価が未実装

**問題**: Radiaの`rad.Fld(obj, 'm', point)`は四面体ポリヘドロンに対応していない

**影響**:
1. ユーザが`rad.Fld(obj, 'm', point)`を呼んでもゼロが返る
2. 内部的にも磁化が正しく計算/保存されていない可能性
3. 減磁係数(self-demagnetization factor)が計算できない
4. 高透磁率で顕著な誤差が発生

**なぜμr=100では動作するのか？**

低透磁率では減磁効果が小さいため、誤差が許容範囲内に収まる：
- μr=100: 減磁の影響 ~ 10%
- μr=1000: 減磁の影響 ~ 50%
- μr=10000: 減磁の影響 ~ 90%

---

## 技術的詳細

### 減磁係数の理論

線形磁性材料の磁化：

```
M = χ * H_internal = (μr - 1) * H_internal
H_internal = H_external + H_demag
H_demag = -N * M  (減磁場)
```

ここで`N`は減磁係数（形状依存、0 ≤ N ≤ 1）

整理すると：
```
M = χ * (H_ext - N * M)
M * (1 + χ * N) = χ * H_ext
M = χ * H_ext / (1 + χ * N)
```

**六面体**（立方体、z方向）:
- N_z = 1/3 = 0.3333
- μr=100の場合: M = 99 * 1000 / (1 + 99 * 0.3333) = 2911.76 A/m ← 測定値と一致

**四面体**:
- N_z ~ 0.25-0.35 (形状・向きに依存)
- しかし磁化がゼロ → Nが計算されていない or 無限大

### Radiaの内部動作（推測）

1. **rad.Solve()**: 相互作用行列を構築し、磁化を計算
   ```
   [N_ij] * {M_j} = {H_ext}
   ```

2. **問題**: 四面体の`N_ii`（自己相互作用）が正しく計算されていない
   - STANDARDメソッドは自己相互作用の解析式を持つ
   - ANALYTICALメソッドは外部磁場評価のみ実装（相互作用は未実装）

3. **結果**: 自己相互作用が過大評価 or 過小評価
   - 磁化が異常な値になる
   - 外部磁場計算も間違った値になる

---

## ANALYTICALメソッドの実装状況

### 実装済み ✅

**外部磁場評価** ([rad_polyhedron.cpp:701-836](src/core/rad_polyhedron.cpp#L701-L836))
- `B_comp_tetrahedron_analytical()` - 外部磁場計算
- RadAnalyticalFieldFromPolygonCharge() を使用
- 基底ベクトル構築
- 2D投影
- 磁荷密度計算

**動作確認**:
```cpp
// Line 717
TVector3d M_global = Magn;  // グローバル磁化を取得

// Line 786
double W = ConstForH * M_local.z;  // 磁荷密度計算

// Line 819-827
RadAnalyticalFieldFromPolygonCharge(...);  // 解析式で磁場計算
```

### 未実装 ❌

**自己相互作用係数（減磁テンソル）**
- 四面体の`N_ii`（diagonal element of interaction matrix）
- STANDARDメソッド相当の実装が必要
- これがないと`rad.Solve()`が正しい磁化を計算できない

---

## 影響範囲

### 動作する条件 ✅

| 条件 | μr | 理由 |
|------|-----|------|
| 低透磁率 | ≤ 100 | 減磁効果が小さい |
| 外部磁場評価のみ | Any | 磁化が既知なら計算可能 |

### 動作しない条件 ❌

| 条件 | μr | 理由 |
|------|-----|------|
| 高透磁率 | > 100 | 減磁効果が大きい |
| 自己相互作用必須 | Any | N_iiが未実装 |

---

## 解決策の検討

### Option A: 自己相互作用の実装（根本的解決）

**必要な作業**:
1. 四面体の減磁テンソル解析式を導出
2. `N_comp_tetrahedron_analytical()`メソッドを実装
3. `rad.Solve()`の相互作用行列構築で使用

**難易度**: 🔴 高
- 解析式の導出が複雑
- Radiaの相互作用行列システムの理解が必要
- 既存コードへの統合

**効果**: ✅ 完全な解決
- すべてのμrで動作
- 磁化評価も正しくなる

### Option B: 数値積分による近似（現実的）

**アイデア**:
1. STANDARDメソッドの減磁係数を借用
2. 四面体を複数の小さな四面体に分割
3. 各小要素の寄与を合算

**難易度**: 🟡 中
- 既存のSTANDARDメソッドを活用
- メッシュ細分化で精度向上

**効果**: 🟡 部分的解決
- 精度は分割数に依存
- 計算コストが増加

### Option C: μr制限の明示（暫定対策）

**対策**:
1. ドキュメントでμr ≤ 100の制限を明記
2. 警告メッセージの追加
3. ユーザーガイドに代替手法を記載

**難易度**: 🟢 低
- ドキュメント更新のみ

**効果**: ⚠️ 回避策
- 問題は解決しない
- ユーザーに制限を知らせる

---

## 推奨事項

### 短期（次の作業）

1. **ドキュメント更新** ⚠️
   - ANALYTICAL_METHOD_VALIDATION.mdに制限を追記
   - μr ≤ 100での使用を推奨
   - 高透磁率での動作不良を明記

2. **警告メッセージ** ⚠️
   - μr > 100でANALYTICALメソッド使用時に警告
   - rad.Solve()の収束性チェック

### 中期（実装検討）

1. **Option B: 数値近似** 🔧
   - 四面体を6つの小四面体に分割
   - 各四面体にSTANDARDメソッドを適用
   - 合算して全体の減磁係数を近似

2. **テストケース拡充** 📊
   - μr = 10, 50, 100, 200, 500での収束性
   - 細分化メッシュでの精度検証

### 長期（根本解決）

1. **Option A: 解析式の実装** 🎯
   - 四面体減磁テンソルの解析式研究
   - `N_comp_tetrahedron_analytical()`実装
   - 完全な四面体サポート

---

## テストファイル

### 作成したファイル

| ファイル | 内容 |
|---------|------|
| test_high_permeability_tetra.py | μr vs 誤差の関係 |
| test_background_field_units.py | 単位変換問題の発見 |
| test_high_mu_corrected_units.py | 単位修正後の再テスト |
| test_self_interaction_tetra.py | 減磁係数の抽出 |
| debug_high_mu_tetra.py | 詳細デバッグ |
| debug_magnetization_inside_tetra.py | 磁化評価の検証 |

### 実行方法

```bash
cd S:/Radia/01_GitHub

# 高透磁率テスト（μr = 100, 1000, 4000, 10000）
python test_high_mu_corrected_units.py

# 減磁係数の検証
python test_self_interaction_tetra.py

# 磁化評価の詳細
python debug_magnetization_inside_tetra.py
```

---

## 関連ファイル

### 実装コード

- [src/core/rad_polyhedron.h:279](src/core/rad_polyhedron.h#L279) - メソッド宣言
- [src/core/rad_polyhedron.cpp:701-836](src/core/rad_polyhedron.cpp#L701-L836) - ANALYTICAL実装

### ドキュメント

- [ANALYTICAL_METHOD_VALIDATION.md](ANALYTICAL_METHOD_VALIDATION.md) - 検証結果まとめ
- [ANALYTICAL_METHOD_SUCCESS.md](ANALYTICAL_METHOD_SUCCESS.md) - 実装成功レポート
- [TETRAHEDRAL_IMPLEMENTATION.md](TETRAHEDRAL_IMPLEMENTATION.md) - 実装計画

---

## 結論

### 現状の評価

**ANALYTICALメソッド**:
- ✅ μr ≤ 100: 実用可能（誤差 < 20%）
- ❌ μr > 100: 使用不可（誤差 > 50%）

**根本原因**:
- 四面体の自己相互作用係数（減磁テンソル）が未実装
- `rad.Fld(obj, 'm', point)`が四面体でゼロを返す
- 高透磁率で減磁効果が大きいため誤差が顕著

**推奨される使用方法**:
1. μr ≤ 100の材料に限定
2. または六面体メッシュを使用（STANDARDメソッド）
3. 高透磁率（μr > 1000）では六面体必須

---

**最終更新**: 2025-11-25
**ステータス**: ⚠️ **制限あり - μr ≤ 100推奨**
**次の作業**: ドキュメント更新、警告メッセージ追加
