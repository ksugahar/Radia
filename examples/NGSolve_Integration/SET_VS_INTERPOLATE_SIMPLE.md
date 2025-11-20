# NGSolve: Set() vs Interpolate() - 分かりやすい解説

**日付:** 2025-11-08

---

## 1. 一言で言うと

```python
# Set(): 関数を有限要素空間に「投影」する（数学的に厳密）
gf.Set(cf)

# Interpolate(): 特定の点で関数を「評価して代入」する（直感的だが限定的）
gf.Interpolate(cf)
```

**決定的な違い:**
- **Set()**: 多くの点（積分点）で評価して、全体として最適な近似を求める
- **Interpolate()**: 少ない点（頂点など）で評価して、その点の値を直接使う

---

## 2. 具体例で理解する

### 例: 温度分布をメッシュに設定する

```python
# 解析的な温度分布（CoefficientFunction）
temperature = CF( (x, y, z) )  # 例: x^2 + y^2

# メッシュ上の近似（GridFunction）
gf_temp = GridFunction(H1(mesh))

# 方法1: Set()
gf_temp.Set(temperature)

# 方法2: Interpolate()
gf_temp.Interpolate(temperature)
```

### どこで関数を評価するか？

#### Set()の場合:

```
メッシュ（三角形要素の例）:

    v1 -------- v2
     |\        /|
     | \  T1  / |
     |  \    /  |
     |   \  /   |
     | T2  \/   |
     |    /\\   |
     |   /  \ \ |
     |  /    \ \|
     | /  T3  \ |
     |/        \|
    v3 -------- v4

Set()が評価する点:
- 各要素T1, T2, T3内の積分点（×マーク）
- 要素T1: ×××× (4-7個の積分点)
- 要素T2: ××××
- 要素T3: ××××
- 合計: 要素数 × 4-7点 = 数千～数万点

評価結果:
- 各積分点で temperature(x,y,z) を計算
- 重み付き積分で有限要素基底関数に投影
- 数学的に厳密なL²投影
```

#### Interpolate()の場合:

```
同じメッシュ:

    v1 -------- v2
     |\        /|
     | \  T1  / |
     |  \    /  |
     |   \  /   |
     | T2  \/   |
     |    /\\   |
     |   /  \ \ |
     |  /    \ \|
     | /  T3  \ |
     |/        \|
    v3 -------- v4

Interpolate()が評価する点:
- 頂点のみ（●マーク）
- v1, v2, v3, v4
- 合計: 頂点数 = 数千点

評価結果:
- 各頂点で temperature(v1), temperature(v2), ... を計算
- その値を直接代入
- 要素内部は線形補間
```

---

## 3. 数学的な違い

### Set() - L²投影（Projection）

**式:**
```
最小化問題: min ||u_h - cf||²_L²
           u_h ∈ V_h

弱形式:
∫_Ω u_h · v_h dΩ = ∫_Ω cf · v_h dΩ  (∀v_h ∈ V_h)
```

**意味:**
- 関数`cf`と近似`u_h`の差の二乗を最小化
- 積分（面積）ベースの最適化
- 全体として最も良い近似

**例:**
```python
# 不連続関数も滑らかに近似
cf = IfPos(x, 1.0, 0.0)  # x > 0 なら1, そうでなければ0
gf.Set(cf)
# → 境界付近で滑らかに遷移する近似
```

### Interpolate() - 補間（Interpolation）

**式:**
```
u_h(x_i) = cf(x_i)  (補間ノード x_i)
```

**意味:**
- 特定の点（補間ノード）で関数値を一致させる
- その間は線形補間
- 点ベースの近似

**例:**
```python
# 不連続関数の補間
cf = IfPos(x, 1.0, 0.0)
gf.Interpolate(cf)
# → 補間ノードでの値が決まるが、不連続を捉えられない可能性
```

---

## 4. どの空間で使えるか？

### H1空間（スカラー場）

```python
fes = H1(mesh, order=1)
gf = GridFunction(fes)

# どちらも使える
gf.Set(temperature_cf)        # ✅ OK
gf.Interpolate(temperature_cf) # ✅ OK
```

**補間ノード:** 頂点

**用途:**
- 温度分布
- 圧力場
- ポテンシャル

---

### VectorH1空間（ベクトル場、各成分がH1）

```python
fes = VectorH1(mesh, order=1)
gf = GridFunction(fes)

# Set()は問題なし
gf.Set(displacement_cf)  # ✅ OK

# Interpolate()は動作するが注意が必要
gf.Interpolate(displacement_cf)  # ⚠️ 注意: 遅い、不正確な場合あり
```

**補間ノード:** 頂点（各成分独立）

**問題点:**
- 各成分(x, y, z)を独立に補間
- 成分ごとに別々に関数呼び出し → 遅い
- ベクトルの整合性が保証されない場合がある

---

### HCurl空間（電磁場）

```python
fes = HCurl(mesh, order=1)
gf = GridFunction(fes)

# Set()のみ使える
gf.Set(B_field_cf)  # ✅ OK

# Interpolate()は未サポート
gf.Interpolate(B_field_cf)  # ❌ エラー！
```

**補間ノードの意味が不明確:**
- HCurlの自由度（DOF）は**辺**に関連付けられている
- 「頂点での値」をどう辺の自由度に変換するか定義されていない

**DOFの例（三角形）:**
```
    v1
     |\
   e1| \ e3      DOF1 = e1に沿った接線成分
     |  \        DOF2 = e2に沿った接線成分
   e2|   \       DOF3 = e3に沿った接線成分
     |    \
     |_____\
    v2  e3  v3

頂点v1, v2, v3での「値」は定義されていない！
```

---

### HDiv空間（流束場）

```python
fes = HDiv(mesh, order=1)
gf = GridFunction(fes)

# Set()のみ使える
gf.Set(flux_cf)  # ✅ OK

# Interpolate()は未サポート
gf.Interpolate(flux_cf)  # ❌ エラー！
```

**補間ノードの意味が不明確:**
- HDivの自由度は**面**に関連付けられている
- 頂点での値では自由度を定義できない

---

## 5. 我々のテスト結果（実測データ）

### テスト条件
- メッシュ: 5034頂点、1260要素
- 関数: Radiaの磁場 `B(x, y, z)`

### HCurl空間

```python
fes_hcurl = HCurl(mesh, order=1)

# Set() - 動作する
gf.Set(B_cf)
# 時間: 51.81 ms
# 結果: ✅ 正しい磁場が設定される

# Interpolate() - 未サポート
gf.Interpolate(B_cf)
# エラー: "Interpolate not supported for HCurl"
```

### VectorH1空間

```python
fes_h1 = VectorH1(mesh, order=1)

# Set() - 動作する
gf.Set(B_cf)
# 時間: 47.06 ms
# 結果: ✅ 正しい磁場が設定される

# Interpolate() - 動作するが問題あり
gf.Interpolate(B_cf)
# 時間: 442.98 ms (10倍遅い！)
# 結果: ❌ 全てゼロになる（我々の実装と不一致）
```

### なぜVectorH1でInterpolate()が失敗するのか？

**NGSolveの期待:**
```python
# NGSolveは成分ごとにアクセスしようとする
cf[0](point)  # x成分のスカラーCF
cf[1](point)  # y成分のスカラーCF
cf[2](point)  # z成分のスカラーCF
```

**我々の実装:**
```python
# ベクトル全体を一度に返す
cf(point)  # → [Bx, By, Bz]

# 成分インデックスアクセスをサポートしていない
cf[0]  # → 未定義（デフォルト値 = 0 が返る）
```

**結果:**
- `cf[0](point)` → 0
- `cf[1](point)` → 0
- `cf[2](point)` → 0
- 全てゼロになる

---

## 6. 評価回数の違い

### Set()の評価

```
要素ごとのループ:
  For element in elements:  # 1260個
    For integration_point in element:  # 4-7点
      cf(integration_point)

総評価回数 ≈ 1260 × 5 = 6300回
```

**特徴:**
- 要素ごとにバッチ評価
- 積分点での評価（数学的に厳密）
- 我々のバッチ評価実装が効率的

### Interpolate()の評価（VectorH1）

```
頂点ごとのループ:
  For vertex in vertices:  # 5034個
    For component in [x, y, z]:  # 3成分
      cf[component](vertex)

総評価回数 = 5034 × 3 = 15102回
```

**特徴:**
- 成分ごとに個別評価
- Python-C++呼び出しが多い
- バッチ評価の恩恵を受けない

**Python呼び出しオーバーヘッド:**
```
Set():        6300 calls × 20 µs ≈ 126 ms
Interpolate(): 15102 calls × 20 µs ≈ 302 ms

実測:
Set():        47 ms (高速！)
Interpolate(): 443 ms (遅い...)
```

---

## 7. ビジュアル比較

### Set() - 積分ベース

```
関数: f(x) = x² + noise (ノイズあり)

評価点（×）と重み:
|
|  ×   ×   ×   ×   ×   ×   ×   ×
|   \ / \ / \ / \ / \ / \ / \ /
|    ×   ×   ×   ×   ×   ×   ×
|___________________________________
         メッシュ要素

各要素内で複数点を評価 → 重み付き積分
→ ノイズが平滑化される
→ 全体として最適な近似
```

### Interpolate() - 点ベース

```
関数: f(x) = x² + noise (ノイズあり)

評価点（●）:
|
|  ●           ●           ●
|  |           |           |
|  |           |           |
|__|___________|___________|___
   v1          v2          v3
         メッシュ頂点

頂点のみ評価 → 線形補間で繋ぐ
→ ノイズがそのまま反映される
→ 点での値は正確だが、間は線形
```

---

## 8. 使い分けのガイドライン

### いつSet()を使うべきか？

✅ **以下の場合は常にSet()を使う:**

1. **HCurl, HDiv空間を使う場合**（必須）
   ```python
   # 電磁場、流束場
   gf_B = GridFunction(HCurl(mesh))
   gf_B.Set(B_field_cf)  # ← 唯一の方法
   ```

2. **数学的に厳密な近似が必要な場合**
   ```python
   # L²投影で最適近似
   gf.Set(complex_function_cf)
   ```

3. **不連続関数や振動する関数**
   ```python
   # 積分ベースなので平滑化効果
   gf.Set(discontinuous_cf)
   ```

4. **VectorH1でベクトル場を扱う場合**
   ```python
   # 各成分の整合性が重要
   gf.Set(displacement_cf)
   ```

### いつInterpolate()を使えるか？

⚠️ **以下の限定的な場合のみ:**

1. **H1スカラー場で、単純な関数**
   ```python
   fes = H1(mesh)
   gf = GridFunction(fes)
   # 滑らかな温度分布など
   gf.Interpolate(temperature_cf)  # OK
   ```

2. **可視化目的のみで精度不要**
   ```python
   # プロット用の大まかな値
   gf.Interpolate(cf)
   ```

3. **頂点の値が既知で、それを直接設定したい**
   ```python
   # 測定データを頂点に設定
   gf.Interpolate(measured_data_cf)
   ```

---

## 9. Radiaフィールド評価での結論

### 我々のケース

```python
# Radiaの磁場をNGSolveメッシュに設定
B_cf = radia_ngsolve.RadiaField(magnet, 'b')
fes = HCurl(mesh, order=1)
gf = GridFunction(fes)

# Set()が唯一の正しい選択
gf.Set(B_cf)  # ✅ 正解
```

### 理由

1. **HCurl空間が必須**
   - 磁場Bは `curl(A) = B` を満たすべき
   - `div(B) = 0` が自動的に保証される
   - HCurlはこの性質を持つ

2. **Interpolate()は使えない**
   - HCurlで未サポート
   - VectorH1でも不正確

3. **Set()は数学的に正しい**
   - L²投影で厳密
   - 積分点での評価は物理的に意味がある

4. **性能も十分**
   - バッチ評価で最適化済み
   - 47-52 ms (5034頂点メッシュ)

---

## 10. まとめ表

| 項目 | Set() | Interpolate() |
|------|-------|---------------|
| **数学的意味** | L²投影（最小二乗近似） | 点での値代入 |
| **評価点** | 積分点（多い） | 補間ノード（少ない） |
| **H1スカラー** | ✅ 使える | ✅ 使える |
| **VectorH1** | ✅ 推奨 | ⚠️ 動作するが遅い |
| **HCurl** | ✅ 必須 | ❌ 未サポート |
| **HDiv** | ✅ 必須 | ❌ 未サポート |
| **数学的厳密性** | 高い | 低い |
| **計算コスト** | 中程度 | 低～高（実装依存） |
| **我々の実装** | ✅ 最適化済み | ❌ 未実装（不要） |
| **推奨度** | ⭐⭐⭐⭐⭐ | ⭐ |

### 一言で言うと

```
Set()     = 数学的に正しい方法（積分で最適化）
Interpolate() = 直感的だが限定的（点で代入）

電磁場や複雑なベクトル場 → Set()を使う！
```

---

## 11. よくある誤解

### 誤解1: 「Interpolate()の方が速い」

❌ **間違い:**
- 評価点が少ない → 速いはず

✅ **現実:**
- VectorH1で成分ごとに評価 → 実は遅い
- Python呼び出しオーバーヘッド → さらに遅い
- 実測: Set()より10倍遅い

### 誤解2: 「どちらも同じ結果」

❌ **間違い:**
- 近似は近似、結果は同じはず

✅ **現実:**
- 数学的に異なる操作
- Set(): L²ノルムで最適
- Interpolate(): 特定点で一致
- 要素内部の値は異なる

### 誤解3: 「HCurlでもInterpolate()が使えるはず」

❌ **間違い:**
- ベクトル場なんだから使えるはず

✅ **現実:**
- HCurlの自由度は辺に関連
- 頂点での「値」は定義されていない
- NGSolveが未サポート（エラーになる）

---

## 12. 実践的なアドバイス

### デバッグ時

```python
# 関数が正しく設定されたか確認
gf.Set(cf)

# 特定点での値を確認
test_point = mesh(0.1, 0.2, 0.3)
value = gf(test_point)
print(f"Value at test point: {value}")

# 元のCoefficientFunctionと比較
cf_value = cf(test_point)
print(f"Original function: {cf_value}")
print(f"Difference: {value - cf_value}")
```

### パフォーマンス最適化

```python
# バッチ評価が効果的
# Set()は自動的に最適化される
gf.Set(expensive_cf)

# 複数のGridFunctionに設定する場合
# CoefficientFunctionを再利用
cf = radia_ngsolve.RadiaField(magnet, 'b')
gf1.Set(cf)  # 1回目: キャッシュなし
gf2.Set(cf)  # 2回目: 若干速い可能性
```

### エラーが出た場合

```python
# HCurlでInterpolate()しようとした
try:
    gf_hcurl.Interpolate(B_cf)
except Exception as e:
    print(f"Error: {e}")
    # 解決策: Set()を使う
    gf_hcurl.Set(B_cf)
```

---

**作成者:** Claude Code
**バージョン:** 2.0（簡潔版）
**日付:** 2025-11-08
