# Wolfram Language to Python Conversion Notes

## 変換日: 2025-01-29

## 変換したファイル

| Wolfram Language | Python | 状態 | 備考 |
|-----------------|--------|------|------|
| case0.wls | case0.py | ✓ 完了 | 正常動作確認済み |
| case1.wls | case1.py | ✓ 完了 | 正常動作確認済み |
| case2.wls | case2.py | ✓ 完了 | g2未定義の問題をコメント化 |
| case3.wls | case3.py | ✓ 完了 | g2未定義の問題をコメント化 |

## 主な変換内容

### 1. 基本構文の変換

```mathematica
(* Wolfram Language *)
<<Radia`;
radUtiDelAll[];
rmin = 100;
g1 = radObjArcCur[{0,0,0}, {rmin, rmax}, {phimin, phimax}, h, nseg, j];
```

```python
# Python
import radia as rad
rad.UtiDelAll()
rmin = 100
g1 = rad.ObjArcCur([0, 0, 0], [rmin, rmax], [phimin, phimax], h, nseg, j)
```

### 2. 数学関数・定数

| Wolfram | Python | 説明 |
|---------|--------|------|
| `Pi` | `math.pi` | 円周率 |
| `Sin[x]` | `math.sin(x)` | 正弦関数 |
| `Cos[x]` | `math.cos(x)` | 余弦関数 |

### 3. Radia関数のマッピング

| Wolfram Language | Python | 機能 |
|-----------------|--------|------|
| `radUtiDelAll[]` | `rad.UtiDelAll()` | 全オブジェクト削除 |
| `radObjArcCur[...]` | `rad.ObjArcCur(...)` | アーク電流作成 |
| `radObjRecMag[...]` | `rad.ObjRecMag(...)` | 矩形磁石作成 |
| `radMatLin[...]` | `rad.MatLin(...)` | 線形材料定義 |
| `radMatSatIso[...]` | `rad.MatSatIso(...)` | 非線形等方性材料 |
| `radMatApl[obj, mat]` | `rad.MatApl(obj, mat)` | 材料適用 |
| `radObjDrwAtr[...]` | `rad.ObjDrwAtr(...)` | 描画属性設定 |
| `radObjCnt[{...}]` | `rad.ObjCnt([...])` | コンテナ作成 |
| `radObjMltExtRtg[...]` | `rad.ObjMltExtRtg(...)` | 多重押出 |
| `radObjDivMag[...]` | `rad.ObjDivMag(...)` | 磁石分割 |
| `radObjPolyhdr[...]` | `rad.ObjPolyhdr(...)` | 多面体作成 |
| `radFld[obj, "bxbybz", pt]` | `rad.Fld(obj, 'b', pt)` | 磁場計算 |

### 4. データ構造の変換

```mathematica
(* Wolfram: リスト *)
points = {{1,2,3}, {4,5,6}};
faces = {{1,2,3,4}, {5,6,7,8}};
```

```python
# Python: リスト
points = [[1, 2, 3], [4, 5, 6]]
faces = [[1, 2, 3, 4], [5, 6, 7, 8]]
```

### 5. 可視化

**Wolfram Language:**
```mathematica
t = Show[Graphics3D[radObjDrw[g]]];
Export["3DPlot.png", t];
```

**Python:**
現時点では直接的な等価機能なし。将来的には以下のライブラリを使用可能：
- matplotlib + mplot3d
- mayavi
- plotly

### 6. ファイル出力

**Wolfram Language:**
```mathematica
t = radFld[g2, "bxbybz", {0,0,0}]
Export["out.dat", t]
```

**Python:**
```python
field = rad.Fld(g2, 'b', [0, 0, 0])
with open('out.dat', 'w') as f:
    f.write(f"{field[0]}\t{field[1]}\t{field[2]}\n")
```

## 発見された問題点

### case2.wls と case3.wls

**問題:** `g2` オブジェクトが定義されていないのに `radFld[g2, ...]` を呼び出している

**対応:** 
- Python版では該当部分をコメントアウト
- README.mdに注意事項を記載

```python
# Note: g2 is not defined in the original script
# field = rad.Fld(g2, 'b', [0, 0, 0])  # This would fail
print("Note: Field calculation requires defining g2 object")
```

## コーディングスタイル

- インデント: タブ（1タブ = 4スペース相当）
- エンコーディング: UTF-8
- 行末: LF (Unix style)
- ドキュメント文字列: 英語で記述
- コメント: 日本語・英語併用

## テスト結果

### case0.py
```
Container object ID: 4
Magnetic field at origin: Bx=0.000000e+00, By=0.000000e+00, Bz=0.000000e+00 T
Calculation complete. Field data saved to out.dat
```
✓ 正常動作

### case1.py
```
Container object ID: 6
Magnetic field at origin: Bx=0.000000e+00, By=0.000000e+00, Bz=0.000000e+00 T
Calculation complete. Field data saved to out.dat
```
✓ 正常動作

### case2.py
```
Object ID: 3
Geometry created and subdivided successfully
Note: Field calculation requires defining g2 object
```
✓ 正常動作（g2未定義の警告付き）

### case3.py
```
Polyhedron object ID: 1
Cube polyhedron created successfully
Vertices: 8
Faces: 6
Note: Field calculation requires defining g2 object
```
✓ 正常動作（g2未定義の警告付き）

## 今後の改善点

1. **3D可視化機能の追加**
   - matplotlib による基本的な3Dプロット
   - インタラクティブな可視化（mayavi, plotly）

2. **エラーハンドリングの強化**
   - より詳細なエラーメッセージ
   - 例外処理の追加

3. **追加機能**
   - フィールドマップの生成
   - 等高線プロット
   - データのCSV出力

4. **ドキュメントの充実**
   - 各関数の詳細な使用例
   - Radiaの概念説明
   - チュートリアルの追加

---

**変換ツール:** 手動変換 + Python スクリプト（タブ変換）
**検証:** Python 3.12 + Radia 4.32
**状態:** 完了
