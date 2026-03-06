# Radia 練習例題 (Practice Examples)

このフォルダには、Wolfram Language (Mathematica) から Python に変換された Radia の練習例題が含まれています。

## ファイル一覧

### オリジナルファイル (Wolfram Language)
- `case0.wls` - アーク電流と矩形磁石
- `case1.wls` - アーク電流と2つの矩形磁石
- `case2.wls` - 面取り付き多重押出形状
- `case3.wls` - 多面体（立方体）

### 変換後ファイル (Python)
- `case0.py` - アーク電流と矩形磁石
- `case1.py` - アーク電流と2つの矩形磁石
- `case2.py` - 面取り付き多重押出形状
- `case3.py` - 多面体（立方体）

## 実行方法

各Pythonファイルは単独で実行できます：

```bash
cd examples/2019_11_29_Radia_練習
python case0.py
python case1.py
python case2.py
python case3.py
```

または、Radiaモジュールのパスを明示的に指定する必要がある場合：

```bash
export PYTHONPATH="../../dist:../../build/lib/Release:$PYTHONPATH"
python case0.py
```

## 例題の説明

### case0.py - Arc Current with Rectangular Magnet
- アーク状の電流要素と矩形磁石を作成
- 線形材料特性を適用
- 原点での磁場を計算

**主要関数:**
- `rad.ObjArcCur()` - アーク電流要素の作成
- `rad.ObjRecMag()` - 矩形磁石の作成
- `rad.MatLin()` - 線形材料の定義
- `rad.MatApl()` - 材料の適用
- `rad.Fld()` - 磁場計算

### case1.py - Arc Current with Two Rectangular Magnets
- アーク電流と2つの異なる位置の矩形磁石
- 複数オブジェクトをコンテナで管理
- 線形材料特性を適用

**主要関数:**
- `rad.ObjCnt()` - オブジェクトコンテナの作成

### case2.py - Multiple Extrusion with Chamfer
- 面取りを含む複雑な押出形状
- 複数の断面を定義して押出
- 要素の分割

**主要関数:**
- `rad.ObjMltExtRtg()` - 多重押出矩形の作成
- `rad.ObjDivMag()` - 磁石の分割

### case3.py - Polyhedron (Cube)
- 頂点と面の定義から多面体を作成
- 立方体の例

**主要関数:**
- `rad.ObjPolyhdr()` - 多面体の作成

## Wolfram Language から Python への主な変換ルール

1. **関数名の変換**
   - Wolfram: `radObjRecMag` → Python: `rad.ObjRecMag`
   - Wolfram: `radUtiDelAll` → Python: `rad.UtiDelAll`

2. **リストの記法**
   - Wolfram: `{1,2,3}` → Python: `[1, 2, 3]`

3. **数学定数**
   - Wolfram: `Pi` → Python: `math.pi`

4. **モジュールのインポート**
   - Wolfram: `<<Radia\`` → Python: `import radia as rad`

5. **3D描画**
   - Wolfram: `Graphics3D[radObjDrw[g]]`
   - Python: 現時点では未実装（matplotlib等の追加ライブラリが必要）

## 注意事項

- case2.py と case3.py では、オリジナルのスクリプトで `g2` が未定義のまま磁場計算を試みているため、その部分はコメントアウトされています。
- 3D可視化機能は現在実装されていません。必要に応じて matplotlib や mayavi 等のライブラリを使用してください。
- 出力ファイル `out.dat` には磁場の計算結果が保存されます。

## 依存関係

- Python 3.12
- Radia モジュール (../../dist/radia.pyd)
- NumPy (配列操作用)
- Math (数学関数用)

## トラブルシューティング

### ModuleNotFoundError: No module named 'radia'

Radiaモジュールのパスが正しく設定されていません。以下を確認してください：

1. `../../dist/radia.pyd` が存在する
2. または `../../build/lib/Release/radia.pyd` が存在する

### 磁場計算結果がゼロ

- 磁石に磁化が設定されているか確認
- 材料特性が正しく適用されているか確認
- オブジェクトIDが正しいか確認

## 参考

- [Radia 公式ドキュメント](https://www.esrf.fr/Accelerators/Groups/InsertionDevices/Software/Radia)
- オリジナルの Mathematica スクリプト: `case*.wls`
