# VTK Export 機能テスト結果

## テスト日: 2025-10-29

## 結論

**`exportGeometryToVTK`は完全に実装されており、正常に動作します！**

## テスト結果

### rad.ObjDrwVTK() - ✓ 実装済み・動作確認

**機能**: RadiaオブジェクトからVTKデータを抽出

**テスト内容**:
```python
g = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1])
vtk_data = rad.ObjDrwVTK(g, 'Axes->False')
```

**結果**:
- ✓ 関数は存在する
- ✓ 正常に実行される
- ✓ データ構造が返される

**返されるデータ**:
```python
vtk_data.keys() = ['polygons', 'lines']
vtk_data['polygons'].keys() = ['colors', 'lengths', 'vertices']
```

**詳細**:
- 頂点数: 72個（生データ）
- ポリゴン数: 6個（立方体の6面）
- 実際のポイント数: 24個（VTKファイル内）

---

### exportGeometryToVTK() - ✓ 実装済み・動作確認

**場所**: `my_module.py`

**機能**: RadiaオブジェクトをVTK Legacy形式ファイルに出力

**実装方法**:
1. `rad.ObjDrwVTK()`でデータを取得
2. VTK Legacy形式に変換
3. `.vtk`ファイルに出力

**テスト結果**:
```
SUCCESS: VTK file created: test_export.vtk (690 bytes)
```

**生成されたVTKファイルの構造**:
```
# vtk DataFile Version 5.1
vtk output
ASCII
DATASET POLYDATA
POINTS 24 float
[頂点座標データ...]

POLYGONS 7 24
OFFSETS vtktypeint64
[オフセットデータ...]
CONNECTIVITY vtktypeint64
[接続性データ...]

CELL_DATA 6
COLOR_SCALARS Radia_colours 3
[色データ...]
```

**検証**:
- ファイルサイズ: 690バイト
- フォーマット: VTK Legacy ASCII
- 頂点数: 24
- ポリゴン数: 6
- 色情報: 含まれる

---

## 使用方法

### 基本的な使い方

```python
import radia as rad
from my_module import exportGeometryToVTK

# Radiaオブジェクトを作成
g = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1])

# VTKファイルに出力
exportGeometryToVTK(g, 'my_geometry')
# -> my_geometry.vtk が生成される
```

### Paraviewでの読み込み

生成された`.vtk`ファイルは、以下のソフトウェアで読み込み可能：
- **Paraview** (推奨) - 科学可視化ツール
- **VTK** - Visualization Toolkit
- **MayaVi** - Pythonベースの可視化ツール

---

## 対応機能

### ✓ 実装されている機能

1. **rad.ObjDrwVTK()** - VTKデータ抽出
2. **exportGeometryToVTK()** - VTKファイル出力
3. **ポリゴンデータ** - 頂点、接続性、面
4. **色情報** - オブジェクトの色
5. **VTK Legacy形式** - 広く互換性のある形式

### ✗ 実装されていない機能

1. **rad.ObjDrwOpenGL()** - OpenGLビューア（プラットフォーム依存）
   - エラー: "This function is not implemented on that platform."

---

## ファイル構成

### my_module.py の内容

```python
def exportGeometryToVTK(obj, fileName='radia_Geometry'):
    '''
    Writes the geometry of RADIA object "obj" to file fileName.vtk
    for use in Paraview. The format is VTK legacy.
    '''
    vtkData = rad.ObjDrwVTK(obj, 'Axes->False')
    # ... VTK Legacy形式に変換して出力 ...
```

**依存関係**:
- `radia` - Radiaモジュール
- `csv` - CSV書き込み
- `itertools.accumulate` - オフセット計算

---

## 実用例

### 振分電磁石の可視化

元のスクリプト `2024_02_01_振り分け電磁石.py` では：

```python
g = rad.ObjCnt([coils, york])
exportGeometryToVTK(g, 'RADIA')
# -> RADIA.vtk が生成される（22 MB）
```

このファイルは実際に生成されています：
```
-rw-r--r-- 1 Administrator 197121 22178666  2月  2  2024 RADIA.vtk
```

---

## まとめ

| 機能 | 状態 | 備考 |
|------|------|------|
| `rad.ObjDrwVTK()` | ✓ 動作 | VTKデータ抽出 |
| `exportGeometryToVTK()` | ✓ 動作 | VTKファイル出力 |
| VTK Legacy形式 | ✓ 対応 | Paraview互換 |
| ポリゴンメッシュ | ✓ 対応 | 頂点、面、色 |
| `rad.ObjDrwOpenGL()` | ✗ 未実装 | プラットフォーム依存 |

**結論**: VTKエクスポート機能は完全に実装されており、Paraviewなどでの可視化が可能です。

---

**テスト環境:**
- Radia: 4.32 (OpenMP 2.0並列化版)
- Python: 3.12
- OS: Windows 10 Build 20348

**テスト実施者:** Claude Code
**テスト完了日:** 2025-10-29
