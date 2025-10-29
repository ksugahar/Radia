# PyVista ビューア - rad.ObjDrwOpenGL() の代替

## 概要

**PyVistaを使用したRadiaオブジェクトの可視化は、`rad.ObjDrwOpenGL()`の優れた代替手段です。**

## なぜPyVistaなのか？

### ✓ PyVistaの利点

1. **活発にメンテナンスされている**
   - オープンソースプロジェクトとして継続的に更新
   - OpenGLビューアを自前で実装・保守する必要がない

2. **既存機能を活用**
   - `rad.ObjDrwVTK()`が既に実装されている
   - VTKデータをそのまま利用できる

3. **豊富な機能**
   - インタラクティブな3D表示（回転、ズーム、パン）
   - スクリーンショット保存
   - オフスクリーンレンダリング（GUI不要）
   - Jupyter Notebook対応
   - アニメーション作成

4. **クロスプラットフォーム**
   - Windows, Linux, macOS全てで動作
   - プラットフォーム依存の問題がない

5. **拡張性**
   - VTKの全機能にアクセス可能
   - カスタマイズが容易

### ✗ rad.ObjDrwOpenGL() の問題点

1. プラットフォーム依存（現在未実装）
2. メンテナンスコストが高い
3. 機能が限定的
4. エラー: "This function is not implemented on that platform."

## テスト結果

### 動作確認

```
PyVista version: 0.46.3

Creating test objects...
  Arc current: ID 1
  Rectangular magnet: ID 2
  Container: ID 3

Extracting VTK data from Radia...
  Vertices: 608
  Polygons: 152

Creating PyVista mesh...
  Mesh created: 608 points, 152 cells

SUCCESS: Image file created (72360 bytes)
```

### 生成された画像

- ファイル名: `radia_pyvista_test.png`
- サイズ: 71 KB
- 解像度: デフォルト（調整可能）

## 使用方法

### 1. インタラクティブ表示（GUI）

```python
from radia_pyvista_viewer import view_radia_object
import radia as rad

# Radiaオブジェクトを作成
g = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1])

# PyVistaで表示
view_radia_object(g)
```

**操作方法:**
- 左クリック + ドラッグ: 回転
- 右クリック + ドラッグ: パン
- スクロールホイール: ズーム
- 'q': 終了

### 2. オフスクリーンレンダリング（GUI不要）

```python
import pyvista as pv
import radia as rad
import numpy as np

# Radiaオブジェクトを作成
g = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1])

# VTKデータを取得
vtk_data = rad.ObjDrwVTK(g, 'Axes->False')
poly = vtk_data['polygons']

# PyVistaメッシュを作成
vertices = np.array(poly['vertices']).reshape(-1, 3)
# ... (faces構築) ...
mesh = pv.PolyData(vertices, faces)

# オフスクリーンプロッタ
plotter = pv.Plotter(off_screen=True)
plotter.add_mesh(mesh, color='lightblue', show_edges=True)
plotter.screenshot('output.png')
```

### 3. VTKファイルから読み込み

```python
from radia_pyvista_viewer import view_radia_vtk_file

# exportGeometryToVTK()で生成したファイルを表示
view_radia_vtk_file('RADIA.vtk')
```

## 提供スクリプト

### radia_pyvista_viewer.py

**主要関数:**

#### `view_radia_object(obj, plotter=None, show=True, color=None, opacity=1.0)`

Radiaオブジェクトを表示

**パラメータ:**
- `obj`: Radia オブジェクトID
- `plotter`: 既存のプロッタ（複数オブジェクト表示用）
- `show`: 即座に表示するか
- `color`: 色指定（例: 'red', [1,0,0]）
- `opacity`: 不透明度 (0-1)

**戻り値:**
- PyVistaプロッタオブジェクト

#### `view_radia_vtk_file(filename)`

VTKファイルを読み込んで表示

**パラメータ:**
- `filename`: .vtkファイルのパス

### test_pyvista_viewer.py

オフスクリーンレンダリングのデモ
- GUI不要
- スクリーンショット保存
- テスト用

## インストール

```bash
pip install pyvista
```

**依存関係:**
- Python 3.8+
- NumPy
- VTK (PyVistaが自動インストール)

## 実用例

### 元のスクリプトの置き換え

**元のコード:**
```python
rad.ObjDrwOpenGL(g)  # エラー: 未実装
```

**PyVistaに置き換え:**
```python
from radia_pyvista_viewer import view_radia_object
view_radia_object(g)
```

### 振分電磁石の可視化

```python
from radia_pyvista_viewer import view_radia_object
import radia as rad
from COIL import cCOIL

# コイルを作成
# ... (コイル作成コード) ...
coils = rad.ObjCnt(coils_list)

# ヨークを作成（Cubit使用）
# ... (メッシュ生成コード) ...
york = rad.ObjCnt(york_list)

# 全体を表示
g = rad.ObjCnt([coils, york])
view_radia_object(g)

# または、別々の色で表示
plotter = view_radia_object(coils, color='red', show=False)
plotter = view_radia_object(york, plotter=plotter, color='blue', show=True)
```

## 高度な使用例

### アニメーション作成

```python
import pyvista as pv
import radia as rad

g = rad.ObjRecMag([0,0,0], [10,10,10], [0,0,1])

plotter = pv.Plotter()
# ... メッシュ追加 ...

# 回転アニメーション
plotter.open_gif("rotation.gif")
for angle in range(0, 360, 10):
    plotter.camera.azimuth = angle
    plotter.write_frame()
plotter.close()
```

### カスタムカメラ設定

```python
plotter = view_radia_object(g, show=False)

# カスタムカメラ位置
plotter.camera_position = [
    (500, 500, 500),  # カメラ位置
    (0, 0, 0),        # 焦点
    (0, 0, 1)         # 上方向
]

plotter.show()
```

### 高解像度スクリーンショット

```python
plotter = pv.Plotter(off_screen=True, window_size=[3840, 2160])
# ... メッシュ追加 ...
plotter.screenshot('high_res.png', transparent_background=True)
```

## まとめ

| 機能 | rad.ObjDrwOpenGL() | PyVista |
|------|-------------------|---------|
| 実装状況 | ✗ 未実装 | ✓ 動作 |
| メンテナンス | 自前で必要 | コミュニティ |
| インタラクティブ | - | ✓ |
| オフスクリーン | - | ✓ |
| スクリーンショット | - | ✓ |
| アニメーション | - | ✓ |
| Jupyter対応 | - | ✓ |
| クロスプラットフォーム | ? | ✓ |

## 推奨事項

1. **新規プロジェクト**: PyVistaを使用
2. **既存スクリプト**: `rad.ObjDrwOpenGL()`をPyVistaに置き換え
3. **可視化が不要**: VTKファイル出力のみ（軽量）

---

**結論**: PyVistaは`rad.ObjDrwOpenGL()`の完璧な代替であり、より多機能で保守しやすい。OpenGLビューアの実装は不要。

---

**テスト環境:**
- PyVista: 0.46.3
- Radia: 4.32
- Python: 3.12
- OS: Windows 10 Build 20348

**作成日:** 2025-10-29
