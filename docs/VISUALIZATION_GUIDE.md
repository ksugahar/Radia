# Visualization Guide: Radia-NGSolve Viewer Selection and Workflows

## 全体フロー

```
┌─────────────────────────────────────────────────────────────────┐
│                    Radia-NGSolve Workflow                        │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  [1] CAD Modeling                                               │
│      Coreform Cubit / FreeCAD / STEP files                      │
│                          ↓                                       │
│  [2] Mesh Generation                                            │
│      Netgen / Cubit → .vol file (表面要素自動生成)               │
│      GMSH → .msh file                                           │
│                          ↓                                       │
│  [3] Visualization (用途別)                                      │
│      ├─ .msh確認: GMSH GUI (ネイティブ、最速)                    │
│      ├─ 形状確認: Netgen GUI (軽量、正確)                        │
│      ├─ 開発確認: PyVista (迅速)                                │
│      ├─ 論文図表: ParaView (高品質)                              │
│      └─ 統合確認: webgui (ブラウザ)                              │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

**重要**: Jupyterノートブック（.ipynb）はClaude Codeとの相性が悪いため、通常のPythonスクリプト（.py）を推奨します。PyVistaもwebguiも.pyスクリプトから完全に使用可能です。

### Jupyter Notebook vs Python Script

**Python Script（.py）推奨の理由**:

| 観点 | Python Script (.py) | Jupyter Notebook (.ipynb) |
|------|-------------------|-------------------------|
| **Claude Code編集** | ✅ 容易（通常のテキスト） | ❌ 困難（JSON構造） |
| **バージョン管理** | ✅ Git差分が明確 | ❌ JSON差分が読みにくい |
| **ファイルサイズ** | ✅ 小さい | ❌ 実行結果で肥大化 |
| **デバッグ** | ✅ 標準デバッガ使用可 | ❌ セル単位のみ |
| **自動化** | ✅ CI/CDで実行容易 | ❌ 追加設定必要 |
| **PyVista** | ✅ インタラクティブウィンドウ | ⚠️ インライン表示（制約あり） |
| **webgui** | ✅ ブラウザ自動起動 | ⚠️ ノートブック内表示 |

**結論**: 特別な理由がない限り、通常のPythonスクリプト（.py）を使用してください。

```python
# demo.py - 通常のPythonスクリプト
import pyvista as pv

# PyVista: インタラクティブウィンドウが開く（Jupyter不要）
grid = pv.read('field.vts')
grid.plot(scalars='B_magnitude', cmap='coolwarm')  # ← ウィンドウが開く
```

```python
# webgui_demo.py - 通常のPythonスクリプト
from ngsolve.webgui import Draw

# NGSolve webgui: ブラウザが自動的に開く（Jupyter不要）
Draw(B_gf, mesh, 'B_field')  # ← ブラウザタブが開く
```

---

## ビューワー比較表

| ビューワー | 品質 | 速度 | 用途 | インストール |
|----------|------|------|------|------------|
| **GMSH GUI** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | **.msh用ベストビューワー（最優先）** | GMSH同梱 |
| **PyVista** | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | Python統合（デフォルト） | `pip install pyvista` |
| **ParaView** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | **論文図表（最高品質）** | 単体アプリ |
| **NGSolve webgui** | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | Radiaフィールド統合 | `pip install ngsolve` |
| **Netgen GUI** | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | .vol形式 / 形状・メッシュ確認 | NGSolve同梱 |

### ビューワー選択フローチャート

```
目的は？
  │
  ├─ .mshメッシュ確認・デバッグ（開発中）
  │    → GMSH GUI ⭐⭐⭐⭐⭐ (最優先)
  │
  ├─ Python統合・自動化
  │    → PyVista ⭐⭐⭐⭐ (推奨)
  │
  ├─ 論文図表作成
  │    → ParaView ⭐⭐⭐⭐⭐ (最高品質)
  │
  ├─ Radiaフィールド探索
  │    → NGSolve webgui ⭐⭐⭐⭐ (Radia統合)
  │
  └─ 形状・メッシュ品質チェック
       → Netgen GUI ⭐⭐⭐ (軽量・正確)
```

### ユースケース別推奨

| ユースケース | 推奨ビューワー | 理由 |
|-------------|--------------|------|
| **.mshメッシュ確認** | **GMSH GUI** | ネイティブ、最速、変換不要 |
| **デバッグ** | **GMSH GUI** | 軽量、即座に確認 |
| **フィールドデータ(.msh内)** | **GMSH GUI** | View統合 |
| **開発中の確認** | PyVista | 迅速、スクリプト統合 |
| **パラメータスタディ** | PyVista | バッチ処理、自動化 |
| **メッシュ品質チェック** | Netgen GUI | 軽量、専用機能 |
| **フィールド分布探索** | NGSolve webgui | インタラクティブ |
| **論文投稿図** | ParaView | 高品質、ベクトル出力 |
| **プレゼン資料** | ParaView | 高解像度、美しい |
| **アニメーション動画** | ParaView | キーフレーム機能 |
| **Jupyter解析** | PyVista + webgui | ノートブック統合 |
| **CI/CD回帰テスト** | PyVista (headless) | 自動化可能 |

---

## Detailed Viewer Evaluations

### 1. GMSH GUI

**評価**: ⭐⭐⭐⭐⭐ (5/5) - `.msh`ファイル用として最高

#### なぜGMSH GUIが.mshに最適か

| 利点 | 詳細 |
|------|------|
| **ネイティブ統合** | `.msh`形式専用設計、変換不要 |
| **軽量・高速** | 起動・読込が非常に速い |
| **フィールドデータ対応** | `.msh`内の解データを直接表示 |
| **Python API連携** | `gmsh.view`, `gmsh.plugin`で自動化可能 |
| **ポストプロセス機能** | 等値面、ベクトル、ストリームライン内蔵 |
| **追加インストール不要** | GMSH同梱、依存関係なし |

#### 他のビューワーとの違い

```
GMSH GUI:     .msh → 直接表示 ✅
PyVista:      .msh → NGSolve → 変換 → 表示 ⚠️
ParaView:     .msh → 読込 → 表示 ✅（が、GMSH固有機能なし）
```

#### 使用例

```bash
# 最も簡単で高速
gmsh geometry.msh

# または Python API経由
python
>>> import gmsh
>>> gmsh.initialize()
>>> gmsh.open('geometry.msh')
>>> gmsh.fltk.run()
```

#### フィールドデータ統合

GMSH `.msh`ファイルにはメッシュだけでなく、**フィールドデータ（View）**も含められます：

```python
import gmsh

gmsh.initialize()
gmsh.open('geometry.msh')

# フィールドデータ追加（例: 磁束密度）
view_tag = gmsh.view.add("B_field")
gmsh.view.addListData(view_tag, "ST", num_elements, data_list)

# 保存（メッシュ+フィールドデータ）
gmsh.write('geometry_with_field.msh')

# GMSH GUIで開く → フィールドが自動表示
gmsh.fltk.run()  # GUIを起動
```

**GMSH GUIでの操作**:
- Tools → Visibility: フィールド表示ON/OFF
- Tools → Options → View: カラーマップ、等値面設定
- Plugins: Streamlines, Isosurface, Cutなど

#### Python API経由の自動化

**スクリプトでの可視化制御**:

```python
import gmsh
import numpy as np

gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)

# メッシュ読込
gmsh.open('coil_surface.msh')

# View追加（例: 電流密度）
view = gmsh.view.add("Current Density")

# データ設定（ST = Scalar Triangle）
# 各三角形要素に対してスカラー値を設定
num_triangles = 100  # 例
data = []
for i in range(num_triangles):
    # Triangle vertices (x1,y1,z1, x2,y2,z2, x3,y3,z3)
    # + scalar value
    data.extend([0, 0, 0,  1, 0, 0,  0, 1, 0,  1.5e6])  # Example

gmsh.view.addListData(view, "ST", num_triangles, data)

# カラーマップ設定
gmsh.view.option.setNumber(view, "ColormapNumber", 2)  # Jet colormap

# GUI起動
gmsh.fltk.run()
gmsh.finalize()
```

#### ポストプロセス機能

**GMSH GUI内蔵機能**:

| 機能 | 説明 | メニュー |
|------|------|---------|
| **Isosurface** | 等値面表示 | Plugins → Isosurface |
| **Streamlines** | 流線表示 | Plugins → Streamlines |
| **CutPlane** | 切断面表示 | Plugins → CutPlane |
| **Skin** | 外側表面のみ表示 | Plugins → Skin |
| **Smooth** | データ平滑化 | Plugins → Smooth |

**使用例**:
```
1. GMSH GUIで.mshを開く
2. Plugins → Isosurface
   - View: B_field を選択
   - Value: 0.5 (等値面の値)
   - Run
3. 新しいViewが生成される（等値面のみ）
```

#### パフォーマンス比較

**起動時間比較** (Windows):

| ビューワー | 起動時間 | メモリ使用量 |
|----------|---------|------------|
| **GMSH GUI** | **<1秒** | **~50MB** |
| ParaView | ~5秒 | ~200MB |
| PyVista | ~2秒 | ~100MB |

**大規模メッシュ読込** (100万要素):

| ビューワー | 読込時間 |
|----------|---------|
| **GMSH GUI** | **5秒** |
| ParaView | 8秒 |
| PyVista | 10秒（変換含む） |

#### キーボードショートカット

| キー | 機能 |
|------|------|
| `0` | メッシュ表示ON/OFF |
| `1-9` | View 1-9の表示切替 |
| `Shift+a` | 軸表示ON/OFF |
| `e` | 要素エッジ表示ON/OFF |
| `v` | Viewパネル表示 |
| `t` | Toolsパネル表示 |

#### マウス操作

| 操作 | 機能 |
|------|------|
| 左ドラッグ | 回転 |
| 中ドラッグ | 平行移動 |
| ホイール | ズーム |
| ダブルクリック | オブジェクト選択 |

#### スクリプト例: メッシュ品質確認

```python
import gmsh

gmsh.initialize()
gmsh.open('coil_surface.msh')

# メッシュ統計表示
gmsh.plugin.setNumber("MeshQuality", "Measure", 1)  # 1=SICN
gmsh.plugin.run("MeshQuality")

# GUI起動（品質がカラーマップで表示される）
gmsh.fltk.run()
gmsh.finalize()
```

#### 実例: メッシュ生成 → 即座に可視化

```python
import gmsh
import numpy as np

gmsh.initialize()
gmsh.model.add("coil")

# コイル形状生成（前述の方法）
# ... (geometry definition)

# 表面メッシュ生成
gmsh.model.mesh.generate(2)

# ファイル保存せずに直接GUI表示
gmsh.fltk.run()  # ← これだけ！

gmsh.finalize()
```

**利点**: ファイルI/O不要、生成と同時に確認

#### 実例: フィールドデータ可視化

```python
import gmsh
import numpy as np

gmsh.initialize()
gmsh.open('coil_surface.msh')

# 電流密度データ生成（例）
elements_2d = gmsh.model.mesh.getElements(2)
num_triangles = len(elements_2d[1][0])

# 三角形毎に電流密度を計算（例: 1.5e6 A/m^2）
view = gmsh.view.add("Current Density [A/m^2]")
data_list = []

for i in range(num_triangles):
    # Get triangle vertices
    elem_nodes = gmsh.model.mesh.getElement(elements_2d[0][0],
                                             elements_2d[1][0][i])
    coords = []
    for node in elem_nodes[1]:
        coord = gmsh.model.mesh.getNode(node)[0]
        coords.extend(coord)

    # Add scalar value (current density)
    coords.append(1.5e6)
    data_list.extend(coords)

gmsh.view.addListData(view, "ST", num_triangles, data_list)

# カラーマップ設定
gmsh.view.option.setNumber(view, "ColormapNumber", 2)  # Jet
gmsh.view.option.setNumber(view, "RangeType", 2)  # Custom range
gmsh.view.option.setNumber(view, "CustomMin", 0)
gmsh.view.option.setNumber(view, "CustomMax", 2e6)

# GUI起動
gmsh.fltk.run()
gmsh.finalize()
```

---

### 2. PyVista

**評価**: ⭐⭐⭐⭐ (4/5) - 開発デフォルト推奨

**長所**:
- ✅ Pythonネイティブ - スクリプト統合が容易
- ✅ Jupyter Notebook/Lab完全対応
- ✅ 迅速な可視化（開発イテレーション高速）
- ✅ VTK全機能へのPythonicアクセス
- ✅ アニメーション・GIF出力が簡単
- ✅ インタラクティブウィジェット（スライダー、チェックボックス）
- ✅ ヘッドレス実行可能（CI/CDパイプライン）
- ✅ 高速レンダリング
- ✅ 十分な表示品質

**短所**:
- ❌ 論文品質の微調整はParaViewより劣る
- ❌ ベクトルグラフィックス出力は限定的
- ❌ 大規模データ（>10M cells）では遅い

**Radia-NGSolve統合での使用例**:

```python
import pyvista as pv
import radia as rad
from radia_ngsolve import RadiaField

rad.FldUnits('m')

# RadiaフィールドをVTS出力
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])
rad.FldVTS(magnet, 'field.vts',
           [-0.1, 0.1], [-0.1, 0.1], [0.02, 0.15],
           41, 41, 27, 1, 0, 1.0)

# PyVistaで可視化
grid = pv.read('field.vts')
plotter = pv.Plotter()
plotter.add_mesh(grid, scalars='B_magnitude', cmap='coolwarm', opacity=0.8)
plotter.add_arrows(grid.points, grid['B_field'], mag=0.01, color='black')
plotter.show()

# Jupyter統合
grid.plot(scalars='B_magnitude', cmap='coolwarm', jupyter_backend='static')
```

**.mshファイル読込（NGSolve経由）**:

```python
from ngsolve import Mesh
import pyvista as pv
import numpy as np

# .mshファイル読込（NGSolve経由）
mesh = Mesh('geometry.msh')

# PyVista形式に変換
points = []
cells = []
for el in mesh.Elements3D():
    vertices = [mesh.vertices[v.nr].point for v in el.vertices]
    points.extend(vertices)
    cell_indices = list(range(len(cells)*4, len(cells)*4 + 4))
    cells.append([4] + cell_indices)  # [n_points, v0, v1, v2, v3]

points_array = np.array(points)
cells_array = np.hstack(cells)

# PyVistaメッシュ作成
grid = pv.UnstructuredGrid(cells_array, np.array([10]*len(cells)), points_array)

# 表示
plotter = pv.Plotter()
plotter.add_mesh(grid, show_edges=True, color='lightblue')
plotter.show()
```

**スライス表示**:

```python
# Z=0平面でスライス
slice_z = grid.slice(normal='z', origin=[0, 0, 0])
slice_z.plot(show_edges=True)
```

**ベストプラクティス**:
- 開発中の可視化確認
- パラメータスタディの自動化
- バッチ処理（複数ケースの比較）
- CI/CDでの回帰テスト可視化

---

### 3. ParaView

**評価**: ⭐⭐⭐⭐⭐ (5/5) - 論文品質最高

**長所**:
- ✅ 最高品質のレンダリング
- ✅ ベクトルグラフィックス出力（SVG, PDF）
- ✅ 高解像度ラスタ画像（300+ DPI）
- ✅ 複雑なフィルタチェーン（Glyph, Contour, StreamTracer）
- ✅ アニメーション・キーフレーム
- ✅ 大規模データ対応（分散並列）
- ✅ カメラ・照明の完全制御

**短所**:
- ❌ GUI操作が必要（スクリプト化は可能だが複雑）
- ❌ 迅速なイテレーションには不向き
- ❌ Jupyter統合は限定的（pvpythonは別プロセス）

**使用方法**:

```bash
# ParaViewで開く
paraview geometry.msh
```

**操作手順**:

1. **File → Open** → `geometry.msh`
2. **Apply** をクリック
3. フィルター適用:
   - **Filters → Slice**: 断面表示
   - **Filters → Clip**: 切断表示
   - **Filters → Glyph**: ベクトル表示
4. **File → Save Screenshot**: 高解像度エクスポート

**Radia-NGSolve統合での使用例**:

```bash
# 1. Radiaフィールド → VTS出力（Python）
python generate_field_vts.py

# 2. ParaViewで開く
paraview field.vts

# 3. ParaView GUI操作:
#    - Glyph filter: ベクトル矢印
#    - Contour: 等値面
#    - Slice: 断面
#    - Camera: アングル調整
#    - Lighting: 照明設定

# 4. 高解像度エクスポート:
#    File > Save Screenshot
#    - Resolution: 3000x2000 (300 DPI at 10x6.67 inch)
#    - Format: PNG (raster) or SVG (vector)
```

**ParaViewスクリプト自動化** (pvpython):

```python
# publication_figure.py
from paraview.simple import *

# Load VTS
reader = XMLStructuredGridReader(FileName=['field.vts'])

# Add glyph
glyph = Glyph(Input=reader, GlyphType='Arrow')
glyph.ScaleFactor = 0.01
glyph.GlyphMode = 'Every Nth Point'
glyph.Stride = 2

# Render
Show(glyph)
view = GetActiveView()
view.CameraPosition = [0.3, 0.2, 0.5]
view.CameraFocalPoint = [0, 0, 0]

# Save high-res image
SaveScreenshot('figure.png', view, ImageResolution=[3000, 2000])
```

**ベストプラクティス**:
- 論文投稿用図表
- プレゼンテーション資料
- 高解像度ポスター
- アニメーション動画（MP4, AVI）

---

### 4. NGSolve webgui

**評価**: ⭐⭐⭐⭐ (4/5) - インタラクティブ探索

**長所**:
- ✅ NGSolve完全統合
- ✅ メッシュ + フィールド同時表示
- ✅ WebGL - ブラウザで動作
- ✅ Jupyter統合（同じノートブック内）
- ✅ Radia CoefficientFunctionと連携
- ✅ リアルタイム更新（パラメータ変更）
- ✅ .pyスクリプトから使用可能
- ✅ OCC形状を正確に表示（近似なし）

**短所**:
- ❌ VTSファイル直接読込不可（NGSolve GridFunctionのみ）
- ❌ 高度なフィルタ機能なし
- ❌ 論文品質の微調整困難
- ❌ エクスポート形式が限定的

**Radia-NGSolve統合での使用例**:

```python
from ngsolve import *
from ngsolve.webgui import Draw
from radia_ngsolve import RadiaField
import radia as rad

rad.FldUnits('m')

# Cubitメッシュ → NGSolve
mesh = Mesh('model.msh')

# Radia磁石
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])

# Radia CoefficientFunction
B_cf = RadiaField(magnet, 'b')

# GridFunctionに投影
fes = HDiv(mesh, order=2)
B_gf = GridFunction(fes)
B_gf.Set(B_cf)

# インタラクティブ表示（Jupyter内）
Draw(B_gf, mesh, name='B_field', vectors={'grid_size': 10})

# メッシュも同時表示
Draw(mesh)
```

**OCC形状の正確な表示**:

```python
from ngsolve.webgui import Draw
from netgen.occ import Box, Pnt

# OCC形状（正確）
occ_magnet = Box(Pnt(-0.02, -0.02, -0.01), Pnt(0.02, 0.02, 0.01))
Draw(occ_magnet, name='Magnet')  # 形状は完全に正確

# フィールド（GridFunction）
Draw(B_gf, mesh, 'B_field')  # ブラウザでインタラクティブ
```

**ベストプラクティス**:
- フィールド分布のクイック確認
- メッシュ品質チェック
- パラメータ最適化時のリアルタイムフィードバック
- 教育・デモンストレーション

---

### 5. Netgen GUI

**評価**: ⭐⭐⭐ (3/5) - .vol形式専用、形状・メッシュ確認に最適

**長所**:
- ✅ Netgen/NGSolveネイティブ（Tcl/Tk GUI）
- ✅ 形状（OCC）を**正確に**表示（近似なし）
- ✅ メッシュ品質可視化（アスペクト比、角度など）
- ✅ 軽量・高速起動（ブラウザ不要）
- ✅ STL/STEP/IGES読込
- ✅ 統合ワークフロー（形状確認→メッシュ生成→品質チェック）
- ✅ 通常のPythonスクリプト（.py）から使用可能
- ✅ 表面要素で外形を正確に表示

**短所**:
- ❌ `.msh`ファイルは直接読込不可
- ❌ フィールドデータ表示は限定的
- ❌ 論文品質レンダリング不可
- ⚠️ 古いGUI（Tcl/Tk）- ただし軽量で安定
- ⚠️ 表面要素が必須（通常は問題なし）

**使用例**:

```python
from netgen.occ import OCCGeometry, Box, Pnt
from netgen.gui import StartGUI

# Radia磁石をOCC形状に変換
occ_magnet = Box(Pnt(-0.02, -0.02, -0.01), Pnt(0.02, 0.02, 0.01))
geo = OCCGeometry(occ_magnet)

# Netgen GUIで確認（ネイティブウィンドウが開く）
StartGUI()
geo.Draw()  # 形状を正確に表示

# メッシュ生成
mesh = geo.GenerateMesh(maxh=0.005)
mesh.Draw()  # メッシュ品質確認
```

**.msh → .vol変換（NGSolve経由）**:

```python
from ngsolve import Mesh

mesh = Mesh('geometry.msh')
mesh.ngmesh.Save('geometry.vol')  # .vol形式で保存

# Netgen GUIで開く
# python utils/netgen_vol_viewer.py geometry.vol
```

**webguiとの使い分け**:

| 用途 | netgen.gui | ngsolve.webgui |
|------|-----------|---------------|
| 形状確認 | **✅ 推奨** | ⚠️ ブラウザオーバーヘッド |
| メッシュ品質 | **✅ 推奨** | ❌ 限定的 |
| フィールド可視化 | ❌ 不可 | **✅ 推奨** |
| 軽量・高速 | **✅ ネイティブGUI** | ⚠️ ブラウザ必要 |

**Windowsファイル関連付け**:

.volファイルをダブルクリックでNetgen GUIで開く設定：

```bash
# 自動設定（管理者権限）
cd S:\Radia\01_GitHub\utils
setup_vol_file_association.bat
```

詳細: [VOL_FILE_ASSOCIATION.md](file://S:/Radia/01_GitHub/utils/VOL_FILE_ASSOCIATION.md)

**ベストプラクティス**:
- メッシュ生成前の形状確認（**webguiより正確**）
- メッシュ品質チェック（アスペクト比、角度）
- 境界条件ラベルの確認
- CADインポート後の形状検証
- Cubit → Netgen ワークフローの確認

---

## Geometry Accuracy: VTS vs OCC

### 問題: VTS Exportによる形状近似

#### rad.FldVTS()の挙動

```python
# Radia解析オブジェクト（完全に正確）
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])

# VTSエクスポート（格子点でのフィールド値計算）
rad.FldVTS(magnet, 'field.vts',
           [-0.1, 0.1], [-0.1, 0.1], [0.02, 0.15],
           41, 41, 27, 1, 0, 1.0)
```

**VTSファイルに含まれる情報**:
- ✅ **フィールド値**: 格子点での B, H, A, Phi（正確）
- ❌ **形状情報**: 失われる（格子の外形で近似）

**結果**:
- ParaViewで`field.vts`を開くと、フィールド分布は正確だが、磁石の形状は格子で近似される
- 長方形磁石が格子状の境界で表示される（完全な長方形ではない）

### ビューワー別の形状精度

| 方法 | 形状精度 | フィールド精度 | 手間 | 品質 |
|------|---------|--------------|------|-----|
| **PyVista + VTS** | ❌ 近似 | ✅ 正確 | 低 | 良 |
| **ParaView + VTS** | ❌ 近似 | ✅ 正確 | 中 | 優 |
| **webgui + OCC** | ✅ 正確 | ✅ 正確 | 中 | 良 |
| **ParaView + STL/STEP + VTS** | ✅ 正確 | ✅ 正確 | 高 | 優 |

### 解決策: ParaView STL+VTS Overlay

形状とフィールドの両方を正確に表示するには：

```python
# 1. 形状をSTL/STEPでエクスポート（正確）
from netgen.occ import Box, Pnt, OCCGeometry

box = Box(Pnt(-0.02, -0.02, -0.01), Pnt(0.02, 0.02, 0.01))
geo = OCCGeometry(box)
mesh = geo.GenerateMesh(maxh=0.002)
mesh.Export('magnet_shape.stl', 'STL Format')

# 2. フィールドをVTSでエクスポート
rad.FldVTS(magnet, 'field.vts', ...)

# 3. ParaViewで両方をoverlayして高品質図表作成
# - magnet_shape.stl: 半透明で形状表示
# - field.vts: カラーマップでフィールド表示
```

### 結論

**ユーザーの懸念「paraview介すと必ずしも形状がそのままではない」は正しい**:
- VTS形式は構造格子なので、形状は格子で近似される
- 解析的な長方形や円柱が、格子の境界で表現される

**解決策**:
1. **形状不要**: PyVista/ParaView + VTS（現在の実装で十分）
2. **形状重要（開発）**: NGSolve webgui + OCC形状（`rad.ExportOCC()` TODO）
3. **形状重要（論文）**: ParaView + STL/STEP + VTS overlay（手動、または自動化TODO）

---

## Recommended Workflows by Use Case

### メッシュ生成パターン

#### パターンA: Netgen直接生成（推奨）

```python
from netgen.occ import Box, Pnt, OCCGeometry
from netgen.gui import StartGUI

# 形状作成
box = Box(Pnt(-0.05, -0.05, -0.05), Pnt(0.05, 0.05, 0.05))
geo = OCCGeometry(box)

# メッシュ生成（表面要素自動）
mesh = geo.GenerateMesh(maxh=0.01)

print(f"Volume elements:  {mesh.ne}")
print(f"Surface elements: {mesh.nse}")  # > 0 (自動生成)

# 保存
mesh.Save('magnet.vol')  # 表面要素を含む

# Netgen GUIで確認
StartGUI()
mesh.Draw()  # ← 問題なく表示される
```

**結果**: 表面要素が自動的に含まれる ✅

#### パターンB: Cubit → Netgen変換

```python
import cubit
import cubit_mesh_export
from ngsolve import Mesh

# Cubitでメッシュ生成
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("import step 'motor_rotor.step'")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")

# Sideset定義（これが表面要素になる）
cubit.cmd("sideset 1 surface all")
cubit.cmd("sideset 1 name 'boundary'")

# Netgenへ変換（表面要素保持）
ngmesh = cubit_mesh_export.export_netgen(cubit)
mesh = Mesh(ngmesh)

# 保存
mesh.ngmesh.Save('motor_rotor.vol')  # 表面要素を含む
```

**結果**: Cubit sidesetが表面要素に変換される ✅

メッシュ生成と表面要素の詳細については [MESH_GUIDE.md](MESH_GUIDE.md) を参照。

#### パターンC: GMSH → NGSolve

```python
from ngsolve import Mesh

# GMSHメッシュ読込（NGSolveが自動変換）
mesh = Mesh('geometry.msh')

# 表面要素は自動的に認識される
print(f"Surface elements: {mesh.ngmesh.nse}")  # > 0

# .volとして保存
mesh.ngmesh.Save('geometry.vol')
```

**結果**: GMSH境界要素がNetgen表面要素に変換される ✅

### 可視化ワークフロー（用途別）

#### フロー1: .mshメッシュ確認（最優先）

```
GMSH mesh generate → geometry.msh → **GMSH GUI**
                                        ↓
                                  即座に確認（変換不要）
```

#### フロー2: 形状・メッシュ確認

```
Netgen / Cubit → .vol → **Netgen GUI** (軽量、形状正確)
```

```python
from netgen.meshing import Mesh
from netgen.gui import StartGUI

mesh = Mesh()
mesh.Load('magnet.vol')

StartGUI()
mesh.Draw()
```

#### フロー3: 開発・デバッグ

```
Cubit → Netgen mesh → Radia solve → VTS export → **PyVista**
                                                     ↓
                                              迅速な確認・修正
```

```python
import radia as rad
import pyvista as pv

rad.FldUnits('m')

# Radia磁石
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])

# VTSエクスポート
rad.FldVTS(magnet, 'field.vts', ...)

# PyVista可視化（迅速）
grid = pv.read('field.vts')
grid.plot(scalars='B_magnitude', cmap='coolwarm')
```

#### フロー4: 論文投稿

```
Cubit → Netgen mesh → Radia solve → VTS export → **ParaView**
                                                     ↓
                                         高解像度PNG/SVG出力
```

形状とフィールドの両方が重要な場合：
```
Radia → ExportOCC() → STL/STEP ──┐
     ↓                            ├→ ParaView overlay
     └→ FldVTS() → VTS ───────────┘
                  ↓
            形状正確 + フィールド正確
```

#### フロー5: インタラクティブ探索

```
Cubit → Netgen mesh → Radia CF → NGSolve GridFunction → **webgui**
                                                            ↓
                                                  ブラウザで即座に確認
```

### 推奨フローチャート（総合）

```
メッシュ生成
    ↓
ファイル形式は？
    ├─ .msh → GMSH GUI（最速確認）
    └─ .vol → 用途で選択
                ├─ 形状確認 → Netgen GUI (最軽量)
                ├─ 開発確認 → PyVista (迅速)
                ├─ 論文図表 → ParaView (高品質)
                └─ 統合探索 → webgui (正確)
```

### よくある質問

#### Q: 表面要素がないと困る？

**A: 実用上は問題ありません。**

理由：
- Netgen生成メッシュ: 自動的に表面要素が含まれる
- Cubit変換メッシュ: sideset定義で表面要素生成
- NGSolveサンプル: すべて表面要素を含む

稀に体積要素のみの場合: ParaView/PyVistaで切断面表示

#### Q: Windowsでダブルクリックで開きたい

**A: 可能です。**

```bash
# 自動設定（管理者権限）
cd S:\Radia\01_GitHub\utils
setup_vol_file_association.bat
```

詳細: [VOL_FILE_ASSOCIATION.md](file://S:/Radia/01_GitHub/utils/VOL_FILE_ASSOCIATION.md)

#### Q: ParaViewで形状が近似される問題は？

**A: 2つの解決策があります。**

1. **形状確認にはNetgen GUI使用**（正確、軽量）
2. **論文図表ではSTL+VTS overlay**（正確、高品質）

詳細は上記「Geometry Accuracy: VTS vs OCC」セクションを参照。

---

## Implementation Status

### Phase 1: 基本可視化（完成済み）
- ✅ `rad.FldVTS()` - VTSエクスポート
- ✅ PyVistaによる基本プロット
- ✅ ParaViewでの手動可視化

### Phase 2: NGSolve統合強化（進行中）
- ✅ `RadiaField` CoefficientFunction
- ⏳ `radia_ngsolve` Python API改善
- ⏳ webgui連携スクリプト集

### Phase 3: 高度な可視化（TODO）
- ⏳ ParaView自動化スクリプト（pvpython）
- ⏳ PyVistaアニメーション生成
- ⏳ `rad.ExportOCC()` - Radia解析オブジェクト → OCC形状変換
  - `ObjRecMag` → `netgen.occ.Box`
  - `ObjCylMag` → `netgen.occ.Cylinder`
  - `ObjSphMag` → `netgen.occ.Sphere`
- ⏳ STL/STEPエクスポート自動化
- ⏳ ParaView overlay自動化スクリプト

### Phase 4: 統合ビューワー（将来）
- ⏳ カスタムPyVistaインターフェース
- ⏳ Jupyter Widgetによるパラメータ制御
- ⏳ WebアプリケーションUI（Dash/Streamlit）

---

## References

### インストールコマンド

```bash
# GMSH（メッシュ生成+ビューワー）
pip install gmsh  # PythonライブラリとGUI両方含む

# PyVista（開発用推奨）
pip install pyvista

# ParaView（論文用）
# https://www.paraview.org/download/ からダウンロード

# NGSolve webgui（フィールド探索）
pip install ngsolve  # webgui同梱
```

### 参考実装

**EMPY_Field** (`S:\NGSolve\EMPY\EMPY_Field`):
- Radia解析オブジェクトのOCC変換実装例
- 正確な形状をOCCで表現

### 関連ドキュメント

- [MESH_GUIDE.md](MESH_GUIDE.md) - メッシュ生成と表面要素の詳細
- [MESH_GUIDE.md](MESH_GUIDE.md) - GMSHメッシュ生成ワークフロー
- [VOL_FILE_ASSOCIATION.md](file://S:/Radia/01_GitHub/utils/VOL_FILE_ASSOCIATION.md) - .volファイル関連付け設定

---

**作成日**: 2026-02-12
**更新日**: 2026-02-22
**対象**: Radia-NGSolve統合フレームワークの可視化ガイド（ビューワー選択、形状精度、ワークフロー）
