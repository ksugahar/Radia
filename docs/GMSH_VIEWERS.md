# GMSH Mesh Viewers Comparison

## Overview

GMSH `.msh` ファイルは複数のビューワーで表示可能です。用途に応じて最適なビューワーを選択してください。

---

## ビューワー比較

| ビューワー | 品質 | 速度 | 用途 | インストール |
|----------|------|------|------|------------|
| **GMSH GUI** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐⭐⭐ | **.msh用ベストビューワー（最優先）** | GMSH同梱 |
| **PyVista** | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | Python統合 | `pip install pyvista` |
| **ParaView** | ⭐⭐⭐⭐⭐ | ⭐⭐⭐ | **論文図表（最高品質）** | 単体アプリ |
| **NGSolve webgui** | ⭐⭐⭐⭐ | ⭐⭐⭐⭐ | Radiaフィールド統合 | `pip install ngsolve` |
| Netgen GUI | ⭐⭐⭐ | ⭐⭐⭐⭐⭐ | .vol形式のみ | NGSolve同梱 |

---

## 1. GMSH GUI（.msh用ベストビューワー）

**評価**: ⭐⭐⭐⭐⭐ (5/5) - **.mshファイルに最適**

### 利点
- ✅ **ネイティブ形式**（変換不要、最速）
- ✅ **軽量・高速**（起動<1秒、メモリ~50MB）
- ✅ **フィールドデータ統合**（`.msh`内のViewを直接表示）
- ✅ **ポストプロセス機能**（等値面、流線、切断面内蔵）
- ✅ **Python API**（スクリプト自動化可能）
- ✅ **追加インストール不要**（GMSH同梱）

### 使用例

```bash
# 最も簡単で高速
gmsh geometry.msh
```

### Python API経由

```python
import gmsh

gmsh.initialize()
gmsh.open('coil_surface.msh')

# フィールドデータ追加（例: 電流密度）
view = gmsh.view.add("Current Density")
# ... データ設定

# GUI起動
gmsh.fltk.run()
gmsh.finalize()
```

### ポストプロセス機能

**GMSH GUI内蔵機能**:
- Plugins → Isosurface（等値面）
- Plugins → Streamlines（流線）
- Plugins → CutPlane（切断面）
- Plugins → Skin（外側表面のみ）

**詳細**: [GMSH_NATIVE_VIEWER.md](GMSH_NATIVE_VIEWER.md)

---

## 2. PyVista（Python統合）

**評価**: ⭐⭐⭐⭐ (4/5)

### 利点
- ✅ Python統合（Jupyter不要の.pyスクリプト）
- ✅ 高速レンダリング
- ✅ インタラクティブ操作
- ✅ スクリプト自動化
- ✅ 十分な表示品質

### 使用例

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

### スライス表示

```python
# Z=0平面でスライス
slice_z = grid.slice(normal='z', origin=[0, 0, 0])
slice_z.plot(show_edges=True)
```

---

## 3. ParaView（論文品質）

**評価**: ⭐⭐⭐⭐⭐ (5/5) - 最高品質

### 利点
- ✅ 最高品質レンダリング
- ✅ ベクトルグラフィックス出力（SVG, PDF）
- ✅ 高度なフィルター（Streamlines, Glyphs）
- ✅ アニメーション

### 使用方法

```bash
# ParaViewで開く
paraview geometry.msh
```

### 操作手順

1. **File → Open** → `geometry.msh`
2. **Apply** をクリック
3. フィルター適用:
   - **Filters → Slice**: 断面表示
   - **Filters → Clip**: 切断表示
   - **Filters → Glyph**: ベクトル表示
4. **File → Save Screenshot**: 高解像度エクスポート

### 論文図表作成ワークフロー

```bash
# 1. GMSHでメッシュ生成
python generate_mesh.py  # → geometry.msh

# 2. ParaViewで開く
paraview geometry.msh

# 3. フィルター適用、レンダリング調整

# 4. 高解像度エクスポート
# File → Save Screenshot
# - Resolution: 300 DPI
# - Format: PNG (高解像度) or SVG (ベクトル)
```

---

## 4. NGSolve webgui（Radiaフィールド統合）

**評価**: ⭐⭐⭐⭐ (4/5)

### 利点
- ✅ ブラウザベース（インストール不要）
- ✅ インタラクティブ操作
- ✅ .pyスクリプトから使用可能
- ✅ フィールドデータ統合

### 使用例

```python
from ngsolve import Mesh
from ngsolve.webgui import Draw

# .msh読込
mesh = Mesh('geometry.msh')

# ブラウザで表示
Draw(mesh)
```

### フィールドデータ表示

```python
from ngsolve import *
from ngsolve.webgui import Draw
from radia_ngsolve import RadiaField
import radia as rad

rad.FldUnits('m')

# Radiaオブジェクト
magnet = rad.ObjRecMag([0, 0, 0], [0.04, 0.04, 0.02], [0, 0, 954930])

# メッシュ読込
mesh = Mesh('geometry.msh')

# RadiaFieldをGridFunctionに投影
B_cf = RadiaField(magnet, 'b')
B_gf = GridFunction(HDiv(mesh, order=2))
B_gf.Set(B_cf)

# ブラウザで表示
Draw(B_gf, mesh, name='B_field')
```

---

## 5. Netgen GUI（.vol形式専用）

**評価**: N/A - GMSH .mshには非対応

### 注意
- ❌ `.msh`ファイルは直接読込不可
- ✅ `.vol`ファイル（Netgen形式）のみ対応

### 代替手段

```python
# .msh → .vol変換（NGSolve経由）
from ngsolve import Mesh

mesh = Mesh('geometry.msh')
mesh.ngmesh.Save('geometry.vol')  # .vol形式で保存

# Netgen GUIで開く
# python utils/netgen_vol_viewer.py geometry.vol
```

---

## 推奨ワークフロー

### メッシュ確認・デバッグ（最優先）

```bash
# GMSH GUI最優先（最速、ネイティブ）
gmsh geometry.msh
```

**理由**:
- ✅ 最速（変換不要、起動<1秒）
- ✅ 軽量（メモリ~50MB）
- ✅ ネイティブ形式
- ✅ フィールドデータ統合

### Python統合（スクリプト自動化）

```python
# PyVista（Pythonスクリプト内で完結）
from ngsolve import Mesh
import pyvista as pv

mesh = Mesh('geometry.msh')
# ... PyVista変換・表示
```

**理由**:
- Pythonスクリプトで完結
- バッチ処理可能
- カスタム解析統合

### 論文図表作成

```bash
# ParaView（最高品質レンダリング）
paraview geometry.msh
```

**理由**:
- 最高品質レンダリング
- 高解像度エクスポート
- ベクトルグラフィックス対応

### Radiaフィールドデータ探索

```python
# NGSolve webgui（Radia統合）
from ngsolve import Mesh
from ngsolve.webgui import Draw
from radia_ngsolve import RadiaField

mesh = Mesh('geometry.msh')
B_gf.Set(RadiaField(magnet, 'b'))
Draw(B_gf, mesh)
```

**理由**:
- Radiaフィールドと統合
- ブラウザベース
- インタラクティブ

---

## ビューワー選択フローチャート

```
目的は？
  │
  ├─ メッシュ確認・デバッグ（開発中）
  │    → GMSH GUI ⭐⭐⭐⭐⭐ (最優先)
  │
  ├─ Python統合・自動化
  │    → PyVista ⭐⭐⭐⭐ (推奨)
  │
  ├─ 論文図表作成
  │    → ParaView ⭐⭐⭐⭐⭐ (最高品質)
  │
  └─ Radiaフィールド探索
       → NGSolve webgui ⭐⭐⭐⭐ (Radia統合)
```

---

## インストールコマンド

```bash
# PyVista（開発用推奨）
pip install pyvista

# ParaView（論文用）
# https://www.paraview.org/download/ からダウンロード

# NGSolve webgui（フィールド探索）
pip install ngsolve  # webgui同梱

# GMSH（メッシュ生成+ビューワー）
pip install gmsh  # PythonライブラリとGUI両方含む
```

---

## まとめ

### 最優先推奨: GMSH GUI（.msh用ベストビューワー）

**理由**:
- ✅ **ネイティブ形式**（変換不要）
- ✅ **最速**（起動<1秒、読込最速）
- ✅ **軽量**（メモリ~50MB）
- ✅ **フィールドデータ統合**
- ✅ **ポストプロセス機能内蔵**
- ✅ **Python API自動化可能**

**使用シーン**: メッシュ確認、デバッグ、クイック可視化、フィールドデータ確認

### Python統合: PyVista

**理由**:
- Pythonスクリプトで完結
- バッチ処理・自動化容易
- カスタム解析統合

**使用シーン**: スクリプト自動化、Jupyter統合、カスタムポストプロセス

### 論文図表: ParaView

**理由**:
- 最高品質レンダリング
- 高解像度出力
- ベクトルグラフィックス

**使用シーン**: 学術論文、プレゼンテーション資料

### Radiaフィールド: NGSolve webgui

**理由**:
- Radiaフィールドと統合
- インタラクティブ探索
- ブラウザベース

**使用シーン**: Radiaフィールドの可視化、パラメータスタディ

---

**作成日**: 2026-02-12
**対象**: GMSH .mshファイルの可視化
**関連ドキュメント**:
- [GMSH_WORKFLOW.md](GMSH_WORKFLOW.md)
- [VIEWER_COMPARISON.md](VIEWER_COMPARISON.md)
- [VISUALIZATION_WORKFLOW.md](VISUALIZATION_WORKFLOW.md)
