# Netgen GUI Limitations and Workarounds

## 問題: 体積要素のみのメッシュは表示されない

### Netgen GUIの表示仕様

Netgen GUIは**表面要素（Surface elements）** を主に表示します：
- ✅ **表面要素あり**: メッシュの外側の境界が表示される
- ❌ **体積要素のみ**: 何も表示されないか、切断面（clipping plane）が必要

### 確認方法

```python
# check_vol_surface_elements.py を使用
python utils/check_vol_surface_elements.py mesh.vol
```

または：

```python
from netgen.meshing import Mesh

mesh = Mesh()
mesh.Load('mesh.vol')

print(f"Volume elements:  {mesh.ne}")   # 体積要素
print(f"Surface elements: {mesh.nse}")  # 表面要素

if mesh.nse == 0:
    print("⚠️  Warning: No surface elements")
    print("   Netgen GUI may not display this mesh")
```

---

## NGSolve標準メッシュの状況

NGSolveの`share/ngsolve/`メッシュは**すべて表面要素を含む**：

| ファイル | 体積要素 | 表面要素 | Netgen GUI |
|---------|---------|---------|-----------|
| cube.vol | 756 | 338 (Triangle) | ✅ OK |
| coil.vol | 1709 | あり | ✅ OK |
| coilshield.vol | 1798 | 376 (Tri+Quad) | ✅ OK |
| chip.vol | 0 | あり (Surface-only) | ✅ OK |
| beam.vol | 31 | あり | ✅ OK |
| shaft.vol | 1622 | あり | ✅ OK |

**結論**: NGSolveの標準.volファイルはNetgen GUIで表示可能。

---

## 体積要素のみのメッシュの対処法

### 方法1: ParaViewで切断面表示（推奨）

```python
# 1. VTSにエクスポート
import radia as rad
rad.FldUnits('m')
rad.FldVTS(magnet, 'field.vts', ...)

# 2. ParaViewで開く
paraview field.vts

# 3. Filters > Slice
#    - Origin: [0, 0, 0.05]
#    - Normal: [0, 0, 1]
#    - Apply

# 4. Filters > Clip
#    - Clip Type: Plane
#    - Normal: [0, 0, 1]
#    - Apply
```

### 方法2: PyVistaで切断面表示

```python
import pyvista as pv

# Load VTS
grid = pv.read('field.vts')

# Create slice at z=0.05m
slice_z = grid.slice(normal='z', origin=[0, 0, 0.05])
slice_z.plot(scalars='B_magnitude', cmap='coolwarm')

# Or clip half of the domain
clipped = grid.clip(normal='z', origin=[0, 0, 0])
clipped.plot(scalars='B_magnitude', cmap='viridis')
```

### 方法3: NGSolveで表面要素を生成

```python
from ngsolve import *
from netgen.occ import Box, Pnt, OCCGeometry

# Create geometry
box = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))
geo = OCCGeometry(box)

# Generate mesh with surface elements
mesh = Mesh(geo.GenerateMesh(maxh=0.1))

# Export to .vol (includes surface elements)
mesh.ngmesh.Save('mesh_with_surface.vol')

# Now can open in Netgen GUI
from netgen.gui import StartGUI
StartGUI()
mesh.ngmesh.Draw()
```

---

## 推奨ワークフロー

### ケース1: 表面要素あり（標準）

```
.vol file (with surface) → Netgen GUI → 3D surface display ✅
```

**適用**: NGSolveサンプルメッシュ、Netgen生成メッシュ

### ケース2: 体積要素のみ

```
.vol file (volume only) → VTS export → ParaView/PyVista → Slice/Clip ✅
```

**適用**: カスタムメッシュ、内部構造の可視化

---

## ビューワー選択ガイド

| メッシュタイプ | Netgen GUI | ParaView | PyVista | webgui |
|--------------|-----------|----------|---------|--------|
| **表面要素あり** | **✅ 推奨** | ⚠️ オーバースペック | ⚠️ オーバースペック | ❌ 形状のみ |
| **体積要素のみ** | ❌ 表示困難 | **✅ 推奨** | **✅ 推奨** | ❌ GridFunction必要 |
| **形状確認** | **✅ 最適** | ⚠️ メッシュ化必要 | ⚠️ メッシュ化必要 | ✅ OCC直接 |
| **フィールド可視化** | ❌ 不可 | **✅ 最高品質** | ✅ 高速 | ✅ インタラクティブ |

---

## まとめ

### Netgen GUIを使うべき場合

- ✅ メッシュに表面要素が含まれている
- ✅ 形状・メッシュ品質の確認が目的
- ✅ 軽量・高速な表示が必要
- ✅ 統合ワークフロー（形状→メッシュ→確認）

### ParaView/PyVistaを使うべき場合

- ✅ 体積要素のみのメッシュ
- ✅ 内部構造の可視化（切断面）
- ✅ フィールドデータの可視化
- ✅ 論文品質の図表が必要

### 確認コマンド

```bash
# 表面要素の有無を確認
python utils/check_vol_surface_elements.py mesh.vol
```

---

**作成日**: 2026-02-12
**適用対象**: Netgen GUI, NGSolve .vol files
