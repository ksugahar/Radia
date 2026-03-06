# Radia-NGSolve Visualization Workflow (Complete Guide)

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
│                          ↓                                       │
│  [3] Visualization (用途別)                                      │
│      ├─ 形状確認: Netgen GUI (軽量、正確)                        │
│      ├─ 開発確認: PyVista (迅速)                                │
│      ├─ 論文図表: ParaView (高品質)                              │
│      └─ 統合確認: webgui (ブラウザ)                              │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

---

## フェーズ1: メッシュ生成

### パターンA: Netgen直接生成（推奨）

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

---

### パターンB: Cubit → Netgen変換

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

---

### パターンC: GMSH → NGSolve

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

---

## フェーズ2: 可視化（用途別）

### ケース1: 形状・メッシュ確認（開発初期）

**推奨**: Netgen GUI

```python
from netgen.meshing import Mesh
from netgen.gui import StartGUI

mesh = Mesh()
mesh.Load('magnet.vol')

StartGUI()
mesh.Draw()
```

**利点**:
- ✅ 軽量・高速
- ✅ 形状が正確（OCC）
- ✅ メッシュ品質確認可能
- ✅ 表面要素で外形表示

**制限**:
- ❌ フィールドデータは表示不可
- ⚠️ 表面要素が必須（通常は問題なし）

---

### ケース2: フィールド分布確認（開発中）

**推奨**: PyVista

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

**利点**:
- ✅ フィールド可視化が迅速
- ✅ スクリプト自動化容易
- ✅ 切断面表示が簡単

**制限**:
- ⚠️ 形状は格子で近似（VTS制限）
- ⚠️ 論文品質はParaViewに劣る

---

### ケース3: 論文投稿図表

**推奨**: ParaView + 形状overlay

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

**利点**:
- ✅ 形状が正確（STL/STEP）
- ✅ フィールドも正確（VTS）
- ✅ 最高品質レンダリング
- ✅ ベクトルグラフィックス出力

**手間**:
- ⚠️ 手動操作が必要（または pvpython自動化）

---

### ケース4: インタラクティブ探索

**推奨**: NGSolve webgui + OCC

```python
from ngsolve.webgui import Draw
from netgen.occ import Box, Pnt

# OCC形状（正確）
occ_magnet = Box(Pnt(-0.02, -0.02, -0.01), Pnt(0.02, 0.02, 0.01))
Draw(occ_magnet, name='Magnet')  # 形状は完全に正確

# フィールド（GridFunction）
Draw(B_gf, mesh, 'B_field')  # ブラウザでインタラクティブ
```

**利点**:
- ✅ 形状が正確（OCC直接表示）
- ✅ フィールドも正確（GridFunction）
- ✅ ブラウザベース（インストール不要）
- ✅ .pyスクリプトから使用可能

**制限**:
- ⚠️ レンダリング品質はParaViewに劣る

---

## 推奨フローチャート

```
メッシュ生成
    ↓
表面要素あり？ ← 通常は自動的に「YES」
    ↓
  YES → 用途で選択
          ├─ 形状確認 → Netgen GUI (最軽量)
          ├─ 開発確認 → PyVista (迅速)
          ├─ 論文図表 → ParaView (高品質)
          └─ 統合探索 → webgui (正確)
    ↓
   NO（稀）→ ParaView/PyVista (切断面)
```

---

## よくある質問

### Q: 表面要素がないと困る？

**A: 実用上は問題ありません。**

理由：
- Netgen生成メッシュ: 自動的に表面要素が含まれる
- Cubit変換メッシュ: sideset定義で表面要素生成
- NGSolveサンプル: すべて表面要素を含む

稀に体積要素のみの場合: ParaView/PyVistaで切断面表示

---

### Q: Windowsでダブルクリックで開きたい

**A: 可能です。**

```bash
# 自動設定（管理者権限）
cd S:\Radia\01_GitHub\utils
setup_vol_file_association.bat
```

詳細: [VOL_FILE_ASSOCIATION.md](file://S:/Radia/01_GitHub/utils/VOL_FILE_ASSOCIATION.md)

---

### Q: ParaViewで形状が近似される問題は？

**A: 2つの解決策があります。**

1. **形状確認にはNetgen GUI使用**（正確、軽量）
2. **論文図表ではSTL+VTS overlay**（正確、高品質）

詳細: [GEOMETRY_ACCURACY_COMPARISON.md](file://S:/Radia/01_GitHub/docs/GEOMETRY_ACCURACY_COMPARISON.md)

---

## まとめ

### 結論

| 用途 | ビューワー | 形状精度 | フィールド | 手間 |
|------|----------|---------|-----------|------|
| **形状確認** | **Netgen GUI** | ✅ 正確 | ❌ なし | 低 |
| **開発確認** | **PyVista** | ⚠️ 近似 | ✅ あり | 低 |
| **論文図表** | **ParaView** | ✅ 正確* | ✅ あり | 中 |
| **統合探索** | **webgui** | ✅ 正確 | ✅ あり | 低 |

*ParaViewで形状正確 = STL/STEP overlay使用時

### 実用上の推奨

1. **メッシュ生成**: Netgen/Cubit（表面要素自動）
2. **形状確認**: Netgen GUI（軽量、正確）
3. **開発中**: PyVista（迅速、.pyスクリプト）
4. **論文**: ParaView（高品質、STL overlay）

**表面要素の制限は実用上問題になりません。**

---

**作成日**: 2026-02-12
**対象**: Radia-NGSolve統合ワークフロー
**関連ドキュメント**:
- [VIEWER_COMPARISON.md](file://S:/Radia/01_GitHub/docs/VIEWER_COMPARISON.md)
- [GEOMETRY_ACCURACY_COMPARISON.md](file://S:/Radia/01_GitHub/docs/GEOMETRY_ACCURACY_COMPARISON.md)
- [SURFACE_ELEMENTS_FAQ.md](file://S:/Radia/01_GitHub/docs/SURFACE_ELEMENTS_FAQ.md)
- [NETGEN_GUI_LIMITATIONS.md](file://S:/Radia/01_GitHub/docs/NETGEN_GUI_LIMITATIONS.md)
