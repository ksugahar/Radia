# Surface Elements FAQ

## Q1: 表面要素が必須なのは困る - どうすればいい？

**A: 実用上は問題ありません。** 理由：

### 自動生成される場合

以下のワークフローでは、表面要素は**自動的に含まれます**：

| ワークフロー | 表面要素 | 理由 |
|------------|---------|------|
| **Netgen直接生成** | ✅ 自動 | `geo.GenerateMesh()` が境界メッシュを自動生成 |
| **NGSolve Mesh()** | ✅ 自動 | STEP/OCC読込時に境界を認識 |
| **Cubit → export_netgen()** | ✅ 自動 | Cubit sidesetが境界要素に変換される |
| **GMSH → NGSolve** | ✅ 自動 | .mshファイルに境界要素が含まれる |

**つまり、通常のメッシュ生成では何も気にする必要がありません。**

---

## Q2: 体積要素のみのメッシュを使いたい場合は？

**A: 意図的に体積要素のみにする必要はありません。** しかし、どうしても必要な場合：

### 対処法1: ParaView/PyVistaを使う（推奨）

```python
import pyvista as pv

# 体積メッシュを切断面で表示
grid = pv.read('volume_only.vts')
slice_z = grid.slice(normal='z', origin=[0, 0, 0])
slice_z.plot(scalars='B_magnitude')
```

### 対処法2: 表面要素を後から生成

```python
from ngsolve import *

# 体積メッシュを読込
mesh = Mesh('volume_only.vol')

# NGSolveが自動的に境界を認識
# （ただし、元のメッシュに境界情報が必要）

# 再エクスポートすると表面要素が含まれる
mesh.ngmesh.Save('with_surface.vol')
```

---

## Q3: Netgen GUIで何も表示されない - なぜ？

**A: 以下の原因が考えられます：**

### 原因1: 体積要素のみ（表面要素なし）

```python
# 確認方法
from netgen.meshing import Mesh
mesh = Mesh()
mesh.Load('problem.vol')

print(f"Volume elements:  {mesh.ne}")
print(f"Surface elements: {mesh.nse}")  # これが0だと表示されない

if mesh.nse == 0:
    print("原因: 表面要素がありません")
```

**解決策**: メッシュ生成時に境界を含める（通常は自動）

### 原因2: メッシュサイズが小さすぎる/大きすぎる

Netgen GUIでズームが必要な場合があります：
- マウスホイールでズーム
- Viewメニュー → Center

### 原因3: ファイルが壊れている

```python
# 再生成を試す
from netgen.occ import Box, Pnt, OCCGeometry

geo = OCCGeometry(Box(Pnt(-1,-1,-1), Pnt(1,1,1)))
mesh = geo.GenerateMesh(maxh=0.2)
mesh.Save('test.vol')
```

---

## Q4: Cubitメッシュでも表面要素は含まれる？

**A: はい、含まれます。**

### Cubit → Netgen ワークフロー

```python
import cubit
import cubit_mesh_export
from ngsolve import Mesh
from netgen.gui import StartGUI

# Cubitでメッシュ生成
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("import step 'model.step'")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")

# Sidesetを定義（これが表面要素になる）
cubit.cmd("sideset 1 surface all")
cubit.cmd("sideset 1 name 'boundary'")

# Netgenへエクスポート（表面要素が含まれる）
ngmesh = cubit_mesh_export.export_netgen(cubit)
mesh = Mesh(ngmesh)

# 確認
print(f"Volume elements:  {mesh.ngmesh.ne}")
print(f"Surface elements: {mesh.ngmesh.nse}")  # > 0 のはず

# Netgen GUIで表示
StartGUI()
mesh.ngmesh.Draw()
```

**重要**: Cubitで**sideset**を定義すると、それが表面要素として変換されます。

---

## Q5: NGSolveサンプルメッシュは大丈夫？

**A: すべて表面要素を含みます。**

### 検証結果（TEST_RESULTS.mdより）

| ファイル | 体積要素 | 表面要素 | 状態 |
|---------|---------|---------|------|
| cube.vol | 756 | 338 (Triangle) | ✅ OK |
| coil.vol | 1709 | あり | ✅ OK |
| coilshield.vol | 1798 | 376 (Tri+Quad) | ✅ OK |
| beam.vol | 31 | あり | ✅ OK |
| shaft.vol | 1622 | あり | ✅ OK |
| chip.vol | 0 | あり (Surface-only) | ✅ OK |
| doubleglazing.vol | 0 | あり (Surface-only) | ✅ OK |
| square.vol | 0 | あり (Surface-only) | ✅ OK |

**結論**: NGSolveの全サンプルメッシュはNetgen GUIで表示可能。

---

## Q6: .volファイルが表示できるか事前確認する方法は？

**A: チェックスクリプトを使用：**

```bash
python utils/check_vol_surface_elements.py mesh.vol
```

**出力例**:

```
Analyzing: mesh.vol
============================================================

Mesh Statistics:
  Vertices:        228
  Volume elements: 756
  Surface elements: 338

  Volume element types:
    Tet: 756

  Surface element types:
    Triangle: 338

============================================================
Display Compatibility:
============================================================

✅ Netgen GUI: COMPATIBLE
   - Surface elements present: 338
   - Mesh will be displayed as surface
   - Recommended viewer: Netgen GUI
```

---

## まとめ

### 実用上の結論

| 状況 | 表面要素 | 対応 |
|------|---------|------|
| **Netgen生成メッシュ** | ✅ 自動生成 | 何もしなくてOK |
| **Cubit → Netgen** | ✅ 自動変換 | Sideset定義でOK |
| **GMSH → NGSolve** | ✅ 自動変換 | 何もしなくてOK |
| **NGSolveサンプル** | ✅ すべて含む | 何もしなくてOK |
| **体積要素のみ（稀）** | ❌ なし | ParaView/PyVista使用 |

### 推奨ワークフロー

```
CAD (STEP) → Netgen/Cubit → Mesh generation → .vol file
                                                 ↓
                                         表面要素自動生成
                                                 ↓
                                           Netgen GUI ✅
```

**結論**: 通常のメッシュ生成ワークフローでは、表面要素は自動的に含まれるため、Netgen GUIの制限は実用上問題になりません。

---

## 困った時のチェックリスト

- [ ] `check_vol_surface_elements.py` で表面要素を確認
- [ ] メッシュ生成元を確認（Netgen/Cubit/GMSH？）
- [ ] Cubitの場合、sideset定義を確認
- [ ] NGSolveサンプルメッシュで動作確認
- [ ] それでもダメなら ParaView/PyVista を使用

---

**作成日**: 2026-02-12
**対象**: Netgen GUI, .vol files, NGSolve mesh workflow
