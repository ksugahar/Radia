# Mesh Guide: GMSH, Netgen, and Cubit Workflows for Radia

## Overview

This guide consolidates the mesh generation workflows for Radia, covering GMSH, Netgen, and Coreform Cubit. It explains mesh types, surface elements, the SetGeomInfo API for high-order curving, and troubleshooting procedures.

```
+-----------------------------------------------------------------+
|                    CAD -> Mesh -> Radia Workflow                 |
+-----------------------------------------------------------------+
|                                                                  |
|  CAD (STEP/IGES) -> GMSH / Netgen / Cubit -> Mesh -> Radia      |
|                                                                  |
|  Mesh types by application:                                      |
|    - Magnetic materials (MMM/MSC): Volume mesh (Tet4, Hex8)      |
|    - Conductors (PEEC):            Surface mesh only (Tri3, Quad4)|
|                                                                  |
|  Mesh file formats:                                              |
|    - GMSH:   .msh -> NGSolve -> Radia                            |
|    - Netgen:  .vol -> NGSolve -> Radia                            |
|    - Cubit:  export_netgen() -> NGSolve -> Radia                  |
|                                                                  |
+-----------------------------------------------------------------+
```

---

## 1. Mesh Types (メッシュタイプの使い分け)

### 1.1 Volume Elements (体積要素)

Volume elements -- tetrahedra, hexahedra, wedges -- fill the interior of a 3-D domain. They are required for magnetic material modelling (permanent magnets, soft magnetic materials).

**要求**: 体積メッシュ（Volume Mesh）

| 要素タイプ | GMSH要素 | Radia API | 用途 |
|----------|---------|----------|------|
| Tetrahedron | Tet4 | `ObjTetrahedron()` | 複雑形状 |
| Hexahedron | Hex8 | `ObjHexahedron()` | 構造格子 |
| Wedge/Prism | Wedge6 | `ObjWedge()` | 遷移要素 |

**GMSH生成**:
```python
gmsh.model.mesh.generate(3)  # 3D体積メッシュ
```

### 1.2 Surface Elements (表面要素)

Surface elements -- triangles or quadrilaterals -- cover the **boundary** of a 3-D domain. They are the outer skin of a mesh and serve two purposes:

1. **PEEC conductors**: Surface mesh is all that is needed for surface-current modelling.
2. **Netgen GUI display**: The Netgen GUI renders surface elements, not volume elements directly.

| 要素タイプ | GMSH要素 | 用途 |
|----------|---------|------|
| Triangle | Tri3 | 表面電流分布 / 境界表示 |
| Quadrilateral | Quad4 | 表面電流分布 / 境界表示 |

**PEEC用GMSH生成**:
```python
gmsh.model.mesh.generate(2)  # 2D表面メッシュのみ
```

**重要**: PEECは表面電流モデルのため、**体積メッシュは不要**

**理由**:
- 表皮効果: SIBC (Surface Impedance Boundary Condition) で処理
- 導体内部: 電流密度は指数減衰（表面インピーダンスで表現）
- 計算効率: 表面のみで十分な精度

### 1.3 Auto-Generation in Standard Workflows

In every standard mesh-generation workflow surface elements are created automatically:

| Workflow | Surface elements | Reason |
|----------|-----------------|--------|
| **Netgen direct** (`geo.GenerateMesh()`) | Auto | Boundary mesh generated automatically |
| **NGSolve `Mesh()`** | Auto | STEP/OCC import recognises boundaries |
| **Cubit -> `export_netgen()`** | Auto | Cubit sidesets are converted to boundary elements |
| **GMSH -> NGSolve** | Auto | `.msh` files include boundary elements |

**In short, normal mesh generation requires no extra steps.**

```
CAD (STEP) -> Netgen / Cubit / GMSH -> Mesh generation -> .vol / .msh file
                                                             |
                                                     Surface elements
                                                   generated automatically
                                                             |
                                                       Netgen GUI OK
```

### 1.4 NGSolve Sample Meshes

Every mesh shipped under `share/ngsolve/` already contains surface elements:

| File | Volume elements | Surface elements | Netgen GUI |
|------|----------------|-----------------|-----------|
| cube.vol | 756 | 338 (Triangle) | OK |
| coil.vol | 1709 | Present | OK |
| coilshield.vol | 1798 | 376 (Tri+Quad) | OK |
| beam.vol | 31 | Present | OK |
| shaft.vol | 1622 | Present | OK |
| chip.vol | 0 | Present (Surface-only) | OK |
| doubleglazing.vol | 0 | Present (Surface-only) | OK |
| square.vol | 0 | Present (Surface-only) | OK |

All NGSolve sample `.vol` files display correctly in the Netgen GUI.

---

## 2. GMSH Workflows

### 2.1 Workflow 1: 磁性体（体積メッシュ）

#### CADファイルからの読込

```python
import gmsh
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia
import radia as rad

# GMSH初期化
gmsh.initialize()
gmsh.option.setNumber("General.Terminal", 1)
gmsh.model.add("magnetic_core")

# CADファイル読込（STEP, IGES, BREP, STL対応）
gmsh.merge("core.step")
gmsh.model.geo.synchronize()

# メッシュサイズ設定
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.002)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.005)

# 物理グループ定義（重要！）
volumes = gmsh.model.getEntities(3)
if volumes:
    volume_tags = [v[1] for v in volumes]
    gmsh.model.addPhysicalGroup(3, volume_tags, 1)
    gmsh.model.setPhysicalName(3, 1, "core")

# 体積メッシュ生成
gmsh.model.mesh.generate(3)  # 3D volume mesh

# エクスポート
gmsh.write('core.msh')
gmsh.finalize()

# NGSolve経由でRadia変換
mesh = Mesh('core.msh')
mag_obj = netgen_mesh_to_radia(mesh,
                                material={'magnetization': [0, 0, 0]},
                                units='m',
                                material_filter='core')

# 材料適用
mat = rad.MatLin(1000)  # mu_r = 1000
rad.MatApl(mag_obj, mat)

# 解く
rad.Solve(mag_obj, 0.0001, 1000, 1)
```

### 2.2 Workflow 2: 導体（表面メッシュ / PEEC）

#### コイル表面メッシュ生成

```python
import gmsh
import numpy as np

gmsh.initialize()
gmsh.model.add("coil_surface")

# コイル断面（矩形）をXZ平面で定義
r_inner = 0.048  # m
r_outer = 0.052  # m
z_bottom = -0.001  # m
z_top = 0.001  # m

p1 = gmsh.model.geo.addPoint(r_inner, 0, z_bottom)
p2 = gmsh.model.geo.addPoint(r_outer, 0, z_bottom)
p3 = gmsh.model.geo.addPoint(r_outer, 0, z_top)
p4 = gmsh.model.geo.addPoint(r_inner, 0, z_top)

l1 = gmsh.model.geo.addLine(p1, p2)
l2 = gmsh.model.geo.addLine(p2, p3)
l3 = gmsh.model.geo.addLine(p3, p4)
l4 = gmsh.model.geo.addLine(p4, p1)

loop = gmsh.model.geo.addCurveLoop([l1, l2, l3, l4])
surf = gmsh.model.geo.addPlaneSurface([loop])

# Z軸周りに回転（完全なコイル表面生成）
gmsh.model.geo.revolve(
    [(2, surf)],
    0, 0, 0,  # 回転軸原点
    0, 0, 1,  # 回転軸方向（Z）
    2 * np.pi  # 角度（全周）
)

gmsh.model.geo.synchronize()

# メッシュサイズ
gmsh.option.setNumber("Mesh.CharacteristicLengthMin", 0.0005)
gmsh.option.setNumber("Mesh.CharacteristicLengthMax", 0.001)

# 物理グループ（表面のみ！）
surfaces = gmsh.model.getEntities(2)
surface_tags = [s[1] for s in surfaces]
gmsh.model.addPhysicalGroup(2, surface_tags, 1)
gmsh.model.setPhysicalName(2, 1, "conductor")

# 表面メッシュのみ生成（dim=2）
gmsh.model.mesh.generate(2)  # Surface mesh ONLY

# 確認: 体積要素がないことを確認
vol_elements = gmsh.model.mesh.getElements(3)
if vol_elements[1] and any(len(e) > 0 for e in vol_elements[1]):
    print("WARNING: Volume elements found - PEEC only needs surface!")

gmsh.write('coil_surface.msh')
gmsh.finalize()

# PEEC変換（将来のAPI）
# from peec_mesh_import import surface_mesh_to_peec
# conductor = surface_mesh_to_peec(mesh, sigma=5.8e7)
```

**現状の代替手段**:

```python
# 単純形状の場合は CndLoop を使用
coil = rad.CndLoop([0, 0, 0], 0.05, [0, 0, 1], 'r',
                   0.002, 0.002, 5.8e7, 8, 36)
```

### 2.3 Workflow 3: 磁性体+導体の統合モデル

#### 例: 電磁石（鉄心+コイル）

```python
import gmsh
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia
import radia as rad

# ===============================
# 1. 鉄心（体積メッシュ）
# ===============================
gmsh.initialize()
gmsh.model.add("core")
gmsh.merge("core.step")  # CAD読込

# 体積メッシュ
gmsh.model.mesh.generate(3)
gmsh.write('core.msh')
gmsh.finalize()

mesh_core = Mesh('core.msh')
core_obj = netgen_mesh_to_radia(mesh_core,
                                 material={'magnetization': [0, 0, 0]},
                                 units='m')
mat_iron = rad.MatLin(1000)
rad.MatApl(core_obj, mat_iron)

# ===============================
# 2. コイル（表面メッシュまたは解析形状）
# ===============================
# 現状: CndLoopを使用（簡易コイル）
coil_obj = rad.CndLoop([0, 0, 0], 0.05, [0, 0, 1], 'r',
                       0.002, 0.002, 5.8e7, 8, 36)

# 将来: GMSH表面メッシュからPEEC変換
# gmsh.initialize()
# ... (coil surface mesh generation)
# coil_obj = surface_mesh_to_peec(mesh_coil, sigma=5.8e7)

# ===============================
# 3. 統合して解く
# ===============================
container = rad.ObjCnt([core_obj, coil_obj])
rad.Solve(container, 0.0001, 1000, 1)

# フィールド計算
B = rad.Fld(container, 'b', [0, 0, 0.1])
print(f"Field at (0, 0, 0.1): {B} T")
```

---

## 3. SetGeomInfo API (High-Order Curving)

> **Source**: [ksugahar/ngsolve](https://github.com/ksugahar/ngsolve) fork with SetGeomInfo API (netgen PR [#232](https://github.com/NGSolve/netgen/pull/232)).

### 3.1 Problem Statement

When meshes are imported from external mesh generators (Gmsh, Cubit, etc.) without geometry, `mesh.Curve(order)` fails because the UV parametric coordinates (geominfo) are not set. The `SetGeomInfo` API enables setting geominfo programmatically.

### 3.2 API

```python
Element2d.SetGeomInfo(vertex_index, u, v, trignum=0)
```

**Parameters:**
- `vertex_index`: 0-based index of the vertex within the element
- `u`, `v`: Surface parametric coordinates from the OCC geometry
- `trignum`: Triangle number for STL meshing (default: 0)

### 3.3 Recommended Workflow: Coreform Cubit + Name-based Mapping

```
1. OCC: Create geometry and name faces (name_occ_faces)
2. OCC: Export to STEP (face names preserved)
3. Cubit: Import STEP, generate mesh
4. Export: Use name-based face mapping to Netgen mesh
5. SetGeomInfo: Compute UV parameters analytically
6. mesh.Curve(order): High-order curving works correctly!
```

### 3.4 Code Example

```python
import cubit_mesh_export
from netgen.occ import OCCGeometry, Box, Cylinder, gp_Pnt, gp_Ax2, gp_Dir
from ngsolve import Mesh

# 1. Create geometry in OCC
brick = Box(gp_Pnt(-1,-1,-1), gp_Pnt(1,1,1))
cyl = Cylinder(gp_Ax2(gp_Pnt(0,0,-2), gp_Dir(0,0,1)), 0.3, 4)
shape = brick - cyl

# 2. Name faces (critical for correct mapping!)
cubit_mesh_export.name_occ_faces(shape)

# 3. Export STEP
shape.WriteStep("geometry.step")

# 4. Load geometry reference
geo = OCCGeometry("geometry.step")

# 5. Import into Cubit and mesh (Cubit commands)
# cubit.cmd('import step "geometry.step" noheal')
# cubit.cmd("volume all scheme tetmesh")
# cubit.cmd("volume all size 0.15")
# cubit.cmd("mesh volume all")

# 6. Export with name-based mapping
ngmesh = cubit_mesh_export.export_netgen_with_names(cubit, geo)

# 7. Set UV for curved surfaces
cubit_mesh_export.set_cylinder_geominfo(ngmesh, radius=0.3, height=4.0,
                                         center=(0,0,-2), axis='z')

# 8. High-order curving
mesh = Mesh(ngmesh)
mesh.Curve(2)  # Now works correctly!
```

### 3.5 Available SetGeomInfo Helper Functions

From `cubit_mesh_export` ([PyPI](https://pypi.org/project/Coreform-Cubit-Mesh-Export/)):

```python
set_cylinder_geominfo(ngmesh, radius, height, center=(0,0,0), axis='z')
set_sphere_geominfo(ngmesh, radius, center=(0,0,0))
set_torus_geominfo(ngmesh, major_radius, minor_radius, center=(0,0,0), axis='z')
set_cone_geominfo(ngmesh, base_radius, height, center=(0,0,0), axis='z')
```

### 3.6 Accuracy Results

| Geometry | Curve(2) Error | Curve(3) Error |
|----------|----------------|----------------|
| Complex (Boolean ops) | **0.0021%** | **0.0004%** |
| Cylinder | 0.0027% | 0.0006% |
| Sphere | 0.0027% | 0.0004% |
| Torus | 0.0010% | 0.0003% |

All results achieve **Netgen-native accuracy** (<0.003%).

### 3.7 Requirements

- NGSolve: Build from `ksugahar/ngsolve` branch `feature/setgeominfo`
- Coreform Cubit 2025.3+
- `pip install coreform-cubit-mesh-export`

### 3.8 Examples and Links

- Full working examples: [ksugahar/Coreform_Cubit_Mesh_Export/examples/netgen/](https://github.com/ksugahar/Coreform_Cubit_Mesh_Export/tree/main/examples/netgen)
- Netgen PR: [NGSolve/netgen#232](https://github.com/NGSolve/netgen/pull/232)
- Forum: [Feature Request - SetGeomInfo API](https://forum.ngsolve.org/t/feature-request-python-api-for-high-order-curving-of-externally-imported-meshes/3810)
- PyPI: [coreform-cubit-mesh-export](https://pypi.org/project/Coreform-Cubit-Mesh-Export/)

---

## 4. Surface Elements: Display and Workarounds

### 4.1 Netgen GUI Display Behaviour

The Netgen GUI renders **surface elements**, not volume elements directly:

- **Surface elements present** -- the mesh boundary is displayed normally.
- **Volume elements only (no surface elements)** -- nothing is displayed, or a clipping plane is required to see anything.

### 4.2 Viewer Selection Guide

| Mesh type | Netgen GUI | ParaView | PyVista | webgui |
|-----------|-----------|----------|---------|--------|
| **Surface elements present** | **Recommended** | Overkill | Overkill | Geometry only |
| **Volume elements only** | Cannot display | **Recommended** | **Recommended** | Needs GridFunction |
| **Geometry check** | **Best** | Requires meshing | Requires meshing | OCC direct |
| **Field visualisation** | Not supported | **Best quality** | Fast | Interactive |

### 4.3 When Netgen GUI is the right tool

- The mesh contains surface elements.
- You need to inspect geometry or mesh quality.
- You want a lightweight, fast viewer.
- You are following an integrated shape -> mesh -> review workflow.

### 4.4 When to use ParaView / PyVista instead

- The mesh contains only volume elements.
- You need to visualise internal structure via slicing or clipping.
- You need to visualise field data (B, H, etc.).
- You need publication-quality figures.

### 4.5 Workaround: ParaView slice / clip (recommended for volume-only meshes)

```python
# 1. Export to VTS
import radia as rad
rad.FldVTS(magnet, 'field.vts', ...)

# 2. Open in ParaView
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

### 4.6 Workaround: PyVista slice / clip

```python
import pyvista as pv

# Load VTS
grid = pv.read('field.vts')

# Create slice at z=0.05 m
slice_z = grid.slice(normal='z', origin=[0, 0, 0.05])
slice_z.plot(scalars='B_magnitude', cmap='coolwarm')

# Or clip half of the domain
clipped = grid.clip(normal='z', origin=[0, 0, 0])
clipped.plot(scalars='B_magnitude', cmap='viridis')
```

### 4.7 Workaround: Regenerate mesh with surface elements via NGSolve

```python
from ngsolve import *
from netgen.occ import Box, Pnt, OCCGeometry

# Create geometry
box = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))
geo = OCCGeometry(box)

# Generate mesh (includes surface elements)
mesh = Mesh(geo.GenerateMesh(maxh=0.1))

# Export to .vol (surface elements included)
mesh.ngmesh.Save('mesh_with_surface.vol')

# Open in Netgen GUI
from netgen.gui import StartGUI
StartGUI()
mesh.ngmesh.Draw()
```

### 4.8 Workaround: Add surface elements to an existing volume-only mesh

```python
from ngsolve import *

# Load volume-only mesh
mesh = Mesh('volume_only.vol')

# NGSolve recognises boundaries automatically
# (the original mesh must still contain boundary information)

# Re-export -- surface elements will be included
mesh.ngmesh.Save('with_surface.vol')
```

### 4.9 Cubit Meshes and Surface Elements

When you define a **sideset** in Cubit and export via `export_netgen()`, the sideset surfaces become surface elements:

```python
import cubit
import cubit_mesh_export
from ngsolve import Mesh
from netgen.gui import StartGUI

# Cubit mesh generation
cubit.init(['cubit', '-nojournal', '-batch'])
cubit.cmd("import step 'model.step'")
cubit.cmd("volume all scheme tetmesh")
cubit.cmd("mesh volume all")

# Define sidesets (these become surface elements)
cubit.cmd("sideset 1 surface all")
cubit.cmd("sideset 1 name 'boundary'")

# Export to Netgen (surface elements included)
ngmesh = cubit_mesh_export.export_netgen(cubit)
mesh = Mesh(ngmesh)

# Verify
print(f"Volume elements:  {mesh.ngmesh.ne}")
print(f"Surface elements: {mesh.ngmesh.nse}")  # should be > 0

# Display in Netgen GUI
StartGUI()
mesh.ngmesh.Draw()
```

**Key point**: defining a sideset in Cubit is what creates the corresponding surface elements on export.

---

## 5. Troubleshooting

### 5.1 How to check whether surface elements are present

```bash
python utils/check_vol_surface_elements.py mesh.vol
```

Or programmatically:

```python
from netgen.meshing import Mesh

mesh = Mesh()
mesh.Load('mesh.vol')

print(f"Volume elements:  {mesh.ne}")
print(f"Surface elements: {mesh.nse}")

if mesh.nse == 0:
    print("Warning: No surface elements")
    print("Netgen GUI may not display this mesh")
```

**Example output** of the check script:

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

  Netgen GUI: COMPATIBLE
   - Surface elements present: 338
   - Mesh will be displayed as surface
   - Recommended viewer: Netgen GUI
```

### 5.2 Netgen GUI shows nothing -- possible causes

| Cause | Diagnosis | Fix |
|-------|-----------|-----|
| **No surface elements** | `mesh.nse == 0` | Regenerate mesh or use ParaView/PyVista (see Section 4) |
| **Mesh too small / too large for viewport** | Elements exist but view is empty | Mouse-wheel zoom, or *View > Center* |
| **Corrupt `.vol` file** | Load raises an error or counts are unexpected | Regenerate the mesh from geometry |

To regenerate a simple test mesh and confirm the GUI works:

```python
from netgen.occ import Box, Pnt, OCCGeometry

geo = OCCGeometry(Box(Pnt(-1, -1, -1), Pnt(1, 1, 1)))
mesh = geo.GenerateMesh(maxh=0.2)
mesh.Save('test.vol')
```

### 5.3 Quick-reference checklist

- [ ] Run `check_vol_surface_elements.py` to confirm surface element count.
- [ ] Identify the mesh source (Netgen / Cubit / GMSH).
- [ ] If Cubit, verify that sidesets are defined before export.
- [ ] Test with an NGSolve sample mesh to rule out environment issues.
- [ ] If none of the above helps, switch to ParaView or PyVista.

### 5.4 Practical summary

| Situation | Surface elements | Action |
|-----------|-----------------|--------|
| **Netgen-generated mesh** | Auto-generated | Nothing to do |
| **Cubit -> Netgen** | Auto-converted | Define sidesets |
| **GMSH -> NGSolve** | Auto-converted | Nothing to do |
| **NGSolve samples** | All included | Nothing to do |
| **Volume-only (rare)** | None | Use ParaView / PyVista |

---

## 6. Tool Comparison (GMSH vs Netgen vs Cubit)

| 観点 | GMSH | Netgen | Coreform Cubit |
|------|------|--------|----------------|
| **CAD読込** | STEP/IGES直接 | STEP/OCC | STEP/IGES直接 |
| **ライセンス** | オープンソース | オープンソース | 商用 |
| **NGSolve連携** | .msh直接読込 | ネイティブ | `export_netgen()` |
| **2D/軸対称** | 対応 | 3Dのみ推奨 | 対応 |
| **表面メッシュ** | `generate(2)` | 自動生成 | sidesetで自動 |
| **体積メッシュ** | Tet/Hex対応 | Tet（Hexは外部） | Tet/Hex対応 |
| **六面体メッシュ** | 構造格子のみ | 非対応（外部ツール） | 高品質（推奨） |
| **High-order curving** | 未対応（外部で処理） | ネイティブ | SetGeomInfo API経由 |
| **可視化** | GMSH GUI | Netgen GUI | Cubit GUI |

**推奨**:
- **標準ワークフロー**: GMSH（CAD読込、表面メッシュ、NGSolve統合）
- **単純形状**: Netgen OCC（コード生成、自動メッシュ）
- **高品質Hex / High-order curving**: Coreform Cubit + SetGeomInfo API

### GMSHとNetgenの使い分け

| 用途 | ツール |
|------|--------|
| **CADファイル読込** | GMSH（より対応形式が多い） |
| **単純形状（OCC）** | Netgen（コード生成が簡潔） |
| **表面メッシュのみ** | GMSH（`generate(2)`で明示的） |
| **高品質Tetメッシュ** | Netgen（メッシュ品質が良い） |

---

## 7. FAQ (よくある質問)

### Q1: PEECに体積メッシュは必要ないのか？

**A: 不要です。** PEECは表面電流近似を使用します。

**理由**:
1. **表皮効果**: 高周波では電流は表面に集中
2. **SIBC**: 表面インピーダンスで導体内部の電流分布を表現
3. **計算効率**: 表面メッシュのみで十分な精度

**適用範囲**: 周波数 x サイズ が表皮深さより大きい場合

### Q2: GMSHで六面体メッシュは生成できるか？

**A: 限定的です。**

- **Tet（四面体）**: 完全自動生成
- **Hex（六面体）**: 構造格子のみ
- **複雑形状のHex**: Coreform Cubit推奨

**GMSH Hexメッシュ生成方法**:
```python
# 構造格子（ブロック形状のみ）
gmsh.model.mesh.setTransfiniteSurface(surf_tag)
gmsh.model.mesh.setTransfiniteVolume(vol_tag)
gmsh.model.mesh.setRecombine(3, vol_tag)
```

### Q3: Surface elements are always required?

**A: In practice, no extra steps are needed.** In every standard mesh-generation workflow (Netgen, Cubit, GMSH), surface elements are created automatically. The only edge case is a volume-only mesh, which is rare. See Section 5 for troubleshooting if this occurs.

### Q4: What if `mesh.Curve(order)` fails on an imported mesh?

**A:** This happens because UV parametric coordinates (geominfo) are not set for externally imported meshes. Use the `SetGeomInfo` API (Section 3) to set these coordinates programmatically, or use the `cubit_mesh_export` helper functions for standard geometric shapes.

---

## 8. References

### サンプルスクリプト

| ファイル | 説明 |
|---------|------|
| `examples/visualization/demo_gmsh_cad_import.py` | CAD読込と体積メッシュ |
| `examples/peec_integration/demo_gmsh_surface_mesh.py` | 表面メッシュ（PEEC導体） |
| `examples/visualization/demo_gmsh_workflow.py` | 基本的なGMSH Python API |

### Radiaでの推奨ワークフロー (まとめ)

```
磁性体（永久磁石・鉄心）:
  CAD -> GMSH -> 体積メッシュ(.msh) -> NGSolve -> Radia (MMM/MSC)

導体（コイル・シールド）:
  CAD -> GMSH -> 表面メッシュ(.msh) -> (将来: PEEC変換)
  現状: rad.CndLoop() で代替

統合モデル（電磁石等）:
  磁性体 + 導体 -> rad.ObjCnt() -> rad.Solve()

High-order curving (Cubit):
  OCC -> STEP -> Cubit -> export_netgen() -> SetGeomInfo -> mesh.Curve(order)
```

### キーポイント

1. **GMSH標準**: CAD読込、オープンソース、NGSolve統合
2. **メッシュタイプ**: 磁性体=体積、導体=表面
3. **表面メッシュのみ**: PEECは体積メッシュ不要
4. **NGSolve経由**: `.msh`ファイルをシームレスに読込
5. **High-order curving**: SetGeomInfo APIで外部メッシュもCurve()対応

### External Links

- Netgen PR #232: [NGSolve/netgen#232](https://github.com/NGSolve/netgen/pull/232)
- SetGeomInfo Forum: [Feature Request - SetGeomInfo API](https://forum.ngsolve.org/t/feature-request-python-api-for-high-order-curving-of-externally-imported-meshes/3810)
- cubit_mesh_export PyPI: [coreform-cubit-mesh-export](https://pypi.org/project/Coreform-Cubit-Mesh-Export/)
- [PEEC_INTEGRATION.md](PEEC_INTEGRATION.md)（将来作成予定）

---

**作成日**: 2026-02-22
**対象**: Radia メッシュ生成ワークフロー（GMSH, Netgen, Cubit）
