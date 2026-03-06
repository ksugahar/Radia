# GMSH Workflow for Radia

## Overview

GMSH serves as the standard mesh generator for Radia, supporting both magnetic materials (volume mesh) and conductors (surface mesh).

```
┌─────────────────────────────────────────────────────────────────┐
│                    GMSH → Radia Workflow                         │
├─────────────────────────────────────────────────────────────────┤
│                                                                  │
│  CAD (STEP/IGES) → GMSH → .msh → NGSolve → Radia               │
│                                                                  │
│  用途別メッシュタイプ:                                            │
│    - 磁性体 (MMM/MSC): 体積メッシュ (Tet4, Hex8)                 │
│    - 導体 (PEEC):     表面メッシュのみ (Tri3, Quad4)             │
│                                                                  │
└─────────────────────────────────────────────────────────────────┘
```

---

## メッシュタイプの使い分け

### 磁性体（永久磁石・軟磁性材料）

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

### 導体（PEEC）

**要求**: 表面メッシュのみ（Surface Mesh）

| 要素タイプ | GMSH要素 | 理由 |
|----------|---------|------|
| Triangle | Tri3 | 表面電流分布 |
| Quadrilateral | Quad4 | 表面電流分布 |

**重要**: PEECは表面電流モデルのため、**体積メッシュは不要**

**GMSH生成**:
```python
gmsh.model.mesh.generate(2)  # 2D表面メッシュのみ
```

**理由**:
- 表皮効果: SIBC (Surface Impedance Boundary Condition) で処理
- 導体内部: 電流密度は指数減衰（表面インピーダンスで表現）
- 計算効率: 表面のみで十分な精度

---

## ワークフロー1: 磁性体（体積メッシュ）

### CADファイルからの読込

```python
import gmsh
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia
import radia as rad

rad.FldUnits('m')

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

---

## ワークフロー2: 導体（表面メッシュ）

### コイル表面メッシュ生成

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
    print("⚠️  WARNING: Volume elements found - PEEC only needs surface!")

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

---

## ワークフロー3: 磁性体+導体の統合モデル

### 例: 電磁石（鉄心+コイル）

```python
import gmsh
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia
import radia as rad

rad.FldUnits('m')

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

## GMSH vs MEG vs Netgen 比較

| 観点 | GMSH | MEG | Netgen |
|------|------|-----|--------|
| **CAD読込** | ✅ STEP/IGES直接 | ❌ Cubit経由のみ | ✅ STEP/OCC |
| **ライセンス** | ✅ オープンソース | ❌ ELF固有 | ✅ オープンソース |
| **NGSolve連携** | ✅ .msh直接読込 | ⚠️ 変換必要 | ✅ ネイティブ |
| **2D/軸対称** | ✅ 対応 | ✅ 対応（DIM='K'/'R'） | ⚠️ 3Dのみ推奨 |
| **表面メッシュ** | ✅ `generate(2)` | ⚠️ 体積メッシュ前提 | ✅ 自動生成 |
| **体積メッシュ** | ✅ Tet/Hex対応 | ✅ Tet/Hex対応 | ✅ Tet（Hexは外部） |
| **可視化** | ✅ GMSH GUI | ⚠️ 専用ツール | ✅ Netgen GUI |

**推奨**:
- **標準**: GMSH（CAD読込、表面メッシュ、NGSolve統合）
- **単純形状**: Netgen OCC（コード生成、自動メッシュ）
- **ELF互換**: MEG（ELF/MAGICとの相互運用が必要な場合）

---

## よくある質問

### Q1: PEECに体積メッシュは必要ないのか？

**A: 不要です。** PEECは表面電流近似を使用します。

**理由**:
1. **表皮効果**: 高周波では電流は表面に集中
2. **SIBC**: 表面インピーダンスで導体内部の電流分布を表現
3. **計算効率**: 表面メッシュのみで十分な精度

**適用範囲**: 周波数 × サイズ が表皮深さより大きい場合

### Q2: GMSHで六面体メッシュは生成できるか？

**A: 限定的です。**

- **Tet（四面体）**: 完全自動生成 ✅
- **Hex（六面体）**: 構造格子のみ ⚠️
- **複雑形状のHex**: Coreform Cubit推奨

**GMSH Hexメッシュ生成方法**:
```python
# 構造格子（ブロック形状のみ）
gmsh.model.mesh.setTransfiniteSurface(surf_tag)
gmsh.model.mesh.setTransfiniteVolume(vol_tag)
gmsh.model.mesh.setRecombine(3, vol_tag)
```

### Q3: GMSHとNetgenの使い分けは？

| 用途 | ツール |
|------|--------|
| **CADファイル読込** | GMSH（より対応形式が多い） |
| **単純形状（OCC）** | Netgen（コード生成が簡潔） |
| **表面メッシュのみ** | GMSH（`generate(2)`で明示的） |
| **高品質Tetメッシュ** | Netgen（メッシュ品質が良い） |

---

## サンプルスクリプト

| ファイル | 説明 |
|---------|------|
| `examples/visualization/demo_gmsh_cad_import.py` | CAD読込と体積メッシュ |
| `examples/peec_integration/demo_gmsh_surface_mesh.py` | 表面メッシュ（PEEC導体） |
| `examples/visualization/demo_gmsh_workflow.py` | 基本的なGMSH Python API |

---

## まとめ

### Radiaでの推奨ワークフロー

```
磁性体（永久磁石・鉄心）:
  CAD → GMSH → 体積メッシュ(.msh) → NGSolve → Radia (MMM/MSC)

導体（コイル・シールド）:
  CAD → GMSH → 表面メッシュ(.msh) → (将来: PEEC変換)
  現状: rad.CndLoop() で代替

統合モデル（電磁石等）:
  磁性体 + 導体 → rad.ObjCnt() → rad.Solve()
```

### キーポイント

1. ✅ **GMSH標準**: CAD読込、オープンソース、NGSolve統合
2. ✅ **メッシュタイプ**: 磁性体=体積、導体=表面
3. ✅ **表面メッシュのみ**: PEECは体積メッシュ不要
4. ✅ **NGSolve経由**: `.msh`ファイルをシームレスに読込

---

**作成日**: 2026-02-12
**対象**: Radia-GMSH統合ワークフロー
**関連ドキュメント**:
- [VISUALIZATION_WORKFLOW.md](VISUALIZATION_WORKFLOW.md)
- [PEEC_INTEGRATION.md](PEEC_INTEGRATION.md)（将来作成予定）
