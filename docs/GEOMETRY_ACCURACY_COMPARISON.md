# Geometry Accuracy Comparison for Radia Visualization

## 問題: VTS Exportによる形状近似

### rad.FldVTS()の挙動

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

---

## ビューワー別の形状精度

### 1. PyVista + VTS（形状近似）

**精度**:
- ❌ 形状: 格子で近似（VTSの制限）
- ✅ フィールド: 正確（格子点での値）

**使用例**:
```python
import pyvista as pv

grid = pv.read('field.vts')
grid.plot(scalars='B_magnitude')  # フィールド正確、形状近似
```

**適用場面**: フィールド分布の確認（形状の正確性は不要）

---

### 2. ParaView + VTS（形状近似）

**精度**:
- ❌ 形状: 格子で近似（VTSの制限）
- ✅ フィールド: 正確（格子点での値）
- ✅ レンダリング: 最高品質

**使用例**:
```bash
paraview field.vts
```

**適用場面**: 論文図表（フィールド分布）、形状は別途表示

---

### 3. ParaView + STL/STEP + VTS（形状正確）

**精度**:
- ✅ 形状: 正確（STL/STEPで表現）
- ✅ フィールド: 正確（VTSで表現）
- ✅ レンダリング: 最高品質

**使用例**:
```python
# 1. 形状をSTL/STEPでエクスポート
export_radia_geometry_to_stl(magnet, 'magnet.stl')

# 2. フィールドをVTSでエクスポート
rad.FldVTS(magnet, 'field.vts', ...)

# 3. ParaViewで両方を重ね合わせ
# - magnet.stl: 形状（半透明）
# - field.vts: フィールド（カラーマップ）
```

**適用場面**: 論文図表（形状とフィールドの両方が重要）

---

### 4. NGSolve webgui + OCC形状（形状正確）

**精度**:
- ✅ 形状: 正確（OCC CAD表現）
- ✅ フィールド: 正確（GridFunction投影）
- ⚠️ レンダリング: ParaViewより劣るが十分

**使用例**:
```python
from ngsolve.webgui import Draw
from netgen.occ import Box, Pnt

# OCC形状（正確）
occ_magnet = Box(Pnt(-0.02, -0.02, -0.01), Pnt(0.02, 0.02, 0.01))
Draw(occ_magnet, name='Magnet')

# フィールド（GridFunction）
Draw(B_gf, mesh, 'B_field')
```

**適用場面**: 開発・デバッグ、インタラクティブ探索、形状の正確性が必要

---

## 推奨ワークフロー

### ケース1: フィールド分布の確認（形状精度不要）

```
Radia → FldVTS() → PyVista
                  ↓
            迅速な確認（形状近似でOK）
```

**理由**: 開発中はフィールド分布の傾向が分かれば十分

---

### ケース2: 論文図表（フィールドのみ）

```
Radia → FldVTS() → ParaView
                  ↓
            高品質レンダリング（形状近似でOK）
```

**理由**: フィールド分布図では形状の境界は重要でない

---

### ケース3: 論文図表（形状とフィールド両方重要）

```
Radia → ExportOCC() → STL/STEP ──┐
     ↓                            ├→ ParaView overlay
     └→ FldVTS() → VTS ───────────┘
                  ↓
            形状正確 + フィールド正確
```

**手順**:
1. 形状をSTL/STEPでエクスポート（正確）
2. フィールドをVTSでエクスポート（格子）
3. ParaViewで重ね合わせ表示

---

### ケース4: インタラクティブ探索（形状精度必要）

```
Radia → ExportOCC() → OCC shape ──┐
     ↓                             ├→ NGSolve webgui
     └→ RadiaField CF → GridFunction ┘
                  ↓
            ブラウザベース、形状正確
```

**理由**: webguiはOCC形状をそのまま表示（近似なし）

---

## 実装状況

### 完成済み

- ✅ `rad.FldVTS()` - フィールドVTSエクスポート
- ✅ PyVista可視化
- ✅ ParaView可視化（VTSのみ）
- ✅ NGSolve webgui基本機能

### 開発中（TODO）

- ⏳ `rad.ExportOCC()` - Radia解析オブジェクト → OCC形状変換
  - `ObjRecMag` → `netgen.occ.Box`
  - `ObjCylMag` → `netgen.occ.Cylinder`
  - `ObjSphMag` → `netgen.occ.Sphere`
- ⏳ STL/STEPエクスポート自動化
- ⏳ ParaView overlay自動化スクリプト

### 参考実装

**EMPY_Field** (`S:\NGSolve\EMPY\EMPY_Field`):
- Radia解析オブジェクトのOCC変換実装例
- 正確な形状をOCCで表現

---

## まとめ

### 形状精度が重要な場合

| 方法 | 形状精度 | フィールド精度 | 手間 | 品質 |
|------|---------|--------------|------|-----|
| **webgui + OCC** | ✅ 正確 | ✅ 正確 | 中 | 良 |
| **ParaView + STL/STEP + VTS** | ✅ 正確 | ✅ 正確 | 高 | 優 |

### 形状精度が不要な場合

| 方法 | 形状精度 | フィールド精度 | 手間 | 品質 |
|------|---------|--------------|------|-----|
| **PyVista + VTS** | ❌ 近似 | ✅ 正確 | 低 | 良 |
| **ParaView + VTS** | ❌ 近似 | ✅ 正確 | 中 | 優 |

---

## 結論

**ユーザーの懸念「paraview介すと必ずしも形状がそのままではない」は正しい**:
- VTS形式は構造格子なので、形状は格子で近似される
- 解析的な長方形や円柱が、格子の境界で表現される

**解決策**:
1. **形状不要**: PyVista/ParaView + VTS（現在の実装で十分）
2. **形状重要（開発）**: NGSolve webgui + OCC形状（`rad.ExportOCC()` TODO）
3. **形状重要（論文）**: ParaView + STL/STEP + VTS overlay（手動、または自動化TODO）

**推奨開発順序**:
1. Phase 1（完成済み）: PyVista/ParaView + VTS（フィールド重視）
2. Phase 2（TODO）: `rad.ExportOCC()` 実装（形状精度対応）
3. Phase 3（TODO）: ParaView overlay自動化（論文品質）

---

**作成日**: 2026-02-12
**参考**: EMPY_Field (`S:\NGSolve\EMPY\EMPY_Field`) for OCC conversion examples
