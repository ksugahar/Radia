# GMSH Models for PEEC Analysis

This directory contains Coreform Cubit scripts to generate PEEC conductor meshes in GMSH format.

## Workflow

```
Coreform Cubit → .msh ファイル → Python (NGSolve) → PEEC 解析
```

## 1. Generate Mesh (Cubit → GMSH)

```bash
# Run with Coreform Cubit's Python
"${CUBIT_PATH:-<Coreform Cubit 2025.8+>/bin}/python3/python.exe" generate_coil_cubit.py
```

**Output**: `circular_coil.msh` (GMSH v2.2 format, surface mesh only)

## 2. View Mesh in GMSH GUI

```bash
gmsh circular_coil.msh
```

**GMSH GUI 操作**:
- 回転: 左ドラッグ
- ズーム: マウスホイール
- メッシュ表示: Press '0'
- エッジ表示: Press 'e'

## 3. Analyze with PEEC

```bash
cd ..
python demo_gmsh_to_peec.py
```

## File Description

| File | Description |
|------|-------------|
| `generate_coil_cubit.py` | Cubit script to generate toroidal coil surface mesh (basic) |
| `generate_coil_with_ports.py` | Cubit script with PORT definitions for PEEC |
| `circular_coil.msh` | Generated GMSH v2.2 mesh file (surface only) |
| `circular_coil_with_ports.msh` | GMSH mesh with port physical groups |
| `PORT_DEFINITION.md` | Guide for PEEC port specification |
| `README.md` | This file |

## Mesh Requirements for PEEC

**CRITICAL**: PEEC requires **SURFACE MESH ONLY**, not volume mesh.

| Mesh Type | PEEC | Reason |
|-----------|------|--------|
| Surface (Tri3, Quad4) | OK Required | Surface current distribution |
| Volume (Tet4, Hex8) | X Not needed | SIBC handles skin effect |

## Coil Parameters

| Parameter | Value | Notes |
|-----------|-------|-------|
| Mean radius | 50 mm | Coil center radius |
| Wire width | 4 mm | Radial direction |
| Wire height | 4 mm | Axial direction |
| Mesh size | 1 mm | Target element size |

## Format: GMSH v2.2

**Why v2.2**:
- GMSH GUI visualization
- Maximum compatibility
- Simpler structure than v4.1

**Note**: For NGSolve FEM computation, use `export netgen "mesh.vol"` (.vol format) instead of .msh. The .msh format is used here for GMSH visualization and PEEC surface mesh input only.

## Next Steps

1. OK Generate mesh: `generate_coil_cubit.py` (COMPLETED)
2. OK View mesh: `gmsh circular_coil.msh` (COMPLETED)
3. OK Load mesh: `demo_gmsh_to_peec.py` with GMSH Python API (COMPLETED)
4. TODO: Implement PEEC matrix calculation
5. TODO: Port definition for impedance calculation

---

**Created**: 2026-02-12
**Purpose**: GMSH-based PEEC workflow for Radia
