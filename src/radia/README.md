# Radia Python Utilities

Python utility modules for Radia visualization, data export, and NGSolve integration.

## Files

### Mesh Import

#### nastran_mesh_import.py

**Unified Nastran mesh import for Radia.**

Imports Nastran mesh files (.bdf, .nas, .dat) and creates Radia geometry.

**Supported Element Types:**
- **CHEXA** - 8-node hexahedron
- **CPENTA** - 6-node wedge/prism
- **CPYRAM** - 5-node pyramid
- **CTETRA** - 4-node tetrahedron
- **CTRIA3** - 3-node triangle (surface mesh, grouped by material ID)
- **GRID/GRID*** - Node definitions (fixed-width format)

**Usage:**
```python
from nastran_mesh_import import create_radia_from_nastran
import radia as rad

# Radia always uses meters

# Import mesh and create Radia objects
mag_obj = create_radia_from_nastran(
    'York.bdf',
    material={'magnetization': [0, 0, 0]},
    units='mm',
    combine=True
)

# Apply material
mat = rad.MatSatIsoFrm([20000, 2], [0.1, 2], [0.1, 2])
rad.MatApl(mag_obj, mat)
```

**CTRIA3 Surface Meshes:**
Surface triangles are automatically grouped by material ID (property ID).
Each material ID creates one polyhedron from all its triangles.

```python
from nastran_mesh_import import import_nastran_mesh

mesh_data = import_nastran_mesh('sphere.bdf', units='mm')
tria_groups = mesh_data['tria_groups']
# Format: {material_id: {'faces': [[n1,n2,n3], ...], 'node_ids': set(...)}}
```

**Migration from nastran_reader.py (removed 2025-11-23):**
- `nastran_reader.py` has been removed and unified into `nastran_mesh_import.py`
- Use `create_radia_from_nastran()` for direct Radia object creation
- Face topologies (TETRA_FACES, WEDGE_FACES, etc.) are in `netgen_mesh_import.py`

#### netgen_mesh_import.py

**Import Netgen/NGSolve meshes to Radia.**

Provides face topology constants and mesh conversion utilities.

**Face Topology Constants:**
```python
from netgen_mesh_import import (
    TETRA_FACES,    # Tetrahedron (4 faces, 1-indexed for Radia)
    HEX_FACES,      # Hexahedron (6 faces)
    WEDGE_FACES,    # Wedge/Prism (5 faces)
    PYRAMID_FACES   # Pyramid (5 faces)
)
```

**ObjTetrahedron API (推奨):**

v1.4.0以降、`ObjTetrahedron` を使用すると面定義を省略できます。

```python
import radia as rad
import numpy as np

rad.UtiDelAll()
# Radia always uses meters (メートル単位)

# 四面体の頂点座標 (メートル単位)
vertices = [
    [0.0, 0.0, 0.0],   # v0
    [1.0, 0.0, 0.0],   # v1
    [0.5, 0.866, 0.0], # v2
    [0.5, 0.289, 0.816] # v3
]

# 四面体オブジェクトを作成 (面は自動生成)
obj = rad.ObjTetrahedron(vertices, [0, 0, 0])

# 材料を適用
chi = 999  # mu_r = 1000
mat = rad.MatLin(chi)
rad.MatApl(obj, mat)
```

**Note:** 内部的には `TETRA_FACES = [[1,3,2], [1,2,4], [2,3,4], [3,1,4]]` が自動的に適用されます。

**使用例 (netgen_mesh_to_radia で自動変換):**

```python
from netgen.occ import Box, Pnt, OCCGeometry
from ngsolve import Mesh
from netgen_mesh_import import netgen_mesh_to_radia
import radia as rad

rad.UtiDelAll()

# Netgen でメッシュを生成
cube = Box(Pnt(-0.5, -0.5, -0.5), Pnt(0.5, 0.5, 0.5))
cube.mat('magnetic')
geo = OCCGeometry(cube)
mesh = Mesh(geo.GenerateMesh(maxh=0.3))

# Radia オブジェクトに変換
radia_obj = netgen_mesh_to_radia(
    mesh,
    material={'magnetization': [0, 0, 0]},
    units='m',
    material_filter='magnetic'
)

# 磁化を取得
all_M = rad.ObjM(radia_obj)
M_list = [m[1] for m in all_M]  # [[center, M], ...] から M を抽出
M_avg_z = np.mean([m[2] for m in M_list])
```

**重要な注意事項:**

1. **1-indexed**: Radia の内部 `ObjPolyhdr` API は **1-indexed** の面定義を要求します (Python では `ObjHexahedron`/`ObjTetrahedron` を使用)
2. **単位**: Radia always uses meters - 座標はメートル単位にする
3. **ObjM の戻り値**: コンテナに対する `ObjM` は `[[center1, M1], [center2, M2], ...]` を返す

### NGSolve Integration

#### RadiaField (integrated into _radia_pybind.pyd since v2.5.0)

**C++ CoefficientFunction integration for NGSolve.**

High-performance NGSolve CoefficientFunction wrappers for Radia magnetostatics fields. Since v2.5.0, `RadiaField` is integrated into the main `_radia_pybind.pyd` module and accessed as `rad.RadiaField()`. No separate module is needed.

**Features:**
- Exact Radia field evaluation in NGSolve
- Both NGSolve and Radia use meters (no unit conversion needed)
- Thread-safe Python GIL handling
- Compatible with both 2D and 3D meshes
- Three field types: B-field, H-field, A-field

**Usage:**
```python
from ngsolve import *
import radia as rad

# Create Radia hexahedral magnet
# ObjHexahedron auto-generates face topology from 8 vertices
vertices = [[-10,-10,-15], [10,-10,-15], [10,10,-15], [-10,10,-15],
            [-10,-10,15], [10,-10,15], [10,10,15], [-10,10,15]]
Mr = 1.2 / (4 * 3.14159 * 1e-7)  # Br=1.2T -> A/m
magnet = rad.ObjHexahedron(vertices, [0, 0, Mr])

# Create mesh (m)
mesh = Mesh(...)

# CoefficientFunction (exact)
B_cf = rad.RadiaField(magnet, 'b')
B = B_cf(mesh(0, 0, 0.02))
```

**See also:** `examples/ngsolve_integration/` for complete examples.

### coil_builder.py

**Modern fluent interface for constructing complex coil geometries.**

Elegant object-oriented design for multi-segment coil paths with automatic state tracking.

**Features:**
- Fluent method chaining
- Automatic position/orientation tracking
- Type-safe with abstract base classes
- Automatic cross-section transformation
- Direct conversion to Radia objects

**Usage:**
```python
from coil_builder import CoilBuilder

coil = (CoilBuilder(current=1000)
	.set_start([0, 0, 0])
	.set_cross_section(width=20, height=20)
	.add_straight(length=100)
	.add_arc(radius=50, arc_angle=180, tilt=90)
	.add_straight(length=100)
	.add_arc(radius=50, arc_angle=180, tilt=90)
	.to_radia())
```


**VTS export for ParaView visualization.**

Export Radia magnetic field to VTS (VTK XML Structured Grid) format using C++ implementation.

**Usage:**
```python
import radia as rad

# Radia always uses meters

# Create hexahedral magnet (30x30x10 mm, magnetization 1.2T in z)
vertices = [[-0.015,-0.015,-0.005], [0.015,-0.015,-0.005], [0.015,0.015,-0.005], [-0.015,0.015,-0.005],
            [-0.015,-0.015,0.005], [0.015,-0.015,0.005], [0.015,0.015,0.005], [-0.015,0.015,0.005]]
mag = rad.ObjHexahedron(vertices, [0, 0, 954930])

# Export to VTS
           [-0.05, 0.05], [-0.05, 0.05], [0.01, 0.05],
           21, 21, 11)
```

## Visualization Workflow

### Option 1: ParaView (VTS Export)
Best for publication-quality figures, batch processing.
```python
import radia as rad
# Radia always uses meters
# ... create magnet ...
           [-0.1, 0.1], [-0.1, 0.1], [0.0, 0.2],
           21, 21, 21)
```

### Option 2: NGSolve Integration
For coupled FEM simulations.
```python
import radia as rad
B_cf = rad.RadiaField(magnet, 'b')
```

## References

- **NGSolve Integration**: `examples/ngsolve_integration/`
- **Coil Builder Examples**: `docs/complex_coil_geometry/`
- **ParaView**: https://www.paraview.org/
- **PyVista**: https://docs.pyvista.org/

Last Updated: 2026-01-09
