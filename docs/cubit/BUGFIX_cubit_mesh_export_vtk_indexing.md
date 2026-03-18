# Bug Fix: cubit_mesh_export VTK Node Index Mapping (2025-11-22)

## Problem

`cubit_mesh_export.export_vtk()` generated VTK files that ParaView could not read due to incorrect node index mapping, especially when Cubit node IDs have gaps (non-contiguous numbering).

**Symptoms**:
- ParaView silently fails to open VTK files (no error message, just exits)
- VTK file contains node indices that exceed the number of points
- Example: 569 points, but indices up to 606 in CELLS section
- Occurs when `compress` command is not used in Cubit journal files, resulting in non-contiguous node IDs

**Root Cause**:
- Cubit node IDs are 1-indexed and **may be non-contiguous** (e.g., 1, 2, 5, 10, ...)
- VTK requires 0-indexed **contiguous** indices (0, 1, 2, 3, ...)
- Original code used `range(cubit.get_node_count()+1)` to iterate through nodes
  - This is inefficient and may miss nodes if node IDs have large gaps
  - Used `node_id - 1` for CELLS indices, which is incorrect for non-contiguous IDs
  - **Mismatch**: If Cubit node 606 exists but is the 300th node written, VTK index should be 299, not 605

## Example

**Cubit node IDs**: 1, 2, 5, 10, 100, 200, 606 (7 nodes total)

**Original code (WRONG)**:
```
POINTS 7 float
<coord for node 1>   # VTK index 0
<coord for node 2>   # VTK index 1
<coord for node 5>   # VTK index 2
<coord for node 10>  # VTK index 3
<coord for node 100> # VTK index 4
<coord for node 200> # VTK index 5
<coord for node 606> # VTK index 6

CELLS 1 9
8 0 1 4 9 99 199 605 ...  # WRONG! Indices 605 > 6 (out of range)
```

**Fixed code (CORRECT)**:
```
POINTS 7 float
<coord for node 1>   # VTK index 0
<coord for node 2>   # VTK index 1
<coord for node 5>   # VTK index 2
<coord for node 10>  # VTK index 3
<coord for node 100> # VTK index 4
<coord for node 200> # VTK index 5
<coord for node 606> # VTK index 6

CELLS 1 9
8 0 1 2 3 4 5 6 ...  # CORRECT! All indices in range [0, 6]
```

## Solution

**Improved approach (2025-11-22 update)**: Collect actual node IDs from elements instead of iterating through all possible node IDs.

This approach is:
1. **More efficient**: Only processes nodes actually used in the mesh
2. **More robust**: Works correctly even with large gaps in node numbering
3. **Consistent with other export functions**: Uses the same pattern as `export_Gmsh_ver2()`

```python
# First, collect all unique node IDs from all elements
node_list = set()
for block_id in cubit.get_block_id_list():
    elem_types = ["hex", "tet", "wedge", "pyramid", "tri", "face", "edge", "node"]
    for elem_type in elem_types:
        if elem_type == "hex":
            func = getattr(cubit, f"get_block_{elem_type}es")
        else:
            func = getattr(cubit, f"get_block_{elem_type}s")
        for element_id in func(block_id):
            node_ids = cubit.get_expanded_connectivity(elem_type, element_id)
            node_list.update(node_ids)

# Write points in sorted order and create mapping
fid.write(f'POINTS {len(node_list)} float\n')
node_id_to_vtk_index = {}
vtk_index = 0
for node_id in sorted(node_list):
    coord = cubit.get_nodal_coordinates(node_id)
    fid.write(f'{coord[0]} {coord[1]} {coord[2]}\n')
    node_id_to_vtk_index[node_id] = vtk_index
    vtk_index += 1
```

Then use mapping for all cell connectivity:

```python
# BEFORE (WRONG)
fid.write(f'8 {node_list[0]-1} {node_list[1]-1} ...\n')

# AFTER (CORRECT)
fid.write(f'8 {node_id_to_vtk_index[node_list[0]]} {node_id_to_vtk_index[node_list[1]]} ...\n')
```

## Files Modified

**File**: `src/radia/cubit_mesh_export.py`

**Function**: `export_vtk(cubit, FileName: str)` (v1.5.1+: ORDER parameter removed, auto-detection)

**Changes** (Updated 2025-11-22):
1. Lines 500-544: Completely refactored node collection and export:
   - Changed from `range(cubit.get_node_count()+1)` iteration to collecting actual node IDs from elements
   - Added element type iteration similar to `export_Gmsh_ver2()`
   - Sort node IDs before writing for consistency
   - Create `node_id_to_vtk_index` mapping dictionary
2. Lines 546+: All cell connectivity now uses `node_id_to_vtk_index[node_id]` instead of `node_id-1`

**Elements affected**:
- Tetrahedra (tet4, tet10)
- Hexahedra (hex8, hex20)
- Wedges (wedge6, wedge15)
- Pyramids (pyramid5, pyramid13)
- Triangles (tri3, tri6)
- Quadrilaterals (quad4, quad8)
- Edges (edge2, edge3)
- Nodes (point)

## Testing

### Before Fix
```bash
$ paraview York_cubit_mesh.vtk
# ParaView exits silently (file unreadable)
```

Verification:
```bash
$ awk 'NR>=483 && NR<=770 {for(i=2;i<=NF;i++) if($i>568) print "Line " NR ": index " $i " > 568"}' York_cubit_mesh.vtk | wc -l
# Output: 1152 (many out-of-range indices)
```

### After Fix
```bash
$ paraview York_cubit_mesh.vtk
# ParaView successfully opens and displays mesh
```

Verification:
```bash
$ awk 'NR>=483 && NR<=770 {for(i=2;i<=NF;i++) if($i>568) print "Line " NR ": index " $i " > 568"}' York_cubit_mesh.vtk | wc -l
# Output: 0 (no out-of-range indices)
```

## Impact

- **All VTK files generated by `cubit_mesh_export.export_vtk()` were affected**
- **This bug prevented ParaView from reading Cubit-generated meshes**
- **Fix enables proper visualization of FEM meshes in ParaView**

## Regenerating Fixed VTK Files

For existing projects using `cubit_mesh_export.export_vtk()`:

### Option 1: With Cubit installed

1. Update `cubit_mesh_export.py` with the fix
2. Re-run mesh generation scripts or use regeneration scripts:

   **For Cubit repository examples**:
   ```bash
   cd S:\Radia\01_GitHub
   python examples/cubit/vtk/regenerate_example_vtk_files.py
   ```

   **For Radia electromagnet example**:
   ```bash
   cd S:\Radia\01_GitHub\examples\electromagnet
   python regenerate_cubit_mesh_vtk.py
   ```

### Option 2: Without Cubit (direct Nastran to VTK conversion)

If you have a `.bdf` (Nastran) file but not Cubit:

```bash
cd S:\Radia\01_GitHub\examples\electromagnet
python convert_nastran_to_vtk.py
```

This reads `York.bdf` directly and creates `York_cubit_mesh.vtk` with correct indexing.

## Related Issues

- Radia project: `examples/electromagnet/York_cubit_mesh.vtk` regenerated
- Any other projects using `cubit_mesh_export.export_vtk()` should regenerate VTK files

---

**Fixed By**: Claude Code AI Assistant
**Date**: 2025-11-22
**Repository**: S:\Radia\01_GitHub
**Severity**: Critical (breaks ParaView visualization)
**Type**: Index mapping bug
