# Cubit Inductance Panel - Development Plan

## Goal

Cubit GUI パネル: **Tools > Radia-NGSolve > Inductance**
Surface mesh + ngsolve.bem LaplaceSL でインダクタンスを抽出する。

## Architecture

```
Cubit GUI (PyQt5)                    External Python (NGSolve + Cubit)
┌─────────────────────┐              ┌─────────────────────────────────┐
│ InductanceDialog     │   QProcess   │ calc_inductance.py              │
│ - Source block combo │──────────────│ 1. NGSolve import FIRST         │
│ - Sink block combo   │   cub5 file  │ 2. cubit import (-noinit)       │
│ - Conductivity       │              │ 3. open cub5                    │
│ - Curve order        │   JSON       │ 4. export_netgen(cubit)         │
│ - Frequency          │◄─────────────│ 5. HDivSurface + LaplaceSL      │
│ - Result table       │              │ 6. Port vector + L extraction   │
└─────────────────────┘              │ 7. print(json.dumps(result))    │
                                     └─────────────────────────────────┘
```

## Completed Steps

### Step 1: HDivSurface on Cubit mesh
- `cubit_mesh_export.export_netgen(cubit)` -> NGSolve Mesh
- `HDivSurface(mesh, order=0)` -> RT0 surface current DOFs
- Verified: 888 tris -> 1332 DOFs

### Step 2: LaplaceSL matrix assembly
- **CORRECT pattern** (from ngbem_peec.py):
  ```python
  L_op = LaplaceSL(u.Trace() * ds) * v.Trace() * ds
  L_dense = extract_dense(L_op.mat, ndof)
  L_matrix = MU_0 * L_dense
  ```
- **WRONG** (hangs): `BilinearForm(LaplaceSL(...)).Assemble()`
- Verified: 870 DOFs, 32.8s extraction time

### Step 3: Uniform current -> L extraction (prototype)
- `e = ones(ndof) / ndof`
- `L_total = 1 / (e @ solve(L_matrix, e))`
- Result: pipeline works end-to-end

### Error Investigation (3154% -> 2470% error)

**Phase 1: FreeDofs mismatch (SOLVED)**
- Volume mesh (dim=3): HDivSurface ndof includes interior edges -> FreeDofs < 100%
- Surface-only mesh (El3D=0, El2D=N): FreeDofs = 100%
- Fix: Build Netgen mesh with only Element2D from Cubit surface tris

**Phase 2: Seam edge issue (INVESTIGATED)**
- Cubit normals are ALL correct (verified on brick, cylinder, sphere, torus)
- Bad DOF is on ACIS seam line (y~0.03, z~0, 2 tris sharing seam edge)
- Surface-only mesh: 1 negative diagonal DOF at seam
- STEP reimport changes seam position but doesn't eliminate it (2 neg DOFs)
- OCC-generated mesh has NO seam issues (verified: 0 neg DOFs)

**Phase 3: Surface Area SUCCESS with volume mesh + SetGeomInfo**
- STEP reimport + tetmesh + export_NetgenMesh(geometry=OCC) + SetGeomInfo
- Cylinder (R=0.5, H=1): Order 1=-0.30%, Order 2=-0.0006%, Order 3=+0.0001%
- This workflow is correct and validated for Integrate(CF(1), BND)

**Phase 4: BEM LaplaceSL with volume mesh (IN PROGRESS)**
- Same workflow as Phase 3 but with LaplaceSL instead of Integrate
- LaplaceSL operator construction + MatVec is slow (~minutes for 1323 FreeDofs)
- Need to verify neg diagonal count with STEP reimport + volume mesh
- If neg=0: BEM works on this foundation
- If neg>0: seam still affects HDivSurface edge orientation for BEM

## Fundamental Workflow (ALL tools must use this)

```
Cubit create geometry
  → STEP export
  → reset + STEP reimport heal (ACIS→OCC seam fix)
  → Cubit mesh (tet/hex/surface)
  → export_netgen(cubit, geometry=OCCGeometry(step_file))
  → SetGeomInfo (cylinder/torus/sphere)
  → mesh.Curve(order)
  → Integrate / LaplaceSL / etc.
```

**Without this workflow**:
- mesh.Curve(order) has no effect (no geometry reference)
- BEM LaplaceSL gets seam artifacts (negative L diagonal)
- High-order elements don't improve accuracy (defeats CEFC presentation purpose)

**Current violations to fix**:
- calc_volume.py Order>1: creates OCC mesh instead of using Cubit mesh (WRONG)
- calc_surface.py Order>1: same (WRONG)
- calc_inductance.py: no STEP reimport (seam breaks BEM)

## Next Steps

### Step 4: Cubit surface mesh -> dim=2 Netgen mesh [CRITICAL]

**Problem**: `export_NetgenMesh()` creates dim=3 volume mesh.
`HDivSurface` on dim=3 mesh gives ndof for ALL edges (interior + boundary).
BEM only needs boundary edges -> interior DOFs have zero L entries -> rank-deficient.

**Solution**: Create `export_surface_netgen()` in `cubit_mesh_export.py`:
1. Read surface triangles/quads from Cubit blocks
2. Build `Netgen Mesh(dim=2)` with `Element2D` only (no `Element3D`)
3. Map Cubit block names to Netgen boundary/material labels

```python
# Target API
ngmesh = cubit_mesh_export.export_surface_netgen(cubit)  # dim=2
mesh = Mesh(ngmesh)
# HDivSurface(mesh, order=0) -> ndof = surface edges only
```

**Reference**: `verify_inductance.py` creates dim=2 mesh via OCC Revolve.
Netgen `Mesh(dim=2)` + `Element2D` is the correct pattern.

**Verification**: After export, check:
- `mesh.dim == 2`
- `fes.ndof == n_edges` (no interior DOFs)
- `np.min(np.diag(L_matrix)) > 0` (positive definite)

### Step 5: Validate with closed torus (no gap)

Use `export_surface_netgen()` to export Cubit's closed torus surface mesh.
Compare L_BEM vs Neumann formula. Should match `verify_inductance.py` results.

Parameters: R=0.05, a=0.005, R/a=10 (Neumann accurate).

### Step 6: Source/Sink port vector (gap model)

- Cubit blocks "source"/"sink" -> boundary labels in dim=2 mesh
- Port vector: +1 on source edges, -1 on sink edges
- L = 1 / (e^T @ L^{-1} @ e)
- Reference: `ngbem_peec.py` `_build_source_sink_vector()`

### Step 7: calc_inductance.py 書き直し

- `export_surface_netgen(cubit)` for dim=2 mesh
- Direct LaplaceSL (not BEMExtractor)
- stdout redirect, JSON output
- Same QProcess pattern as Volume/Surface Area

### Step 8: GUI テスト

- Cubit torus+gap model
- Tools > Radia-NGSolve > Inductance
- Source/Sink block selection -> Extract -> result table

### Step 9: Curve order + quad elements

- Cubit quad mesh (pave/map scheme) + SetGeomInfo
- mesh.Curve(order) for high-order BEM
- Compare L vs order (convergence study for CEFC paper)

## Key Lessons Learned

### ngsolve.bem API
- `LaplaceSL(u.Trace() * ds)` returns PotentialOperator
- `PotentialOperator * v.Trace() * ds` returns IntegralOperator
- `.mat` で BaseMatrix を取得 (dense)
- `BilinearForm().Assemble()` は使わない（ハングする）

### Cubit メッシュサイズ制御
- `surface all size X` は曲率で自動リファインされる
- `curve all interval N` が直接的
- Surface mesh only (`mesh surface all`) で BEM 用メッシュ生成

### DOF 目安
| DOFs | Time | Status |
|------|------|--------|
| 300 | ~5s | 快適 |
| 1000 | ~30s | 許容 |
| 3000 | ~5min | 限界 |
| 10000+ | --- | H-matrix 必要 |

### Cubit Panel Two-Process Model
- NGSolve before cubit (numpy DLL 衝突回避)
- QProcess (非ブロッキング、subprocess.run はフリーズ)
- -noinit flag (.cubit 再生防止)
- JSON stdout (cubit C-level 出力と混在 -> `{` 行検索)
- _find_main_window: menu count 最大のウィンドウを選択

## Next: Current Density Visualization

**cubit.set_nodal_variable(node_ids, "J_magnitude", values)** で Cubit に直接設定。
1. calc_inductance.py が per-node 電流密度を JSON で返す
2. パネル側（Cubit Python）で set_nodal_variable を呼ぶ
3. Cubit の結果表示ボタン（虹色アイコン）で可視化
4. Exodus export にも nodal variable が含まれる

GMSH post-processing export も並行サポート（GmshPostExport 既存）。

## Presentation Target

CEFC 2026 Thessaloniki: Cubit hex mesh + BEM inductance extraction + Curve order comparison

### BEM neg diag: Deep Investigation (2026-03-19)

**Confirmed facts:**
- Cubit tri winding is correct (verified via tet 4th node outward check)
- STEP reimport reverses `surface.normal_at()` but NOT tri winding
- NormalizeNumbering (rotate to smallest vertex first) causes segfault
- brick (flat surfaces): neg=0 always — BEM works perfectly
- torus (curved): neg=1 consistently with Cubit mesh
- OCC native torus: neg=0 always — BEM works
- Save/Load .vol doesn't fix neg
- Compress() doesn't fix neg
- geometry=OCC vs no geometry: same neg count
- SetGeomInfo doesn't fix neg (it sets UV for Curve, not edge orientation)

**Remaining hypothesis:**
- Netgen's OCC mesh generator handles degenerated/seam edges specially
  (libsrc/occ/occgeom.cpp: seam detection, BRep_Tool::Degenerated)
- export_NetgenMesh doesn't replicate this special handling
- The neg DOF corresponds to an edge near the OCC seam line
- May need Netgen C++ level fix or a Python-accessible API

**Workarounds:**
1. Use brick/flat geometries (neg=0, fully working)
2. For curved: exclude neg DOFs from L matrix (1-2 DOFs only)
3. Request Netgen patch for Element2D edge orientation normalization

### Root Cause FOUND: ds(label) mismatch (2026-03-19)

- `export_NetgenMesh(geometry=OCC)` sets boundary names from OCC face names
  (e.g. 'face_0', 'face_1'), NOT from Cubit block names
- `ds('conductor')` finds no matching boundary → LaplaceSL HANGS
- `ds('face_0')` works; `ds` (no label, all boundaries) works
- This is the same issue for ALL labeled ds operations on geometry-mapped meshes
- **Fix**: Use `ds` (all boundaries) or query `mesh.GetBoundaries()` for actual names
