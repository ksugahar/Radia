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

**Phase 2: Surface normal orientation (CURRENT ISSUE)**
- Surface-only mesh: ndof=1224, FreeDofs=100%, 1 bad DOF (negative diagonal)
- Even excluding bad DOF: L_BEM=1608nH vs L_Neumann=63nH (+2470% error)
- L matrix has negative eigenvalues -> not positive definite
- **Root cause**: Cubit `get_connectivity('tri', tid)` returns node order that
  does NOT guarantee outward-pointing normal. Netgen OCC automatically ensures
  consistent outward normals; manual Element2D construction does not.
- **Solution**: Use Cubit's surface normal to check/fix triangle orientation:
  1. `cubit.surface(sid).normal_at(centroid)` -> expected outward normal
  2. Compare with cross product of tri edges -> flip if opposite
  3. OR use `domin`/`domout` FaceDescriptor to let Netgen know orientation

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

## Presentation Target

CEFC 2026 Thessaloniki: Cubit hex mesh + BEM inductance extraction + Curve order comparison
