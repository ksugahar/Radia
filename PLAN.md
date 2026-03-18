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

### Error Investigation (3154% error)
- **Root cause**: Cubit volume mesh -> `export_NetgenMesh` creates dim=3 mesh
  with `HDivSurface` having ndof=2448 (all edges) but FreeDofs=870 (boundary only)
- The full matrix includes interior DOFs with zero interaction -> L matrix is rank-deficient
- **Working code** (`verify_inductance.py`): uses Netgen OCC surface mesh (dim=2),
  `HDivSurface` ndof = surface edges only, no interior DOFs
- **Solution**: Cubit -> surface mesh export as dim=2 Netgen mesh (not dim=3 volume)
  OR: use `definedon=mesh.Boundaries(label)` + filter FreeDofs correctly

## Remaining Steps

### Step 4: Source/Sink port vector
- Gap のある導体: source 面に +1, sink 面に -1 の電流を注入
- Port vector `e[i]` は source/sink 面の edge DOFs にのみ非ゼロ
- 内部 edge の電流は BEM が自動的に決定
- 参考: `ngbem_peec.py` の `_build_source_sink_vector()`

### Step 5: Cubit block -> boundary label mapping
- Cubit block "source"/"sink" -> NGSolve boundary label
- `cubit_mesh_export.export_netgen()` でブロック名が boundary label になる
- `mesh.Boundaries("source")` で source 面の DOFs を取得

### Step 6: calc_inductance.py 書き直し
- BEMExtractor 不使用（直接 LaplaceSL を使う、シンプル）
- Surface mesh only (volume mesh 不要)
- stdout redirect で print 抑制
- JSON 出力: L, R, port info

### Step 7: GUI テスト
- Cubit で torus+gap モデル
- Tools > Radia-NGSolve > Inductance
- Source/Sink 選択 -> Extract -> 結果表示

### Step 8: Curve order 対応
- Cubit quad mesh + SetGeomInfo -> mesh.Curve(order)
- 高次要素で BEM 精度向上

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
