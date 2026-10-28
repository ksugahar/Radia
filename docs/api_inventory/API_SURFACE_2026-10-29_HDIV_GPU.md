# Targeted HDiv GPU API ownership audit — 2026-10-29

| Family | Classification | Boundary |
| --- | --- | --- |
| `PreparedHDivNonlinear` | user-intent | One explicit prepared material solve owner; no implicit fallback. |
| Shared HDiv/HMatrix/BH/Newton implementation | method-keep | Radia-owned private native implementation, shared with MEX. |
| `prepare_hdiv_host_setup` | method-keep | Radia charge/material preparation; NGSolve still owns FE/mesh operations. |
| Former Simulink-owned Python solver/binary staging | deprecated-drop | Removed; callers use Radia directly, without a compatibility alias. |

No new mesh, CAD or generic solver implementation is added. The MATLAB
publisher remains a thin consumer. During development, first move ownership and
preserve kernel bytes, then rebuild and revalidate the exact active development
revision before adopting its results. This focused audit does not reclassify
unrelated package APIs or certify the experimental solver.

Scope: ownership/build preparation against main2f1dbcc0, not numerical acceptance. Next phases: exact admitted build/basic smoke, then rebuild/revalidation against the evolving development revision. Performance remains outside scope.
