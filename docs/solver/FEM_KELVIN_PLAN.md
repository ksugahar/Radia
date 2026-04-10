# FEM Kelvin Panel - Implementation Plan

## Current Status (2026-04-09)

### SIBC Implementation (Fixed 2026-04-09)

**Interface approach (correct formulation)**:
- Workpiece volume meshed as air (nu = nu0)
- Robin BC `+jw/Z_s * u.Trace() * v.Trace() * ds("wp_surface")` on internal interface
- H_t extracted via BND integral with **tangential projection** (`|A_t|^2 = |A|^2 - (A.n)^2`)
- Validated: H_t -2.9%, P -2.3% vs BEM (copper cylinder, 7 kHz)

**Key fixes applied**:
- Robin sign: `-jw/Z_s` -> `+jw/Z_s`
- Tangential projection via `specialcf.normal(3)` (without it, H_t ~ 0)
- Boundary name: `sibc` (hole) + `wp_surface` (interface) both supported
- Kelvin detection: skipped for non-Kelvin meshes (was causing spurious center offset)

**Hole approach is wrong**: PEC baseline + Robin perturbation. Systematic +6% error.

### Cubit Path (Primary, Implemented)
- User creates model via .jou (journal) in Cubit panel
- Webcut sphere + `volume copy move nomesh` + `copy mesh surface`
- `IdentifyPeriodicBoundaries` + `Mesh(ngmesh)` re-wrap
- BDDC+BVP solver (fes_order=1: PARDISO, fes_order=2: BDDC+CG)
- Verified: L=88.91 nH (+0.4% vs analytical 88.55 nH)
- High-order: curve_order=2 + fes_order=2 verified

### OCC/STEP Path (Future, for organizations without Cubit license)
- OCC fallback via `build_occ_ih_mesh_3d()` in `calc_common.py`
- Uses `netgen.occ.Identify()` for periodic (mesher generates matching mesh)
- No Cubit dependency
- Target: standalone Python script (STEP file input, not panel subprocess)
- Use case: universities/companies without Coreform Cubit license

## Remaining Tasks

### Panel (Cubit)
- [ ] Cubit GUI integration test (actual Solve button)
- [ ] Source/sink auto-detection from coil gap geometry
- [ ] Heat coupling: P_density -> Q -> T(x,t)
- [x] SIBC fix: tangential projection + Robin sign (2026-04-09)

### OCC/STEP Path
- [ ] Standalone script with argparse (STEP file input)
- [ ] Netgen mesh generation from STEP
- [ ] No Cubit dependency
- [ ] Documentation for non-Cubit users

### Solver
- [ ] `Periodic(HCurl(order>1))` investigation (currently works via BDDC+BVP)
- [ ] CompactAMS for order=1 (requires Periodic-aware gradient matrix)
