# FEM Kelvin Panel - Implementation Plan

## Current Status (2026-03-29)

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

### OCC/STEP Path
- [ ] Standalone script with argparse (STEP file input)
- [ ] Netgen mesh generation from STEP
- [ ] No Cubit dependency
- [ ] Documentation for non-Cubit users

### Solver
- [ ] `Periodic(HCurl(order>1))` investigation (currently works via BDDC+BVP)
- [ ] CompactAMS for order=1 (requires Periodic-aware gradient matrix)
