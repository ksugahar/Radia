# FEM Kelvin Panel - Implementation Plan

## Current Status (2026-04-11)

### HCurl A-Formulation (calc_fem_kelvin.py)

The primary solver uses **HCurl** (vector potential A) with Kelvin transformation
for open boundary induction heating problems.

**Kelvin modulation (HCurl)**:
- Bilinear: `nu_cf * curl(u) * curl(v) * dx(bonus_intorder=4)`
- Kelvin factor: `nu' = (r'/R)^2 * nu0` (reluctivity, NOT permeability)
- Source: `J * v * dx("coil")` -- no modulation needed (J=0 in Kelvin domain)
- GND: optional for HCurl (gauge regularization provides uniqueness)

### SIBC Implementation (Fixed 2026-04-09)

**SIBC = Robin BC on conductor surface (hole approach)**:
- Workpiece subtracted from mesh (hole); conductor interior NOT solved
- Robin BC `+jw/Z_s * u.Trace() * v.Trace() * ds("sibc")` on hole boundary
- H_t extracted via BND integral with **tangential projection** (`|A_t|^2 = |A|^2 - (A.n)^2`)
- Z_s for solid cylinder: `rho * gamma * I1(ga)/I0(ga)` (cylindrical Bessel)

**Validated (2026-04-14, 2D axisymmetric Kelvin)**:
Full-resolution vs SIBC (hole + Robin): L < 1%, P < 2% for Cu/Steel/Al,
R/delta = 3 to 160.  Script: `validation_test/eddy_current_analytical_validation/reference_2d_axisym.py`

**Key fixes applied**:
- Robin sign: `+jw/Z_s` (positive)
- Tangential projection via `specialcf.normal(3)` (without it, H_t ~ 0)
- Kelvin detection: skipped for non-Kelvin meshes

### Cubit Path (Primary, Implemented)

- User creates model via .jou (journal) in Cubit panel
- 2-sphere architecture: interior (coil+air) + exterior (Kelvin), same radius R
- Webcut both spheres + `copy mesh surface` for 1:1 node correspondence
- `export netgen` writes periodic identification as translation
- PARDISO (direct) or BDDC+BVP (iterative) solver
- Verified: L=88.91 nH (+0.4% vs analytical 88.55 nH)
- High-order: curve_order=2 + fes_order=2 verified

### GND Vertex (Updated 2026-04-11)

| Formulation | GND Required? | Cubit Command |
|-------------|---------------|---------------|
| H1 (phi, Omega) | **Essential** | `create vertex; nodeset N name "GND"` |
| HCurl (A) | Optional | Gauge reg `reg*nu0*u*v*dx` suffices |

GND vertex at exterior sphere center (maps to physical infinity).
calc_fem_kelvin.py auto-detects "GND" in boundary labels.

### OCC/STEP Path (Future, for organizations without Cubit license)

- OCC fallback via `build_occ_ih_mesh_3d()` in `calc_common.py`
- Uses `netgen.occ.Identify()` for periodic (mesher generates matching mesh)
- No Cubit dependency
- Target: standalone Python script (STEP file input, not panel subprocess)
- Use case: universities/companies without Coreform Cubit license

## Remaining Tasks

### Panel (Cubit)
- [ ] FES order 2+ / mesh refinement convergence study (BEM reference)
- [ ] GND vertex in sample journal files
- [ ] Cubit GUI integration test (actual Solve button)
- [ ] Source/sink auto-detection from coil gap geometry
- [ ] Heat coupling: P_density -> Q -> T(x,t)
- [x] SIBC fix: tangential projection + Robin sign (2026-04-09)
- [x] GND vertex documentation in MCP knowledge (2026-04-11)
- [x] HCurl nu/mu documentation in MCP knowledge (2026-04-11)

### OCC/STEP Path
- [ ] Standalone script with argparse (STEP file input)
- [ ] Netgen mesh generation from STEP
- [ ] No Cubit dependency
- [ ] Documentation for non-Cubit users

### Solver
- [ ] `Periodic(HCurl(order>1))` investigation (currently works via BDDC+BVP)
- [ ] CompactAMS for order=1 (requires Periodic-aware gradient matrix)

### MCP Knowledge
- [x] kelvin_knowledge.py: KELVIN_HCURL_3D topic added (2026-04-11)
- [x] cubit_scripting_knowledge.py: 2-sphere workflow with GND (2026-04-11)
- [x] server.py: hcurl_3d topic registered (2026-04-11)
- [ ] MCP nu/mu description update for HCurl (verify consistency)
- [ ] 100-go deploy
