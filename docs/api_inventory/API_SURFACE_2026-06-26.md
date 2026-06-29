# Radia Proprietary API Surface -- Inventory (Stock-Take) and Reduction Plan

Policy basis: CLAUDE.md "Reduce Proprietary API Surface -- Plumbing to netgen/ngsolve, Methods Stay".
Decision rule: *"Does netgen/ngsolve (or MKL/OCC/GMSH/Cubit) already provide this?"* YES -> plumbing-delete/delegate; NO (a method NGSolve lacks) -> keep (and demote the user-facing pybind surface over time).

This is a READ-ONLY inventory. No files were edited. Counts are per API FAMILY, not per function.
Generated 2026-06-26 by the `api-inventory` skill (`.agents/skills/api-inventory/inventory_workflow.js`).

## Update 2026-06-26: caller verification + Phase 0 executed

Caller analysis (grep across src / examples / tests) was run before any deletion
-- the inventory is advisory; callers decide. Corrections + actions:

- **`Trf*` (TrfTrsl/TrfRot/TrfInv/TrfCmbL/TrfOrnt) -> RECLASSIFIED to method-keep,
  NOT plumbing-delete.** `TrfOrnt(obj, trf)` adds a SYMMETRY COPY that contributes
  to the field (not mere CAD placement), and the family is covered by
  `tests/test_transformations.py` + `test_group_operations.py` (23 tests pass).
  OCC/.vol does not replace field-contributing symmetry replication. The original
  "delete Trf*" row below is SUPERSEDED -- keep them. (`TrfMlt` was already removed
  earlier in favor of IMA symmetry.)
- **`create_hex_mesh_grid` -> REMOVED.** Dead CplMag-era helper (CplMag itself is
  removed), zero callers. (`netgen_mesh_import.py`)
- **`MatMagFixed` / `MatMagLinear` / `MatMagCurve` -> REMOVED** (pybind surface;
  rebuilt + verified). Unused skeleton trio (all behaved as fixed M); use `MatPM`
  / direct `ObjHexahedron(..., M)`. The radentry C-API `RadMatMag*` is left in
  place. (`src/lib/radia_pybind.cpp`)
- **`create_sphere_mesh` / `create_box_mesh` (dielectric_solver) -> DEFERRED.**
  Real caller `docs/peec_integration/demos/applications/demo_dielectric.py`; migrate
  to `ngsolve.occ` before deletion.
- **`esim_vtk_export` -> DEFERRED.** Real caller `tests/test_esim_integration.py`;
  migrate the test off `ESIMVTKOutput` (to NGSolve VTKOutput) before deletion.

The body tables below are the original synthesis; this section is the
authoritative delta. Re-run the `api-inventory` skill to refresh.

## Executive summary

The radia surface is overwhelmingly **method-keep**: it exists to provide what NGSolve cannot -- analytic open-boundary field (rad.Fld), MMM/MSC + yano-MSC, HDiv-VIM (the sole VIM), axifem Henrotte basis, DtN/FEM-Kelvin open boundary, PEEC, BEM, sparsesolv preconditioners, CLN/PRIMA MOR, the analytical_formulas reference layer, stream-function coil design, and the maglev/ECB application. Mass deletion is the wrong frame.

The real reduction is **demotion, not deletion**. The 2-layer API target is already partly built: `SoftIron(...).solve(...).field(...)` over `.vol -> soft_iron_from_mesh` hides `ObjHexahedron` exactly as the policy wants. The work is to keep the element/coil primitives + containers + mesh->element bridges as INTERNAL C++/representation and un-pybind their hand-built-mesh user surface behind SoftIron / Magnet / CoilBuilder.

Genuine **plumbing-delete** is small and clean: the `Trf*` transform algebra (delegate to OCC/.vol), a few mesh-generation helpers (`create_hex_mesh_grid`, `create_sphere/box_mesh`), and visualization/CAD-kernel glue (`GmshPostExport`, `vol_sol_viewer`, `kelvin_geometry`, `step_mesh_builder`, `_b3d_shim`). **deprecated-drop** is also small and safe: the `MatMagFixed/Linear/Curve` no-op trio, `esim_vtk_export` (already de-exported), and the already-removed `CndLoop/CplMag/FldVTS/Rwg/beam_tracking` shells.

One hard constraint: **HACApK is under active development by another agent -- do-not-touch.** Categorize it (method-keep) but make no edits to `src/ext/HACApK/`, `src/core/rad_hacapk.*`, any `HLU*` entry, `HACApKBEMManager`/`HACApKPEECManager`, or `BuildHMatrix`/`MatVec`.

## Counts by bucket (families)

| bucket | family count |
|---|---|
| plumbing-delete | 11 |
| method-keep | 47 |
| method-demote | 12 |
| user-intent | 14 |
| deprecated-drop | 5 |

## plumbing-delete -- netgen/ngsolve/OCC/GMSH/Cubit already provide this

| group | members (sample) | count | rationale | replacement / action |
|---|---|---|---|---|
| Transformations Trf* | TrfTrsl, TrfRot, TrfInv, TrfCmbL/R, TrfOrnt | 6 | generic rigid-body/affine transforms | apply placement in OCC/.vol before meshing; delete transform algebra |
| Structured grid gen | netgen_mesh_import.create_hex_mesh_grid | 1 | mesh GENERATION (the named example) | NGSolve MakeStructured3DMesh / OCC / Cubit -> .vol |
| Primitive surface meshers | dielectric_solver.create_sphere_mesh / create_box_mesh | 2 | ad-hoc mesh gen for primitives | ngsolve.occ Sphere/Box -> GenerateMesh |
| Kelvin exterior geom builder | kelvin_geometry.add_kelvin_exterior_domain | ~3 | OCC offset-sphere + face naming = CAD plumbing | OCC (ngsolve.occ) / .vol with kelvin labels |
| STEP->yoke+air+Kelvin mesher | step_mesh_builder.build_mesh_from_step | ~8 | OCC geom + symmetry cut + meshing | OCC/Netgen or Cubit -> .vol -- BUT MEMORY flags KEEP (accel route); owner sign-off before delete |
| OCC coil profile/sweep glue | coil_profile_occ.coil_from_build123d_sweep | ~6 | build123d/OCC face->profile + sweep | build123d/OCC -> feed CoilBuilder; keep only thin glue |
| Round-body faceted magnets | round_bodies.cyl_mag / sphere_mag | 3 | hand-assembles ObjHexahedron geometry | Magnet(...) intent + .vol; faceting -> OCC/netgen |
| GmshPostExport viz writer | GmshPostExport, vol2msh, convert_sol_to_msh | ~40 | visualization & mesh export | GMSH .msh v4.1 owns format -- DEMOTE/keep as sole curved-export glue until NGSolve/GMSH covers it |
| .vol/.sol viewer | radia.tools.vol_sol_viewer | 2 | mesh/solution visualization | netgen native .vol viewer / NGSolve Draw / GMSH |
| OCP CAD-inspection shim | _b3d_shim.import_step/section/Face/Solid | ~15 | re-implements build123d/OCC kernel API | build123d/OCC; retire shim when import cost solved upstream |
| Cubit plugin install tooling | install_panels.py, setup_cubit.py | ~30 | deploy plumbing, not a Radia API | cubit-mesh-export's cubit-plugin-install (Tier-2 policy) |

## method-keep -- methods NGSolve lacks (KEEP)

| group | members (sample) | count | rationale | note |
|---|---|---|---|---|
| Analytic field eval | rad.Fld, FldLst, FldInt, FldFrc | ~10 | crown-jewel open-boundary analytic field/force | also user-intent entry |
| Coil-source primitives | ObjRecCur, ObjArcCur, ObjRaceTrk, ObjFlmCur | 4 | mesh-free Biot-Savart current sources | leaves under CoilBuilder |
| Background field | ObjBckg (callable) | 1 | external H_ext into the demag solve | callable form only |
| Hysteresis materials | MatPlayHysteresis, MatEnergyHysteresis, MatHys* | 7 | B-input Play/Energy hysteresis + state stepping | Egger/Hane lineage |
| Interaction-matrix probes | GetInteractMatrix, BuildMomentSystem, MomentSystemDenseRaw | ~7 | collocation MMMM / classic MMM matrix introspection | DEMOTE to internal |
| H-LU debug + self-test | HLUSelfTest*, HLUSetParallel | ~25 | H-LU diagnostics for maintained HACApK routes | keep internal-only |
| HDiv-VIM kernels | _hdiv_vim_assemble, _HDivVimHMatrix, solve_nonlinear_picard | ~25 | FEEC H(div) RT demag (the sole VIM) | underscore-internal; via radia.vim |
| BEM / equivalence-source | _AssembleSLDL_Galerkin, _EquivalenceSource*, HACApKBEMManager | ~12 | HACApK Laplace Galerkin BEM | manager build/matvec **do-not-touch** |
| Cuboid analytic helpers | _average_B_in_box, _average_demag_tensor | 2 | closed-form cuboid avg-B / demag tensor | analytical reference |
| Stream-function ACA-TSVD | _stream_aca_tsvd, StreamTSVD, RegularizedTSVD, aca_tsvd | ~8 | SF coil-design method | radia.stream_function |
| PEEC C++ + circuit solvers | PEECBuilder, PEECCircuitSolver, CoupledPEECSolver, PEECHACApKSolver | ~10 | FastImp-style L,R,C,M extraction | HACApK backend **do-not-touch** |
| PEEC MOR (PRIMA/Lanczos) | PRIMAHACApKModel, LanczosReducer | ~8 | model-order reduction | radia.lanczos_reduction |
| PEEC proximity/shield/import | solve_proximity_iterative, add_shield_mesh, surface_mesh_to_peec | ~8 | PEEC applications | mesh src via .vol already |
| FastHenry parser | FastHenryParser | 1 | PEEC interop format | -- |
| PEEC AC field bundles | FilamentBundleAC, WorkpieceSurfaceBundle, delta_L_telegen | ~5 | complex AC Biot-Savart + Telegen DeltaL | C++ kernels |
| Verilog-A/SPICE gen | VerilogAGenerator | ~4 | circuit-extraction deliverable | standard Verilog-A, not a wrapper |
| Scalar BEM-SIBC | ScalarBIESIBCSolver, scalar_bie_sibc | ~3 | scalar BIE + SIBC eddy current | uses ngsolve.bem operators |
| ESIM cell + table | ESIMCellProblemSolver, ESITable, generate_esi_table_from_bh_curve | ~6 | 1D nonlinear SIBC cell problem | re-exported (IH) |
| ESIM workpiece/coupled | ESIMWorkpiece, ESIMCoupledSolver, WPTCoupledSolver | ~10 | IH/WPT coupling | re-exported |
| Kelvin A-potential solvers | kelvin_solver.solve_*_A_kelvin, kelvin_material | ~12 | open-boundary FEM (DtN/Kelvin core) | built on ngsolve, KEEP |
| Kelvin pullback support | kelvin_source.kelvin_pullback_vector, *_nu_factor_*_cf | ~35 | differential-forms machinery of FEM-Kelvin | internal method support |
| Kelvin validate / identify | kelvin_validate.*, kelvin_identify_ngsolve.add_kelvin_identification | ~6 | cross-check + post-hoc Periodic Identify | kelvin-identify skill |
| Reduced potential solvers | ScalarPotentialSolver, VectorPotentialSolver | 2 | Radia+NGSolve weak coupling | the complement-NGSolve method |
| Dielectric/unified IE | DielectricSolver, UnifiedSurfaceSolver | ~6 | surface IE (loop-star/dielectric BEM) | (mesh helpers split to plumbing) |
| Equivalence source | NearFieldSource, reconstruct_static_H | ~3 | Huygens surface source | C++ rad_equivalence_source |
| Cohomology / de Rham | cohomology.betti_numbers, cohomology_basis, CohomologyCutSolver | ~8 | gmsh-free T-Omega cohomology | no NGSolve equivalent |
| Infinite-element / DtN | infinite_element.add_exterior_ie, dtn_surface_matrix | ~5 | static IE/DtN open boundary | DtN/Kelvin core |
| Biot-Savart filament field | biot_savart.h_filament, a_segments_batch | ~8 | mesh-free analytic coil source | coil_builder family |
| Analytic magnet closed-forms | SphericalMagnet, CuboidMagnet, CylindricalMagnet, RingMagnet | ~12 | exact magnet fields | analytical reference layer |
| IMA image-method | ima_field.parse_image_spec, add_ima_images | ~4 | Image Method of Analysis | backs Solve(image=) |
| Clebsch/SF potential | clebsch_potential.AxisymStreamFunctionSolver | 2 | SF/Clebsch coil design | stream-function family |
| Energy/Play hys models | energy_play_model.EnergyBasedPlayModel | 2 | nonlinear-material method | mirrors C++ Mat*Hysteresis |
| Hysteresis I/O + fit | hysteresis_io.build_shape_functions, convert_play_to_energy | ~9 | play-operator identification | material-method support |
| Coil-from-CAD extraction | coil_from_cad.extract_centerline, build_peec_from_path | ~20 | STEP->centerline/filament method | OCC kernel delegated; STEP-only policy |
| Coil topology classifier | coil_topology.extract_coil_topology, detect_cap_faces | ~6 | cap-aware spine reasoning | -- |
| radia.vim package | build_demag, hdiv_demag_solve, DemagOperator, soft_iron_from_mesh/_from_vol | ~40 | FEEC HDiv-VIM, the sole VIM | C++ rad_hdiv_vim.cpp |
| radia.bem package | sibc_hacapk, coil_inductance_ngsolve | ~8 | HACApK Galerkin BEM (+ ngsolve.bem wrap) | sibc_hacapk HACApK = do-not-touch |
| radia.open_boundary | eddy_dtn, cauer_ladder, kelvin_fem_radial_dtn, steklov_spectrum | ~18 | exact DtN open boundary + Cauer ladder | DtN/FEM-Kelvin core |
| radia.maglev | mixed_galerkin.alpha, ecb.lorentz, simulink.export | ~20 | levitation/ECB application method | radia.<domain> |
| radia.analytical_formulas | ellipsoid, shielding, conductor_impedance, cuboid_average_field | ~13 | closed-form reference layer | validation ground-truth |
| axifem.pyd | H1Henrotte, AxiHenrotteStiffnessBFI, AxiHenrotteSigmaMassBFI | ~8 | Henrotte axisymmetric magnetic FE | named keep example |
| cln_core.pyd | CLN transient ROM | ~6 | Cauer Ladder Network ROM | consumed by maglev |
| peec_matrices.pyd | PyPEECBuilder, PEECMatrices, MutualInductance | ~12 | PEEC L,R,C,M assembly | C++ rad_peec_matrices |
| sparsesolv_ngsolve.pyd | CompactAMSPreconditioner, COCRSolver, SparseSolvSolver | ~6 | HYPRE-free TaskManager AMS/AMG/COCR | src/ext/sparsesolv |
| radia_ngsolve.create_voxel_cf | create_voxel_cf | 1 | Radia-field -> ngsolve.VoxelCoefficient | correct delegation, keeps Radia field |

## method-demote -- KEEP the method, un-pybind / hide behind the intent layer

| group | members (sample) | count | rationale | action |
|---|---|---|---|---|
| Geometry primitives | ObjHexahedron, ObjTetrahedron, ObjWedge, ObjRecMag, ObjMltExt*, ObjCylMag | ~11 | internal element repr behind SoftIron/Magnet | un-pybind hand-built-mesh surface; keep as C++ repr |
| Container/introspection | ObjCnt, ObjM, ObjSetM, ObjDpl, ObjGeoVol, ObjDegFre | ~10 | assembly glue tied to the proprietary repr | demote once intent objects own assembly |
| mesh->element bridge | netgen_mesh_to_radia, create_radia_hexahedron/tetrahedron/wedge | ~8 | internal .vol->element conversion | keep behind soft_iron_from_mesh; un-expose create_radia_* |
| Round-body builders | round_bodies.cyl_mag/sphere_mag | 3 | faceted-primitive shape helper | prefer analytical_magnet (exact) or .vol |
| Coil parallel geometry | CoilBlock, CoilArc, CoilLoop, CoilRacetrack, CoilAssembly | 5 | overlaps CoilBuilder | demote behind CoilBuilder/CoilSpec |
| Coil profiles | RectProfile, CircleProfile, AnnularProfile, PolygonProfile | 6 | supporting representation | demote to CoilBuilder-facing |
| Interaction-matrix probes | GetInteractMatrix, MomentSystemDenseRaw, HMatrixDensify, BuildMatrix | ~8 | raw matrix dumps, not user API | internal/research only |
| ESIM solver internals | AndersonAccelerator, HantilaBIESolver, ESIMMultiportSolver | ~8 | accelerators/variants not in __init__ | reach via ESIMCoupledSolver/IH intent |
| Equivalence-source class | NearFieldSource | ~3 | keep as solver class | shrink direct pybind later |
| _heat_panel .vol parsers | _parse_vol_bcnames, _parse_vol_materials | ~4 | duplicate netgen label reading | use Mesh.GetMaterials/GetBoundaries; panel stays |
| analysis.py facade | UnifiedAnalysis, PEECAnalysisSolver, MMMAnalysisSolver | ~25 | orchestration over PEEC/MMM, no own kernel | candidate internal orchestration layer |
| HLU/PEEC self-tests | HLUSelfTest* (~25), _TestPEECHACApKSanity | ~26 | self-tests shipped on pybind surface | move to tests/ -- but **owned by HACApK agent, do-not-touch** |

## user-intent -- the intended USER layer (KEEP / PROMOTE)

| group | members (sample) | count | note |
|---|---|---|---|
| SoftIron intent object | SoftIron(geometry, mu_r=/bh_table=).solve().field() | 1 | THE canonical soft-iron user entry; primitives hide behind it |
| 2-layer selector wrappers | Solve(demag_backend=), ObjCnt, SolverConfig, set/get_demag_backend, UtiDelAll | 6 | backend = a switch, not two APIs (IMA via image=) |
| Field eval | Fld, FldLst, FldInt, FldFrc (+ tuning FldCmpPrc/FldLenTol) | ~10 | primary user field-eval entry |
| RadiaField CF bridge | RadiaField, as_voxel_cf, PrepareCache | ~5 | Radia field -> NGSolve coupling |
| Materials (constitutive) | MatLin, MatSatIsoTab, MatSatAniso, MatSatLamTab, MatPM, MatMvsH, MatApl + EMMaterial | ~12 | the materials part of the intent layer |
| Solve/config | Solve, SolveNonl, SolverConfig, GetSolveStats | 6 | demag solve user entry |
| Utility | UtiDelAll, UtiVer | 3 | mandatory cleanup + version |
| CoilBuilder | CoilBuilder, StraightSegment, ArcSegment, Loft* + CoilSpec/build_coil_from_spec | ~8 | mesh-free coil source intent API |
| Coil convenience | circular_loop, helmholtz_pair, solenoid | ~8 | simple-shape constructors |
| stream_function.py | StreamTSVD, aca_tsvd, solve | ~7 | SF method user wrapper |
| streamfunction_volume.py | design_volume_coil, extract_wires | ~12 | SF volume driver |
| PySide6 panels | radia_ih, radia_em, radia_motor, radia_pcb, radia_streamfunction | 6 | Layer-3 GUI user layer |
| Panel base + notebooks | ModePanel, AnalysisWindow, CommandWorkbench, IHWorkbench | ~15 | GUI/notebook infrastructure |
| Design specs + IH optimizer | EMDesignSpec, IHDesignSpec, IHOptimizer, IHWorkpieceContext | ~11 | application-domain config/workflow |

## deprecated-drop -- dead / no-op / already-removed

| group | members (sample) | count | rationale | action |
|---|---|---|---|---|
| PM no-op material trio | MatMagFixed, MatMagLinear, MatMagCurve | 3 | all behave as fixed M today (no real demag) | fold into Magnet(M=)/MatPM until full PM demag lands |
| ESIM VTK export | ESIMVTKOutput, export_esim_*_vtk | ~4 | viz export, already NOT re-exported | NGSolve VTKOutput + GmshPostExport |
| Removed conductor/viz shells | CndLoop, CndRecBlock, CplMag*, Rwg*, FldVTS, beam_tracking | ~8 | already gone (NOTE comments only) | peec_topology/coupled; GmshPostExport; Xsuite |
| OCP CAD shim (cold) | _b3d_shim (import_step/Face/Solid/...) | ~15 | OCC-kernel re-impl as perf shim | build123d/OCC; retire when upstream import is fast |
| Cubit install tooling | install_panels.py, setup_cubit.py | ~30 | deploy glue, not an API | cubit-mesh-export cubit-plugin-install |

## Phased reduction plan

Constraint for ALL phases: **do NOT touch HACApK** -- `src/ext/HACApK/`, `src/core/rad_hacapk.*`, every `HLU*`, `HACApKBEMManager`/`HACApKPEECManager`, `BuildHMatrix`/`MatVec`, and the HACApK backends of `peec_hacapk_solver.py` / `bem.sibc_hacapk` / `_ChargeGramHMatrix`. These stay method-keep, owned by another agent.

### Phase 0 -- free, safe deletions (low risk, clear replacement)
- Delete the `Trf*` family (radia_pybind.cpp:4364-4406); document OCC/.vol placement as the path.
- Delete `create_hex_mesh_grid`, `dielectric_solver.create_sphere_mesh/create_box_mesh` -> point to MakeStructured3DMesh / ngsolve.occ.
- Drop `MatMagFixed/Linear/Curve` (no-op trio) -> route to MatPM / Magnet(M=).
- Remove the dead `Cnd*/CplMag*/Rwg*/FldVTS/beam_tracking` NOTE comments and any residual references.
- Stop re-exporting / retire `esim_vtk_export` callers (already de-exported).

### Phase 1 -- delegate CAD/viz plumbing
- Replace `kelvin_geometry`, `coil_profile_occ` geometry with OCC/build123d glue feeding the intent layer.
- Replace `_heat_panel` `.vol` text parsers with `Mesh.GetMaterials()/GetBoundaries()`.
- Demote `_b3d_shim` to a documented internal perf shim; plan removal when build123d cold-start is fast.
- Keep `GmshPostExport` and `vol_sol_viewer` for now (sole curved-export / viewer glue); revisit when NGSolve/GMSH cover curved export. Get owner sign-off before removing `step_mesh_builder` (MEMORY: KEEP for accel route).

### Phase 2 -- consolidate the coil/geometry user API onto CoilBuilder
- Make `CoilBuilder` + `CoilSpec` the single coil entry; demote `CoilBlock/Arc/Loop/Racetrack/Assembly` and the `*Profile` classes behind it.
- Demote `round_bodies` faceting; steer field use to `analytical_magnet`/`cylindrical_magnet` (exact) or `.vol`.

### Phase 3 -- un-pybind the primitive/container layer behind intent objects
- Introduce/confirm `Magnet(...)` as the PM intent twin of `SoftIron`.
- Make `SoftIron` / `Magnet` / `CoilBuilder` + `.vol -> soft_iron_from_mesh/_from_vol` the ONLY documented geometry path; demote `ObjHexahedron/Tetrahedron/Wedge/RecMag/...`, `ObjCnt/ObjM/ObjSetM`, and `netgen_mesh_to_radia`/`create_radia_*` to internal representation (un-pybind gradually -- CoilBuilder/panels/examples migrate first, then remove).

### Phase 4 -- tidy debug/self-test surface (coordinate with owners)
- Move `HMatrixDensify`, `GetInteractMatrix`, moment-system probes, and the `*SelfTest*`/`_TestPEECHACApKSanity` entries off the shipped pybind surface into tests/ -- but the HLU/HACApK ones are owned by the active-dev agent; hand off rather than edit.

### Target end-state user layer
`SoftIron` / `Magnet` / `CoilBuilder` (+ `CoilSpec`) / `rad.Fld` / `rad.Solve` (+ `SolverConfig`, `set_demag_backend`, `image=`) / materials (`MatLin`, `MatSatIsoTab`, `MatPM`, `MatPlayHysteresis`, `EMMaterial`) / `RadiaField` CF / `UtiDelAll` / `UtiVer`, plus the named method modules `radia.vim | bem | open_boundary | axifem | stream_function | maglev | analytical_formulas | peec | sparsesolv_ngsolve | cln`. Everything else becomes internal kernel/representation.

---
*Read-only inventory. Grounded in src/lib/radia_pybind.cpp, src/radia/\*.py, and the shipped .pyd set (axifem, cln_core, peec_matrices, sparsesolv_ngsolve, _radia_pybind; note mmm_core.pyd is listed in CLAUDE.md but the MMM/MSC kernel actually ships inside _radia_pybind.pyd). No edits made.*
