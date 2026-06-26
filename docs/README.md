# Radia Documentation

**Version:** 4.55.0

For installation, see the top-level [README.md](../README.md) Quick Start (covers
the pinned production install of `radia[cubit,gui]` + `radia-mcp` + `cubit-mesh-export`,
multi-user lab deploy, verify, and troubleshooting).

For release-by-release changes, see [CHANGELOG.md](../CHANGELOG.md).

> **Where is the canonical PEEC / FEM / Cubit knowledge?**
> Per CLAUDE.md "MCP Knowledge Placement Policy", the **single source of
> truth** for PEEC / FEM / Cubit / build123d / GMSH knowledge is the
> `radia-mcp` package's knowledge modules at
> [packages/radia-mcp/src/radia_mcp/](../packages/radia-mcp/src/radia_mcp/),
> not `docs/`.  `docs/` is the **academic / historical reference layer** --
> useful for citations and architectural overviews, but not the place to
> look for "how do I do X today".  For runnable how-tos and the current
> recipe for any task, query `radia-mcp` via Claude / your MCP client
> (e.g. `peec_inductance(topic="step_authoring")` returns Cubit + build123d
> recipes for auto-detect-friendly STEPs).

## API Reference

- [API_REFERENCE.md](api/API_REFERENCE.md) - Complete Python API reference: quick start, elements, solvers, NGSolve integration
- [UNIFIED_FIELD_API_DESIGN.md](api/UNIFIED_FIELD_API_DESIGN.md) - Unified `rad.Fld()` dispatch design for static and PEEC solvers

## Solver Architecture

- [SOLVER_ARCHITECTURE.md](solver/SOLVER_ARCHITECTURE.md) - Solver design philosophy and architecture overview
- [EDDY_CURRENT_METHODS.md](solver/EDDY_CURRENT_METHODS.md) - Conductor eddy current modeling: method comparison (NGSolve + ngbem)
- [IMA_SYMMETRY_DESIGN.md](solver/IMA_SYMMETRY_DESIGN.md) - Image symmetry implementation for MSC hexahedra
- [NGBEM_INTEGRATION_DESIGN.md](solver/NGBEM_INTEGRATION_DESIGN.md) - Unified PEEC Loop-Star + MMM + MSC architecture with ngbem

## Multipole-Moment MMM

- [multipole_moment_mmm/ACA_MOMENT_DESIGN.ipynb](multipole_moment_mmm/ACA_MOMENT_DESIGN.ipynb) - The production moment formulation for MMM/MSC. It should not be described with the old Yano-centered label: the contribution is the symbolic multipole-moment derivation that closes 3-DOF MMM and 5/6-DOF surface-charge elements by monopole, dipole, and residual-quadrupole conditions. This keeps matrix entries local and cheap compared with the HDiv Galerkin charge-Gram route, while retaining the open-boundary MMM workflow.

## FEEC / HDiv-type VIM (multipole-moment MMM complement)

- [hdiv_vim/README.md](hdiv_vim/README.md) - The **HDiv-type Volume Integral Method**: a symmetric FEEC H(div) demag operator `N = BᵀGB` whose loop modes are **field-null by construction** (de Rham). It complements the canonical multipole-moment MMM MSC backend with curved/high-order geometry, general FEEC loops, and symmetry-model machinery. Validated (feec 85/85): linear demag (sphere/spheroid/triaxial exact vs analytic), nonlinear (damped Newton; cube/C-yoke `<1%` vs shipped Radia; `analytic_gram` required for `div M ≠ 0`), distorted-mesh μr-independence, **curved + high-order** (`~10-30×` accuracy-per-DOF vs flat Radia), and **symmetry models** (1/2, 1/4, 1/8). The runnable layer is the radia-mcp `hdiv_vim(topic=...)` tool.
- [loop_star_breakdown.md](loop_star_breakdown.md) - The *problem* the HDiv-type VIM solves: the high-μ magnetostatic loop-mode breakdown ↔ the low-frequency EFIE/MoM breakdown (same cause, same Loop-Star remedy).

## Kelvin Transformation

- [KELVIN_TRANSFORMATION.md](kelvin/KELVIN_TRANSFORMATION.md) - Theory, API, workflow, and references for open boundary magnetostatics
- [DTN_SPECTRUM_COARSE_MESH.md](kelvin/DTN_SPECTRUM_COARSE_MESH.md) - **Coarse-mesh accuracy as a DtN-spectrum property** (Kameari's coarse-mesh demonstration, reframed). The exterior Dirichlet-to-Neumann operator `Λ_ext` has the closed-form eigenvalue ladder `−(n+1)/R` (3D) / `−n/R` (2D); the discrete `Λ_h` already lands the low multipoles on that ladder on the coarsest mesh (dipole 0.07%), and the **isolated** Kelvin open-boundary error (~0.1%) sits ~45× below the interior FEM error — readable off the operator before any solve, and separated from the interior discretisation a field-refinement study conflates. Includes the BEM `Λ_h` spectrum, the Kelvin polynomial-image / order-threshold mechanism (mode `n` exact iff FEM order `≥ n`), the 2D static-apparatus / rotating-machine cross-section, and the real two-sphere periodic-Kelvin validation. Runnable layer: `dtn_coarse_mesh(topic=...)`.

## PEEC (Partial Element Equivalent Circuit)

- [PEEC_PANEL_IMPLEMENTATION.md](peec/PEEC_PANEL_IMPLEMENTATION.md) - Panel-based 2D surface integration implementation
- [PEEC_CONDUCTOR_MODELING_GUIDE.md](peec/PEEC_CONDUCTOR_MODELING_GUIDE.md) - Conductor modeling via `coil_from_cad.py` (5-predicate classification dispatch, RMF, adaptive resampling, cap-centroid endpoint anchoring -- updated for v4.55.0)
- [PEEC_MSC_COUPLING.md](peec/PEEC_MSC_COUPLING.md) - PEEC conductor + MSC magnetic material coupled solver
- [PEEC_SHIELD_CONDUCTOR.md](peec/PEEC_SHIELD_CONDUCTOR.md) - Shield conductor modeling (`peec_shield.py`)
- [PEEC_SURFACE_IMPEDANCE.md](peec/PEEC_SURFACE_IMPEDANCE.md) - Surface impedance formulation and SPICE export
- [PEEC_VALIDATION_PLAN.md](peec/PEEC_VALIDATION_PLAN.md) - Systematic validation plan for PEEC solver and PEEC-MSC coupling
- [NPORT_BLOCK_LANCZOS_SPICE.md](peec/NPORT_BLOCK_LANCZOS_SPICE.md) - N-port Block Lanczos algorithm for SPICE-compatible circuit extraction

> **For "how do I author a STEP that PEEC can solve?"**: query
> `radia-mcp` `peec_inductance(topic="step_authoring")` for Cubit + build123d
> recipes, anti-patterns, and a 10-line build123d probe script that
> verifies a STEP is auto-detect-friendly BEFORE running the panel.
> The MCP knowledge is the runnable layer; `PEEC_CONDUCTOR_MODELING_GUIDE.md`
> is the architectural / theoretical overview.

## Cauer Ladder Network (CLN)

- [CAUER_LADDER_NETWORK.md](cln/CAUER_LADDER_NETWORK.md) - Foundational CLN (Tanimoto-Kameari method): iterative orthogonalization, Cauer-II ladder synthesis, three formulations (A-T, T-Ω, A-Φ), 2D / 3D variants, gauge / constraint options, Kelvin transformation coupling
- **Mixed Galerkin (CLN + HOIBC)** — bulk CLN Krylov modes + HOIBC surface envelope, coupled via the Schur complement. Single conductor admittance Y(s) with **no `d` parameter** and wall-band error 0.001–0.33% (geometry dependent). Scripts: `examples/mixed_galerkin/`. **Superseded the Warburg-Schur termination as of 2026-06-12** (see `memory/project_warburg_schur_deprecated_2026_06_12.md` for the history).
- [BEM_CLN.md](cln/BEM_CLN.md) - Multi-conductor BEM-CLN: per-element polarizability + integral-equation coupling for N-conductor clusters (Paper 2, IH workpiece + coils, paired transformer windings)
- [CLN_3D_CUBOID.md](cln/CLN_3D_CUBOID.md) - 3D Cu cuboid benchmark: HCurl FEM + BEM Cauer 3-way validation

## Analytical Reference Formulas

- [analytical_formulas.md](analytical_formulas.md) - Closed-form formulas covering Wakao-Igarashi-Fujiwara-Kameari Part 1-9 (IEE Japan SA / RM technical meetings, 2002-2007). Group B + C: ellipsoid demag/torque, AC vector locus, magnetic shielding, 2D rectangular magnet, thin-plate eddy current, Fabri solenoid axial field, three-phase line (triangle / planar / helical), K(k) / E(k) Hastings approximations, Gauss-Legendre. Group D (Part 6/8/9 extensions): plate Joule dissipation, AC thin-shell shielding, magnetic-shell interior fields, planar surface impedance, full Bessel cylindrical-conductor AC impedance, Gauss-Patterson nested quadrature, cuboid average B. Source: [src/radia/analytical_formulas/](../src/radia/analytical_formulas/), tests: [tests/analytical_formulas/](../tests/analytical_formulas/), notebook: [docs/analytical_formulas/analytical_formulas.ipynb](analytical_formulas/analytical_formulas.ipynb).

## Coil Design / Inverse Source

- [stream_function.md](stream_function.md) - **(ACA+)+TSVD least-norm solver** (stream function method, generalised). Kernel-agnostic field-synthesis / inverse-source solver `A phi = B` (M field points x N basis sources, M < N): TSVD-regularised pseudo-inverse accelerated by ACA+ low-rank recompression. ACA+ delegated to HACApK (`cHACApK_acaplus`); the matrix entry `A(i,j)` is a caller callback built from Radia's existing field (Biot-Savart for coils, MMM/MSC for magnets) via `radia_field_kernel`. Methods 2/3 (IEEJ SA-25-020). Source: [src/radia/stream_function.py](../src/radia/stream_function.py) + `src/core/rad_stream_function.cpp`, tests: [tests/test_stream_function.py](../tests/test_stream_function.py), examples: [examples/stream_function/](../examples/stream_function/).

## Visualization

- [MESH_GUIDE.md](visualization/MESH_GUIDE.md) - Mesh generation workflows (Cubit + Netgen)

## Cubit Mesh Export

- [Function_Reference.md](cubit/Function_Reference.md) - All plugin commands and Python API reference
- [Vol_vs_Step_Labels.md](cubit/Vol_vs_Step_Labels.md) - Mesh format routing: .vol carries labels; .step is geometry-only. Why cubit-mesh-export is the single chokepoint for labeled .vol.
- [Cubit_Element_Order.md](cubit/Cubit_Element_Order.md) - Element order control (1st/2nd order) in Coreform Cubit
- [export_NetgenMesh.md](cubit/export_NetgenMesh.md) - Netgen .vol export (order 1-5, recommended)
- [export_Gmsh.md](cubit/export_Gmsh.md) - Gmsh v4.1 export (order 1-3)
- [export_vtk.md](cubit/export_vtk.md) - VTK Legacy export (order 1-2)
- [export_Nastran.md](cubit/export_Nastran.md) - Nastran BDF export (order 1-2)
- [export_meg.md](cubit/export_meg.md) - MEG export (Gifu Univ. FEM mesh format)
- [export_femeem.md](cubit/export_femeem.md) - FEMEEM format export (Gifu Univ.)
- [export_exodus.md](cubit/export_exodus.md) - Exodus II export (Cubit native)

## Removed (historical reference)

The following features were removed; their docs were deleted to prevent confusion:

- **Nastran BDF input** (Radia-side `.bdf` mesh import): removed.  Use Netgen `.vol` (preferred) or GMSH v4.1 as the mesh interchange format into NGSolve.  Note: **Cubit-side Nastran BDF export** (`export jmag_nastran`, [export_Nastran.md](cubit/export_Nastran.md)) is unaffected and still ships.
- **Scattered-field Robin RHS** (`docs/FEM_SCATTERED_FIELD.md`): removed 2026-04-24; the formulation could not be made stable for MSC coupling.  The total-field formulation in `calc_fem_kelvin.py` is the shipped path.

See [CHANGELOG.md](../CHANGELOG.md) for the full removal history per release.
