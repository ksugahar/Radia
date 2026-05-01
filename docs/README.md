# Radia Documentation

**Version:** 4.6.0

## API Reference

- [API_REFERENCE.md](api/API_REFERENCE.md) - Complete Python API reference: quick start, elements, solvers, NGSolve integration
- [UNIFIED_FIELD_API_DESIGN.md](api/UNIFIED_FIELD_API_DESIGN.md) - Unified `rad.Fld()` dispatch design for static and PEEC solvers

## Solver Architecture

- [SOLVER_ARCHITECTURE.md](solver/SOLVER_ARCHITECTURE.md) - Solver design philosophy and architecture overview
- [EDDY_CURRENT_METHODS.md](solver/EDDY_CURRENT_METHODS.md) - Conductor eddy current modeling: method comparison (NGSolve + ngbem)
- [ELF_COMPATIBILITY_GUIDE.md](solver/ELF_COMPATIBILITY_GUIDE.md) - Compatibility guide between Radia and ELF/MAGIC solvers
- [IMA_SYMMETRY_DESIGN.md](solver/IMA_SYMMETRY_DESIGN.md) - Image symmetry implementation for MSC hexahedra
- [NGBEM_INTEGRATION_DESIGN.md](solver/NGBEM_INTEGRATION_DESIGN.md) - Unified PEEC Loop-Star + MMM + MSC architecture with ngbem

## Kelvin Transformation

- [KELVIN_TRANSFORMATION.md](kelvin/KELVIN_TRANSFORMATION.md) - Theory, API, workflow, and references for open boundary magnetostatics

## PEEC (Partial Element Equivalent Circuit)

- [PEEC_PANEL_IMPLEMENTATION.md](peec/PEEC_PANEL_IMPLEMENTATION.md) - Panel-based 2D surface integration implementation
- [PEEC_CONDUCTOR_MODELING_GUIDE.md](peec/PEEC_CONDUCTOR_MODELING_GUIDE.md) - Conductor modeling: filament model via `coil_from_cad.py`
- [PEEC_MSC_COUPLING.md](peec/PEEC_MSC_COUPLING.md) - PEEC conductor + MSC magnetic material coupled solver
- [PEEC_SHIELD_CONDUCTOR.md](peec/PEEC_SHIELD_CONDUCTOR.md) - Shield conductor modeling (`peec_shield.py`)
- [PEEC_SURFACE_IMPEDANCE.md](peec/PEEC_SURFACE_IMPEDANCE.md) - Surface impedance formulation and SPICE export
- [PEEC_VALIDATION_PLAN.md](peec/PEEC_VALIDATION_PLAN.md) - Systematic validation plan for PEEC solver and PEEC-MSC coupling
- [NPORT_BLOCK_LANCZOS_SPICE.md](peec/NPORT_BLOCK_LANCZOS_SPICE.md) - N-port Block Lanczos algorithm for SPICE-compatible circuit extraction

## Analytical Reference Formulas

- [analytical_formulas.md](analytical_formulas.md) - Closed-form formulas (demag factor of rotational ellipsoid, AC vector locus, magnetic shielding, 2D rectangular magnet, thin plate eddy current, ...). PDF-equation -> code cross-reference for the Wakao-Igarashi-Fujiwara-Kameari series. Source: [src/radia/analytical_formulas/](../src/radia/analytical_formulas/), tests: [tests/analytical_formulas/](../tests/analytical_formulas/), examples: [examples/analytical_formulas/](../examples/analytical_formulas/).

## Visualization

- [MESH_GUIDE.md](visualization/MESH_GUIDE.md) - Mesh generation workflows (Cubit + Netgen)

## Cubit Mesh Export

- [Function_Reference.md](cubit/Function_Reference.md) - All plugin commands and Python API reference
- [Cubit_Element_Order.md](cubit/Cubit_Element_Order.md) - Element order control (1st/2nd order) in Coreform Cubit
- [export_NetgenMesh.md](cubit/export_NetgenMesh.md) - Netgen .vol export (order 1-5, recommended)
- [export_Gmsh.md](cubit/export_Gmsh.md) - Gmsh v4.1 export (order 1-3)
- [export_Nastran.md](cubit/export_Nastran.md) - Nastran BDF export (order 1-2)
- [export_vtk.md](cubit/export_vtk.md) - VTK Legacy export (order 1-2)
- [export_meg.md](cubit/export_meg.md) - ELF/MAGIC MEG export
- [export_femeem.md](cubit/export_femeem.md) - FEMEEM format export (Gifu Univ.)
- [export_exodus.md](cubit/export_exodus.md) - Exodus II export (Cubit native)
