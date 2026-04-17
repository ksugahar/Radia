# Radia Documentation

**Version:** 2.3.0

## API Reference

- [API_REFERENCE.md](api/API_REFERENCE.md) - Complete Python API reference: quick start, elements, solvers, NGSolve integration
- [UNIFIED_FIELD_API_DESIGN.md](api/UNIFIED_FIELD_API_DESIGN.md) - Unified `rad.Fld()` dispatch design for static and PEEC solvers

## Solver Architecture

- [SOLVER_ARCHITECTURE.md](solver/SOLVER_ARCHITECTURE.md) - Solver design philosophy and architecture overview
- [EDDY_CURRENT_METHODS.md](solver/EDDY_CURRENT_METHODS.md) - Conductor eddy current modeling: method comparison (NGSolve + ngbem)
- [ELF_COMPATIBILITY_GUIDE.md](solver/ELF_COMPATIBILITY_GUIDE.md) - Compatibility guide between Radia and ELF/MAGIC solvers
- [IMA_SYMMETRY_DESIGN.md](solver/IMA_SYMMETRY_DESIGN.md) - Image symmetry implementation for MSC hexahedra
- [NGBEM_INTEGRATION_DESIGN.md](solver/NGBEM_INTEGRATION_DESIGN.md) - Unified PEEC Loop-Star + MMM + MSC architecture with ngbem

## PEEC (Partial Element Equivalent Circuit)

- [PEEC_PANEL_IMPLEMENTATION.md](peec/PEEC_PANEL_IMPLEMENTATION.md) - Panel-based 2D surface integration implementation
- [PEEC_CONDUCTOR_MODELING_GUIDE.md](peec/PEEC_CONDUCTOR_MODELING_GUIDE.md) - Conductor modeling guide and implementation plan
- [PEEC_MSC_COUPLING.md](peec/PEEC_MSC_COUPLING.md) - PEEC conductor + MSC magnetic material coupled solver
- [PEEC_SHIELD_CONDUCTOR.md](peec/PEEC_SHIELD_CONDUCTOR.md) - Shield conductor modeling (`peec_shield.py`)
- [PEEC_SURFACE_IMPEDANCE.md](peec/PEEC_SURFACE_IMPEDANCE.md) - Surface impedance formulation and SPICE export
- [PEEC_VALIDATION_PLAN.md](peec/PEEC_VALIDATION_PLAN.md) - Systematic validation plan for PEEC solver and PEEC-MSC coupling
- [NPORT_BLOCK_LANCZOS_SPICE.md](peec/NPORT_BLOCK_LANCZOS_SPICE.md) - N-port Block Lanczos algorithm for SPICE-compatible circuit extraction

## Visualization

- [MESH_GUIDE.md](visualization/MESH_GUIDE.md) - Mesh generation workflows for GMSH, Netgen, and Cubit

## Cubit Mesh Export

- [Function_Reference.md](cubit/Function_Reference.md) - Function reference index for the radia Cubit plugin and `radia_cubit_mesh` module
- [Cubit_Element_Order.md](cubit/Cubit_Element_Order.md) - Element order control (1st/2nd order) in Coreform Cubit
- [export_Gmsh_ver2.md](cubit/export_Gmsh_ver2.md) - Export to Gmsh format version 2.2
- [export_Gmsh_ver4.md](cubit/export_Gmsh_ver4.md) - Export to Gmsh format version 4.1
- [export_NetgenMesh.md](cubit/export_NetgenMesh.md) - Export Cubit mesh to Netgen with high-order curving support
- [export_Nastran.md](cubit/export_Nastran.md) - Export to NX Nastran bulk data format
- [export_exodus.md](cubit/export_exodus.md) - Export to Exodus II format
