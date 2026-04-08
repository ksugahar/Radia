"""
cubit-mesh-export: Cubit plugin binaries and mesh consistency checking.

The Cubit plugin (radia_cubit.ccm) provides mesh export commands:
  radia_export netgen "model.vol" order N    -- High-order curved .vol
  radia_export gmsh "model.msh" version 2   -- GMSH v2.2
  radia_export gmsh "model.msh" version 4   -- GMSH v4.1

Consistency checking (does NOT require Cubit):
  check-vol model.vol              # CLI
  from cubit_mesh_export.check import check_consistency  # API

cubit-plugin-install deploys plugin binaries to Cubit.
"""

__version__ = "0.5.0"

__all__ = ["__version__"]
