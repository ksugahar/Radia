"""
cubit-mesh-export: Cubit plugin binaries and mesh consistency checking.

The Cubit plugin (cubit_mesh_export.ccm) provides mesh export commands:
  export netgen "model.vol" order N   -- High-order curved .vol
  export gmsh "model.msh"              -- GMSH v4.1 raw data
                                         + model.geo launch companion

Consistency checking (does NOT require Cubit):
  check-vol model.vol --strict-labels
  from cubit_mesh_export.check import check_consistency  # API

cubit-plugin-install deploys plugin binaries to Cubit.
"""

__version__ = "0.14.12"

# Compatibility window with the main radia package. The Cubit plugin
# binaries bundled here (.ccm/.pyd; .ccl was removed in radia 4.80.0)
# are rebuilt alongside radia's Cubit toolbar, calc_*.py scripts, and notebook
# workbench layer. A mismatch is reported by cubit-plugin-install verification.
COMPAT_RADIA_MIN = "4.5.0"
COMPAT_RADIA_MAX = "4.999.999"  # bumped when we cut the next radia minor

__all__ = ["__version__", "COMPAT_RADIA_MIN", "COMPAT_RADIA_MAX"]
