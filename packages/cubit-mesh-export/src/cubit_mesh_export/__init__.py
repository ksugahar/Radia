"""
cubit-mesh-export: Cubit plugin binaries and mesh consistency checking.

The Cubit plugin (cubit_mesh_export.ccm) provides mesh export commands:
  export netgen "model.vol" order N   -- High-order curved .vol
  export gmsh "model.msh"              -- GMSH v4.1 (lab standard)

Consistency checking (does NOT require Cubit):
  check-vol model.vol              # CLI
  from cubit_mesh_export.check import check_consistency  # API

cubit-plugin-install deploys plugin binaries to Cubit.
"""

__version__ = "0.14.0"

# Compatibility window with the main radia package. The Cubit plugin
# binaries bundled here (.ccm/.pyd; .ccl was removed in radia 4.80.0)
# are rebuilt alongside radia's
# register_toolbar.py and calc_*.py at each release; the Python side
# expects a matching minor series on the plugin. A mismatch triggers a
# loud warning from radia's _check_plugin_freshness at Cubit startup.
COMPAT_RADIA_MIN = "4.5.0"
COMPAT_RADIA_MAX = "4.999.999"  # bumped when we cut the next radia minor

__all__ = ["__version__", "COMPAT_RADIA_MIN", "COMPAT_RADIA_MAX"]
