"""
cubit-mesh-export: Cubit plugin binaries and mesh consistency checking.

The Cubit plugin (cubit_mesh_export.ccm) provides mesh export commands:
  export netgen "model.vol" order N   -- High-order curved .vol
  export gmsh "model.msh"             -- GMSH v4.1 raw data
                                        + model.geo launch companion
  export nastran_bdf "model.bdf"      -- Nastran mesh interchange

Consistency checking (does NOT require Cubit):
  check-vol model.vol --strict-labels
  from cubit_mesh_export.check import check_consistency  # API

cubit-plugin-install deploys plugin binaries to Cubit.
"""

# Import Netgen before the bundled curver on Windows. Netgen registers the
# directory containing nglib.dll with the process DLL search path; without
# that initialization, importing cubit_mesh_curver in a fresh Python process
# fails even though the exact pinned Netgen wheel is installed.
import netgen as _netgen  # noqa: F401

__version__ = "0.14.14"

# Compatibility window with the main radia package. The Cubit plugin
# binaries bundled here (.ccm/.pyd; .ccl was removed in radia 4.80.0)
# are rebuilt alongside Radia's Cubit toolbar and mesh-validation layer. A
# mismatch is reported by cubit-plugin-install verification.
COMPAT_RADIA_MIN = "4.5.0"
COMPAT_RADIA_MAX = "4.999.999"  # bumped when we cut the next radia minor

__all__ = ["__version__", "COMPAT_RADIA_MIN", "COMPAT_RADIA_MAX"]
