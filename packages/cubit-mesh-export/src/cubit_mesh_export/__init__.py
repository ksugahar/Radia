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

__version__ = "2.0.0"

__all__ = ["__version__"]
