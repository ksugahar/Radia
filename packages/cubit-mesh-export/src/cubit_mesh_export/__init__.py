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

__version__ = "1.0.4"

# Optional interoperability window for installations that also use Radia.
# The standalone Cubit plugin and MCP do not import Radia; a mismatch is
# reported only when explicitly requested by cubit-plugin-install.
COMPAT_RADIA_MIN = "4.5.0"
COMPAT_RADIA_MAX = "5.0.0"  # exact Radia 5 release candidate; no open-ended major range

__all__ = ["__version__", "COMPAT_RADIA_MIN", "COMPAT_RADIA_MAX"]
