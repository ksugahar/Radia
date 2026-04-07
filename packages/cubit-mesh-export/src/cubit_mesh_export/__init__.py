"""
cubit-mesh-export: High-order curved mesh export from Coreform Cubit.

Exports Cubit meshes to Netgen .vol format with:
  - Arbitrary-order curving (order 1-5) via ACIS geometry projection
  - Material labels (block/entity name)
  - Boundary labels (sideset/block/entity name)
  - Edge labels (BBND, sideset on curves)
  - Companion JSON with CAD reference values for consistency checking

Usage:
    # Export (requires Cubit)
    from cubit_mesh_export import extract_curved_mesh
    ng_mesh = extract_curved_mesh(cubit, order=3)

    # Check (does NOT require Cubit)
    check-vol model.vol              # CLI
    from cubit_mesh_export.check import check_consistency  # API
"""

__version__ = "0.4.0"

from .bridge import extract_curved_mesh

__all__ = ["extract_curved_mesh", "__version__"]
