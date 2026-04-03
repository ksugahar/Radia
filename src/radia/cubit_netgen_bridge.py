"""
Cubit -> Netgen bridge: extract mesh via cubit Python API, curve via C++.

This module re-exports from cubit_mesh_export (canonical source).
Install: pip install cubit-mesh-export  (or pip install radia[cubit])

Usage:
    import netgen  # MUST import before cubit to avoid DLL conflicts
    import cubit
    from cubit_netgen_bridge import extract_curved_mesh

    cubit.init(['cubit', '-nojournal', '-batch', '-nographics'])
    cubit.cmd('open "model.cub5"')
    ng_mesh = extract_curved_mesh(cubit, order=3)
    mesh = ngsolve.Mesh(ng_mesh)
"""

try:
    from cubit_mesh_export.bridge import extract_curved_mesh, _set_labels
except ImportError:
    raise ImportError(
        "cubit_mesh_export is required for Cubit mesh export.\n"
        "Install with: pip install cubit-mesh-export  (or pip install radia[cubit])"
    ) from None

__all__ = ["extract_curved_mesh"]
