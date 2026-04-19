"""gmsh_post — Gmsh .msh post-processing tools (MSH v4.1 first-class).

Complementary to the radia-coupled `mcp-server-gmsh` (authoring help):
this module OPERATES on existing .msh files — inspect, validate,
convert between versions, write $NodeData / $ElementData view blocks
onto an existing mesh.

Primary pipeline served:
    Cubit `export mesh` (v2.2, flat) → gmsh_post_convert → v4.1 entity-
    grouped → ngsolve / radia / gmsh GUI consume cleanly.
"""
