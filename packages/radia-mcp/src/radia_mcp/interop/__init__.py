"""radia_mcp.interop — universal CAD-MCP mesh backend.

Any CAD tool that can produce a STEP file can route its output
through this server to get a Cubit hex mesh via `cubit_mesh_auto`.
Adapters for OpenSCAD (CLI) and FreeCAD (FreeCADCmd subprocess) are
shipped; users of other CAD MCPs (Blender, Onshape, AutoCAD, …) feed
their exported STEP into `any_step_to_cubit_hex` directly.
"""
