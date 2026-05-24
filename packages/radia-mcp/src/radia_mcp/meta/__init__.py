"""radia_mcp.meta — cross-server catalog + health check.

Single entry-point for LLMs to discover the 29-server radia_mcp
ecosystem. Adopted 2026-05-24 inspired by wjc9011/COMSOL_Multiphysics_MCP
where each server has a *_status() tool for self-introspection.

For radia_mcp the problem is harder: 29 distinct servers, each with
its own focus. The meta tool returns:
  - server catalog (name, description, primary tools, related servers)
  - per-server reachability via importlib (can the subpackage
    even be imported? — surfaces broken installs early)
  - cross-server map for "which server should I call for X?"
"""
