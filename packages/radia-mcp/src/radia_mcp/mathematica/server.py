"""
Mathematica MCP Server (radia_mcp.mathematica)

A standalone FastMCP server exposing Wolfram Mathematica as MCP tools via
the `wolframscript` command-line driver.

Usage:
    mcp-server-mathematica              # Start MCP server (stdio transport)
    mcp-server-mathematica --selftest   # Run self-test (no Mathematica needed)

Tools:
    mathematica_evaluate(code, timeout)   — low-level Wolfram eval bridge
    mathematica_status()                   — diagnostic: is wolframscript OK?
    mathematica_simplify(expr, assumptions)
    mathematica_to_tex(expr)               — TeXForm for paper writing
    mathematica_check_identity(lhs, rhs)   — verify LHS == RHS
    mathematica_vector_calc(expr, op)      — curl / div / grad / laplacian / ...
    mathematica_unit_convert(quantity, target)
    mathematica_solve(equations, unknowns)
    mathematica_integrate(integrand, var, lower=None, upper=None)
    mathematica_differentiate(expr, var, order=1)

Without Mathematica installed, all tools return {"ok": False, "error": ...}
cleanly; the server itself still loads and can be probed via --selftest.

History: promoted on 2026-05-20 from
  legacy private source tree
into radia-mcp as a standalone subpackage.  The original in
mcp_server_document is left in place for backward compatibility.
"""

import sys

from mcp.server.fastmcp import FastMCP
from ..common import register_status_tool

from . import tools as _tools

mcp = FastMCP("mcp-server-mathematica")


# ============================================================
# Tool registration
# ============================================================
# Auto-register all mathematica_* callables from tools.py and helpers.py.
# helpers.py functions are re-exported through tools.py, so iterating
# tools alone is enough.
_REGISTERED: list[str] = []
for _name in dir(_tools):
    if _name.startswith("mathematica_"):
        _fn = getattr(_tools, _name)
        if callable(_fn):
            mcp.tool()(_fn)
            _REGISTERED.append(_name)


# ============================================================
# Entry point
# ============================================================



register_status_tool(
    mcp,
    server_name='mcp-server-mathematica',
    description='Mathematica recipes: vector calc, Kelvin transform, symbolic Maxwell, evaluation pipeline',
    subpackage='radia_mcp.mathematica',
    related_servers=["differential-forms", "radia-ngsolve"],
    # avoid collision with the existing `mathematica_status` tool
    # (wolframscript health probe — different from server introspection)
    tool_name="mathematica_server_status",
)


def main():
    if "--selftest" in sys.argv:
        print("Mathematica MCP server self-test:")
        print(f"  Registered tools ({len(_REGISTERED)}):")
        for name in sorted(_REGISTERED):
            print(f"    - {name}")
        # Probe wolframscript availability without requiring it.
        st = _tools.mathematica_status()
        if st["available"]:
            print(f"  wolframscript: OK at {st['path']}")
            print(f"  Mathematica version: {st['version']}")
        else:
            print(f"  wolframscript: NOT AVAILABLE ({st.get('license_msg', 'n/a')})")
            print("  (this is fine; server still loads -- tools just return ok=False)")
        print("  PASSED")
        return

    mcp.run(transport="stdio")


if __name__ == "__main__":
    main()
