"""MCP server for NGSolve weak-form variational formulations (ngsolve-mcp.weak-form)."""

from __future__ import annotations

import argparse
import asyncio

from mcp.server.fastmcp import FastMCP

from radia_mcp.common import register_status_tool

mcp = FastMCP("ngsolve-mcp.weak-form")

_TOPICS = {
    "overview": (
        "ngsolve-mcp.weak-form: NGSolve variational weak-form patterns.\n"
        "Topics: magnetostatics_a, eddy_current, hdiv_flux, symbolic_bfi."
    ),
    "magnetostatics_a": (
        "Vector-potential A-formulation in HCurl:\n"
        "  a += SymbolicBFI(nu * curl(u) * curl(v))\n"
        "  f += SymbolicLFI(J * v)"
    ),
    "eddy_current": (
        "Frequency-domain A-phi / A* eddy-current weak form:\n"
        "  a += SymbolicBFI(nu * curl(u) * curl(v) + 1j * omega * sigma * u * v)"
    ),
    "hdiv_flux": (
        "Mixed HDiv-L2 flux-potential formulation:\n"
        "  Enforces normal flux continuity across material interfaces."
    ),
    "symbolic_bfi": (
        "Use BFI(grad/curl/div) with CoefficientFunction material laws\n"
        "and assemble with BilinearForm(fes, symmetric=True/False)."
    ),
}


@mcp.tool()
def weak_form_usage(topic: str = "overview") -> str:
    """Return NGSolve weak-form formulation guidance for the given topic."""
    key = topic.strip().lower()
    if key in _TOPICS:
        return _TOPICS[key]
    return f"Unknown topic '{topic}'. Available: {', '.join(sorted(_TOPICS))}"


# Registered after every tool: the status registration applies the
# annotation contract to the tools that exist at this point.
register_status_tool(
    mcp,
    "ngsolve-mcp.weak-form",
    "NGSolve weak-form variational formulations (HCurl/HDiv/H1, BilinearForm, LinearForm, SymbolicBFI).",
    "ngsolve_mcp.weak_form",
    related_servers=["radia-mcp.ngsolve", "fem", "ngsolve-mcp.gauging"],
)


def _selftest() -> None:
    tools = asyncio.run(mcp.list_tools())
    if "weak_form_usage" not in {t.name for t in tools}:
        raise SystemExit("weak_form_usage tool missing")
    print("ngsolve-mcp.weak-form selftest passed")


def main() -> None:
    parser = argparse.ArgumentParser(description="ngsolve-mcp.weak-form MCP server")
    parser.add_argument("--selftest", action="store_true", help="Verify tool registration")
    if parser.parse_args().selftest:
        from radia_mcp.common.utf8_stdout import use_utf8_stdout

        use_utf8_stdout()
        _selftest()
        return
    from radia_mcp.common.windows_stdio import run_stdio

    run_stdio(mcp)


if __name__ == "__main__":
    main()
