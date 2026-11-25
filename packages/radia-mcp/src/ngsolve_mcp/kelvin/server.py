"""MCP server for Kelvin transformation open-boundary methods (ngsolve-mcp.kelvin)."""

from __future__ import annotations

import argparse
import asyncio

from mcp.server.fastmcp import FastMCP

from radia_mcp.common import register_status_tool
from radia_mcp.meta.catalog import CATALOG

mcp = FastMCP("ngsolve-mcp.kelvin")

_TOPICS = {
    "overview": (
        "ngsolve-mcp.kelvin: Kelvin transformation for unbounded magnetostatic/eddy-current domains.\n"
        "Topics: inversion, jacobian, interface_coupling, verify_first."
    ),
    "inversion": (
        "Spherical/ellipsoidal Kelvin inversion maps exterior r > R0 to interior r' = R0^2 / r,\n"
        "preserving harmonicity with scaled material weight (R0 / r')^2."
    ),
    "jacobian": (
        "In NGSolve, apply the Jacobian-weighted bilinear form on the inverted region\n"
        "so that exterior energy and far-field decay match exact open-boundary conditions."
    ),
    "verify_first": """# Verify-first checks for a periodic Kelvin pair

2. Periodic slave DOFs: build the underlying HCurl space with the intended
   essential constraints, then wrap it with Periodic. Compare FreeDofs().NumSet()
   before and after. A nonempty Kelvin identification should remove free slave
   DOFs in this setup. Confirm expected labels, pairing, orientation and phase
   before proceeding; a lower count alone does not certify the pairing.

3. Interior/exterior ratio: set a nonzero tangential test field on kelvin_int in
   the periodic space and evaluate its squared trace norms on kelvin_int and
   kelvin_ext with boundary integration (BND). For a complex trace use the
   conjugate in the squared magnitude. Reject zero interior norm. Compare
   norm_ext / norm_int to one for the matching unit-phase Kelvin interface.
   The existing IH verify-first example uses abs(ratio - 1.0) < 1e-3; do not
   generalize this criterion to interfaces with different weighting or phase.""",
    "interface_coupling": (
        "Match tangential trace (HCurl) or normal flux (HDiv) continuously across r = R0\n"
        "between the physical interior mesh and the Kelvin-inverted exterior mesh."
    ),
}


@mcp.tool()
def kelvin_usage(topic: str = "overview") -> str:
    """Return Kelvin transformation guidance for the given topic."""
    key = topic.strip().lower()
    if key in _TOPICS:
        return _TOPICS[key]
    return f"Unknown topic '{topic}'. Available: {', '.join(sorted(_TOPICS))}"


# Registered after every tool: the status registration applies the
# annotation contract to the tools that exist at this point.
register_status_tool(
    mcp,
    "ngsolve-mcp.kelvin",
    "Kelvin transformation open-boundary exterior inversion for NGSolve/Radia.",
    "ngsolve_mcp.kelvin",
    related_servers=CATALOG["ngsolve-mcp.kelvin"]["related"],
)


def _selftest() -> None:
    tools = asyncio.run(mcp.list_tools())
    if "kelvin_usage" not in {t.name for t in tools}:
        raise SystemExit("kelvin_usage tool missing")
    print("ngsolve-mcp.kelvin selftest passed")


def main() -> None:
    parser = argparse.ArgumentParser(description="ngsolve-mcp.kelvin MCP server")
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
