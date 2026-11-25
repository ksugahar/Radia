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
        "Topics: magnetostatics_a, eddy_current, sibc_robin, ih_reference, hdiv_flux, symbolic_bfi."
    ),
    "magnetostatics_a": (
        "Vector-potential A-formulation in HCurl:\n"
        "  a += SymbolicBFI(nu * curl(u) * curl(v))\n"
        "  f += SymbolicLFI(J * v)"
    ),
    "eddy_current": """# Frequency-domain A-phi weak form

Use exp(+i omega t) and scaled phi = electric potential / (i omega), so
E = -i omega (A + grad(phi)). Use HCurl for A and an H1 scalar space on
conductors for phi. Fix the scalar gauge on each independent conductor and
handle the HCurl gradient nullspace consistently with the port constraints.

With (A, phi) and (v, q) from the compound trial and test spaces, use the
current NGSolve integrator syntax (nu = inverse permeability, sigma = conductivity):

    a = BilinearForm(fes)
    a += nu * InnerProduct(curl(A), curl(v)) * dx
    a += 1j * omega * sigma * InnerProduct(A + grad(phi), v + grad(q)) * dx("conductor")

Use the supplied conductor region selector in place of the example label.
Use one consistent test-space convention throughout. The q test enforces
charge conservation; specify terminal and insulating boundary conditions
before choosing source terms. A-only variants require separate assumptions.
""",
    "sibc_robin": """# SIBC as a Robin term in FEM

In a frequency-domain FEM weak form, integrate the impedance contribution
over the selected boundary using tangential traces. For an A-only form:

    a += alpha * InnerProduct(A.Trace(), v.Trace()) * ds("sibc")

Here alpha is the coefficient after eliminating the exterior tangential
magnetic field with the chosen SIBC. Derive its sign and frequency factor
from the time convention, normal orientation and E-to-A relation; do not
substitute Zs itself. For A-phi, include the tangential scalar gradient where
the impedance law uses the total tangential electric field. Verify that the
chosen convention gives nonnegative dissipated power for a passive surface.""",
    "ih_reference": """# Independent NGSolve FEM reference for induction heating

For an axisymmetric reference, use the azimuthal A formulation with the
axisymmetric integration weight and the correct axis condition. For a fully
three-dimensional coil use HCurl A-phi with conductor and terminal constraints.
Resolve the skin layer for a volume reference; document the assumptions and
impedance coefficient if using the independent SIBC alternative.

Match geometry, frequency, material properties and excitation conventions
before comparing. Report complex coil impedance from V/I and time-average
workpiece power from the conductor Joule integral or the consistent surface
loss integral. Record whether amplitudes are peak or RMS. Check power balance,
gauge and mesh convergence independently. This topic is a reference workflow
description, not a numerical qualification or a production Radia route.""",
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
