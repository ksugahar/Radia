"""MCP server for spice-circuit-lab.

Exposes conversion, linting, topology comparison, and small circuit-design
seed helpers as MCP tools so AI agents (Claude Code, Cursor, etc.) can
author, convert, validate, and bootstrap LTspice .asc / SPICE .cir /
schemdraw circuits on demand.

Run via the console script ``mcp-spice-circuit-lab``.  The legacy
``mcp-ltspice`` command remains available.
"""
from __future__ import annotations

import sys
from typing import List, Optional

from mcp.server.fastmcp import FastMCP

from . import conversion
from . import cli as _cli
from .knowledge import buck_seed as _buck_seed
from .knowledge import circuit_knowledge as _circuit_knowledge
from .patentability import patentability_search_plan as _patentability_search_plan


mcp = FastMCP("mcp-spice-circuit-lab")


@mcp.tool()
def netlist_to_schemdraw(netlist: str, name: str = "circuit") -> str:
    """Convert a SPICE netlist to a runnable schemdraw Python script.

    Args:
        netlist: SPICE netlist text (with .end). E.g.
            'V1 in 0 AC 1\\nR1 in out 1k\\nC1 out 0 1u\\n.ac dec 20 1 100k\\n.end'
        name: Circuit name for the output file (default 'circuit').

    Returns:
        Runnable Python script that uses schemdraw to draw the circuit.
        Supported elements: R, C, L, V, I, D, BJT (NPN/PNP), MOSFET,
        JFET, opamp.
    """
    return conversion.netlist_to_schemdraw(netlist, name)


@mcp.tool()
def schemdraw_to_netlist(script: str, title: str = "circuit") -> str:
    """Convert a schemdraw Python script to a SPICE netlist.

    Args:
        script: schemdraw Python script text (must create a Drawing).
        title: Title for the netlist (default 'circuit').

    Returns:
        SPICE netlist (.cir) text ready for LTspice simulation.
    """
    return conversion.schemdraw_to_netlist(script, title)


@mcp.tool()
def netlist_to_asc(netlist: str,
                   asy_search_dirs: Optional[List[str]] = None) -> str:
    """Convert a SPICE netlist (.cir) to an LTspice schematic (.asc).

    Args:
        netlist: SPICE netlist text.
        asy_search_dirs: optional list of directory paths to search for
            vendor `.asy` symbol files (e.g. LTspiceControlLibrary,
            LTspicePowerSim). Combined with the ``LTSPICE_ASY_SEARCH_PATH``
            env var.

    Returns:
        LTspice .asc schematic text. Can be saved as a .asc file and
        opened in LTspice.
    """
    return conversion.netlist_to_asc(netlist, asy_search_dirs=asy_search_dirs)


@mcp.tool()
def asc_to_netlist(asc_text: str,
                   use_ltspice: Optional[bool] = None,
                   asy_search_dirs: Optional[List[str]] = None) -> str:
    """Convert an LTspice schematic (.asc) to a SPICE netlist.

    Args:
        asc_text: LTspice .asc schematic text.
        use_ltspice: backend selection.
            - ``None`` (default): **auto** — use LTspice's own
              ``-netlist`` when LTspice.exe is installed (canonical,
              ground-truth topology), else fall back to the pure-Python
              extractor. Best fidelity where LTspice exists; portable
              everywhere.
            - ``True``: force LTspice (falls back on error).
            - ``False``: force pure-Python (deterministic, no LTspice
              dependency).
        asy_search_dirs: optional list of vendor `.asy` search dirs.

    Returns:
        SPICE netlist (.cir) text.
    """
    return conversion.asc_to_netlist(
        asc_text, use_ltspice=use_ltspice, asy_search_dirs=asy_search_dirs,
    )


@mcp.tool()
def check_circuit(text: str, fmt: str,
                  asy_search_dirs: Optional[List[str]] = None,
                  use_ltspice: bool = False) -> dict:
    """Lint a circuit: round-trip drift + static netlist checks.

    Same logic as the ``ltspice-convert --check`` CLI command, exposed
    so AI agents can validate their own generated SPICE without
    shelling out.

    Args:
        text: file content (.asc text for ``fmt='asc'``, SPICE netlist
            for ``fmt='cir'``, Python script for ``fmt='py'``).
        fmt: one of ``'asc'``, ``'cir'``, ``'py'``.
        asy_search_dirs: optional list of vendor `.asy` search dirs.
        use_ltspice: backend for the asc round-trip extraction.
            ``False`` (default) = pure-Python on both ends, so the check
            is deterministic and measures the converter's own
            self-consistency. Pass ``True`` to validate against LTspice's
            canonical netlister instead (requires LTspice installed).

    Returns:
        Dict with keys:

        - ``ok`` (bool): True iff no warnings.
        - ``info`` (list[str]): informational messages
          (round-trip component counts, topology verdict, etc.).
        - ``warnings`` (list[str]): things the user should fix —
          component-count drift, topology drift, unparsed lines,
          orphan/undefined `.model` references, duplicate instance
          names, floating nodes, undefined ``{PARAM}`` references, etc.

    Example agent workflow: after generating a netlist, call
    ``check_circuit(netlist, 'cir')`` and refuse to ship the netlist
    if ``warnings`` is non-empty.
    """
    try:
        info, warn = _cli.check_text(text, fmt, asy_search_dirs,
                                     use_ltspice=use_ltspice)
    except Exception as e:
        return {'ok': False, 'info': [], 'warnings': [f'{type(e).__name__}: {e}']}
    return {'ok': not warn, 'info': info, 'warnings': warn}


@mcp.tool()
def info_circuit(text: str, fmt: str,
                 asy_search_dirs: Optional[List[str]] = None) -> dict:
    """Summarise a circuit: component-type counts, symbol kinds,
    `.subckt` block count, `.asy` resolution rate.

    Same logic as ``ltspice-convert --info --json``.

    Args:
        text: file content (.asc, .cir, or .py).
        fmt: one of ``'asc'``, ``'cir'``, ``'py'``.
        asy_search_dirs: optional vendor `.asy` search dirs.

    Returns:
        Dict containing (depending on fmt):

        - ``format``, ``size_bytes``
        - ``component_count``, ``component_types`` (e.g. ``{'R': 4, 'C': 2}``)
        - ``symbol_kinds`` (.asc only)
        - ``symbols_total``, ``symbols_asy_resolved`` (.asc only)
        - ``subckt_blocks``
    """
    return _cli.info_text(text, fmt, asy_search_dirs)


@mcp.tool()
def compare_topology(netlist_a: str, netlist_b: str) -> dict:
    """Check whether two SPICE netlists have the same connectivity.

    Node-rename-invariant: anonymous node renumbering (``N001`` vs
    ``net3``) and benign R/C/L pin swaps do NOT count as a difference;
    only genuine rewiring does. Use this to confirm an edit changed
    *only* what you intended -- e.g. after changing a resistor value,
    ``compare_topology(before, after)`` should return ``equivalent:
    true``. If you moved a wire, it returns ``false``.

    Args:
        netlist_a: first SPICE netlist text.
        netlist_b: second SPICE netlist text.

    Returns:
        Dict with keys:

        - ``equivalent`` (bool): True iff the two circuits are the same
          up to node renaming.
        - ``components_a`` / ``components_b`` (int): component counts.
        - ``pin_incidences_a`` / ``pin_incidences_b`` (int): total
          pin-to-node connections on each side (a quick tell for added
          or dropped pins).
    """
    from .topology import topology_equivalent
    try:
        equal, info = topology_equivalent(netlist_a, netlist_b)
    except Exception as e:
        return {'equivalent': False, 'error': f'{type(e).__name__}: {e}'}
    return {'equivalent': equal, **info}


@mcp.tool()
def circuit_knowledge(topic: str = "") -> dict:
    """Return compact circuit-design rules by topic.

    Args:
        topic: Topic hint such as ``"buck"``, ``"switching"``,
            ``"asc conversion"``, or ``"opamp"``.

    Returns:
        Dict with ``topic`` and a list of public design/checking rules.
    """
    return _circuit_knowledge(topic)


@mcp.tool()
def buck_seed(
    vin_v: float,
    vout_v: float,
    iout_a: float,
    fsw_hz: float = 100_000.0,
    ripple_fraction: float = 0.25,
) -> dict:
    """Create a first-pass asynchronous buck-converter simulation seed.

    Args:
        vin_v: Input voltage.
        vout_v: Target output voltage.
        iout_a: Target output/load current.
        fsw_hz: PWM switching frequency.
        ripple_fraction: Target inductor peak-to-peak ripple fraction
            relative to load current.

    Returns:
        Dict containing sizing calculations and an LTspice-ready SPICE
        netlist.  This is an open-loop seed, not a finished supply.
    """
    seed = _buck_seed(
        vin_v=vin_v,
        vout_v=vout_v,
        iout_a=iout_a,
        fsw_hz=fsw_hz,
        ripple_fraction=ripple_fraction,
    )
    return {"calculations": seed.to_dict(), "netlist": seed.to_netlist()}


@mcp.tool()
def patentability_search_plan(
    title: str,
    features: List[str],
    effects: Optional[List[str]] = None,
    domains: Optional[List[str]] = None,
    include_japanese: bool = True,
) -> dict:
    """Create prior-art search queries for patentability triage.

    This is a non-legal search aid for circuit or engineering ideas.  It
    prepares Google Scholar, Google Patents, J-PlatPat, and web query
    strings plus report questions.  A human or agent must still inspect
    the references before judging novelty or inventive step.
    """
    return _patentability_search_plan(
        title=title,
        features=features,
        effects=effects or [],
        domains=domains or [],
        include_japanese=include_japanese,
    )


def main() -> int:
    """Entry point for the ``mcp-spice-circuit-lab`` console script."""
    try:
        mcp.run()
    except KeyboardInterrupt:
        return 0
    return 0


if __name__ == "__main__":
    sys.exit(main())
