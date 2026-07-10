"""spice-circuit-lab — circuit-aware SPICE/LTspice conversion tools.

The project was originally published as ``ltspice-converter``.  The old
package import remains available for compatibility.

Public API
----------
- netlist_to_schemdraw(netlist, name) -> str
- schemdraw_to_netlist(script, title) -> str
- netlist_to_asc(netlist) -> str
- asc_to_netlist(asc_text, use_ltspice=None) -> str   # None = auto (LTspice if installed)
- topology_signature(netlist) -> str            (node-rename-invariant)
- topology_equivalent(netlist_a, netlist_b) -> (bool, info)
- circuit_knowledge(topic) -> dict
- buck_seed(vin_v, vout_v, iout_a, fsw_hz=...) -> BuckSeed
- patentability_search_plan(title, features, effects, domains) -> dict

CLI / MCP server
----------------
- mcp-ltspice                (FastMCP stdio server)
"""
from __future__ import annotations

from .conversion import (
    netlist_to_schemdraw,
    schemdraw_to_netlist,
    netlist_to_asc,
    asc_to_netlist,
)
from .topology import topology_signature, topology_equivalent
from .knowledge import circuit_knowledge, buck_seed, BuckSeed
from .patentability import patentability_search_plan

__all__ = [
    "netlist_to_schemdraw",
    "schemdraw_to_netlist",
    "netlist_to_asc",
    "asc_to_netlist",
    "topology_signature",
    "topology_equivalent",
    "circuit_knowledge",
    "buck_seed",
    "BuckSeed",
    "patentability_search_plan",
]

__version__ = "0.4.0"
