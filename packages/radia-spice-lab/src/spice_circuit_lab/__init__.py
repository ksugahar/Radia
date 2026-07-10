"""spice-circuit-lab public API.

The project used to be named ``ltspice-converter``.  The legacy
``ltspice_converter`` import remains supported; this package is the new
circuit-aware public name.
"""
from __future__ import annotations

from ltspice_converter import *  # noqa: F401,F403
from ltspice_converter import __version__  # noqa: F401
from ltspice_converter.knowledge import (  # noqa: F401
    buck_seed,
    circuit_knowledge,
)

__all__ = [
    "netlist_to_schemdraw",
    "schemdraw_to_netlist",
    "netlist_to_asc",
    "asc_to_netlist",
    "topology_signature",
    "topology_equivalent",
    "circuit_knowledge",
    "buck_seed",
]
