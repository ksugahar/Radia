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
from ltspice_converter.measure import (  # noqa: F401
    parse_ltspice_measure_lines,
    parse_ltspice_step_lines,
    parse_spice_scalar,
    summarize_measure_log,
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
    "parse_spice_scalar",
    "parse_ltspice_measure_lines",
    "parse_ltspice_step_lines",
    "summarize_measure_log",
]
