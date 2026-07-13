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
- summarize_measure_log(log_text) -> dict
- summarize_stepped_measure_log(log_text) -> dict

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
from .measure import (
    parse_ltspice_measure_lines,
    parse_ltspice_stepped_measure_tables,
    parse_ltspice_step_lines,
    parse_spice_scalar,
    summarize_measure_log,
    summarize_stepped_measure_log,
)
from .voltage_multiplier_gate import cockcroft_walton_stage_gate
from .monte_carlo_gate import monte_carlo_tolerance_family_gate
from .bipolar_efficiency_gate import bipolar_converter_efficiency_gate
from .bipolar_rail_gate import bipolar_rail_power_quality_gate
from .series_rlc_gate import series_rlc_complex_impedance_gate
from .noise_gate import rc_thermal_noise_psd_gate
from .three_phase_delta_gate import balanced_three_phase_delta_rl_gate

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
    "parse_spice_scalar",
    "parse_ltspice_measure_lines",
    "parse_ltspice_stepped_measure_tables",
    "parse_ltspice_step_lines",
    "summarize_measure_log",
    "summarize_stepped_measure_log",
    "cockcroft_walton_stage_gate",
    "monte_carlo_tolerance_family_gate",
    "bipolar_converter_efficiency_gate",
    "bipolar_rail_power_quality_gate",
    "series_rlc_complex_impedance_gate",
    "rc_thermal_noise_psd_gate",
    "balanced_three_phase_delta_rl_gate",
]

__version__ = "0.4.1"
