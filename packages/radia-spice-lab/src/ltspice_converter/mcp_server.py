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
from .learning_quality import build_balanced_learning_profile as _build_balanced_learning_profile
from .measure import summarize_measure_log as _summarize_measure_log
from .measure import summarize_stepped_measure_log as _summarize_stepped_measure_log
from .patentability import patentability_search_plan as _patentability_search_plan
from .filter_gates import sallen_key_filter_family_gate as _sallen_key_filter_family_gate
from .hysteresis_gate import hysteretic_inductor_cycle_gate as _hysteretic_inductor_cycle_gate
from .rectifier_gate import half_wave_rectifier_gate as _half_wave_rectifier_gate


mcp = FastMCP("mcp-spice-circuit-lab")

@mcp.tool()
def half_wave_rectifier_gate(vin_peak_v: float, frequency_hz: float, capacitance_f: float, load_ohm: float, vout_avg_v: float, vout_pp_v: float, diode_avg_a: float) -> dict:
    """Gate capacitor-input half-wave rectifier DC, ripple, and current balance."""
    try: return _half_wave_rectifier_gate(vin_peak_v,frequency_hz,capacitance_f,load_ohm,vout_avg_v,vout_pp_v,diode_avg_a)
    except ValueError as exc: return {"policy":"half_wave_rectifier_gate_v1","status":"invalid_input","error":str(exc)}


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
def balanced_learning_profile() -> dict:
    """Return the ten-stage equal public/source MCP learning contract."""

    return _build_balanced_learning_profile()


@mcp.tool()
def parse_measure_log(log_text: str) -> dict:
    """Parse LTspice `.measure` results from log text into a stable schema.

    Args:
        log_text: LTspice `.log` text, or the relevant scalar-result lines.

    Returns:
        Dict with schema ``radia-spice-lab.measure-log.v1``, parsed measure rows,
        `.step` rows, duplicate-name warnings, and an ``ok`` gate.  This tool
        intentionally does not read local files, so public MCP callers must pass
        the log text they intend to summarize.
    """
    return _summarize_measure_log(log_text)


@mcp.tool()
def parse_stepped_measure_log(log_text: str) -> dict:
    """Pair LTspice stepped ``Measurement`` table rows with `.step` values.

    The caller supplies log text explicitly; this public tool does not read a
    local file.  The result gate rejects missing steps, incomplete tables,
    duplicate measurement names, and non-finite values.
    """

    return _summarize_stepped_measure_log(log_text)


@mcp.tool()
def sallen_key_filter_family_gate(rows: list[dict]) -> dict:
    """Gate multiple unity-gain low-pass variants against two-pole theory."""

    try:
        return _sallen_key_filter_family_gate(rows)
    except (TypeError, ValueError) as exc:
        return {"schema": "radia-spice-lab.sallen-key-filter-family.v1",
                "status": "invalid_input", "ok": False, "error": str(exc)}


@mcp.tool()
def hysteretic_inductor_cycle_gate(
    cycle_rows: list[dict],
    expected_current_peak_a: float,
    expected_copper_energy_j: float,
    voltage_thd: float,
) -> dict:
    """Gate settled hysteresis cycles by energy, flux closure, and harmonic evidence."""

    try:
        return _hysteretic_inductor_cycle_gate(
            cycle_rows,
            expected_current_peak_a=expected_current_peak_a,
            expected_copper_energy_j=expected_copper_energy_j,
            voltage_thd=voltage_thd,
        )
    except (TypeError, ValueError) as exc:
        return {"schema": "radia-spice-lab.hysteretic-inductor-cycle.v1",
                "status": "invalid_input", "ok": False, "error": str(exc)}


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
