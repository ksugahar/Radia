"""Public conversion entry points: netlist <-> schemdraw <-> ASC."""
from __future__ import annotations


def netlist_to_schemdraw(netlist: str, name: str = "circuit") -> str:
    """Convert a SPICE netlist to a runnable schemdraw Python script.

    Generates a Python script that draws the circuit using schemdraw.
    Supports R, C, L, V, I, D, BJT (NPN/PNP), MOSFET, JFET, opamp.

    Args:
        netlist: SPICE netlist text (with .end). E.g.
            'V1 in 0 AC 1\\nR1 in out 1k\\nC1 out 0 1u\\n.ac dec 20 1 100k\\n.end'
        name: Circuit name for the output file (default 'circuit').

    Returns:
        Runnable Python script that generates a PDF schematic.
    """
    from .parser.cir_to_schemdraw import CirToSchemdraw
    converter = CirToSchemdraw()
    return converter.convert_string(netlist, name)


def schemdraw_to_netlist(script: str, title: str = "circuit") -> str:
    """Convert a schemdraw Python script to a SPICE netlist.

    Executes the schemdraw script, extracts circuit topology from element
    anchors, and generates a SPICE netlist.

    Supports R, C, L, V, I, D, BJT, MOSFET, JFET, Opamp.

    Args:
        script: schemdraw Python script text (must create a Drawing).
        title: Title for the netlist (default 'circuit').

    Returns:
        SPICE netlist text (.cir format) ready for LTspice simulation.
    """
    from .parser.schemdraw_to_cir import schemdraw_script_to_cir
    return schemdraw_script_to_cir(script, title)


def netlist_to_asc(netlist: str,
                   asy_search_dirs=None) -> str:
    """Convert a SPICE netlist (.cir) to an LTspice schematic (.asc) text.

    Args:
        netlist: SPICE netlist text.
        asy_search_dirs: optional list of directory paths to search for
            vendor `.asy` symbol files. Combined with the
            ``LTSPICE_ASY_SEARCH_PATH`` env var.

    Returns:
        LTspice .asc schematic text.
    """
    from .parser.netlist_to_asc import NetlistToAsc
    converter = NetlistToAsc(asy_search_dirs=asy_search_dirs)
    return converter.convert_string(netlist)


def asc_to_netlist(asc_text: str, use_ltspice=None,
                   asy_search_dirs=None) -> str:
    """Convert an LTspice schematic (.asc) text to a SPICE netlist.

    Args:
        asc_text: LTspice .asc schematic text.
        use_ltspice: backend selection.

            - ``None`` (default): **auto** — use LTspice's own
              ``-netlist`` when ``LTspice.exe`` is installed (canonical,
              ground-truth topology and anonymous-node numbering), and
              transparently fall back to the pure-Python extractor when
              it is not. This is the most faithful option on a machine
              with LTspice and still works everywhere.
            - ``True``: force the LTspice backend (errors fall back to
              pure-Python).
            - ``False``: force the pure-Python NetlistExtractor (no
              external dependency, byte-for-byte reproducible across
              machines — use this in CI or when you need determinism).
        asy_search_dirs: optional list of directory paths for vendor
            `.asy` symbol files.

    Returns:
        SPICE netlist (.cir) text.
    """
    import tempfile
    from pathlib import Path as _Path
    from .parser.asc_parser import asc_to_netlist as _impl

    with tempfile.NamedTemporaryFile(
        mode="w", suffix=".asc", delete=False, encoding="utf-8"
    ) as tmp:
        tmp.write(asc_text)
        tmp_path = tmp.name
    try:
        _dirs = (
            [_Path(d) for d in asy_search_dirs]
            if asy_search_dirs else None
        )
        # use_ltspice=None propagates to prefer_ltspice=None, which the
        # low-level impl resolves to "auto" (LTspice if available).
        return _impl(tmp_path, prefer_ltspice=use_ltspice,
                     asy_search_dirs=_dirs)
    finally:
        _Path(tmp_path).unlink(missing_ok=True)


# Backwards-compatible aliases (mcp-server-document used these names)
circuit_netlist_to_schemdraw = netlist_to_schemdraw
circuit_schemdraw_to_netlist = schemdraw_to_netlist
