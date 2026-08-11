"""Compatibility checks for names used before LTspice joined Radia."""
from __future__ import annotations


def test_legacy_top_level_imports_forward_to_radia() -> None:
    import ltspice_converter
    import radia.ltspice as canonical
    import spice_circuit_lab

    assert ltspice_converter.netlist_to_asc is canonical.netlist_to_asc
    assert spice_circuit_lab.netlist_to_asc is canonical.netlist_to_asc
    assert ltspice_converter.__version__ == canonical.__version__
    assert spice_circuit_lab.__version__ == canonical.__version__


def test_legacy_submodule_import_uses_canonical_source() -> None:
    from ltspice_converter.topology import topology_signature as legacy_signature
    from radia.ltspice.topology import topology_signature

    netlist = "R1 in out 1k\nC1 out 0 1u\n.end\n"
    assert legacy_signature(netlist) == topology_signature(netlist)
