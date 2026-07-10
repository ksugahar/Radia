#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Regression tests for asc-extractor fidelity fixes (oracle-driven,
validated against LTspice's own netlister).

The jumper / A-device fixes resolve symbol pins from LTspice's
``lib.zip``; the relevant tests are skipped when it is not installed.
"""
from __future__ import annotations

import pytest

from ltspice_converter.parser.asc_parser import (
    AscParser, NetlistExtractor, AsyParser,
)

_has_jumper = AsyParser.get_terminal_offsets("jumper", "R0") is not None
needs_libzip = pytest.mark.skipif(
    not _has_jumper, reason="LTspice lib.zip (jumper.asy) not available"
)


JUMPER_ASC = """Version 4
SHEET 1 880 680
SYMBOL jumper 0 0 R0
SYMATTR InstName J1
FLAG -32 64 nodeA
FLAG 32 64 nodeB
"""


@needs_libzip
def test_jumper_unions_its_two_pins():
    """A jumper symbol is a zero-ohm short: its two pins must land in the
    same net group (so the nodes they touch collapse into one)."""
    p = AscParser()
    p.parse_string(JUMPER_ASC)
    ext = NetlistExtractor(p)
    groups = ext._build_net_groups()
    a, b = (-32, 64), (32, 64)
    group_with_a = next((g for g in groups if a in g), None)
    assert group_with_a is not None
    assert b in group_with_a, "jumper pins were not merged into one net"


@needs_libzip
def test_jumper_emits_no_component():
    """The jumper itself must not appear in the netlist (LTspice drops it)."""
    p = AscParser()
    p.parse_string(JUMPER_ASC)
    netlist = NetlistExtractor(p).extract()
    assert "J1" not in netlist
    assert "jumper" not in netlist.lower()


@needs_libzip
def test_symbol_class_read_from_asy_prefix():
    """A symbol unmapped/generically-X adopts the SPICE class its own .asy
    declares via SYMATTR Prefix. A crystal (xtal, Prefix C) must netlist
    as a capacitor (C…), not a generic subcircuit (X§…)."""
    from ltspice_converter.parser.asc_parser import AsyParser
    # Confirm the fixture data the fix relies on.
    assert AsyParser.get_symbol_attr("xtal", "Prefix").upper() == "C"
    asc = """Version 4
SHEET 1 880 680
SYMBOL xtal 100 100 R0
SYMATTR InstName Y1
SYMATTR Value {Cs}
FLAG 116 100 a
FLAG 116 164 b
"""
    p = AscParser()
    p.parse_string(asc)
    netlist = NetlistExtractor(p).extract()
    # Y1 must be emitted as a capacitor class, not X.
    assert "X§Y1" not in netlist
    assert "C§Y1" in netlist or "\nCY1 " in netlist


def test_a_named_instance_keeps_a_prefix_not_x():
    """An InstName starting with 'A' (LTspice special-function device)
    must not be rewritten to 'X§A...'. Terminals resolve from nearby wire
    endpoints, so this actually reaches the prefixing guard (no .asy
    needed)."""
    # Unknown symbol whose pins are estimated from the wires around it.
    asc = """Version 4
SHEET 1 880 680
SYMBOL SpecialFunctions\\\\sample 100 100 R0
SYMATTR InstName A1
WIRE 100 100 240 100
WIRE 100 132 240 132
FLAG 240 100 in
FLAG 240 132 out
"""
    p = AscParser()
    p.parse_string(asc)
    sym = next(s for s in p.symbols if s.inst_name == "A1")
    ext = NetlistExtractor(p)
    ext.extract()  # populate node map
    line = ext._symbol_to_spice(sym)
    assert line is not None and line.lstrip().startswith("A1 "), line
    assert "X§A1" not in line
    assert "A1" not in ext._name_remap, "A-device InstName was wrongly X-prefixed"
