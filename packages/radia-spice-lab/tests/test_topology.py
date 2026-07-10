#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Tests for the node-rename-invariant topology (connectivity) check.

All circuits here are author-authored minimal netlists -- no textbook,
LTspice-bundled, or third-party derivative content (see README
distribution policy).
"""
from __future__ import annotations

from ltspice_converter import topology_signature, topology_equivalent
from ltspice_converter import cli


RC = """* RC
V1 in 0 AC 1
R1 in mid 1k
C1 mid 0 1u
.ac dec 20 1 100k
.end"""


def test_node_rename_is_invariant():
    """Renaming an internal node must not change the signature."""
    renamed = RC.replace("mid", "whatever_name")
    assert topology_signature(RC) == topology_signature(renamed)
    assert topology_equivalent(RC, renamed)[0]


def test_symmetric_passive_pin_swap_is_benign():
    """Flipping a resistor end-for-end is the same circuit."""
    swapped = RC.replace("R1 in mid 1k", "R1 mid in 1k")
    assert topology_equivalent(RC, swapped)[0]


def test_real_rewire_is_detected():
    """Moving the capacitor from `mid` to `in` is a different circuit."""
    rewired = RC.replace("C1 mid 0 1u", "C1 in 0 1u")
    assert not topology_equivalent(RC, rewired)[0]


def test_ground_is_anchored():
    """Connecting to ground is not the same as connecting to a node."""
    grounded = "V1 a 0 1\nR1 a 0 1k\n.end"
    floating = "V1 a b 1\nR1 a b 1k\n.end"
    assert not topology_equivalent(grounded, floating)[0]


def test_diode_polarity_is_significant():
    """Asymmetric elements keep pin order -- a flipped diode differs."""
    fwd = "V1 a 0 1\nD1 a b 1N4148\nR1 b 0 1k\n.end"
    rev = "V1 a 0 1\nD1 b a 1N4148\nR1 b 0 1k\n.end"
    assert not topology_equivalent(fwd, rev)[0]


def test_value_change_preserves_topology():
    """Editing only a component value must keep topology equivalent --
    this is the agent's 'I changed a value, not the wiring' guarantee.
    """
    changed = RC.replace("R1 in mid 1k", "R1 in mid 4k7")
    assert topology_equivalent(RC, changed)[0]


def test_pin_incidence_reported_on_drift():
    """When pin counts differ, the info dict surfaces both totals."""
    a = "V1 a 0 1\nR1 a 0 1k\n.end"
    b = "V1 a 0 1\nR1 a 0 1k\nR2 a 0 2k\n.end"
    equal, info = topology_equivalent(a, b)
    assert not equal
    assert info["pin_incidences_a"] == 4
    assert info["pin_incidences_b"] == 6


def test_check_text_clean_circuit_has_no_topology_warning():
    """A standard R/C/L/V circuit round-trips with topology preserved."""
    info, warn = cli.check_text(RC, "cir")
    assert any("topology" in i for i in info)
    assert not any("topology drift" in w for w in warn)


def test_empty_netlist_signature_is_stable():
    assert topology_signature(".end") == topology_signature("* just a comment\n.end")


def test_subckt_inline_params_are_not_pins():
    """X<name> <nodes> <subckt> <param=value...> -- inline parameters
    (LTspice opamps emit them) must not be mistaken for pins."""
    from ltspice_converter.parser.netlist_to_asc import NetlistParser
    np = NetlistParser().parse_string(
        "X1 inp inn vp vm out level2 Avol=1Meg GBW=10Meg Rin=500Meg\n.end")
    c = np.components[0]
    nodes = [c.node_pos, c.node_neg, c.node_ctrl, c.node_ctrl2, c.node_out]
    nodes = [n for n in nodes if n] + list(c.extra_nodes)
    assert nodes == ["inp", "inn", "vp", "vm", "out"]
    assert c.value == "level2"          # subckt name, not "Rin=500Meg"


def test_subckt_no_params_unchanged():
    """The legacy no-param form still parses last-token-as-subckt-name."""
    from ltspice_converter.parser.netlist_to_asc import NetlistParser
    np = NetlistParser().parse_string("X1 a b c mysub\n.end")
    c = np.components[0]
    assert c.node_pos == "a" and c.node_neg == "b"
    assert list(c.extra_nodes) == ["c"]
    assert c.value == "mysub"
