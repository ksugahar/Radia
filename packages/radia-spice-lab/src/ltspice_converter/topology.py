"""Node-rename-invariant topology (connectivity) signatures for SPICE
netlists.

Why this exists
---------------
The round-trip metrics this package historically reported -- *component
count* and *GND-pin position* -- are necessary but not sufficient.  A
circuit can keep every component and every ground connection yet have
its internal wiring scrambled: a multi-pin vendor symbol whose ``.asy``
is not on the search path comes back from ``.asc -> .cir -> .asc`` with
a different (often inflated) pin list.  Count says "100 % clean"; the
circuit is silently rewired.

``topology_signature`` collapses a netlist to a hash that depends only
on *how things are connected*, not on what the nodes are named.  Two
netlists with equal signatures are the same circuit up to node
renaming; unequal signatures mean the connectivity genuinely changed.

How it works
------------
We run Weisfeiler-Leman (1-WL) colour refinement on the bipartite
component<->node incidence graph:

* Ground (``0`` / ``gnd``) is anchored to a fixed colour -- it must not
  be relabelled, because "connected to ground" is not the same as
  "connected to some internal node".
* Component colours start from the SPICE class letter (R/C/L/D/Q/...),
  so a resistor is never confused with a capacitor.
* Two-terminal *passives* (R, C, L) treat their two pins as unordered
  (a resistor flipped end-for-end is the same resistor).  Every other
  element keeps pin order, because polarity / pin role matters
  (diode anode vs cathode, BJT C/B/E, op-amp +/-, ...).

1-WL never assigns different colours to two genuinely isomorphic
graphs, so ``topology_equivalent`` has *no* false alarms: it will never
report drift on a circuit that actually round-tripped correctly.  It
can in principle miss drift between two WL-indistinguishable graphs
(extremely rare for real circuits), which only makes it conservative.
"""
from __future__ import annotations

import hashlib
from typing import Dict, List, Tuple

from .parser.netlist_to_asc import NetlistParser, ComponentType, Component

# Element classes whose two terminals are electrically interchangeable.
# A flipped resistor/capacitor/inductor is the same circuit, so we treat
# their pins as an unordered set during refinement.
_SYMMETRIC_TYPES = {
    ComponentType.RESISTOR,
    ComponentType.CAPACITOR,
    ComponentType.INDUCTOR,
}


def _is_ground(node: str) -> bool:
    return node == "0" or node.lower() == "gnd"


def _stable_hash(obj) -> str:
    """Process-independent short hash (builtin hash() is salted per run)."""
    return hashlib.sha1(repr(obj).encode("utf-8")).hexdigest()[:16]


def component_pins(comp: Component) -> List[str]:
    """Ordered pin -> node list for one parsed :class:`Component`.

    Mirrors the field layout the parser fills: positive, negative, then
    the optional control / output terminals, then any extra subcircuit
    pins.  Blank fields (unused terminals) are dropped.
    """
    nodes = [comp.node_pos, comp.node_neg]
    for n in (comp.node_ctrl, comp.node_ctrl2, comp.node_out):
        if n:
            nodes.append(n)
    nodes.extend(comp.extra_nodes)
    return [n for n in nodes if n]


def _parse(netlist: str) -> List[Component]:
    return NetlistParser().parse_string(netlist).components


def topology_signature(netlist: str, rounds: int = 8) -> str:
    """Return a node-rename-invariant hash of the netlist's connectivity.

    Args:
        netlist: SPICE netlist text.
        rounds: max WL refinement rounds (refinement also stops early
            once colours are stable).

    Returns:
        A short hex string.  Equal strings  => same circuit topology up
        to node renaming (and pin-swap for R/C/L).  Different strings
        => the wiring genuinely differs.
    """
    comps = _parse(netlist)
    if not comps:
        return _stable_hash(("empty",))

    # Incidence: node name -> list of (component index, pin index)
    node_inc: Dict[str, List[Tuple[int, int]]] = {}
    comp_pins: List[List[str]] = []
    for ci, c in enumerate(comps):
        pins = component_pins(c)
        comp_pins.append(pins)
        for pi, node in enumerate(pins):
            node_inc.setdefault(node, []).append((ci, pi))

    is_sym = [c.comp_type in _SYMMETRIC_TYPES for c in comps]
    # Anchor ground; everything else starts identical.
    node_color: Dict[str, str] = {
        n: ("G" if _is_ground(n) else "N") for n in node_inc
    }
    comp_color: List[str] = [c.comp_type.value for c in comps]

    for _ in range(rounds):
        new_comp = []
        for ci, pins in enumerate(comp_pins):
            ncols = [node_color[n] for n in pins]
            pin_sig = tuple(sorted(ncols)) if is_sym[ci] else tuple(ncols)
            new_comp.append(_stable_hash((comp_color[ci], pin_sig)))
        new_node = {}
        for n, inc in node_inc.items():
            sig = sorted(
                (comp_color[ci], "S" if is_sym[ci] else pi)
                for ci, pi in inc
            )
            new_node[n] = _stable_hash((node_color[n], tuple(sig)))
        if new_comp == comp_color and new_node == node_color:
            break  # converged
        comp_color, node_color = new_comp, new_node

    return _stable_hash(tuple(sorted(comp_color)))


def _pin_incidence_count(netlist: str) -> int:
    return sum(len(component_pins(c)) for c in _parse(netlist))


def topology_equivalent(netlist_a: str, netlist_b: str) -> Tuple[bool, dict]:
    """Compare two netlists for node-rename-invariant connectivity.

    Returns ``(equivalent, info)``.  ``info`` carries a small,
    human-readable diff (component counts and total pin-incidence on
    each side) so callers can produce a useful message without
    re-parsing.
    """
    comps_a = _parse(netlist_a)
    comps_b = _parse(netlist_b)
    sig_a = topology_signature(netlist_a)
    sig_b = topology_signature(netlist_b)
    pins_a = sum(len(component_pins(c)) for c in comps_a)
    pins_b = sum(len(component_pins(c)) for c in comps_b)
    info = {
        "components_a": len(comps_a),
        "components_b": len(comps_b),
        "pin_incidences_a": pins_a,
        "pin_incidences_b": pins_b,
    }
    return sig_a == sig_b, info
