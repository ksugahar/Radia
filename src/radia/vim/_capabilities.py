"""Single source of truth for production HDiv-VIM element/geometry pairs.

NGSolve's ``order`` names the finite-element space; ``mesh.GetCurveOrder()``
names the independent geometry map.  Their useful pairings differ between
2D and 3D because the Piola map is dimension dependent, so callers must not
derive one order from the other with an ad-hoc arithmetic rule.
"""
from __future__ import annotations

from dataclasses import dataclass


@dataclass(frozen=True)
class HDivCapability:
    dimension: int
    topology: str
    hdiv_order: int
    geometry_orders: tuple[int, ...]
    recommended_geometry_order: int


_CAPABILITIES = (
    HDivCapability(2, "tri", 1, (1, 2), 2),
    HDivCapability(2, "tri", 2, (1, 2, 3), 3),
    HDivCapability(2, "quad", 1, (1, 2), 2),
    HDivCapability(2, "quad", 2, (1, 2, 3), 3),
    HDivCapability(2, "tri-quad", 1, (1, 2), 2),
    HDivCapability(2, "tri-quad", 2, (1, 2, 3), 3),
    HDivCapability(3, "tet", 1, (1, 2), 2),
    HDivCapability(3, "tet", 2, (1, 2), 2),
    HDivCapability(3, "hex", 1, (1, 2), 2),
    HDivCapability(3, "hex", 2, (1, 2), 2),
    HDivCapability(3, "wedge", 1, (1, 2), 2),
    HDivCapability(3, "wedge", 2, (1, 2), 2),
)
_BY_KEY = {(c.dimension, c.topology, c.hdiv_order): c for c in _CAPABILITIES}


def classify_hdiv_topology(dimension: int, vertex_counts) -> str:
    """Return the production topology name from element vertex counts."""
    counts = frozenset(int(v) for v in vertex_counts)
    if int(dimension) == 2:
        if counts == {3}:
            return "tri"
        if counts == {4}:
            return "quad"
        if counts == {3, 4}:
            return "tri-quad"
    elif int(dimension) == 3:
        if counts == {4}:
            return "tet"
        if counts == {6}:
            return "wedge"
        if counts == {8}:
            return "hex"
    raise ValueError(
        "HDiv-VIM supports 2D TRI/QUAD mixtures or one pure 3D topology "
        "(TET, HEX, WEDGE); got dimension=%r, vertex_counts=%s."
        % (dimension, sorted(counts)))


def validate_hdiv_configuration(dimension: int, vertex_counts, hdiv_order: int,
                                geometry_order: int) -> HDivCapability:
    """Validate and return the production capability for one mesh/space pair."""
    topology = classify_hdiv_topology(dimension, vertex_counts)
    key = (int(dimension), topology, int(hdiv_order))
    capability = _BY_KEY.get(key)
    if capability is None:
        raise ValueError(
            "HDiv-VIM supports HDiv order in {1,2}; got dimension=%d, topology=%s, order=%r."
            % (int(dimension), topology, hdiv_order))
    geometry_order = max(1, int(geometry_order))
    if geometry_order not in capability.geometry_orders:
        raise ValueError(
            "HDiv-VIM does not support geometry order %d for %dD %s BDM%d; supported geometry orders are %s. "
            "Geometry and field orders are independent Piola-FEM choices, not an automatic p+1 rule."
            % (geometry_order, capability.dimension, capability.topology,
               capability.hdiv_order, capability.geometry_orders))
    return capability


def hdiv_capabilities() -> tuple[HDivCapability, ...]:
    """Return the immutable production capability table for docs/lint/tests."""
    return _CAPABILITIES


__all__ = [
    "HDivCapability", "classify_hdiv_topology", "validate_hdiv_configuration",
    "hdiv_capabilities",
]
