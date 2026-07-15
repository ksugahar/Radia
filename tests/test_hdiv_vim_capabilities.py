"""Production HDiv-VIM field/geometry-order capability contract."""

import pytest

pytest.importorskip("ngsolve")

from radia.vim import hdiv_capabilities  # noqa: E402
from radia.vim._capabilities import validate_hdiv_configuration  # noqa: E402


def test_capability_table_is_dimension_and_topology_explicit():
    table = {
        (c.dimension, c.topology, c.hdiv_order):
            (c.geometry_orders, c.recommended_geometry_order)
        for c in hdiv_capabilities()
    }
    assert table[(2, "quad", 1)] == ((1, 2), 2)
    assert table[(2, "quad", 2)] == ((1, 2, 3), 3)
    assert table[(3, "tet", 2)] == ((1, 2), 2)
    assert table[(3, "hex", 2)] == ((1, 2), 2)
    assert table[(3, "wedge", 2)] == ((1, 2), 2)


def test_configuration_validation_has_no_global_p_plus_one_rule():
    validate_hdiv_configuration(2, {4}, 2, 3)
    validate_hdiv_configuration(3, {8}, 2, 2)
    validate_hdiv_configuration(3, {6}, 2, 2)

    with pytest.raises(ValueError, match="does not support geometry order 3"):
        validate_hdiv_configuration(3, {8}, 2, 3)
    with pytest.raises(ValueError, match="does not support geometry order 3"):
        validate_hdiv_configuration(3, {6}, 2, 3)
