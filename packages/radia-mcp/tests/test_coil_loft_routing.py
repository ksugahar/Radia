from radia_mcp.electromagnet.server import electromagnet_usage


def test_loft_routing_is_explicit_and_bounded():
    text = electromagnet_usage('coilbuilder')
    assert 'to_radia_loft_filaments(nw, nh, n_arc=64)' in text
    assert 'thin-filament approximation' in text
    assert 'Require this method to' in text
    assert 'Negative-angle arc' in text
    assert 'self-energy' in text
    assert 'return path' in text


def test_loft_cad_does_not_imply_native_solid_field_support():
    text = electromagnet_usage('coilbuilder')
    assert '0 < angle < 360' in text
    assert 'n_sub >= 4' in text
    assert 'STEP round trips' in text
    assert 'CAD export does not make `to_radia()` support' in text
