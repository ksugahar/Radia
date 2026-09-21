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


def test_circular_loft_routing_requires_explicit_approximation():
    text = electromagnet_usage('coilbuilder')
    assert 'matching rectangular or circular' in text
    assert 'equal-area radial/angular grid' in text
    assert 'cross-type transitions' in text
    assert 'independent integral' in text


def test_closed_loop_contract_is_not_automatic_geometry_certification():
    text = electromagnet_usage('coilbuilder')
    assert 'require_closed=True' in text
    assert 'No return wire is added' in text
    assert 'self-intersections' in text


def test_mixed_omega_cost_and_volume_model_limits():
    text=electromagnet_usage('coilbuilder')
    assert 'timeout means incomplete execution, not numerical disagreement' in text
    assert 'interface projection' in text
    assert 'conductivity and terminal conditions' in text
    assert 'does not claim that coupled volume-conductor support' in text
    assert 'full-turn rectangular/circular loft CAD uses exact revolution' in text
    assert 'audit_segment_volume_overlaps()' in text
    assert 'not clearance or intra-segment validity' in text
    assert 'source RHS' in text
