from radia_mcp.radia_ngsolve.knowledge.radia import get_radia_documentation


def test_arc_manual_separates_native_accuracy_routes():
    text = get_radia_documentation('geometry')
    assert 'analytic rectangular-section integration' in text
    assert 'Non-finite or' in text and 'raises an error' in text
    assert 'Full circles use a regular axis expansion' in text
    assert 'older installed extension' in text
    assert 'analytic rectangular moments and a tail bound' in text
    assert 'require the expected shape and finite values' in text
    assert "0.001, 100, 'man', 'z', 1e6" in text
