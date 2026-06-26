from radia_mcp.radia_ngsolve.knowledge.axifem import get_axifem_documentation


def test_axifem_knowledge_lists_curved_q2_as_shipping_opt_in():
    support = get_axifem_documentation("support_matrix")
    curved = get_axifem_documentation("curved_geometry")
    combined = support + "\n" + curved

    assert "curvedquad=True" in combined
    assert "AxiHenrotteFE_Q2_Curved" in combined
    assert "shipping opt-in" in support
    assert "production C++" in curved

    stale_phrases = [
        "prototype-only",
        "not wired",
        "not production C++",
        "not part of the production C++ implementation",
        "does not dispatch",
        "true curved Q2 quads remain",
    ]
    lowered = combined.lower()
    for phrase in stale_phrases:
        assert phrase.lower() not in lowered
