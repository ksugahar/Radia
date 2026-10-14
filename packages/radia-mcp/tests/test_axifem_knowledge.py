from pathlib import Path

from radia_mcp.radia_ngsolve.knowledge.axifem import get_axifem_documentation


STALE_CURVED_Q2_PHRASES = [
    "prototype-only",
    "not wired",
    "not production C++",
    "not part of the production C++ implementation",
    "does not dispatch",
    "true curved Q2 quads remain",
]


def test_axifem_knowledge_lists_curved_q2_as_shipping_opt_in():
    support = get_axifem_documentation("support_matrix")
    curved = get_axifem_documentation("curved_geometry")
    combined = support + "\n" + curved

    assert "curvedquad=True" in combined
    assert "AxiHenrotteFE_Q2_Curved" in combined
    assert "shipping opt-in" in support
    assert "production C++" in curved

    lowered = combined.lower()
    for phrase in STALE_CURVED_Q2_PHRASES:
        assert phrase.lower() not in lowered


def test_axifem_tool_docstring_does_not_reintroduce_curved_q2_stale_status():
    server_py = (
        Path(__file__).resolve().parents[1]
        / "src"
        / "radia_mcp"
        / "radia_ngsolve"
        / "server.py"
    )
    text = server_py.read_text(encoding="utf-8")

    assert "curvedquad=True" in text
    for phrase in STALE_CURVED_Q2_PHRASES:
        assert phrase.lower() not in text.lower()


def test_axifem_knowledge_explains_henrotte_as_hodge_geometry():
    text = get_axifem_documentation("hodge_geometry")

    assert "Hodge" in text
    assert "exterior derivative" in text
    assert "psi(r,z)" in text
    assert "2*pi*r*A_phi" in text
    assert "s = r^2" in text
    assert "metric/Hodge" in text
    assert "weights" in text
