from radia_mcp.radia_ngsolve.knowledge.basis_functions import (
    get_basis_functions_documentation,
)


def test_basis_functions_knowledge_exposes_both_hdiv_families():
    doc = " ".join(get_basis_functions_documentation("hdiv_rt_bdm").split())

    assert "HDivTetBDM[p]" in doc
    assert "HDivTetRT[p]" in doc
    assert "HDivTrigBDM[p]" in doc
    assert "HDivTrigRT[p]" in doc
    assert "same `ker(div)` dimension" in doc
    assert "HDivHexBDM[p] === HDivHexRT[p]" in doc
