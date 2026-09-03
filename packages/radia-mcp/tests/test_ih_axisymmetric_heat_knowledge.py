from radia_mcp.ih.ih_knowledge import get_induction_heating_documentation
from radia_mcp.radia_ngsolve.knowledge.axifem import get_axifem_documentation


def test_ih_thermal_knowledge_keeps_em_and_heat_spaces_separate():
    thermal = get_induction_heating_documentation("thermal")

    assert "Henrotte for EM, NGSolve H1 for heat" in thermal
    assert "H1Henrotte" in thermal
    assert "ngsolve.H1(mesh, order=2)" in thermal
    assert "Q2 on quadrilateral meshes" in thermal
    assert "P2 on triangular meshes" in thermal
    assert "AxiHenrotteHeatStiffnessBFI" in thermal
    assert "AxiHenrotteHeatMassBFI" in thermal
    assert "fail fast on axis-touching Q2" in thermal


def test_ih_pitfalls_reject_henrotte_temperature_reuse():
    pitfalls = get_induction_heating_documentation("pitfalls")

    assert "Do Not Reuse the Henrotte A_phi Space for Temperature" in pitfalls
    assert "WRONG on an axis-touching Q2 mesh" in pitfalls
    assert "standard scalar H1" in pitfalls
    assert "Do not catch the" in pitfalls
    assert "exception and fall back" in pitfalls


def test_axifem_knowledge_routes_axisymmetric_heat_to_standard_h1():
    support = get_axifem_documentation("support_matrix")
    api = get_axifem_documentation("api")
    combined = support + "\n" + api

    assert "Axisymmetric heat" in combined
    assert "ngsolve.H1(mesh, order=2)" in combined
    assert "Q2 on quadrilateral meshes" in combined
    assert "P2 on triangular meshes" in combined
    assert "fail fast" in combined
