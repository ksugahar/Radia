from radia_mcp.ih.ih_knowledge import get_induction_heating_documentation
from radia_mcp.radia_ngsolve.knowledge.axifem import get_axifem_documentation


def test_ih_thermal_knowledge_keeps_em_and_heat_spaces_separate():
    thermal = get_induction_heating_documentation("thermal")
    normalized = " ".join(thermal.split())

    assert "Henrotte for EM, NGSolve H1 for heat" in thermal
    assert "H1Henrotte" in thermal
    assert "ngsolve.H1(mesh, order=2)" in thermal
    assert "Q2 on quadrilateral meshes" in thermal
    assert "P2 on triangular meshes" in thermal
    assert "AxiHenrotteHeatStiffnessBFI" in thermal
    assert "AxiHenrotteHeatMassBFI" in thermal
    assert "fail fast on axis-touching Q2" in thermal
    assert "--heat-flux-boundaries" in thermal
    assert "--convection-boundaries" in thermal
    assert "--radiation-boundaries" in thermal
    assert "boundary_audit" in thermal
    assert "--surface-label`` fails" in thermal
    assert "true ``r=0`` center axis" in normalized
    assert "zero because the weak form uses ``2*pi*r*ds``" in normalized
    assert "solver fails fast" in normalized
    assert "inner cylindrical surface at ``r>0``" in normalized
    assert "Thermal `.vol` geometry contract: preserve after load" in thermal
    assert "do NOT call wp_mesh.Curve(2)" in thermal
    assert "preserve-input-vol-geometry" in thermal
    assert "post_load_curve_applied=false" in thermal
    assert "297 million times" in thermal


def test_ih_curve_guidance_distinguishes_live_geometry_from_loaded_vol():
    combined = get_induction_heating_documentation("gmsh_mesh")
    normalized = " ".join(combined.split())

    assert "Curve a mesh only while its generating geometry is live" in combined
    assert "Never call `mesh.Curve()`" in combined
    assert "after loading an arbitrary `.vol`" in normalized
    assert "ngsolve-curve-after-vol-import" in combined


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


def test_ih_thermal_knowledge_fixes_cross_mesh_qsurf_handoff_to_p1():
    thermal = get_induction_heating_documentation("thermal")
    rotating = get_induction_heating_documentation("rotating")

    assert "fixed P1 cross-mesh handoff" in thermal
    assert "electromagnetic solve itself uses a higher order" in thermal
    assert "always produces this handoff at order 1" in rotating


def test_ih_thermal_knowledge_prefers_axisym_and_audits_3d_handoff():
    thermal = get_induction_heating_documentation("thermal")
    normalized = " ".join(thermal.split())

    assert "3D EM does not require 3D heat" in thermal
    assert "Do not select 3D heat" in normalized
    assert "experiment constraint, not a formulation requirement" in normalized
    assert "separately generated 2D ``(r,z)`` workpiece mesh" in normalized
    assert "``Curve()`` changes geometry order, not dimension" in normalized
    assert "boundary point locator (``BND``)" in normalized
    assert "default is 128 azimuth samples" in normalized
    assert "``qsurf_projection``" in thermal
    assert "no zero-flux fallback" in normalized
