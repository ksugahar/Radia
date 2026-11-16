// axi_henrotte_fespace.hpp — FESpace returning AxiHenrotteFE elements.
//
// order=1 dispatches to P1 triangles or Q1 axis-aligned quads.
// order=2 dispatches to P2 triangles or Q2 axis-aligned quads.  The P2
// triangle path is curved-mesh aware through NGSolve's element
// transformation; the Q2 quad path is the straight axis-aligned closed-form
// element.

#ifndef RADIA_AXIFEM_AXI_HENROTTE_FESPACE_HPP
#define RADIA_AXIFEM_AXI_HENROTTE_FESPACE_HPP

#include <comp.hpp>

namespace axifem {

using namespace ngcomp;

class AxiHenrotteFESpace : public FESpace {
public:
    int axi_order = 1;   // 1 = P1/Q1, 2 = P2/Q2
    bool curved_quad = false;  // order=2 quad: opt into AxiHenrotteFE_Q2_Curved

    AxiHenrotteFESpace(shared_ptr<MeshAccess> ma, const Flags & flags);

    static DocInfo GetDocu() {
        DocInfo doc;
        doc.short_docu = "Axisymmetric H1 with Henrotte basis (linear in r^2, z).";
        doc.long_docu =
            "FESpace whose nodal shape functions are linear in {1, r^2, z, ...} "
            "on physical (r, z) coordinates instead of the reference element.\n"
            "  order=1 (default): P1 triangle (3 DOFs) or Q1 quad (4 DOFs),\n"
            "                     vertex DOFs only\n"
            "  order=2          : P2 triangle (6 DOFs; curved-mesh aware) or\n"
            "                     Q2 axis-aligned quad (9 DOFs: vertex + edge\n"
            "                     midnode + face center).  Pass curvedquad=True\n"
            "                     to use the opt-in curved Q2 quad element.";
        return doc;
    }

    void Update() override;
    // Mesh generators place axis nodes at |r| ~ 1e-15 of the mesh extent.
    // Every coordinate handed to an element is snapped to r = 0 within
    // axis_tol (1e-12 of the mesh extent), the same tolerance everywhere, so
    // the axis nodes of all elements and boundary traces agree.  Their basis
    // is identically zero (N_i = r_i psi_i / r), so their DOFs are UNUSED.
    double axis_tol = 0.0;
    double SnapAxis(double r) const { return std::fabs(r) <= axis_tol ? 0.0 : r; }
    // Marks every DOF that no volume element of the current level uses as
    // UNUSED: the reserved centre slot of a P2 triangle and, after
    // mesh.Refine(), the coarse-level edges.  Left free, each is a zero
    // row/column that makes K/M singular.
    void UpdateCouplingDofArray() override;
    // The base FinalizeUpdate builds free_dofs BEFORE any ctofdof is set and does
    // not invoke UpdateCouplingDofArray, so we override it to (a) populate ctofdof
    // and (b) drop the UNUSED DOFs from free_dofs.
    void FinalizeUpdate() override;
    FiniteElement & GetFE(ElementId ei, Allocator & lh) const override;
    void GetDofNrs(ElementId ei, Array<DofId> & dnums) const override;
    // NodeId overload — needed so ngsolve.Periodic(H1Henrotte(...)) can
    // discover which DOFs sit on each vertex / edge that participates in
    // a PERIODIC identification, and couple them. The base FESpace's
    // default is empty (returns 0 DOFs), which silently disables Periodic
    // coupling. Phase B3 (2026-05-12).
    void GetDofNrs(NodeId ni, Array<DofId> & dnums) const override;
    // Per-node-type variants — base FESpace::GetDofNrs(NodeId) default
    // implementation actually dispatches to these by node-type, so we
    // override each one explicitly in addition to the unified NodeId
    // variant.
    void GetVertexDofNrs(int vnr, Array<DofId> & dnums) const override;
    void GetEdgeDofNrs(int ednr, Array<DofId> & dnums) const override;
    void GetFaceDofNrs(int fanr, Array<DofId> & dnums) const override;

    string GetClassName() const override { return "AxiHenrotteFESpace"; }
};

}  // namespace axifem

#endif
