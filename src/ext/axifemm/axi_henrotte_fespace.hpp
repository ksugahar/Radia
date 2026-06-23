// axi_henrotte_fespace.hpp — FESpace returning AxiHenrotteFE elements.
//
// order=1 dispatches to P1 triangles or Q1 axis-aligned quads.
// order=2 dispatches to P2 triangles or Q2 axis-aligned quads.  The P2
// triangle path is curved-mesh aware through NGSolve's element
// transformation; the Q2 quad path is the straight axis-aligned closed-form
// element.

#ifndef RADIA_AXIFEMM_AXI_HENROTTE_FESPACE_HPP
#define RADIA_AXIFEMM_AXI_HENROTTE_FESPACE_HPP

#include <comp.hpp>

namespace axifem {

using namespace ngcomp;

class AxiHenrotteFESpace : public FESpace {
public:
    int axi_order = 1;   // 1 = P1/Q1, 2 = P2/Q2

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
            "                     midnode + face center).";
        return doc;
    }

    void Update() override;
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
