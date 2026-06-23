// axi_henrotte_fespace.hpp — FESpace returning AxiHenrotteFE elements.
//
// Phase 2-A: stub (registration, basic Update/GetDofNrs). Real GetFE that
// dispatches between Q1 axis-aligned and P1 triangle based on the mesh
// element type lands in Phase 2-C.

#ifndef RADIA_AXIFEMM_AXI_HENROTTE_FESPACE_HPP
#define RADIA_AXIFEMM_AXI_HENROTTE_FESPACE_HPP

#include <comp.hpp>

namespace axifem {

using namespace ngcomp;

class AxiHenrotteFESpace : public FESpace {
public:
    int axi_order = 1;   // 1 = Q1/P1 (vertex DOFs), 2 = Q2 (vertex+edge+face)
    bool curved_quad = false;  // order=2 quads: sample 9 curved nodes (Q2_Curved) vs axis-aligned

    AxiHenrotteFESpace(shared_ptr<MeshAccess> ma, const Flags & flags);

    static DocInfo GetDocu() {
        DocInfo doc;
        doc.short_docu = "Axisymmetric H1 with Henrotte basis (linear in r^2, z).";
        doc.long_docu =
            "FESpace whose nodal shape functions are linear in {1, r^2, z, ...} "
            "on physical (r, z) coordinates instead of the reference element.\n"
            "  order=1 (default): 3 DOFs/triangle, 4 DOFs/quad (vertex DOFs only)\n"
            "  order=2          : 9 DOFs/quad (vertex + edge midnode + face center;\n"
            "                     all-quad meshes only).";
        return doc;
    }

    void Update() override;
    // order=2 reserves a face-center DOF per VOL element (nv+ne+ei.Nr()), but
    // P2 TRIANGLES never use theirs -> a dead (zero row/col) DOF that leaves the
    // assembled K/M rank-deficient (breaks the eddy generalized eigenvalue solve
    // and any iterative solve).  Mark those trig-center slots UNUSED so they drop
    // out of FreeDofs; quad centers stay LOCAL (Q2 uses them).
    void UpdateCouplingDofArray() override;
    // The base FinalizeUpdate builds free_dofs BEFORE any ctofdof is set and does
    // not invoke UpdateCouplingDofArray, so we override it to (a) populate ctofdof
    // and (b) directly drop the dead P2-triangle center slots from free_dofs.
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
