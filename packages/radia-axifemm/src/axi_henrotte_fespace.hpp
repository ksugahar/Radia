// axi_henrotte_fespace.hpp — FESpace returning AxiHenrotteFE elements.
//
// Phase 2-A: stub (registration, basic Update/GetDofNrs). Real GetFE that
// dispatches between Q1 axis-aligned and P1 triangle based on the mesh
// element type lands in Phase 2-C.

#ifndef RADIA_AXIFEMM_AXI_HENROTTE_FESPACE_HPP
#define RADIA_AXIFEMM_AXI_HENROTTE_FESPACE_HPP

#include <comp.hpp>

namespace radia_axifemm {

using namespace ngcomp;

class AxiHenrotteFESpace : public FESpace {
public:
    AxiHenrotteFESpace(shared_ptr<MeshAccess> ma, const Flags & flags);

    static DocInfo GetDocu() {
        DocInfo doc;
        doc.short_docu = "Axisymmetric H1 with Henrotte basis (linear in r^2, z).";
        doc.long_docu =
            "FESpace whose nodal shape functions are linear in {1, r^2, z, ...} "
            "on physical (r, z) coordinates instead of the reference element. "
            "Currently supports order=2 (3 nodes/triangle, 4 nodes/quad).";
        return doc;
    }

    void Update() override;
    FiniteElement & GetFE(ElementId ei, Allocator & lh) const override;
    void GetDofNrs(ElementId ei, Array<DofId> & dnums) const override;

    string GetClassName() const override { return "AxiHenrotteFESpace"; }
};

}  // namespace radia_axifemm

#endif
