// axi_henrotte_fespace.cpp — minimal FESpace stub. Full impl in Phase 2-C.

#include "axi_henrotte_fespace.hpp"
#include "axi_henrotte_fe.hpp"
#include <python_comp.hpp>

namespace radia_axifemm {

using namespace ngcomp;

AxiHenrotteFESpace::AxiHenrotteFESpace(shared_ptr<MeshAccess> ma, const Flags & flags)
  : FESpace(ma, flags)
{
    // Phase 2-A: vertex-DOF only (one DOF per mesh vertex).
    type = "axihenrotte";
    needs_transform_vec = false;

    // For a real implementation we'd register evaluators here; stubbed for now.
}

void AxiHenrotteFESpace::Update() {
    FESpace::Update();
    // One DOF per vertex.
    SetNDof(ma->GetNV());
}

FiniteElement & AxiHenrotteFESpace::GetFE(ElementId ei, Allocator & lh) const {
    Ngs_Element ngel = ma->GetElement(ei);
    auto vertices = ngel.Vertices();

    if (ngel.GetType() == ET_QUAD && vertices.Size() == 4) {
        // Get vertex coords; assume axis-aligned for Phase 2-A.
        Vec<3> p[4];
        for (int i = 0; i < 4; ++i) p[i] = ma->GetPoint<3>(vertices[i]);
        // Q1 expects (r_a, z_a), (r_b, z_a), (r_b, z_b), (r_a, z_b)
        double r_a = p[0](0), r_b = p[1](0);
        double z_a = p[0](1), z_b = p[2](1);
        return *new (lh) AxiHenrotteFE_Q1_AxisAligned(r_a, r_b, z_a, z_b);
    }
    if (ngel.GetType() == ET_TRIG && vertices.Size() == 3) {
        double rs[3], zs[3];
        for (int i = 0; i < 3; ++i) {
            auto pt = ma->GetPoint<3>(vertices[i]);
            rs[i] = pt(0);
            zs[i] = pt(1);
        }
        return *new (lh) AxiHenrotteFE_P1_Triangle(rs, zs);
    }
    throw Exception("AxiHenrotteFESpace: unsupported element type for "
                    + ToString(ei) + " (only triangles and quads in 2D)");
}

void AxiHenrotteFESpace::GetDofNrs(ElementId ei, Array<DofId> & dnums) const {
    Ngs_Element ngel = ma->GetElement(ei);
    auto vertices = ngel.Vertices();
    dnums.SetSize(vertices.Size());
    for (size_t i = 0; i < vertices.Size(); ++i)
        dnums[i] = vertices[i];
}

// Free creator function (so we can pass it as a function pointer to AddFESpace).
static shared_ptr<FESpace> CreateAxiHenrotteFESpace(
    shared_ptr<MeshAccess> ma, const Flags & flags)
{
    return make_shared<AxiHenrotteFESpace>(ma, flags);
}

void ExportAxiHenrotteFESpace(pybind11::module & m) {
    namespace py = pybind11;

    static bool registered = []() {
        ngcomp::GetFESpaceClasses().AddFESpace(
            "axihenrotte",
            CreateAxiHenrotteFESpace,
            &AxiHenrotteFESpace::GetDocu);
        return true;
    }();
    (void)registered;

    py::class_<AxiHenrotteFESpace, FESpace, shared_ptr<AxiHenrotteFESpace>>(
        m, "AxiHenrotteFESpace",
        "FESpace with FEMM/Henrotte axisymmetric basis on physical (r, z).");

    m.def("H1Henrotte",
          [](shared_ptr<MeshAccess> ma, py::kwargs kwargs)
              -> shared_ptr<FESpace> {
              Flags flags;
              for (auto item : kwargs)
                  flags.SetFlag(py::cast<std::string>(item.first),
                                py::cast<std::string>(py::str(item.second)));
              auto fes = make_shared<AxiHenrotteFESpace>(ma, flags);
              fes->Update();
              fes->FinalizeUpdate();
              return fes;
          },
          py::arg("mesh"),
          "Construct an AxiHenrotteFESpace for the given mesh.");
}

}  // namespace radia_axifemm
