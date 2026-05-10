// axi_henrotte_fespace.cpp — FESpace returning AxiHenrotteFE elements +
// custom DiffOps wired into evaluator[VOL] (Phase 2-C).

#include "axi_henrotte_fespace.hpp"
#include "axi_henrotte_fe.hpp"
#include "axi_henrotte_diffop.hpp"
#include <python_comp.hpp>

namespace radia_axifemm {

using namespace ngcomp;

AxiHenrotteFESpace::AxiHenrotteFESpace(shared_ptr<MeshAccess> ma, const Flags & flags)
  : FESpace(ma, flags)
{
    type = "axihenrotte";
    needs_transform_vec = false;
    axi_order = int(flags.GetNumFlag("order", 1.0));
    if (axi_order != 1 && axi_order != 2)
        throw Exception("AxiHenrotteFESpace: only order=1 (Q1/P1) and order=2 "
                        "(Q2 quad-only) are supported.");
    order = axi_order;  // base-class field; 4 = poly-degree-in-r for Q2

    // Wire our custom DiffOps so that SymbolicBilinearForm can integrate
    // u, grad(u) for u : H1Henrotte.
    //   evaluator[VOL]: invoked when ProxyFunction is sampled with no derivative
    //   additional_evaluators.Set("grad", ...): exposes grad(u) syntax in Python
    //   flux_evaluator[VOL]: used by visualization / flux post-processing
    evaluator[VOL]      = make_shared<AxiHenrotteDiffOpId>();
    flux_evaluator[VOL] = make_shared<AxiHenrotteDiffOpGradient>();
    additional_evaluators.Set("grad", make_shared<AxiHenrotteDiffOpGradient>());
    // Boundary evaluator: needed for `LinearForm += q * v * ds(label)`
    // patterns used by axisymmetric heat / magnetic Neumann RHS.  The
    // basis polynomial restricted to a boundary edge of an axis-aligned
    // quad is a 1D polynomial in (s, z); the DiffOpId evaluation at an
    // edge MIP gives the correct trace because the quad's (s, z) coords
    // are linear functions of the edge parameter.
    // (2026-05-10: added to enable full migration of axisymmetric
    // solvers to Henrotte per the NGSolve-H1-not-for-axisym policy.)
    evaluator[BND]      = make_shared<AxiHenrotteDiffOpIdBnd>();
}

void AxiHenrotteFESpace::Update() {
    FESpace::Update();
    if (axi_order == 1) {
        SetNDof(ma->GetNV());
    } else {
        // Q2 (quad only): vertex DOFs + edge midnode DOFs + face center DOFs.
        SetNDof(ma->GetNV() + ma->GetNEdges() + ma->GetNE());
    }
}

FiniteElement & AxiHenrotteFESpace::GetFE(ElementId ei, Allocator & lh) const {
    Ngs_Element ngel = ma->GetElement(ei);
    auto vertices = ngel.Vertices();

    // Boundary (1D segment) elements: needed for `LinearForm += q*v*ds(label)`
    // patterns in axisymmetric Neumann RHS assembly.  Return the Edge FE.
    if (ei.VB() != VOL && ngel.GetType() == ET_SEGM && vertices.Size() == 2) {
        Vec<3> p0 = ma->GetPoint<3>(vertices[0]);
        Vec<3> p1 = ma->GetPoint<3>(vertices[1]);
        if (axi_order == 1)
            return *new (lh) AxiHenrotteFE_Edge_Q1(p0(0), p0(1), p1(0), p1(1));
        else
            return *new (lh) AxiHenrotteFE_Edge_Q2(p0(0), p0(1), p1(0), p1(1));
    }

    if (axi_order == 2) {
        if (ngel.GetType() != ET_QUAD || vertices.Size() != 4)
            throw Exception("AxiHenrotteFESpace order=2 requires an all-quad "
                            "axis-aligned mesh; element " + ToString(ei) +
                            " is not a quad.");
        Vec<3> p[4];
        for (int i = 0; i < 4; ++i) p[i] = ma->GetPoint<3>(vertices[i]);
        double r_a = p[0](0), r_b = p[1](0);
        double z_a = p[0](1), z_b = p[2](1);
        return *new (lh) AxiHenrotteFE_Q2_AxisAligned(r_a, r_b, z_a, z_b);
    }
    // axi_order == 1 path
    if (ngel.GetType() == ET_QUAD && vertices.Size() == 4) {
        Vec<3> p[4];
        for (int i = 0; i < 4; ++i) p[i] = ma->GetPoint<3>(vertices[i]);
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

    if (axi_order == 1) {
        dnums.SetSize(vertices.Size());
        for (size_t i = 0; i < vertices.Size(); ++i)
            dnums[i] = vertices[i];
        return;
    }
    // order == 2 — distinguish VOL (2D quad: 9 DOFs) from BND (1D segment:
    // 2 vertex DOFs + 1 edge midnode DOF).
    int nv = ma->GetNV();
    int ne = ma->GetNEdges();

    if (ei.VB() != VOL) {
        // Boundary segment (or vertex-only element). Expose vertex + edge DOFs
        // so the user can Dirichlet axis edges.
        auto edges = ngel.Edges();
        dnums.SetSize(vertices.Size() + edges.Size());
        for (size_t i = 0; i < vertices.Size(); ++i)
            dnums[i] = vertices[i];
        for (size_t i = 0; i < edges.Size(); ++i)
            dnums[vertices.Size() + i] = nv + edges[i];
        return;
    }
    // VOL element — must be quad.
    if (ngel.GetType() != ET_QUAD)
        throw Exception("AxiHenrotteFESpace order=2: VOL element must be quad");
    auto edges = ngel.Edges();
    dnums.SetSize(9);
    for (int i = 0; i < 4; ++i) dnums[i] = vertices[i];
    // NGSolve's QUAD local edge order is the tensor-product convention:
    //   edges[0]: bottom (v0 -> v1), at (sm, za)
    //   edges[1]: top    (v2 -> v3), at (sm, zb)
    //   edges[2]: left   (v0 -> v2), at (sa, zm)
    //   edges[3]: right  (v1 -> v3), at (sb, zm)
    // Our Q2 FE local node order is the cyclic-CCW convention:
    //   local 4: bottom, local 5: right, local 6: top, local 7: left.
    // Permute accordingly so dnums[k] matches the FE's local-DOF k.
    dnums[4] = nv + edges[0];   // bottom
    dnums[5] = nv + edges[3];   // right
    dnums[6] = nv + edges[1];   // top
    dnums[7] = nv + edges[2];   // left
    dnums[8] = nv + ne + ei.Nr();
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
          [](shared_ptr<MeshAccess> ma, int order, py::kwargs kwargs)
              -> shared_ptr<FESpace> {
              Flags flags;
              flags.SetFlag("order", double(order));
              for (auto item : kwargs)
                  flags.SetFlag(py::cast<std::string>(item.first),
                                py::cast<std::string>(py::str(item.second)));
              auto fes = make_shared<AxiHenrotteFESpace>(ma, flags);
              fes->Update();
              fes->FinalizeUpdate();
              return fes;
          },
          py::arg("mesh"), py::arg("order") = 1,
          "Construct an AxiHenrotteFESpace for the given mesh.\n"
          "  order=1 (default): Q1 quad / P1 triangle, vertex DOFs only.\n"
          "  order=2          : Q2 quad-only (9 DOFs: 4 vertex + 4 edge + 1 face).");
}

}  // namespace radia_axifemm
