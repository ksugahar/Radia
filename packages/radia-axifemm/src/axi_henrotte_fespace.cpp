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
    if (axi_order != 1 && axi_order != 2 && axi_order != 3)
        throw Exception("AxiHenrotteFESpace: only order=1 (Q1/P1), order=2 "
                        "(Q2 quad-only), and order=3 (Q3 quad-only) are supported.");
    order = axi_order;

    // Wire our custom DiffOps so that SymbolicBilinearForm can integrate
    // u, grad(u) for u : H1Henrotte.
    //   evaluator[VOL]: invoked when ProxyFunction is sampled with no derivative
    //   additional_evaluators.Set("grad", ...): exposes grad(u) syntax in Python
    //   flux_evaluator[VOL]: used by visualization / flux post-processing
    evaluator[VOL]      = make_shared<AxiHenrotteDiffOpId>();
    flux_evaluator[VOL] = make_shared<AxiHenrotteDiffOpGradient>();
    additional_evaluators.Set("grad", make_shared<AxiHenrotteDiffOpGradient>());
}

void AxiHenrotteFESpace::Update() {
    FESpace::Update();
    int nv = ma->GetNV(), ne = ma->GetNEdges(), nf = ma->GetNE();
    if (axi_order == 1) {
        SetNDof(nv);
    } else if (axi_order == 2) {
        SetNDof(nv + ne + nf);              // Q2: 1 vertex + 1 edge + 1 face per
    } else /* axi_order == 3 */ {
        SetNDof(nv + 2 * ne + 4 * nf);      // Q3: 1 vertex + 2 edge + 4 face per
    }
}

FiniteElement & AxiHenrotteFESpace::GetFE(ElementId ei, Allocator & lh) const {
    Ngs_Element ngel = ma->GetElement(ei);
    auto vertices = ngel.Vertices();

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
    if (axi_order == 3) {
        if (ngel.GetType() != ET_QUAD || vertices.Size() != 4)
            throw Exception("AxiHenrotteFESpace order=3 requires an all-quad "
                            "axis-aligned mesh; element " + ToString(ei) +
                            " is not a quad.");
        Vec<3> p[4];
        for (int i = 0; i < 4; ++i) p[i] = ma->GetPoint<3>(vertices[i]);
        double r_a = p[0](0), r_b = p[1](0);
        double z_a = p[0](1), z_b = p[2](1);
        return *new (lh) AxiHenrotteFE_Q3_AxisAligned(r_a, r_b, z_a, z_b);
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
    int nv = ma->GetNV();
    int ne = ma->GetNEdges();

    if (axi_order == 2) {
        // order == 2 - distinguish VOL (2D quad: 9 DOFs) from BND (1D segment:
        // 2 vertex DOFs + 1 edge midnode DOF).
        if (ei.VB() != VOL) {
            auto edges = ngel.Edges();
            dnums.SetSize(vertices.Size() + edges.Size());
            for (size_t i = 0; i < vertices.Size(); ++i)
                dnums[i] = vertices[i];
            for (size_t i = 0; i < edges.Size(); ++i)
                dnums[vertices.Size() + i] = nv + edges[i];
            return;
        }
        if (ngel.GetType() != ET_QUAD)
            throw Exception("AxiHenrotteFESpace order=2: VOL element must be quad");
        auto edges = ngel.Edges();
        dnums.SetSize(9);
        for (int i = 0; i < 4; ++i) dnums[i] = vertices[i];
        // NGSolve's QUAD local edge order is the tensor-product convention:
        //   edges[0]: bottom, edges[1]: top, edges[2]: left, edges[3]: right
        // Our Q2 FE local node order is cyclic CCW (bottom, right, top, left).
        dnums[4] = nv + edges[0];
        dnums[5] = nv + edges[3];
        dnums[6] = nv + edges[1];
        dnums[7] = nv + edges[2];
        dnums[8] = nv + ne + ei.Nr();
        return;
    }
    // order == 3: 16 DOFs/quad (4 vertex + 8 edge + 4 face). Q3 FE uses
    // NGSolve's tensor-product local edge order DIRECTLY (no permutation).
    if (ei.VB() != VOL) {
        // BND segment: 2 vertex + 2 edge midnode DOFs.
        auto edges = ngel.Edges();
        dnums.SetSize(vertices.Size() + 2 * edges.Size());
        for (size_t i = 0; i < vertices.Size(); ++i)
            dnums[i] = vertices[i];
        for (size_t i = 0; i < edges.Size(); ++i) {
            dnums[vertices.Size() + 2*i + 0] = nv + 2 * edges[i] + 0;
            dnums[vertices.Size() + 2*i + 1] = nv + 2 * edges[i] + 1;
        }
        return;
    }
    if (ngel.GetType() != ET_QUAD)
        throw Exception("AxiHenrotteFESpace order=3: VOL element must be quad");
    auto edges = ngel.Edges();
    dnums.SetSize(16);
    for (int i = 0; i < 4; ++i) dnums[i] = vertices[i];
    for (int e = 0; e < 4; ++e) {
        // edges[0]: bottom, edges[1]: top, edges[2]: left, edges[3]: right
        // Q3 FE local indices: 4-5 bottom, 6-7 top, 8-9 left, 10-11 right
        dnums[4 + 2*e + 0] = nv + 2 * edges[e] + 0;
        dnums[4 + 2*e + 1] = nv + 2 * edges[e] + 1;
    }
    int face_base = nv + 2 * ne + 4 * ei.Nr();
    for (int k = 0; k < 4; ++k) dnums[12 + k] = face_base + k;
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
