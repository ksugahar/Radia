// axi_henrotte_fespace.cpp — FESpace returning AxiHenrotteFE elements +
// custom DiffOps wired into evaluator[VOL] (Phase 2-C).

#include "axi_henrotte_fespace.hpp"
#include "axi_henrotte_fe.hpp"
#include "axi_henrotte_diffop.hpp"
#include <python_comp.hpp>
#include <algorithm>
#include <cmath>

namespace axifem {

using namespace ngcomp;

namespace {

struct AxisAlignedQuadBounds {
    double r_min;
    double r_max;
    double z_min;
    double z_max;
};

double QuadCoordTol(const Vec<3> p[4])
{
    double span = 1.0;
    for (int i = 0; i < 4; ++i) {
        span = std::max(span, std::abs(p[i](0)));
        span = std::max(span, std::abs(p[i](1)));
    }
    return 1.0e-12 * span;
}

bool Near(double a, double b, double tol)
{
    return std::abs(a - b) <= tol;
}

AxisAlignedQuadBounds RequireAxisAlignedQuad(const Vec<3> p[4], ElementId ei, int order)
{
    double tol = QuadCoordTol(p);
    AxisAlignedQuadBounds b{p[0](0), p[0](0), p[0](1), p[0](1)};
    for (int i = 1; i < 4; ++i) {
        b.r_min = std::min(b.r_min, p[i](0));
        b.r_max = std::max(b.r_max, p[i](0));
        b.z_min = std::min(b.z_min, p[i](1));
        b.z_max = std::max(b.z_max, p[i](1));
    }
    bool has_ll = false, has_lr = false, has_ul = false, has_ur = false;
    bool axis_aligned = (b.r_max > b.r_min + tol) && (b.z_max > b.z_min + tol);
    for (int i = 0; i < 4 && axis_aligned; ++i) {
        const bool r_lo = Near(p[i](0), b.r_min, tol);
        const bool r_hi = Near(p[i](0), b.r_max, tol);
        const bool z_lo = Near(p[i](1), b.z_min, tol);
        const bool z_hi = Near(p[i](1), b.z_max, tol);
        axis_aligned = (r_lo || r_hi) && (z_lo || z_hi);
        has_ll = has_ll || (r_lo && z_lo);
        has_lr = has_lr || (r_hi && z_lo);
        has_ul = has_ul || (r_lo && z_hi);
        has_ur = has_ur || (r_hi && z_hi);
    }
    axis_aligned = axis_aligned && has_ll && has_lr && has_ul && has_ur;
    if (!axis_aligned) {
        string hint = (order == 2)
            ? "use H1Henrotte(..., order=2, curvedquad=True) for skewed/curved quads"
            : "use structured axis-aligned quads or triangle elements";
        throw Exception("AxiHenrotteFESpace: non-axis-aligned quad at "
                        + ToString(ei) + " is not valid for the closed-form Q"
                        + ToString(order) + " path; " + hint);
    }
    return b;
}

}  // namespace

AxiHenrotteFESpace::AxiHenrotteFESpace(shared_ptr<MeshAccess> ma, const Flags & flags)
  : FESpace(ma, flags)
{
    type = "axihenrotte";
    needs_transform_vec = false;
    // Honor complex=True (time-harmonic A_phi): the base 2-arg FESpace ctor does
    // not auto-set iscomplex for this custom space, so reflect the flag here.
    // Accept either the boolean DEFINE flag (pybind now passes complex=True as
    // such) or a legacy "True" string flag.  With iscomplex set, GridFunctions
    // store COMPLEX DOFs (needed for solve_axi_eddy_harmonic post-processing);
    // the complex CalcMatrix overloads in the DiffOps then evaluate them.
    iscomplex = flags.GetDefineFlag("complex") ||
                (flags.GetStringFlag("complex", "") == string("True"));
    axi_order = int(flags.GetNumFlag("order", 1.0));
    if (axi_order != 1 && axi_order != 2)
        throw Exception("AxiHenrotteFESpace: only order=1 (Q1/P1) and order=2 "
                        "(Q2 quad or P2 triangle) are supported.");
    order = axi_order;  // base-class field; 4 = poly-degree-in-r for Q2

    // curved=True: order-2 QUADS sample 9 curved node positions from the mesh's
    // (Curve(2)-d) ElementTransformation and use the AxiHenrotteFE_Q2_Curved
    // isoparametric element -- follows curved boundaries to O(curve^2), the
    // proper curved element for eddy-current / rounded-conductor problems.
    // Default False keeps the axis-aligned closed-form Q2 (and the heat path).
    curved_quad = flags.GetDefineFlag("curvedquad") ||
                  (flags.GetStringFlag("curvedquad", "") == string("True"));

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
        // order==2: vertex DOFs + edge midnode DOFs + face-center DOFs.
        // Q2 (quad) uses all three; P2 (trig) uses only the first two. We
        // allocate NV + NEdges + NE so each quad gets its face-center via
        // ei.Nr(); trig elements have an unused slot at ei.Nr() (harmless).
        SetNDof(ma->GetNV() + ma->GetNEdges() + ma->GetNE());
    }
}

void AxiHenrotteFESpace::UpdateCouplingDofArray() {
    // Vertices/edges are shared between elements (INTERFACE); element face-center
    // slots are LOCAL for quads (Q2 uses them) but UNUSED for triangles (P2 has
    // no center DOF -> the reserved nv+ne+ei.Nr() slot would be a dead zero
    // row/col that makes K/M singular).  Marking it UNUSED removes it from
    // FreeDofs so the generalized eddy eigenvalue solve (and direct solves) see a
    // full-rank system.
    ctofdof.SetSize(GetNDof());
    ctofdof = INTERFACE_DOF;
    if (axi_order == 2) {
        int nv = ma->GetNV();
        int ne = ma->GetNEdges();
        for (auto i : Range(ma->GetNE(VOL))) {
            DofId center = nv + ne + i;
            if (ma->GetElement(ElementId(VOL, i)).GetType() == ET_TRIG)
                ctofdof[center] = UNUSED_DOF;   // P2 triangle: no center DOF
            else
                ctofdof[center] = LOCAL_DOF;     // Q2 quad: element-local center
        }
    }
}

void AxiHenrotteFESpace::FinalizeUpdate() {
    FESpace::FinalizeUpdate();          // builds dirichlet_dofs, free_dofs, etc.
    if (axi_order != 2) return;
    UpdateCouplingDofArray();           // expose UNUSED trig-center coupling type
    // The base built free_dofs before ctofdof existed, so drop the dead
    // P2-triangle center slots (nv+ne+ei.Nr()) from the free-dof sets directly.
    int nv = ma->GetNV();
    int ne = ma->GetNEdges();
    for (auto i : Range(ma->GetNE(VOL))) {
        if (ma->GetElement(ElementId(VOL, i)).GetType() != ET_TRIG) continue;
        DofId c = nv + ne + i;
        if (free_dofs && c < free_dofs->Size()) free_dofs->Clear(c);
        if (external_free_dofs && c < external_free_dofs->Size())
            external_free_dofs->Clear(c);
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
        if (ngel.GetType() == ET_QUAD && vertices.Size() == 4) {
            if (curved_quad) {
                // Sample the 9 node positions from the (possibly Curve(2)-d)
                // element transformation, in our Q2 node order: 4 corners (v0-v3),
                // then bottom/right/top/left edge midnodes, then the face center.
                const double xi_refs[9]  = {0,1,1,0, 0.5,1,0.5,0, 0.5};
                const double eta_refs[9] = {0,0,1,1, 0,0.5,1,0.5, 0.5};
                auto & eltrans = ma->GetTrafo(ei, lh);
                double rs[9], zs[9];
                for (int k = 0; k < 9; ++k) {
                    IntegrationPoint ip(xi_refs[k], eta_refs[k], 0.0);
                    auto & mip = eltrans(ip, lh);
                    auto pt = mip.GetPoint();
                    rs[k] = pt(0); zs[k] = pt(1);
                }
                return *new (lh) AxiHenrotteFE_Q2_Curved(rs, zs);
            }
            Vec<3> p[4];
            for (int i = 0; i < 4; ++i) p[i] = ma->GetPoint<3>(vertices[i]);
            auto b = RequireAxisAlignedQuad(p, ei, 2);
            double r_a = b.r_min, r_b = b.r_max;
            double z_a = b.z_min, z_b = b.z_max;
            return *new (lh) AxiHenrotteFE_Q2_AxisAligned(r_a, r_b, z_a, z_b);
        }
        if (ngel.GetType() == ET_TRIG && vertices.Size() == 3) {
            // P2 triangle: 6 DOFs at 3 vertices + 3 edge midpoints.
            // Node order: [v0, v1, v2, m01, m12, m20].
            //
            // Edge midpoint positions are obtained via the mesh's element
            // transformation. For a straight-edge mesh this reduces to the
            // geometric midpoint of the chord; for a mesh that has been
            // `.Curve(p)`-d the transformation returns the curved-geometry
            // position of the mid-edge node so that the P2 FE can follow
            // curved boundaries (e.g. a sphere) exactly to the order of the
            // curve.
            //
            // NGSolve's reference triangle (verified 2026-05-12 in
            // examples/CLN/scripts/axifem/test_ngsolve_ref_tri_vertices.py):
            //   ref (1, 0) <-> mesh vertex 0 (V0)
            //   ref (0, 1) <-> mesh vertex 1 (V1)
            //   ref (0, 0) <-> mesh vertex 2 (V2)
            // So mid-edges are at:
            //   m01 = midpoint V0-V1 -> ref (0.5, 0.5)
            //   m12 = midpoint V1-V2 -> ref (0, 0.5)
            //   m20 = midpoint V2-V0 -> ref (0.5, 0)
            const double xi_refs[6]  = { 1.0, 0.0, 0.0, 0.5, 0.0, 0.5 };
            const double eta_refs[6] = { 0.0, 1.0, 0.0, 0.5, 0.5, 0.0 };
            auto & eltrans = ma->GetTrafo(ei, lh);
            double rs[6], zs[6];
            for (int k = 0; k < 6; ++k) {
                IntegrationPoint ip(xi_refs[k], eta_refs[k], 0.0);
                auto & mip = eltrans(ip, lh);
                auto pt = mip.GetPoint();
                rs[k] = pt(0);
                zs[k] = pt(1);
            }
            return *new (lh) AxiHenrotteFE_P2_Triangle(rs, zs);
        }
        throw Exception("AxiHenrotteFESpace order=2: element " + ToString(ei) +
                        " is neither a quad (Q2) nor a triangle (P2).");
    }
    // axi_order == 1 path
    if (ngel.GetType() == ET_QUAD && vertices.Size() == 4) {
        Vec<3> p[4];
        for (int i = 0; i < 4; ++i) p[i] = ma->GetPoint<3>(vertices[i]);
        auto b = RequireAxisAlignedQuad(p, ei, 1);
        double r_a = b.r_min, r_b = b.r_max;
        double z_a = b.z_min, z_b = b.z_max;
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
    // VOL element — quad (Q2, 9 DOFs) or triangle (P2, 6 DOFs).
    if (ngel.GetType() == ET_TRIG) {
        // P2 triangle: 3 vertex DOFs + 3 edge midpoint DOFs.
        // Local order [v0, v1, v2, m01, m12, m20] matches FE's CalcShape.
        auto edges = ngel.Edges();
        dnums.SetSize(6);
        for (int i = 0; i < 3; ++i) dnums[i] = vertices[i];
        // NGSolve's TRIG local edge order (verified empirically 2026-05-12):
        //   edges[0] connects (V0, V2)  -> midpoint = m20
        //   edges[1] connects (V1, V2)  -> midpoint = m12
        //   edges[2] connects (V0, V1)  -> midpoint = m01
        // Permute to our local-DOF order [m01, m12, m20]:
        dnums[3] = nv + edges[2];   // m01 = edge V0-V1
        dnums[4] = nv + edges[1];   // m12 = edge V1-V2
        dnums[5] = nv + edges[0];   // m20 = edge V0-V2
        return;
    }
    if (ngel.GetType() != ET_QUAD)
        throw Exception("AxiHenrotteFESpace order=2: VOL element must be quad or triangle");
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

void AxiHenrotteFESpace::GetVertexDofNrs(int vnr, Array<DofId> & dnums) const {
    dnums.SetSize(0);
    dnums.Append(vnr);
}

void AxiHenrotteFESpace::GetEdgeDofNrs(int ednr, Array<DofId> & dnums) const {
    dnums.SetSize(0);
    if (axi_order == 1) return;  // order==1 has no edge DOF
    const int nv = ma->GetNV();
    dnums.Append(nv + ednr);
}

void AxiHenrotteFESpace::GetFaceDofNrs(int fanr, Array<DofId> & dnums) const {
    dnums.SetSize(0);
    if (axi_order == 1) return;
    // 2D mesh: face-DOF slot per element (P2 trig: harmless unused slot,
    // Q2 quad: face-center DOF). Periodic 2D normally doesn't query NT_FACE
    // but we expose it for symmetry.
    const int nv = ma->GetNV();
    const int ne = ma->GetNEdges();
    dnums.Append(nv + ne + fanr);
}

void AxiHenrotteFESpace::GetDofNrs(NodeId ni, Array<DofId> & dnums) const {
    // Needed for ngsolve.Periodic(H1Henrotte(...)) to discover which DOFs
    // sit on each vertex / edge participating in a PERIODIC identification.
    // PeriodicFESpace::Update walks identifications and calls
    // space->GetDofNrs(NodeId(NT_VERTEX|NT_EDGE|NT_FACE, idx), dnums) on
    // the wrapped space; the base FESpace default returns an empty array
    // and silently disables coupling.
    //
    // DOF layout mirrors AxiHenrotteFESpace::Update():
    //   order==1: dnums[i] = vertex i,                  range [0,  NV)
    //   order==2: dnums[i] = vertex i,                  range [0,  NV)
    //             dnums[NV + j] = edge midnode j,       range [NV, NV+NE)
    //             dnums[NV + NE + k] = face slot k      range [NV+NE, NV+NE+NEl)
    //                                                   (P2 trig: unused slot,
    //                                                    Q2 quad: face-center DOF)
    // Periodic identification in 2D only walks NT_VERTEX and NT_EDGE, so the
    // face slot is never queried for 2D periodicity; we expose it anyway in
    // case a future 3D-axisymmetric extension wires Q2 face DOFs to NT_FACE.
    dnums.SetSize(0);
    if (axi_order == 1) {
        if (ni.GetType() == NT_VERTEX)
            dnums.Append(ni.GetNr());
        return;
    }
    // axi_order == 2
    const int nv = ma->GetNV();
    const int ne = ma->GetNEdges();
    if (ni.GetType() == NT_VERTEX) {
        dnums.Append(ni.GetNr());
    } else if (ni.GetType() == NT_EDGE) {
        dnums.Append(nv + ni.GetNr());
    } else if (ni.GetType() == NT_FACE) {
        // 2D-axisymmetric meshes treat each surface element as a "face" node;
        // P2 trig leaves the slot unused (harmless), Q2 quad uses it for the
        // face-center DOF.
        dnums.Append(nv + ne + ni.GetNr());
    }
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
              for (auto item : kwargs) {
                  std::string key = py::cast<std::string>(item.first);
                  py::handle val = item.second;
                  // Type-correct flag conversion.  The previous code stringified
                  // EVERY value, so complex=True became the STRING flag
                  // "complex"="True" -- which GetDefineFlag("complex") (a boolean
                  // DEFINE flag) does not see, leaving the space REAL.  bool must
                  // become a define flag, numbers a num flag (note: Python bool is
                  // a subclass of int, so test bool FIRST).
                  if (py::isinstance<py::bool_>(val)) {
                      if (py::cast<bool>(val))
                          flags.SetFlag(key);                 // boolean define flag
                  } else if (py::isinstance<py::int_>(val) ||
                             py::isinstance<py::float_>(val)) {
                      flags.SetFlag(key, py::cast<double>(val));
                  } else {
                      flags.SetFlag(key, py::cast<std::string>(py::str(val)));
                  }
              }
              auto fes = make_shared<AxiHenrotteFESpace>(ma, flags);
              fes->Update();
              fes->FinalizeUpdate();
              return fes;
          },
          py::arg("mesh"), py::arg("order") = 1,
          "Construct an AxiHenrotteFESpace for the given mesh.\n"
          "  order=1 (default): Q1 quad / P1 triangle, vertex DOFs only.\n"
          "  order=2          : Q2 axis-aligned quad (9 DOFs) or P2 triangle\n"
          "                     (6 DOFs; curved-mesh aware via mesh.Curve(2)).");
}

}  // namespace axifem
