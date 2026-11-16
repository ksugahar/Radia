// axi_henrotte_diffop.hpp — DiffOps wiring AxiHenrotteFE into NGSolve.
//
// NGSolve's BilinearForm/LinearForm/GridFunction machinery evaluates the space
// through these DifferentialOperators:
//
//   AxiHenrotteDiffOpId       — A_phi basis N_i = r_i psi_i / r (VOL)
//   AxiHenrotteDiffOpIdBnd    — the same basis on boundary segments (BND)
//   AxiHenrotteDiffOpGradient — physical (dN_i/dr, dN_i/dz) (VOL)
//
// All values come from AxiHenrotteBaseFE::CalcVdof, the basis of the native
// stiffness and sigma-mass integrators, so symbolic loads, traces and
// post-processing use the field convention of the assembled matrices.  The
// gradient is returned in physical coordinates by the element itself (the
// s-isoparametric map is not NGSolve's), so no Jacobian of the mesh map is
// applied here.  The complex overloads serve complex (time-harmonic)
// GridFunctions: the basis is real, only the coefficients are complex.
//
// Because the element evaluates its own map, every entry point first requires
// that the integration point was mapped by the mesh's own transformation
// (RequireMeshMap).  An integration-local deformation (dx(deformation=gf))
// maps the points and the measure but cannot change the element's basis, so
// it is refused instead of returning the undeformed gradient.

#ifndef RADIA_AXIFEM_AXI_HENROTTE_DIFFOP_HPP
#define RADIA_AXIFEM_AXI_HENROTTE_DIFFOP_HPP

#include <comp.hpp>
#include <algorithm>
#include <cmath>
#include "axi_henrotte_fe.hpp"

namespace axifem {

using namespace ngfem;

// Refuse an integration point whose transformation is not the mesh map: the
// point and the Jacobian are compared with those of the mesh's own
// transformation at the same reference point.  The mesh path computes both
// with the same code, so they agree to round-off.
inline void RequireMeshMap(const ngcomp::MeshAccess & ma,
                           const BaseMappedIntegrationPoint & mip, LocalHeap & lh)
{
    const ElementTransformation & trafo = mip.GetTransformation();
    if (!trafo.BelongsToMesh(&ma))
        throw Exception("AxiHenrotte: the integration point belongs to another mesh");
    HeapReset hr(lh);
    const ElementTransformation & mesh_trafo = ma.GetTrafo(trafo.GetElementId(), lh);
    if (&mesh_trafo == &trafo) return;
    const BaseMappedIntegrationPoint & ref = mesh_trafo(mip.IP(), lh);
    if (mip.DimSpace() != 2 || ref.DimSpace() != 2 || mip.DimElement() != ref.DimElement())
        throw Exception("AxiHenrotte: the integration point is not on the 2D (r, z) mesh map");
    double scale = 0.0, dj = 0.0, dx = 0.0;
    auto compare = [&](const auto & a, const auto & b) {
        const auto & J = a.GetJacobian();
        const auto & Jr = b.GetJacobian();
        for (int i = 0; i < J.Height(); ++i)
            for (int j = 0; j < J.Width(); ++j) {
                scale = std::max(scale, std::abs(Jr(i, j)));
                dj = std::max(dj, std::abs(J(i, j) - Jr(i, j)));
            }
        for (int i = 0; i < 2; ++i)
            dx = std::max(dx, std::abs(a.GetPoint()(i) - b.GetPoint()(i)));
    };
    if (mip.DimElement() == 2)
        compare(static_cast<const MappedIntegrationPoint<2, 2> &>(mip),
                static_cast<const MappedIntegrationPoint<2, 2> &>(ref));
    else
        compare(static_cast<const MappedIntegrationPoint<1, 2> &>(mip),
                static_cast<const MappedIntegrationPoint<1, 2> &>(ref));
    const double tol = 1e-12 * scale;
    if (dx > tol || dj > tol)
        throw Exception("AxiHenrotte: the integration point is not mapped by the mesh"
                        " transformation (an integration-local deformation such as"
                        " dx(deformation=...) is not supported); point offset "
                        + ToString(dx) + ", Jacobian offset " + ToString(dj));
}

template <bool kBnd>
class AxiHenrotteDiffOpValue : public DifferentialOperator {
public:
    explicit AxiHenrotteDiffOpValue(shared_ptr<ngcomp::MeshAccess> ma)
        : DifferentialOperator(/*dim=*/1, /*blockdim=*/1, kBnd ? BND : VOL, /*difforder=*/0),
          ma_(std::move(ma))
    {
        SetDimensions(Array<int>{});   // scalar field, as standard H1
    }

    string Name() const override { return kBnd ? "AxiHenrotte:Id:Bnd" : "AxiHenrotte:Id"; }

    template <typename T>
    static void Fill(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
                     BareSliceMatrix<T, ColMajor> mat)
    {
        const auto & afel = static_cast<const AxiHenrotteBaseFE&>(fel);
        double n[AxiHenrotteBaseFE::kMaxDof];
        afel.CalcVdof(mip.IP(), n, nullptr);
        for (size_t i = 0; i < afel.GetNDof(); ++i) mat(0, i) = n[i];
    }

    template <typename T>
    static void Trans(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
                      FlatVector<T> flux, BareSliceVector<T> values)
    {
        const auto & afel = static_cast<const AxiHenrotteBaseFE&>(fel);
        double n[AxiHenrotteBaseFE::kMaxDof];
        afel.CalcVdof(mip.IP(), n, nullptr);
        for (size_t i = 0; i < afel.GetNDof(); ++i) values(i) = n[i] * flux(0);
    }

    void CalcMatrix(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
                    BareSliceMatrix<double, ColMajor> mat, LocalHeap & lh) const override
    { RequireMeshMap(*ma_, mip, lh); Fill(fel, mip, mat); }

    void CalcMatrix(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
                    BareSliceMatrix<Complex, ColMajor> mat, LocalHeap & lh) const override
    { RequireMeshMap(*ma_, mip, lh); Fill(fel, mip, mat); }

    void ApplyTrans(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
                    FlatVector<double> flux, BareSliceVector<double> values,
                    LocalHeap & lh) const override
    { RequireMeshMap(*ma_, mip, lh); Trans(fel, mip, flux, values); }

    void ApplyTrans(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
                    FlatVector<Complex> flux, BareSliceVector<Complex> values,
                    LocalHeap & lh) const override
    { RequireMeshMap(*ma_, mip, lh); Trans(fel, mip, flux, values); }

    template <typename T>
    static void Forward(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
                        BareSliceVector<T> x, FlatVector<T> flux)
    {
        const auto & afel = static_cast<const AxiHenrotteBaseFE&>(fel);
        double n[AxiHenrotteBaseFE::kMaxDof];
        afel.CalcVdof(mip.IP(), n, nullptr);
        T v = 0.0;
        for (size_t i = 0; i < afel.GetNDof(); ++i) v += n[i] * x(i);
        flux(0) = v;
    }

    void Apply(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
               BareSliceVector<double> x, FlatVector<double> flux, LocalHeap & lh) const override
    { RequireMeshMap(*ma_, mip, lh); Forward(fel, mip, x, flux); }

    void Apply(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
               BareSliceVector<Complex> x, FlatVector<Complex> flux, LocalHeap & lh) const override
    { RequireMeshMap(*ma_, mip, lh); Forward(fel, mip, x, flux); }

    // No SIMD evaluation: NGSolve falls back to the point-wise path on
    // ExceptionNOSIMD (the base class would print a warning first).
    void Apply(const FiniteElement &, const SIMD_BaseMappedIntegrationRule &,
               BareSliceVector<double>, BareSliceMatrix<SIMD<double>>) const override
    { throw ExceptionNOSIMD("AxiHenrotte: no SIMD evaluation"); }
    void Apply(const FiniteElement &, const SIMD_BaseMappedIntegrationRule &,
               BareSliceVector<Complex>, BareSliceMatrix<SIMD<Complex>>) const override
    { throw ExceptionNOSIMD("AxiHenrotte: no SIMD evaluation"); }
    void AddTrans(const FiniteElement &, const SIMD_BaseMappedIntegrationRule &,
                  BareSliceMatrix<SIMD<double>>, BareSliceVector<double>) const override
    { throw ExceptionNOSIMD("AxiHenrotte: no SIMD evaluation"); }
    void AddTrans(const FiniteElement &, const SIMD_BaseMappedIntegrationRule &,
                  BareSliceMatrix<SIMD<Complex>>, BareSliceVector<Complex>) const override
    { throw ExceptionNOSIMD("AxiHenrotte: no SIMD evaluation"); }

private:
    shared_ptr<ngcomp::MeshAccess> ma_;
};

using AxiHenrotteDiffOpId = AxiHenrotteDiffOpValue<false>;
using AxiHenrotteDiffOpIdBnd = AxiHenrotteDiffOpValue<true>;

class AxiHenrotteDiffOpGradient : public DifferentialOperator {
public:
    explicit AxiHenrotteDiffOpGradient(shared_ptr<ngcomp::MeshAccess> ma)
        : DifferentialOperator(/*dim=*/2, /*blockdim=*/1, VOL, /*difforder=*/1),
          ma_(std::move(ma)) {}

    string Name() const override { return "AxiHenrotte:Gradient"; }

    template <typename T>
    static void Fill(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
                     BareSliceMatrix<T, ColMajor> mat)
    {
        const auto & afel = static_cast<const AxiHenrotteBaseFE&>(fel);
        double n[AxiHenrotteBaseFE::kMaxDof], dn[AxiHenrotteBaseFE::kMaxDof][2];
        afel.CalcVdof(mip.IP(), n, dn);
        for (size_t i = 0; i < afel.GetNDof(); ++i) {
            mat(0, i) = dn[i][0];
            mat(1, i) = dn[i][1];
        }
    }

    void CalcMatrix(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
                    BareSliceMatrix<double, ColMajor> mat, LocalHeap & lh) const override
    { RequireMeshMap(*ma_, mip, lh); Fill(fel, mip, mat); }

    void CalcMatrix(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
                    BareSliceMatrix<Complex, ColMajor> mat, LocalHeap & lh) const override
    { RequireMeshMap(*ma_, mip, lh); Fill(fel, mip, mat); }

    template <typename T>
    static void Forward(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
                        BareSliceVector<T> x, FlatVector<T> flux)
    {
        const auto & afel = static_cast<const AxiHenrotteBaseFE&>(fel);
        double n[AxiHenrotteBaseFE::kMaxDof], dn[AxiHenrotteBaseFE::kMaxDof][2];
        afel.CalcVdof(mip.IP(), n, dn);
        T gr = 0.0, gz = 0.0;
        for (size_t i = 0; i < afel.GetNDof(); ++i) { gr += dn[i][0] * x(i); gz += dn[i][1] * x(i); }
        flux(0) = gr; flux(1) = gz;
    }

    void Apply(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
               BareSliceVector<double> x, FlatVector<double> flux, LocalHeap & lh) const override
    { RequireMeshMap(*ma_, mip, lh); Forward(fel, mip, x, flux); }

    void Apply(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
               BareSliceVector<Complex> x, FlatVector<Complex> flux, LocalHeap & lh) const override
    { RequireMeshMap(*ma_, mip, lh); Forward(fel, mip, x, flux); }

    // No SIMD evaluation: NGSolve falls back to the point-wise path on
    // ExceptionNOSIMD (the base class would print a warning first).
    void Apply(const FiniteElement &, const SIMD_BaseMappedIntegrationRule &,
               BareSliceVector<double>, BareSliceMatrix<SIMD<double>>) const override
    { throw ExceptionNOSIMD("AxiHenrotte: no SIMD evaluation"); }
    void Apply(const FiniteElement &, const SIMD_BaseMappedIntegrationRule &,
               BareSliceVector<Complex>, BareSliceMatrix<SIMD<Complex>>) const override
    { throw ExceptionNOSIMD("AxiHenrotte: no SIMD evaluation"); }
    void AddTrans(const FiniteElement &, const SIMD_BaseMappedIntegrationRule &,
                  BareSliceMatrix<SIMD<double>>, BareSliceVector<double>) const override
    { throw ExceptionNOSIMD("AxiHenrotte: no SIMD evaluation"); }
    void AddTrans(const FiniteElement &, const SIMD_BaseMappedIntegrationRule &,
                  BareSliceMatrix<SIMD<Complex>>, BareSliceVector<Complex>) const override
    { throw ExceptionNOSIMD("AxiHenrotte: no SIMD evaluation"); }

private:
    shared_ptr<ngcomp::MeshAccess> ma_;
};

}  // namespace axifem

#endif  // RADIA_AXIFEM_AXI_HENROTTE_DIFFOP_HPP
