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

#ifndef RADIA_AXIFEM_AXI_HENROTTE_DIFFOP_HPP
#define RADIA_AXIFEM_AXI_HENROTTE_DIFFOP_HPP

#include <fem.hpp>
#include "axi_henrotte_fe.hpp"

namespace axifem {

using namespace ngfem;

template <bool kBnd>
class AxiHenrotteDiffOpValue : public DifferentialOperator {
public:
    AxiHenrotteDiffOpValue()
        : DifferentialOperator(/*dim=*/1, /*blockdim=*/1, kBnd ? BND : VOL, /*difforder=*/0)
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
                    BareSliceMatrix<double, ColMajor> mat, LocalHeap &) const override
    { Fill(fel, mip, mat); }

    void CalcMatrix(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
                    BareSliceMatrix<Complex, ColMajor> mat, LocalHeap &) const override
    { Fill(fel, mip, mat); }

    void ApplyTrans(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
                    FlatVector<double> flux, BareSliceVector<double> values,
                    LocalHeap &) const override
    { Trans(fel, mip, flux, values); }

    void ApplyTrans(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
                    FlatVector<Complex> flux, BareSliceVector<Complex> values,
                    LocalHeap &) const override
    { Trans(fel, mip, flux, values); }

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
               BareSliceVector<double> x, FlatVector<double> flux, LocalHeap &) const override
    { Forward(fel, mip, x, flux); }

    void Apply(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
               BareSliceVector<Complex> x, FlatVector<Complex> flux, LocalHeap &) const override
    { Forward(fel, mip, x, flux); }

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
};

using AxiHenrotteDiffOpId = AxiHenrotteDiffOpValue<false>;
using AxiHenrotteDiffOpIdBnd = AxiHenrotteDiffOpValue<true>;

class AxiHenrotteDiffOpGradient : public DifferentialOperator {
public:
    AxiHenrotteDiffOpGradient()
        : DifferentialOperator(/*dim=*/2, /*blockdim=*/1, VOL, /*difforder=*/1) {}

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
                    BareSliceMatrix<double, ColMajor> mat, LocalHeap &) const override
    { Fill(fel, mip, mat); }

    void CalcMatrix(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
                    BareSliceMatrix<Complex, ColMajor> mat, LocalHeap &) const override
    { Fill(fel, mip, mat); }

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
               BareSliceVector<double> x, FlatVector<double> flux, LocalHeap &) const override
    { Forward(fel, mip, x, flux); }

    void Apply(const FiniteElement & fel, const BaseMappedIntegrationPoint & mip,
               BareSliceVector<Complex> x, FlatVector<Complex> flux, LocalHeap &) const override
    { Forward(fel, mip, x, flux); }

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
};

}  // namespace axifem

#endif  // RADIA_AXIFEM_AXI_HENROTTE_DIFFOP_HPP
