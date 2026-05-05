// axi_henrotte_diffop.hpp — DiffOps wiring AxiHenrotteFE into NGSolve.
//
// NGSolve's BilinearForm machinery does not call FE methods directly; it
// goes through DifferentialOperator (DiffOp) classes registered as the
// FESpace's evaluator[VOL] / flux_evaluator[VOL]. We provide:
//
//   AxiHenrotteDiffOpId       — evaluates the scalar shape function value
//   AxiHenrotteDiffOpGradient — evaluates the gradient w.r.t. physical (r, z)
//
// Both cast the generic FiniteElement to AxiHenrotteBaseFE and call its
// CalcShape / CalcDShape, then apply the geometric Jacobian transform
// for the gradient.

#ifndef RADIA_AXIFEMM_AXI_HENROTTE_DIFFOP_HPP
#define RADIA_AXIFEMM_AXI_HENROTTE_DIFFOP_HPP

#include <fem.hpp>
#include "axi_henrotte_fe.hpp"

namespace radia_axifemm {

using namespace ngfem;

// ---------------------------------------------------------------------------
// AxiHenrotteDiffOpId — value evaluator
// ---------------------------------------------------------------------------

class AxiHenrotteDiffOpId : public DifferentialOperator {
public:
    AxiHenrotteDiffOpId()
        : DifferentialOperator(/*dim=*/1, /*dim_diffop=*/1,
                               /*vb=*/VOL, /*difforder=*/0) {}

    string Name() const override { return "AxiHenrotte:Id"; }

    void CalcMatrix(const FiniteElement & fel,
                    const BaseMappedIntegrationPoint & mip,
                    BareSliceMatrix<double, ColMajor> mat,
                    LocalHeap & lh) const override
    {
        const auto & afel = static_cast<const AxiHenrotteBaseFE&>(fel);
        IntegrationPoint ip = mip.IP();
        FlatVector<double> shape(afel.GetNDof(), lh);
        afel.CalcShape(ip, shape);
        for (size_t i = 0; i < afel.GetNDof(); ++i)
            mat(0, i) = shape(i);
    }
};

// ---------------------------------------------------------------------------
// AxiHenrotteDiffOpGradient — gradient evaluator (returns 2D vector in r, z)
// ---------------------------------------------------------------------------

class AxiHenrotteDiffOpGradient : public DifferentialOperator {
public:
    AxiHenrotteDiffOpGradient()
        : DifferentialOperator(/*dim=*/2, /*dim_diffop=*/2,
                               /*vb=*/VOL, /*difforder=*/1) {}

    string Name() const override { return "AxiHenrotte:Gradient"; }

    void CalcMatrix(const FiniteElement & fel,
                    const BaseMappedIntegrationPoint & mip,
                    BareSliceMatrix<double, ColMajor> mat,
                    LocalHeap & lh) const override
    {
        const auto & afel = static_cast<const AxiHenrotteBaseFE&>(fel);
        IntegrationPoint ip = mip.IP();
        FlatMatrix<double> dshape_ref(afel.GetNDof(), 2, lh);
        afel.CalcDShape(ip, dshape_ref);

        // Apply the geometric Jacobian transform: dphi/dx = (J^{-T}) * dphi/dxi
        // where x = (r, z) physical, xi = reference coords.
        const auto & mip2d = static_cast<const MappedIntegrationPoint<2,2>&>(mip);
        Mat<2,2> jac_inv = mip2d.GetJacobianInverse();  // J^{-1}

        // dphi/dx_a = sum_b (J^{-1})_{a,b}^T * dphi/dxi_b
        //           = sum_b (J^{-1})_{b,a} * dphi/dxi_b
        for (size_t i = 0; i < afel.GetNDof(); ++i) {
            mat(0, i) = jac_inv(0, 0) * dshape_ref(i, 0)
                      + jac_inv(1, 0) * dshape_ref(i, 1);
            mat(1, i) = jac_inv(0, 1) * dshape_ref(i, 0)
                      + jac_inv(1, 1) * dshape_ref(i, 1);
        }
    }
};

}  // namespace radia_axifemm

#endif  // RADIA_AXIFEMM_AXI_HENROTTE_DIFFOP_HPP
