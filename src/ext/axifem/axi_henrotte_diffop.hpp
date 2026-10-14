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

#ifndef RADIA_AXIFEM_AXI_HENROTTE_DIFFOP_HPP
#define RADIA_AXIFEM_AXI_HENROTTE_DIFFOP_HPP

#include <fem.hpp>
#include "axi_henrotte_fe.hpp"

namespace axifem {

using namespace ngfem;

// ---------------------------------------------------------------------------
// AxiHenrotteDiffOpId — value evaluator
//
// DifferentialOperator(dim, blockdim, vb, difforder):
//   dim      = dimension of range (output): 1 for scalar value
//   blockdim = block multiplicity for BlockDifferentialOperator (1 for scalar FE)
//   vb       = VOL / BND
//   difforder = differentiation order
// ---------------------------------------------------------------------------

class AxiHenrotteDiffOpId : public DifferentialOperator {
public:
    AxiHenrotteDiffOpId()
        : DifferentialOperator(/*dim=*/1, /*blockdim=*/1,
                               /*vb=*/VOL, /*difforder=*/0)
    {
        // Match standard H1: scalar field → empty Dimensions (not [1]).
        SetDimensions(Array<int>{});
    }

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

    // COMPLEX overload: the shape values are real, but a COMPLEX GridFunction
    // (time-harmonic A_phi) evaluates through this overload.  Without it NGSolve
    // falls back to the base-class stub and returns wrong values ("base class
    // apply").  Body mirrors the real CalcMatrix; double -> Complex is implicit.
    void CalcMatrix(const FiniteElement & fel,
                    const BaseMappedIntegrationPoint & mip,
                    BareSliceMatrix<Complex, ColMajor> mat,
                    LocalHeap & lh) const override
    {
        const auto & afel = static_cast<const AxiHenrotteBaseFE&>(fel);
        IntegrationPoint ip = mip.IP();
        FlatVector<double> shape(afel.GetNDof(), lh);
        afel.CalcShape(ip, shape);
        for (size_t i = 0; i < afel.GetNDof(); ++i)
            mat(0, i) = shape(i);
    }

    void ApplyTrans(const FiniteElement & fel,
                    const BaseMappedIntegrationPoint & mip,
                    FlatVector<double> flux,
                    BareSliceVector<double> values,
                    LocalHeap & lh) const override
    {
        const auto & afel = static_cast<const AxiHenrotteBaseFE&>(fel);
        FlatVector<double> shape(afel.GetNDof(), lh);
        afel.CalcShape(mip.IP(), shape);
        for (size_t i = 0; i < afel.GetNDof(); ++i)
            values(i) = shape(i) * flux(0);
    }

    void ApplyTrans(const FiniteElement & fel,
                    const BaseMappedIntegrationPoint & mip,
                    FlatVector<Complex> flux,
                    BareSliceVector<Complex> values,
                    LocalHeap & lh) const override
    {
        const auto & afel = static_cast<const AxiHenrotteBaseFE&>(fel);
        FlatVector<double> shape(afel.GetNDof(), lh);
        afel.CalcShape(mip.IP(), shape);
        for (size_t i = 0; i < afel.GetNDof(); ++i)
            values(i) = shape(i) * flux(0);
    }
};

// ---------------------------------------------------------------------------
// AxiHenrotteDiffOpIdBnd — boundary-trace value evaluator
//
// NGSolve's SymbolicLinearForm checks DifferentialOperator::VB() to decide
// whether a given evaluator can serve as a test function's trace.  We need
// a separate DiffOp with vb=BND for `LinearForm += q * v * ds(label)` to
// work on H1Henrotte.  CalcMatrix is identical to the VOL variant -- our
// AxiHenrotteFE_Edge_{Q1,Q2} CalcShape already returns the appropriate
// 1D trace values.  The 2026-05-10 Henrotte-only-axisym policy made this
// class load-bearing for axisymmetric Neumann RHS assembly.
// ---------------------------------------------------------------------------

class AxiHenrotteDiffOpIdBnd : public DifferentialOperator {
public:
    AxiHenrotteDiffOpIdBnd()
        : DifferentialOperator(/*dim=*/1, /*blockdim=*/1,
                               /*vb=*/BND, /*difforder=*/0)
    {
        SetDimensions(Array<int>{});
    }

    string Name() const override { return "AxiHenrotte:Id:Bnd"; }

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

    // COMPLEX overload (see AxiHenrotteDiffOpId): boundary-trace value eval for
    // a complex GridFunction.  Mirrors the real CalcMatrix above.
    void CalcMatrix(const FiniteElement & fel,
                    const BaseMappedIntegrationPoint & mip,
                    BareSliceMatrix<Complex, ColMajor> mat,
                    LocalHeap & lh) const override
    {
        const auto & afel = static_cast<const AxiHenrotteBaseFE&>(fel);
        IntegrationPoint ip = mip.IP();
        FlatVector<double> shape(afel.GetNDof(), lh);
        afel.CalcShape(ip, shape);
        for (size_t i = 0; i < afel.GetNDof(); ++i)
            mat(0, i) = shape(i);
    }

    void ApplyTrans(const FiniteElement & fel,
                    const BaseMappedIntegrationPoint & mip,
                    FlatVector<double> flux,
                    BareSliceVector<double> values,
                    LocalHeap & lh) const override
    {
        const auto & afel = static_cast<const AxiHenrotteBaseFE&>(fel);
        FlatVector<double> shape(afel.GetNDof(), lh);
        afel.CalcShape(mip.IP(), shape);
        for (size_t i = 0; i < afel.GetNDof(); ++i)
            values(i) = shape(i) * flux(0);
    }

    void ApplyTrans(const FiniteElement & fel,
                    const BaseMappedIntegrationPoint & mip,
                    FlatVector<Complex> flux,
                    BareSliceVector<Complex> values,
                    LocalHeap & lh) const override
    {
        const auto & afel = static_cast<const AxiHenrotteBaseFE&>(fel);
        FlatVector<double> shape(afel.GetNDof(), lh);
        afel.CalcShape(mip.IP(), shape);
        for (size_t i = 0; i < afel.GetNDof(); ++i)
            values(i) = shape(i) * flux(0);
    }
};

// ---------------------------------------------------------------------------
// AxiHenrotteDiffOpGradient — gradient evaluator (returns 2D vector in r, z)
// ---------------------------------------------------------------------------

class AxiHenrotteDiffOpGradient : public DifferentialOperator {
public:
    AxiHenrotteDiffOpGradient()
        : DifferentialOperator(/*dim=*/2, /*blockdim=*/1,
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

    // COMPLEX overload: gradient (-> B) of a complex time-harmonic A_phi.  The
    // shape derivatives and Jacobian are real; only the DOF coefficients are
    // complex, so the matrix entries are identical to the real CalcMatrix.
    // Required so curl/grad of a COMPLEX GridFunction evaluates correctly
    // instead of hitting the base-class stub.
    void CalcMatrix(const FiniteElement & fel,
                    const BaseMappedIntegrationPoint & mip,
                    BareSliceMatrix<Complex, ColMajor> mat,
                    LocalHeap & lh) const override
    {
        const auto & afel = static_cast<const AxiHenrotteBaseFE&>(fel);
        IntegrationPoint ip = mip.IP();
        FlatMatrix<double> dshape_ref(afel.GetNDof(), 2, lh);
        afel.CalcDShape(ip, dshape_ref);

        const auto & mip2d = static_cast<const MappedIntegrationPoint<2,2>&>(mip);
        Mat<2,2> jac_inv = mip2d.GetJacobianInverse();

        for (size_t i = 0; i < afel.GetNDof(); ++i) {
            mat(0, i) = jac_inv(0, 0) * dshape_ref(i, 0)
                      + jac_inv(1, 0) * dshape_ref(i, 1);
            mat(1, i) = jac_inv(0, 1) * dshape_ref(i, 0)
                      + jac_inv(1, 1) * dshape_ref(i, 1);
        }
    }
};

}  // namespace axifem

#endif  // RADIA_AXIFEM_AXI_HENROTTE_DIFFOP_HPP
