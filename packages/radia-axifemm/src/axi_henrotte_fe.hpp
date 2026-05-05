// axi_henrotte_fe.hpp — Henrotte/Meeker axisymmetric finite element.
//
// Shape functions are linear in {1, r^2, z, ...} on physical (r, z) coordinates
// (NOT on the reference element). The basis is determined by the physical
// vertex positions of each element, which the FE stores at construction.
//
// We follow the NGSolve i-tutorial Unit 9.1 ("MyElement") pattern: inherit
// directly from FiniteElement (not BaseScalarFiniteElement) to avoid the
// SIMD/parallel-evaluator pure-virtual obligations. Our shape and derivative
// access goes through a custom AxiHenrotteBaseFE base class with its own
// virtual interface.

#ifndef RADIA_AXIFEMM_AXI_HENROTTE_FE_HPP
#define RADIA_AXIFEMM_AXI_HENROTTE_FE_HPP

#include <fem.hpp>

namespace radia_axifemm {

using namespace ngfem;

// ---------------------------------------------------------------------------
// AxiHenrotteBaseFE — minimal base for our custom shape-function family.
//   - Subclass of ngfem::FiniteElement (not BaseScalarFiniteElement, to avoid
//     the abstract SIMD evaluators).
//   - Adds CalcShape / CalcDShape virtual interface for our DiffOps.
// ---------------------------------------------------------------------------

class AxiHenrotteBaseFE : public FiniteElement {
public:
    AxiHenrotteBaseFE(int ndof, int order) : FiniteElement(ndof, order) {}

    /// Shape values at reference-element coords ip; returns ndof values.
    /// Our implementation maps ref -> phys (r, z) and evaluates {1, r^2, z, ...}.
    virtual void CalcShape(const IntegrationPoint & ip,
                           BareSliceVector<> shape) const = 0;

    /// Shape gradients w.r.t. reference coords, ndof x dim matrix.
    /// (NGSolve geometry transformation handles the chain rule to physical.)
    virtual void CalcDShape(const IntegrationPoint & ip,
                            BareSliceMatrix<> dshape) const = 0;
};

// ---------------------------------------------------------------------------
// AxiHenrotteFE_Q1_AxisAligned
//   - 4 DOFs/element on an axis-aligned rectangle [r_a, r_b] x [z_a, z_b]
//   - shape_i(r, z) = a_i + b_i*r^2 + c_i*z + d_i*r^2*z
//   - (a, b, c, d)_i fixed by shape_i(r_j, z_j) = delta_ij
// ---------------------------------------------------------------------------

class AxiHenrotteFE_Q1_AxisAligned : public AxiHenrotteBaseFE {
public:
    double r_a, r_b, z_a, z_b;

    AxiHenrotteFE_Q1_AxisAligned(double ra, double rb, double za, double zb);

    ELEMENT_TYPE ElementType() const override { return ET_QUAD; }

    void CalcShape(const IntegrationPoint & ip,
                   BareSliceVector<> shape) const override;
    void CalcDShape(const IntegrationPoint & ip,
                    BareSliceMatrix<> dshape) const override;
};

// ---------------------------------------------------------------------------
// AxiHenrotteFE_P1_Triangle
//   - 3 DOFs/element on a general triangle (r_i, z_i) in (r, z) plane
//   - shape_i(r, z) = a_i + b_i*r^2 + c_i*z
// ---------------------------------------------------------------------------

class AxiHenrotteFE_P1_Triangle : public AxiHenrotteBaseFE {
public:
    double r[3], z[3];
    double alpha[3], beta[3], gamma_[3];

    AxiHenrotteFE_P1_Triangle(const double rs[3], const double zs[3]);

    ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    void CalcShape(const IntegrationPoint & ip,
                   BareSliceVector<> shape) const override;
    void CalcDShape(const IntegrationPoint & ip,
                    BareSliceMatrix<> dshape) const override;
};

}  // namespace radia_axifemm

#endif  // RADIA_AXIFEMM_AXI_HENROTTE_FE_HPP
