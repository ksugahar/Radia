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
// AxiHenrotteFE_Q2_AxisAligned
//   - 9 DOFs/element on an axis-aligned rectangle [r_a, r_b] x [z_a, z_b]
//   - shape monomials (in s = r^2, z):
//       {1, s, s^2, z, s*z, s^2*z, z^2, s*z^2, s^2*z^2}
//     (Lagrange-interpolatory at 9 nodes; s-midpoint convention).
//
//   Node order (matches JSON `node_order_general`):
//      0: (sa, za)    1: (sb, za)    2: (sb, zb)    3: (sa, zb)    (corners)
//      4: (sm, za)    5: (sb, zm)    6: (sm, zb)    7: (sa, zm)    (edge mid)
//      8: (sm, zm)                                                  (face center)
//     where sm = (sa+sb)/2 and zm = (za+zb)/2 — NOTE the s-midpoint
//     convention: physical r at edge midpoints is sqrt(sm) = sqrt((ra^2+rb^2)/2),
//     not (ra+rb)/2.
//
//   Axis-touching case (ra < EPS_AXIS):
//     6 monomials {s, s^2, s*z, s^2*z, s*z^2, s^2*z^2} on 6 non-axis nodes.
//     The 3 axis-side nodes (0, 3, 7) get zero shape functions and zero
//     gradient; physically A_phi vanishes on the symmetry axis.
//     Caller MUST Dirichlet the 3 axis-side global DOFs to zero.
// ---------------------------------------------------------------------------

class AxiHenrotteFE_Q2_AxisAligned : public AxiHenrotteBaseFE {
public:
    double r_a, r_b, z_a, z_b;
    bool   is_axis;            // sa < EPS_AXIS triggers the 6-monomial axis basis
    // Cached Vandermonde inverse (column j = coefficients of Lagrange basis L_j
    // in the monomial basis). For interior elements: 9x9 dense. For axis: 6x6
    // padded into the same 9x9 storage (rows/cols at axis-node indices are 0).
    double Vinv[9][9];
    // Local-index mapping for axis case: which 6 of the 9 nodes are non-axis.
    // For axis-touching: {1, 2, 4, 5, 6, 8}. For interior: identity {0..8}.
    int    nz_idx[9];          // first n_nz entries are valid
    int    n_nz;               // 9 (interior) or 6 (axis)

    AxiHenrotteFE_Q2_AxisAligned(double ra, double rb, double za, double zb);

    ELEMENT_TYPE ElementType() const override { return ET_QUAD; }

    void CalcShape(const IntegrationPoint & ip,
                   BareSliceVector<> shape) const override;
    void CalcDShape(const IntegrationPoint & ip,
                    BareSliceMatrix<> dshape) const override;
};

// ---------------------------------------------------------------------------
// AxiHenrotteFE_Q3_AxisAligned
//   - 16 DOFs/element on an axis-aligned rectangle [r_a, r_b] x [z_a, z_b]
//   - Tensor-product Lagrange in (s = r^2, z) on a 4x4 grid; monomial basis
//     ordered with `a` (s-power) major, `b` (z-power) minor:
//        {1, z, z^2, z^3, s, sz, sz^2, sz^3, s^2, s^2z, s^2z^2, s^2z^3,
//         s^3, s^3z, s^3z^2, s^3z^3}
//   - 16 nodes (4 corners + 8 edge midnodes + 4 face interior). Tensor-product
//     local indexing (matches NGSolve QUAD edge order [bottom, top, left, right])
//     -- no permutation in GetDofNrs.
//
//   Local node order:
//     0..3:  4 corners (NGSolve QUAD vertex order: (sa,za),(sb,za),(sb,zb),(sa,zb))
//     4-5:   edge 0 = bottom (z=za)  midnodes at s=s_t1, s_t2  (low->high s)
//     6-7:   edge 1 = top    (z=zb)  midnodes at s=s_t1, s_t2  (low->high s)
//     8-9:   edge 2 = left   (s=sa)  midnodes at z=z_t1, z_t2  (low->high z)
//     10-11: edge 3 = right  (s=sb)  midnodes at z=z_t1, z_t2  (low->high z)
//     12-15: face interior   at (s_t1,z_t1), (s_t2,z_t1), (s_t1,z_t2), (s_t2,z_t2)
//   where s_t1 = (2sa+sb)/3, s_t2 = (sa+2sb)/3, z_t1 = (2za+zb)/3, z_t2 = (za+2zb)/3.
//
//   Axis-touching case (ra < EPS_AXIS): 12-monomial restricted basis
//     {s^a z^b : 1 <= a <= 3, 0 <= b <= 3}; the 4 axis-side nodes (local
//     indices 0, 3, 8, 9) get zero shape functions.
//
//   *** CONDITIONING WARNING ***
//   The 16x16 Vandermonde for the raw {s^a z^b} basis has a condition number
//   that scales as ~(sb/sa)^6 / (sb-sa)^6, easily exceeding 1e30 for refined
//   meshes (where (sb-sa) is small). The matrix InvertNxN<16> via Gauss-Jordan
//   then produces noise that contaminates the assembled K_V and M_V.
//   On a Cu-disk Hiruma test the *coarse* mesh happens to give the right
//   tau_1 because the dominant mode is robust to noise, but mesh refinement
//   makes results worse, not better. To use Q3 in production, switch to an
//   orthogonal basis (shifted Legendre on [sa, sb] x [za, zb]) -- see TODO.
// ---------------------------------------------------------------------------

class AxiHenrotteFE_Q3_AxisAligned : public AxiHenrotteBaseFE {
public:
    double r_a, r_b, z_a, z_b;
    bool   is_axis;
    double Vinv[16][16];        // padded; axis case fills only 12x12 (rows 0..11, cols nz_idx)
    int    nz_idx[16];          // active (non-zero-shape) local DOF indices
    int    n_nz;                // 16 (interior) or 12 (axis)

    AxiHenrotteFE_Q3_AxisAligned(double ra, double rb, double za, double zb);

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
