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

#ifndef RADIA_AXIFEM_AXI_HENROTTE_FE_HPP
#define RADIA_AXIFEM_AXI_HENROTTE_FE_HPP

#include <fem.hpp>

namespace axifem {

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
// AxiHenrotteFE_Edge_Q1 / _Q2
//
// 1-D restriction of the Q1 / Q2 axis-aligned-quad Henrotte basis to a
// single edge (axis-aligned in (r, z)):
//
//   - Horizontal edge (z = const): basis values are Lagrange in s = r^2
//     across the edge endpoints.
//   - Vertical edge (r = const): basis values are Lagrange in z.
//
// This is needed for `LinearForm += q * v * ds(label)` patterns -- the
// axisymmetric Neumann-BC RHS of heat / magnetic Henrotte solvers.
//
// Edge endpoint convention: (r0, z0) and (r1, z1) are the two segment
// vertices in the order returned by `ma->GetElement(BND ei).Vertices()`.
// The order does NOT need to match the parent quad's vertex numbering
// because Lagrange basis values depend only on the edge endpoint
// coordinates, not on the parent quad's other two corners.
//
// Q1 case: 2 DOFs (the two segment vertices).
// Q2 case: 3 DOFs (vertex0, vertex1, edge midnode at the s-midpoint
// (= NGSolve's segment edge index, listed after the vertices in
// AxiHenrotteFESpace::GetDofNrs for BND elements)).
// ---------------------------------------------------------------------------

class AxiHenrotteFE_Edge_Q1 : public AxiHenrotteBaseFE {
public:
    double r0, z0, r1, z1;

    AxiHenrotteFE_Edge_Q1(double ar0, double az0, double ar1, double az1)
        : AxiHenrotteBaseFE(2, 1),
          r0(ar0), z0(az0), r1(ar1), z1(az1) {}

    ELEMENT_TYPE ElementType() const override { return ET_SEGM; }

    void CalcShape(const IntegrationPoint & ip,
                   BareSliceVector<> shape) const override;
    void CalcDShape(const IntegrationPoint & ip,
                    BareSliceMatrix<> dshape) const override;
};

class AxiHenrotteFE_Edge_Q2 : public AxiHenrotteBaseFE {
public:
    double r0, z0, r1, z1;

    AxiHenrotteFE_Edge_Q2(double ar0, double az0, double ar1, double az1)
        : AxiHenrotteBaseFE(3, 2),
          r0(ar0), z0(az0), r1(ar1), z1(az1) {}

    ELEMENT_TYPE ElementType() const override { return ET_SEGM; }

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

// ---------------------------------------------------------------------------
// AxiHenrotteFE_P2_Triangle
//   - 6 DOFs/element on a P2 triangle in (r, z) plane (3 vertex + 3 edge mid)
//   - psi_i(r, z) = a_i + b_i*r^2 + c_i*z + d_i*r^4 + e_i*r^2*z + f_i*z^2
//     (Lagrange-interpolatory at 6 nodes: psi_i(r_k, z_k) = delta_ik)
//   - V-DOF basis: N_A_i(r, z) = r_i * psi_i(r, z) / r
//
//   Node order: [v0, v1, v2, m01, m12, m20] where v_i are vertices and
//   m_ij is the midpoint of edge (v_i, v_j). This matches the standard
//   NGSolve P2 layout (vertex DOFs first, then edge DOFs).
//
//   For axis vertices (r_i = 0), N_A_i = 0 trivially due to the r_i factor;
//   the corresponding global DOF is auto-decoupled (zero row/column in K, M)
//   and is naturally pinned to A_phi = 0 (axis regularity).
//
//   Phase B1c reference: w:/.../axifem/axifem_p2_triangle.py
// ---------------------------------------------------------------------------

class AxiHenrotteFE_P2_Triangle : public AxiHenrotteBaseFE {
public:
    double r[6], z[6];
    // psi_i(r, z) = sum_j Vinv[j, i] * m_j(r, z), where m = (1, s, z, s^2, s*z, z^2)
    // and s = r^2. Vinv is the inverse of the 6x6 Vandermonde matrix on the
    // 6 P2 nodes. Storing as Vinv[6][6] in row-major.
    double Vinv[6][6];

    AxiHenrotteFE_P2_Triangle(const double rs[6], const double zs[6]);

    ELEMENT_TYPE ElementType() const override { return ET_TRIG; }

    void CalcShape(const IntegrationPoint & ip,
                   BareSliceVector<> shape) const override;
    void CalcDShape(const IntegrationPoint & ip,
                    BareSliceMatrix<> dshape) const override;
};

// ---------------------------------------------------------------------------
// AxiHenrotteFE_Q2_Curved
//   - 9-DOF CURVED biquadratic Henrotte quad on a GENERAL 9-node element in the
//     physical (r, z) plane (the proper curved isoparametric element -- Phase B6
//     port of examples/.../axifem/axifem_quad_q2_curved.py).  Unlike
//     AxiHenrotteFE_Q2_AxisAligned (axis-aligned rectangle, closed-form matrices),
//     this follows curved boundaries to O(curve^2): the FESpace samples the 9 node
//     positions from the mesh's (curved) ElementTransformation, exactly as the P2
//     triangle does for curved triangles.
//   - psi_i(r, z) = sum_j coeffs[j][i] * m_j(r', z') with the even-r-power monomials
//     m = {1, s, s^2, z, sz, s^2 z, z^2, sz^2, s^2 z^2}, s = (r/L_r)^2 (scaled for
//     conditioning), z' = (z - z_c)/L_z; coeffs = Vandermonde^{-1} over the 9 nodes.
//   - geometric map: biquadratic Lagrange on the NGSolve quad ref [0,1]^2 through
//     the 9 nodes (corners 0-3, edge mids 4-7, face center 8), matching the
//     AxiHenrotteFESpace Q2 node order (4=bottom,5=right,6=top,7=left,8=center).
//   - axis-touching nodes (r_i = 0) get the r_i factor in the V-DOF matrices, so the
//     caller MUST Dirichlet the axis-side DOFs (same contract as the other elements).
// ---------------------------------------------------------------------------

class AxiHenrotteFE_Q2_Curved : public AxiHenrotteBaseFE {
public:
    double rn[9], zn[9];        // physical node coords (FEMM/our Q2 node order)
    double L_r, z_c, L_z;       // conditioning scales (s' = (r/L_r)^2, z' = (z-z_c)/L_z)
    double coeffs[9][9];        // coeffs[j][i] = monomial-j coefficient of nodal basis i

    AxiHenrotteFE_Q2_Curved(const double rs[9], const double zs[9]);

    ELEMENT_TYPE ElementType() const override { return ET_QUAD; }

    void CalcShape(const IntegrationPoint & ip,
                   BareSliceVector<> shape) const override;
    void CalcDShape(const IntegrationPoint & ip,
                    BareSliceMatrix<> dshape) const override;
};

}  // namespace axifem

#endif  // RADIA_AXIFEM_AXI_HENROTTE_FE_HPP
