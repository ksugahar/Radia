// axi_henrotte_fe.hpp — Henrotte/Meeker axisymmetric finite elements.
//
// Every element of H1Henrotte exposes one field convention: nodal values of
// A_phi with the V-DOF basis N_i = r_i psi_i / r (psi = 2 pi r A_phi), its
// physical (r, z) gradient, and the same basis in the native stiffness and
// sigma-mass integrators.  See axi_henrotte_siso.hpp for the geometry, the
// shared edge midnode rule and the evaluation target.
//
// We follow the NGSolve i-tutorial Unit 9.1 ("MyElement") pattern: inherit
// directly from FiniteElement (not BaseScalarFiniteElement) to avoid the
// SIMD/parallel-evaluator pure-virtual obligations.

#ifndef RADIA_AXIFEM_AXI_HENROTTE_FE_HPP
#define RADIA_AXIFEM_AXI_HENROTTE_FE_HPP

#include <fem.hpp>
#include "axi_henrotte_siso.hpp"

namespace axifem {

using namespace ngfem;

class AxiHenrotteBaseFE : public FiniteElement {
public:
    static constexpr int kMaxDof = 9;
    AxiHenrotteBaseFE(int ndof, int order) : FiniteElement(ndof, order) {}

    /// V-DOF basis values n[i] = N_i and, if dn is not null, the physical
    /// gradient dn[i] = (dN_i/dr, dN_i/dz) at the NGSolve reference point ip.
    virtual void CalcVdof(const IntegrationPoint & ip, double * n,
                          double (*dn)[2]) const = 0;
};

// Q1 on an axis-aligned rectangle [r_a, r_b] x [z_a, z_b]: psi bilinear in
// (s = r^2, z).  Local vertex order (r_a,z_a), (r_b,z_a), (r_b,z_b), (r_a,z_b)
// (the FESpace refuses quads with another vertex order).
class AxiHenrotteFE_Q1_AxisAligned : public AxiHenrotteBaseFE {
public:
    double r_a, r_b, z_a, z_b;

    AxiHenrotteFE_Q1_AxisAligned(double ra, double rb, double za, double zb);

    ELEMENT_TYPE ElementType() const override { return ET_QUAD; }
    void CalcVdof(const IntegrationPoint & ip, double * n, double (*dn)[2]) const override;
};

// Q2 on an axis-aligned rectangle: psi in {1, s, s^2} x {1, z, z^2}, Lagrange
// at the s-midpoint nodes
//      0: (sa, za)    1: (sb, za)    2: (sb, zb)    3: (sa, zb)    (corners)
//      4: (sm, za)    5: (sb, zm)    6: (sm, zb)    7: (sa, zm)    (edge mid)
//      8: (sm, zm)                                                  (centre)
// with sm = (sa+sb)/2 (siso::EdgeMidnode: horizontal edges take the
// s-midpoint).  Axis-touching (r_a on the axis): the 6 monomials
// {s, s^2} x {1, z, z^2} on the 6 non-axis nodes; the axis nodes carry no basis.
class AxiHenrotteFE_Q2_AxisAligned : public AxiHenrotteBaseFE {
public:
    double r_a, r_b, z_a, z_b;
    bool   is_axis;
    // Column j = monomial coefficients of the Lagrange function of node j.
    double Vinv[9][9];
    int    nz_idx[9];          // first n_nz entries: nodes carrying a basis
    int    n_nz;               // 9 (interior) or 6 (axis)

    AxiHenrotteFE_Q2_AxisAligned(double ra, double rb, double za, double zb);

    ELEMENT_TYPE ElementType() const override { return ET_QUAD; }
    void CalcVdof(const IntegrationPoint & ip, double * n, double (*dn)[2]) const override;
};

// Boundary trace on one segment (order 1: 2 DOFs, order 2: vertex0, vertex1,
// edge midnode).  Used by `LinearForm += q * v * ds(label)` and the DtN trace.
class AxiHenrotteFE_Edge : public AxiHenrotteBaseFE {
public:
    double r0, z0, r1, z1;
    double rm, zm;     // segment midpoint from the element transformation

    AxiHenrotteFE_Edge(int order, double ar0, double az0, double ar1, double az1,
                       double arm, double azm)
        : AxiHenrotteBaseFE(order == 1 ? 2 : 3, order),
          r0(ar0), z0(az0), r1(ar1), z1(az1), rm(arm), zm(azm) {}

    ELEMENT_TYPE ElementType() const override { return ET_SEGM; }
    void CalcVdof(const IntegrationPoint & ip, double * n, double (*dn)[2]) const override;
};

// P1 triangle: linear reference basis on the quadratic (s, z) map.  Nodes
// [v0, v1, v2, m01, m12, m20] are the mesh nodes; the DOFs are the vertices.
class AxiHenrotteFE_P1_Triangle : public AxiHenrotteBaseFE {
public:
    siso::Geometry<siso::P2Tri, siso::P1Tri> geo;

    AxiHenrotteFE_P1_Triangle(const double rs[6], const double zs[6]);

    ELEMENT_TYPE ElementType() const override { return ET_TRIG; }
    void CalcVdof(const IntegrationPoint & ip, double * n, double (*dn)[2]) const override;
};

// P2 triangle, node order [v0, v1, v2, m01, m12, m20] on NGSolve's ET_TRIG
// reference.
class AxiHenrotteFE_P2_Triangle : public AxiHenrotteBaseFE {
public:
    siso::Geometry<siso::P2Tri> geo;

    AxiHenrotteFE_P2_Triangle(const double rs[6], const double zs[6]);

    ELEMENT_TYPE ElementType() const override { return ET_TRIG; }
    void CalcVdof(const IntegrationPoint & ip, double * n, double (*dn)[2]) const override;
};

// Q2 general straight quad (curvedquad=True), nodes: corners 0-3, edge mids
// 4 (bottom), 5 (right), 6 (top), 7 (left), centre 8 on [0,1]^2.
class AxiHenrotteFE_Q2_Curved : public AxiHenrotteBaseFE {
public:
    siso::Geometry<siso::Q2Quad> geo;

    AxiHenrotteFE_Q2_Curved(const double rs[9], const double zs[9]);

    ELEMENT_TYPE ElementType() const override { return ET_QUAD; }
    void CalcVdof(const IntegrationPoint & ip, double * n, double (*dn)[2]) const override;
};

}  // namespace axifem

#endif  // RADIA_AXIFEM_AXI_HENROTTE_FE_HPP
