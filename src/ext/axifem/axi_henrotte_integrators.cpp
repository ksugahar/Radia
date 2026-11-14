// axi_henrotte_integrators.cpp — Closed-form BilinearFormIntegrators.
//
// Element-matrix formulas ported from the validated Python fixtures under
// tests/axifem/_reference_python/ (research copies are retained under
// validation_test/axifem/research/prototypes/):
//   axifem_quad.py
//     _element_matrices_quad_closed_form  (Q1 stiffness)
//     element_sigma_mass_quad             (Q1 sigma mass)
//   axifem_core.py
//     element_matrices                    (P1 triangle stiffness, FEMM prob3big)
//   sigma_mass.py
//     element_sigma_mass                  (P1 triangle sigma mass, Hammer 7-point)
//
// DOF convention (V-DOF): u_j = A_phi at vertex j. This matches the Python
// references (validation_test/axifem/research/verification/disk_axisym_reference.py).
// - Q1 closed form: M_V = T M_phi T with T_jj = 2 pi r_j; axis-vertex rows/
//   cols become zero (caller MUST Dirichlet axis edges).
// - P1 triangle stiffness: -pi factor applied to FEMM raw matrices, matches
//   axifem_core.py post-multiplication.
// - P1 triangle sigma-mass: V-DOF Hammer 7-point integral of N_A_i N_A_j with
//   N_A_j = r_j psi_j / r (the A_phi basis at vertex j).
//
// The user MUST apply Dirichlet on every axis edge ("axis|...") for axis-
// touching elements; otherwise the V-DOF zero rows/cols make the system
// singular.

#include <comp.hpp>
#include <python_comp.hpp>
#include <pybind11/numpy.h>
#include <pybind11/stl.h>
#include <cmath>
#include "axi_henrotte_integrators.hpp"
#include "axi_henrotte_fe.hpp"
#include "axi_henrotte_numeric.hpp"
#include "q2_henrotte_generated.hpp"
#include "q_heat_henrotte_generated.hpp"

namespace axifem {

using namespace ngfem;
using std::log;
using std::abs;

namespace {

constexpr double EPS_AXIS = 1.0e-14;
constexpr double PI = 3.14159265358979323846;

// ---------------------------------------------------------------------------
// Sample a CoefficientFunction at the element centroid (constant per element).
// ---------------------------------------------------------------------------
double SampleAtCentroid(const CoefficientFunction & cf,
                        const ElementTransformation & eltrans,
                        ELEMENT_TYPE et, LocalHeap & lh)
{
    IntegrationPoint ip_center;
    if (et == ET_QUAD)      ip_center = IntegrationPoint(0.5, 0.5, 0.0);
    else /* ET_TRIG */      ip_center = IntegrationPoint(1.0/3.0, 1.0/3.0, 0.0);

    HeapReset hr(lh);
    auto & mip = eltrans(ip_center, lh);
    return cf.Evaluate(mip);
}

void CopyMatrix4(const numeric::Matrix4 & source, FlatMatrix<double> target)
{
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            target(i, j) = source[4 * i + j];
}

void CopyMatrix9(const numeric::Matrix9 & source, FlatMatrix<double> target)
{
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            target(i, j) = source[9 * i + j];
}

void Q1InverseVandermonde(double ra, double rb, double za, double zb,
                          Mat<4,4> & inverse)
{
    const auto values =
        numeric::ComputeQ1InverseVandermonde(ra, rb, za, zb);
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            inverse(i, j) = values[4 * i + j];
}

// ---------------------------------------------------------------------------
// Q2 element matrix (interior or axis-touching, 9 DOFs in V-DOF basis).
//
//   M_phi (monomial)  =  KPhiGeneral / MSigmaPhiGeneral / KPhiAxis / ...
//                        (delivered by q2_henrotte_generated.hpp)
//   M_node (Lagrange) =  Vinv^T * M_phi * Vinv
//   M_V (psi-DOF)     =  T * M_node * T,    T_jj = 2 pi r_node_j
//
// For axis elements the axis nodes (local indices 0, 3, 7) get zero rows and
// columns — caller MUST Dirichlet those global DOFs.
// ---------------------------------------------------------------------------

// Returns r-coord at each of the 9 local nodes (sm convention).
inline void Q2NodeRCoords(double ra, double rb, double r_node[9])
{
    double sa = ra * ra, sb = rb * rb;
    double rm = std::sqrt(0.5 * (sa + sb));   // r at s-midpoint = sqrt((ra^2+rb^2)/2)
    r_node[0] = ra;
    r_node[1] = rb;
    r_node[2] = rb;
    r_node[3] = ra;
    r_node[4] = rm;
    r_node[5] = rb;
    r_node[6] = rm;
    r_node[7] = ra;
    r_node[8] = rm;
}

// Generic V-DOF transform driver.
//   For interior:  9x9 monomial matrix M_phi[81] -> 9x9 V-DOF matrix elmat.
//   For axis    :  6x6 monomial matrix M_phi[36] embedded into 9x9 elmat.
//                  Uses fe.Vinv[k=0..5][j_local] for the 6-monomial axis basis.
template <int N_PHI>
void Q2MonomialToVDof(const AxiHenrotteFE_Q2_AxisAligned & fe,
                      const double M_phi_flat[],
                      FlatMatrix<double> elmat)
{
    static_assert(N_PHI == 9 || N_PHI == 6, "Q2MonomialToVDof: N_PHI must be 6 or 9");

    elmat = 0.0;

    double r_node[9];
    Q2NodeRCoords(fe.r_a, fe.r_b, r_node);
    double T[9];
    for (int j = 0; j < 9; ++j) T[j] = 2.0 * PI * r_node[j];

    // M_phi flat -> M_phi[N_PHI][N_PHI]
    double Mphi[N_PHI][N_PHI];
    for (int i = 0; i < N_PHI; ++i)
        for (int j = 0; j < N_PHI; ++j)
            Mphi[i][j] = M_phi_flat[i * N_PHI + j];

    // Step 1: tmp[k, j_local] = sum_l M_phi[k, l] * Vinv[l, j_local]
    // (Vinv stored as Vinv[mono_index, local_node_index]; for axis case only
    // the first 6 rows of Vinv are non-zero, indices into nz_idx.)
    int n_active = fe.n_nz;
    double tmp[N_PHI][9];
    for (int k = 0; k < N_PHI; ++k) {
        for (int q = 0; q < n_active; ++q) {
            int j = fe.nz_idx[q];
            double s = 0.0;
            for (int l = 0; l < N_PHI; ++l) s += Mphi[k][l] * fe.Vinv[l][j];
            tmp[k][j] = s;
        }
    }

    // Step 2: M_node[i_local, j_local] = sum_k Vinv[k, i_local] * tmp[k, j_local]
    // Then M_V[i, j] = T[i] * M_node[i, j] * T[j].
    for (int p = 0; p < n_active; ++p) {
        int i = fe.nz_idx[p];
        for (int q = 0; q < n_active; ++q) {
            int j = fe.nz_idx[q];
            double s = 0.0;
            for (int k = 0; k < N_PHI; ++k) s += fe.Vinv[k][i] * tmp[k][j];
            elmat(i, j) = T[i] * s * T[j];
        }
    }

    // Symmetrise (numerical noise).
    for (int i = 0; i < 9; ++i)
        for (int j = i + 1; j < 9; ++j) {
            double avg = 0.5 * (elmat(i, j) + elmat(j, i));
            elmat(i, j) = avg;
            elmat(j, i) = avg;
        }
}

void Q2StiffnessElement(const AxiHenrotteFE_Q2_AxisAligned & fe, double mu,
                        FlatMatrix<double> elmat)
{
    CopyMatrix9(numeric::ComputeQ2MagneticElementMatrices(
        fe.r_a, fe.r_b, fe.z_a, fe.z_b, mu, 0.0).stiffness, elmat);
}

void Q2SigmaMassElement(const AxiHenrotteFE_Q2_AxisAligned & fe, double sigma,
                        FlatMatrix<double> elmat)
{
    CopyMatrix9(numeric::ComputeQ2MagneticElementMatrices(
        fe.r_a, fe.r_b, fe.z_a, fe.z_b, 1.0, sigma).sigma_mass, elmat);
}

// ---------------------------------------------------------------------------
// P1 triangle stiffness in V-DOF (matches the P2 triangle / the Q2 closed
// form):  K_ij = (2 pi / mu) r_i r_j INT_T grad(psi_i).grad(psi_j)/r dA,
// psi_i = alpha_i + beta_i r^2 + gamma_i z (the {1, r^2, z} Henrotte basis).
//
// This is the CORRECT axisymmetric V-DOF stiffness: a uniform axial B_z has
// nodal A_phi V_i = B0 r_i/2 with  sum_j r_j^2 psi_j == r^2, so V lies in the
// kernel (interior) and order-1 magnetostatics converges (magnetized sphere
// O(h): -1.8% -> -0.6%, FEMM-P1-like).  It REPLACES the earlier FEMM prob3big
// A=psi port, which is NOT a V-DOF operator (it does not annihilate the uniform
// field -- ||(K V_unif)_interior||/scale ~ 0.6 -- so order-1 magnetostatic
// solves gave a wrong, non-uniform interior field).  The matching V-DOF
// sigma-mass (P1TriangleSigmaMass) is unchanged, so the eddy path stays
// consistent.  Axis-touching cells use the same r_q<=EPS_AXIS skip as P2.
// ---------------------------------------------------------------------------
void P1TriangleStiffness(const double rn[3], const double zn[3], double mu,
                          FlatMatrix<double> elmat)
{
    Mat<3,3> M_vand;
    for (int j = 0; j < 3; ++j) {
        M_vand(j, 0) = 1.0;
        M_vand(j, 1) = rn[j] * rn[j];
        M_vand(j, 2) = zn[j];
    }
    Mat<3,3> inv_vand;
    CalcInverse(M_vand, inv_vand);
    double beta[3]  = { inv_vand(1,0), inv_vand(1,1), inv_vand(1,2) };
    double gamma[3] = { inv_vand(2,0), inv_vand(2,1), inv_vand(2,2) };

    double drxi  = rn[1] - rn[0];
    double dreta = rn[2] - rn[0];
    double dzxi  = zn[1] - zn[0];
    double dzeta = zn[2] - zn[0];
    double detJ  = abs(drxi * dzeta - dreta * dzxi);

    // Hammer 7-point (same rule as P1TriangleSigmaMass).
    const double s15 = 3.872983346207417;
    const double w_a = 9.0 / 80.0;
    const double w_b = (155.0 - s15) / 2400.0;
    const double w_c = (155.0 + s15) / 2400.0;
    const double a1  = (6.0 + s15) / 21.0;
    const double a2  = (9.0 - 2.0 * s15) / 21.0;
    const double b1  = (6.0 - s15) / 21.0;
    const double b2  = (9.0 + 2.0 * s15) / 21.0;
    const double pts[7][3] = {
        {1.0/3.0, 1.0/3.0, w_a},
        {a1, a1, w_b}, {a2, a1, w_b}, {a1, a2, w_b},
        {b1, b1, w_c}, {b2, b1, w_c}, {b1, b2, w_c},
    };

    Mat<3,3> Ke;
    Ke = 0.0;
    for (int q = 0; q < 7; ++q) {
        double xi = pts[q][0], eta = pts[q][1], w = pts[q][2];
        double r_q = rn[0] + drxi * xi + dreta * eta;
        if (r_q <= EPS_AXIS) continue;          // integrable 1/r near axis
        double dpsi_dr[3], dpsi_dz[3];
        for (int j = 0; j < 3; ++j) {
            dpsi_dr[j] = 2.0 * beta[j] * r_q;   // d/dr (alpha + beta r^2 + gamma z)
            dpsi_dz[j] = gamma[j];
        }
        double factor = w * detJ / r_q;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                Ke(i, j) += factor * (dpsi_dr[i]*dpsi_dr[j] + dpsi_dz[i]*dpsi_dz[j]);
    }
    double inv_mu = 1.0 / mu;
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            elmat(i, j) = 2.0 * PI * inv_mu * rn[i] * rn[j] * Ke(i, j);
}

// ---------------------------------------------------------------------------
// P1 triangle sigma-mass in V-DOF via Hammer 7-point quadrature.
// Verbatim port of <private-path>/.../axifem/sigma_mass.py:element_sigma_mass.
//
// Basis for A_phi at vertex j: N_A_j(r, z) = r_j * psi_j(r, z) / r
// where psi_j = alpha_j + beta_j r^2 + gamma_j z is the {1, r^2, z} basis
// satisfying psi_j(r_k, z_k) = delta_jk.
//
// M_e[i, j] = sigma * Integrate( N_A_i N_A_j * 2 pi r dr dz )
//           = 2 pi sigma r_i r_j * Integrate( psi_i psi_j / r dr dz )
//
// Returned matrix is in V-DOF (DOF = A_phi at vertex), matching the FEMM-
// style stiffness convention used by P1TriangleStiffness.
// ---------------------------------------------------------------------------
void P1TriangleSigmaMass(const double rn[3], const double zn[3], double sigma,
                         FlatMatrix<double> elmat)
{
    // Vandermonde for {1, r^2, z}.
    Mat<3,3> M_vand;
    for (int j = 0; j < 3; ++j) {
        M_vand(j, 0) = 1.0;
        M_vand(j, 1) = rn[j] * rn[j];
        M_vand(j, 2) = zn[j];
    }
    Mat<3,3> inv_vand;
    CalcInverse(M_vand, inv_vand);
    // Column j of inv_vand = (alpha_j, beta_j, gamma_j) for psi_j.
    double alpha[3] = { inv_vand(0,0), inv_vand(0,1), inv_vand(0,2) };
    double beta[3]  = { inv_vand(1,0), inv_vand(1,1), inv_vand(1,2) };
    double gamma[3] = { inv_vand(2,0), inv_vand(2,1), inv_vand(2,2) };

    double drxi  = rn[1] - rn[0];
    double dreta = rn[2] - rn[0];
    double dzxi  = zn[1] - zn[0];
    double dzeta = zn[2] - zn[0];
    double detJ  = abs(drxi * dzeta - dreta * dzxi);

    const double s15 = 3.872983346207417;
    const double w_a = 9.0 / 80.0;
    const double w_b = (155.0 - s15) / 2400.0;
    const double w_c = (155.0 + s15) / 2400.0;
    const double a1  = (6.0 + s15) / 21.0;
    const double a2  = (9.0 - 2.0 * s15) / 21.0;
    const double b1  = (6.0 - s15) / 21.0;
    const double b2  = (9.0 + 2.0 * s15) / 21.0;
    const double pts[7][3] = {
        {1.0/3.0, 1.0/3.0, w_a},
        {a1, a1, w_b}, {a2, a1, w_b}, {a1, a2, w_b},
        {b1, b1, w_c}, {b2, b1, w_c}, {b1, b2, w_c},
    };

    Mat<3,3> Me;
    Me = 0.0;
    for (int q = 0; q < 7; ++q) {
        double xi = pts[q][0], eta = pts[q][1], w = pts[q][2];
        double r_q = rn[0] + drxi * xi + dreta * eta;
        double z_q = zn[0] + dzxi * xi + dzeta * eta;
        if (r_q <= 0.0) continue;
        double psi[3];
        for (int j = 0; j < 3; ++j)
            psi[j] = alpha[j] + beta[j] * r_q * r_q + gamma[j] * z_q;
        // V-DOF basis: N_A_j = r_j * psi_j / r
        double NA[3];
        for (int j = 0; j < 3; ++j) NA[j] = rn[j] * psi[j] / r_q;
        // Integrand: sigma * N_A_i * N_A_j * 2 pi r * w * detJ
        double factor = sigma * 2.0 * PI * r_q * w * detJ;
        for (int i = 0; i < 3; ++i)
            for (int j = 0; j < 3; ++j)
                Me(i, j) += factor * NA[i] * NA[j];
    }

    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            elmat(i, j) = 0.5 * (Me(i, j) + Me(j, i));
}

}  // anonymous namespace

// ---------------------------------------------------------------------------
// AxiHenrotteStiffnessBFI
// ---------------------------------------------------------------------------

AxiHenrotteStiffnessBFI::AxiHenrotteStiffnessBFI(
    shared_ptr<CoefficientFunction> amu_cf)
  : mu_cf(amu_cf)
{ SetName("AxiHenrotteStiffnessBFI"); }

void AxiHenrotteStiffnessBFI::CalcElementMatrix(
    const FiniteElement & fel,
    const ElementTransformation & eltrans,
    FlatMatrix<double> elmat,
    LocalHeap & lh) const
{
    elmat = 0.0;
    ELEMENT_TYPE et = fel.ElementType();
    double mu = SampleAtCentroid(*mu_cf, eltrans, et, lh);

    if (auto * q1 = dynamic_cast<const AxiHenrotteFE_Q1_AxisAligned*>(&fel)) {
        CopyMatrix4(numeric::ComputeQ1MagneticStiffness(
                        q1->r_a, q1->r_b, q1->z_a, q1->z_b, mu),
                    elmat);
        return;
    }
    if (auto * q2 = dynamic_cast<const AxiHenrotteFE_Q2_AxisAligned*>(&fel)) {
        Q2StiffnessElement(*q2, mu, elmat);
        return;
    }
    if (auto * q2c = dynamic_cast<const AxiHenrotteFE_Q2_Curved*>(&fel)) {
        q2c->geo.Stiffness(mu, elmat);
        return;
    }
    if (auto * p1 = dynamic_cast<const AxiHenrotteFE_P1_Triangle*>(&fel)) {
        P1TriangleStiffness(p1->r, p1->z, mu, elmat);
        return;
    }
    if (auto * p2 = dynamic_cast<const AxiHenrotteFE_P2_Triangle*>(&fel)) {
        p2->geo.Stiffness(mu, elmat);
        return;
    }
    throw Exception(string("AxiHenrotteStiffnessBFI: unsupported FE type ")
                    + typeid(fel).name());
}

// ---------------------------------------------------------------------------
// AxiHenrotteSigmaMassBFI
// ---------------------------------------------------------------------------

AxiHenrotteSigmaMassBFI::AxiHenrotteSigmaMassBFI(
    shared_ptr<CoefficientFunction> asigma_cf)
  : sigma_cf(asigma_cf)
{ SetName("AxiHenrotteSigmaMassBFI"); }

void AxiHenrotteSigmaMassBFI::CalcElementMatrix(
    const FiniteElement & fel,
    const ElementTransformation & eltrans,
    FlatMatrix<double> elmat,
    LocalHeap & lh) const
{
    elmat = 0.0;
    ELEMENT_TYPE et = fel.ElementType();
    double sigma = SampleAtCentroid(*sigma_cf, eltrans, et, lh);
    if (sigma == 0.0) return;  // air / non-conductor: zero contribution

    if (auto * q1 = dynamic_cast<const AxiHenrotteFE_Q1_AxisAligned*>(&fel)) {
        CopyMatrix4(numeric::ComputeQ1SigmaMass(
                        q1->r_a, q1->r_b, q1->z_a, q1->z_b, sigma),
                    elmat);
        return;
    }
    if (auto * q2 = dynamic_cast<const AxiHenrotteFE_Q2_AxisAligned*>(&fel)) {
        Q2SigmaMassElement(*q2, sigma, elmat);
        return;
    }
    if (auto * q2c = dynamic_cast<const AxiHenrotteFE_Q2_Curved*>(&fel)) {
        q2c->geo.SigmaMass(sigma, elmat);
        return;
    }
    if (auto * p1 = dynamic_cast<const AxiHenrotteFE_P1_Triangle*>(&fel)) {
        P1TriangleSigmaMass(p1->r, p1->z, sigma, elmat);
        return;
    }
    if (auto * p2 = dynamic_cast<const AxiHenrotteFE_P2_Triangle*>(&fel)) {
        p2->geo.SigmaMass(sigma, elmat);
        return;
    }
    throw Exception(string("AxiHenrotteSigmaMassBFI: unsupported FE type ")
                    + typeid(fel).name());
}

// ===========================================================================
// HEAT EQUATION BFI implementations (no T = diag(2 pi r) wrap; nodal-T DOFs)
// ===========================================================================
//
// For axisymmetric heat the integrand is polynomial in (s, z) with NO 1/s
// term, so a dedicated scalar space could integrate axis-touching elements
// with the FULL Q2 9-monomial basis.  H1Henrotte is an electromagnetic space,
// however: its axis Q2 FE exposes the magnetic 6-function basis and zeros the
// three axis DOFs.  Mixing that FE with a full 9x9 heat matrix is invalid, so
// the BFIs below fail fast on axis-touching Q2 elements.  Production heat uses
// standard ngsolve.H1(order=2) with the 2*pi*r weak-form weight.
//
// For the still-supported off-axis Q2 research path, compute the full 9x9
// inverse Vandermonde in-place from element corner + s-midpoint edge midnode +
// face center coordinates.
//
// All matrices in q_heat_henrotte_generated.hpp already include the leading
// pi factor (from 2 pi r dr dz = pi ds dz).

namespace {

[[noreturn]] void ThrowAxisQ2HeatBasisMismatch(const char * integrator_name)
{
    throw Exception(string(integrator_name)
                    + ": axis-touching Q2 AxiHenrotteFE uses the magnetic "
                      "6-function basis (A_phi = 0 on the axis), but the heat "
                      "operator requires a full scalar-temperature basis. "
                      "Use ngsolve.H1(mesh, order=2) with the 2*pi*r "
                      "axisymmetric weight for heat; reserve H1Henrotte for "
                      "the electromagnetic A_phi solve.");
}

inline void Q2HeatInverseVandermonde(double ra, double rb, double za, double zb,
                                      Mat<9,9> & inv_V)
{
    // Q2 monomial basis: {1, s, s^2, z, sz, s^2 z, z^2, sz^2, s^2 z^2}
    // 9 nodes (s-midpoint convention; edge midnodes at s_m = (s_a + s_b)/2,
    // i.e. r = sqrt((ra^2 + rb^2)/2)).
    double sa = ra * ra, sb = rb * rb;
    double sm = 0.5 * (sa + sb);
    double zm = 0.5 * (za + zb);
    // Local node order: (s, z) coords
    //   0: (sa, za)   1: (sb, za)   2: (sb, zb)   3: (sa, zb)   (vertices)
    //   4: (sm, za)   5: (sb, zm)   6: (sm, zb)   7: (sa, zm)   (edge mids)
    //   8: (sm, zm)                                              (face center)
    double s_n[9] = { sa, sb, sb, sa, sm, sb, sm, sa, sm };
    double z_n[9] = { za, za, zb, zb, za, zm, zb, zm, zm };
    Mat<9,9> V;
    for (int j = 0; j < 9; ++j) {
        double s = s_n[j], z = z_n[j];
        V(j, 0) = 1.0;
        V(j, 1) = s;
        V(j, 2) = s * s;
        V(j, 3) = z;
        V(j, 4) = s * z;
        V(j, 5) = s * s * z;
        V(j, 6) = z * z;
        V(j, 7) = s * z * z;
        V(j, 8) = s * s * z * z;
    }
    CalcInverse(V, inv_V);
}

// Q1 heat element matrix: elmat = inv_V^T * K_mono * inv_V * coef.
// No T = diag(2 pi r) wrap (heat DOFs are nodal T directly).
template <int N, typename MonoBuilder>
void HeatElementMatrix_Q1(double ra, double rb, double za, double zb,
                          MonoBuilder build_mono, double coef,
                          FlatMatrix<double> elmat)
{
    Mat<4,4> inv_V;
    Q1InverseVandermonde(ra, rb, za, zb, inv_V);
    double M_flat[16];
    build_mono(ra*ra, rb*rb, za, zb, M_flat);
    Mat<4,4> M_mono;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            M_mono(i, j) = M_flat[i * 4 + j];
    Mat<4,4> M_node = Trans(inv_V) * M_mono * inv_V;
    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j)
            elmat(i, j) = coef * M_node(i, j);
    // Symmetrise (numerical noise).
    for (int i = 0; i < 4; ++i)
        for (int j = i + 1; j < 4; ++j) {
            double avg = 0.5 * (elmat(i, j) + elmat(j, i));
            elmat(i, j) = avg;
            elmat(j, i) = avg;
        }
}

// Q2 heat element matrix: elmat = inv_V^T * K_mono * inv_V * coef.
template <typename MonoBuilder>
void HeatElementMatrix_Q2(double ra, double rb, double za, double zb,
                          MonoBuilder build_mono, double coef,
                          FlatMatrix<double> elmat)
{
    Mat<9,9> inv_V;
    Q2HeatInverseVandermonde(ra, rb, za, zb, inv_V);
    double M_flat[81];
    build_mono(ra*ra, rb*rb, za, zb, M_flat);
    Mat<9,9> M_mono;
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            M_mono(i, j) = M_flat[i * 9 + j];
    Mat<9,9> M_node = Trans(inv_V) * M_mono * inv_V;
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            elmat(i, j) = coef * M_node(i, j);
    for (int i = 0; i < 9; ++i)
        for (int j = i + 1; j < 9; ++j) {
            double avg = 0.5 * (elmat(i, j) + elmat(j, i));
            elmat(i, j) = avg;
            elmat(j, i) = avg;
        }
}

}  // namespace

// ---------------------------------------------------------------------------
// AxiHenrotteHeatStiffnessBFI
// ---------------------------------------------------------------------------

AxiHenrotteHeatStiffnessBFI::AxiHenrotteHeatStiffnessBFI(
    shared_ptr<CoefficientFunction> ak_cf)
  : k_cf(ak_cf)
{ SetName("AxiHenrotteHeatStiffnessBFI"); }

void AxiHenrotteHeatStiffnessBFI::CalcElementMatrix(
    const FiniteElement & fel,
    const ElementTransformation & eltrans,
    FlatMatrix<double> elmat,
    LocalHeap & lh) const
{
    elmat = 0.0;
    ELEMENT_TYPE et = fel.ElementType();
    double k_val = SampleAtCentroid(*k_cf, eltrans, et, lh);

    if (auto * q1 = dynamic_cast<const AxiHenrotteFE_Q1_AxisAligned*>(&fel)) {
        HeatElementMatrix_Q1<4>(q1->r_a, q1->r_b, q1->z_a, q1->z_b,
                                 q_heat::KMonomialQ1, k_val, elmat);
        return;
    }
    if (auto * q2 = dynamic_cast<const AxiHenrotteFE_Q2_AxisAligned*>(&fel)) {
        if (q2->is_axis) {
            ThrowAxisQ2HeatBasisMismatch("AxiHenrotteHeatStiffnessBFI");
        }
        HeatElementMatrix_Q2(q2->r_a, q2->r_b, q2->z_a, q2->z_b,
                             q_heat::KMonomialQ2, k_val, elmat);
        return;
    }
    throw Exception(string("AxiHenrotteHeatStiffnessBFI: unsupported FE type ")
                    + typeid(fel).name()
                    + " (P1 triangle heat support not yet implemented; "
                      "use a structured quad mesh)");
}

// ---------------------------------------------------------------------------
// AxiHenrotteHeatMassBFI
// ---------------------------------------------------------------------------

AxiHenrotteHeatMassBFI::AxiHenrotteHeatMassBFI(
    shared_ptr<CoefficientFunction> arho_c_cf)
  : rho_c_cf(arho_c_cf)
{ SetName("AxiHenrotteHeatMassBFI"); }

void AxiHenrotteHeatMassBFI::CalcElementMatrix(
    const FiniteElement & fel,
    const ElementTransformation & eltrans,
    FlatMatrix<double> elmat,
    LocalHeap & lh) const
{
    elmat = 0.0;
    ELEMENT_TYPE et = fel.ElementType();
    double rho_c = SampleAtCentroid(*rho_c_cf, eltrans, et, lh);

    // Validate the FE contract before the zero-coefficient fast path.  An
    // axis-touching Q2 magnetic basis is incompatible with this heat operator
    // regardless of the material value on a particular element.
    if (auto * q2 = dynamic_cast<const AxiHenrotteFE_Q2_AxisAligned*>(&fel)) {
        if (q2->is_axis) {
            ThrowAxisQ2HeatBasisMismatch("AxiHenrotteHeatMassBFI");
        }
    }

    if (rho_c == 0.0) return;

    if (auto * q1 = dynamic_cast<const AxiHenrotteFE_Q1_AxisAligned*>(&fel)) {
        HeatElementMatrix_Q1<4>(q1->r_a, q1->r_b, q1->z_a, q1->z_b,
                                 q_heat::MMonomialQ1, rho_c, elmat);
        return;
    }
    if (auto * q2 = dynamic_cast<const AxiHenrotteFE_Q2_AxisAligned*>(&fel)) {
        HeatElementMatrix_Q2(q2->r_a, q2->r_b, q2->z_a, q2->z_b,
                             q_heat::MMonomialQ2, rho_c, elmat);
        return;
    }
    throw Exception(string("AxiHenrotteHeatMassBFI: unsupported FE type ")
                    + typeid(fel).name()
                    + " (P1 triangle heat support not yet implemented; "
                      "use a structured quad mesh)");
}

// ---------------------------------------------------------------------------
// Python bindings
// ---------------------------------------------------------------------------

void ExportAxiHenrotteIntegrators(pybind11::module & m)
{
    namespace py = pybind11;

    m.def(
        "q1_magnetic_element_matrices",
        [](double ra, double rb, double za, double zb,
           double mu, double sigma) {
            const auto matrices = numeric::ComputeQ1MagneticElementMatrices(
                ra, rb, za, zb, mu, sigma);
            auto as_numpy = [](const numeric::Matrix4 & values) {
                py::array_t<double> result({
                    static_cast<py::ssize_t>(4), static_cast<py::ssize_t>(4)});
                auto view = result.mutable_unchecked<2>();
                for (py::ssize_t i = 0; i < 4; ++i)
                    for (py::ssize_t j = 0; j < 4; ++j)
                        view(i, j) = values[4 * i + j];
                return result;
            };
            py::dict result;
            result["stiffness"] = as_numpy(matrices.stiffness);
            result["sigma_mass"] = as_numpy(matrices.sigma_mass);
            result["backend"] = "native-pybind";
            result["dof_convention"] = "nodal A_phi (V-DOF)";
            result["node_order"] =
                "(ra,za),(rb,za),(rb,zb),(ra,zb)";
            return result;
        },
        py::arg("ra"), py::arg("rb"), py::arg("za"), py::arg("zb"),
        py::arg("mu"), py::arg("sigma"),
        "Return shared-native Q1 Henrotte stiffness and sigma-mass matrices.");

    m.def(
        "q2_magnetic_element_matrices",
        [](double ra, double rb, double za, double zb,
           double mu, double sigma) {
            const auto matrices = numeric::ComputeQ2MagneticElementMatrices(
                ra, rb, za, zb, mu, sigma);
            auto as_numpy = [](const numeric::Matrix9 & values) {
                py::array_t<double> result({
                    static_cast<py::ssize_t>(9), static_cast<py::ssize_t>(9)});
                auto view = result.mutable_unchecked<2>();
                for (py::ssize_t i = 0; i < 9; ++i)
                    for (py::ssize_t j = 0; j < 9; ++j)
                        view(i, j) = values[9 * i + j];
                return result;
            };
            py::dict result;
            result["stiffness"] = as_numpy(matrices.stiffness);
            result["sigma_mass"] = as_numpy(matrices.sigma_mass);
            result["backend"] = "native-pybind";
            result["dof_convention"] = "nodal A_phi (V-DOF)";
            result["node_order"] =
                "(ra,za),(rb,za),(rb,zb),(ra,zb),(rm,za),(rb,zm),(rm,zb),(ra,zm),(rm,zm)";
            result["axis_touching"] = matrices.axis_touching;
            return result;
        },
        py::arg("ra"), py::arg("rb"), py::arg("za"), py::arg("zb"),
        py::arg("mu"), py::arg("sigma"),
        "Return shared-native Q2 Henrotte stiffness and sigma-mass matrices.");

    py::class_<AxiHenrotteStiffnessBFI, BilinearFormIntegrator,
               shared_ptr<AxiHenrotteStiffnessBFI>>(
        m, "AxiHenrotteStiffnessBFI",
        "Closed-form axisymmetric stiffness integrator (Henrotte basis).\n"
        "Bypasses NGSolve quadrature; uses Mathematica-derived element matrices\n"
        "exact for the 1/r-weighted axisymmetric weak form.")
        .def(py::init<shared_ptr<CoefficientFunction>>(),
             py::arg("mu"));

    py::class_<AxiHenrotteSigmaMassBFI, BilinearFormIntegrator,
               shared_ptr<AxiHenrotteSigmaMassBFI>>(
        m, "AxiHenrotteSigmaMassBFI",
        "Closed-form axisymmetric sigma-mass integrator (Henrotte basis).\n"
        "Q1 quad uses fully closed-form integration; P1 triangle uses Hammer\n"
        "7-point quadrature (exact for polynomial degree 5).")
        .def(py::init<shared_ptr<CoefficientFunction>>(),
             py::arg("sigma"));

    py::class_<AxiHenrotteHeatStiffnessBFI, BilinearFormIntegrator,
               shared_ptr<AxiHenrotteHeatStiffnessBFI>>(
        m, "AxiHenrotteHeatStiffnessBFI",
        "Closed-form axisymmetric HEAT stiffness on Q1/Q2 quad with the\n"
        "Henrotte {1, r^2, z, ...} basis.  Weak form (after s = r^2):\n"
        "  a(T, v) = pi * Integrate( k * [ 4 s d_s T d_s v + d_z T d_z v ]\n"
        "                            ds dz )\n"
        "Nodal-T DOFs (no flux-function transformation). Q2 is supported only\n"
        "off axis; axis-touching Q2 fails fast. Production heat should use\n"
        "ngsolve.H1(order=2) with the 2*pi*r weak-form weight.")
        .def(py::init<shared_ptr<CoefficientFunction>>(),
             py::arg("k"));

    py::class_<AxiHenrotteHeatMassBFI, BilinearFormIntegrator,
               shared_ptr<AxiHenrotteHeatMassBFI>>(
        m, "AxiHenrotteHeatMassBFI",
        "Closed-form axisymmetric HEAT capacity (transient term):\n"
        "  m(T, v) = pi * Integrate( rho_c * T * v ) ds dz.\n"
        "Pass rho_c = rho * c_p [J/(m^3 K)] as a CoefficientFunction.\n"
        "Q2 is supported only off axis; axis-touching Q2 fails fast.")
        .def(py::init<shared_ptr<CoefficientFunction>>(),
             py::arg("rho_c"));
}

}  // namespace axifem
