// axi_henrotte_integrators.cpp — Closed-form BilinearFormIntegrators.
//
// Element-matrix formulas ported verbatim from the validated Python prototype:
//   W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/axifem/axifem_quad.py
//     _element_matrices_quad_closed_form  (Q1 stiffness)
//     element_sigma_mass_quad             (Q1 sigma mass)
//   W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/axifem/axifem_core.py
//     element_matrices                    (P1 triangle stiffness, FEMM prob3big)
//   W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/axifem/sigma_mass.py
//     element_sigma_mass                  (P1 triangle sigma mass, Hammer 7-point)
//
// DOF convention (V-DOF): u_j = A_phi at vertex j. This matches the Python
// references (axifem/disk_hiruma_quad.py and disk_hiruma.py).
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
    double sa = fe.r_a * fe.r_a, sb = fe.r_b * fe.r_b;
    double za = fe.z_a,           zb = fe.z_b;
    if (!fe.is_axis) {
        double M[81];
        q2_henrotte::KPhiGeneral(sa, sb, za, zb, mu, mu, M);
        Q2MonomialToVDof<9>(fe, M, elmat);
    } else {
        double M[36];
        q2_henrotte::KPhiAxis(sb, za, zb, mu, mu, M);
        Q2MonomialToVDof<6>(fe, M, elmat);
    }
}

void Q2SigmaMassElement(const AxiHenrotteFE_Q2_AxisAligned & fe, double sigma,
                        FlatMatrix<double> elmat)
{
    double sa = fe.r_a * fe.r_a, sb = fe.r_b * fe.r_b;
    double za = fe.z_a,           zb = fe.z_b;
    if (!fe.is_axis) {
        double M[81];
        q2_henrotte::MSigmaPhiGeneral(sa, sb, za, zb, sigma, M);
        Q2MonomialToVDof<9>(fe, M, elmat);
    } else {
        double M[36];
        q2_henrotte::MSigmaPhiAxis(sb, za, zb, sigma, M);
        Q2MonomialToVDof<6>(fe, M, elmat);
    }
}

// ---------------------------------------------------------------------------
// P1 triangle stiffness in V-DOF (matches P2TriangleStiffness / the Q2 closed
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
// Verbatim port of W:/.../axifem/sigma_mass.py:element_sigma_mass.
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

// ---------------------------------------------------------------------------
// P2 triangle stiffness in V-DOF via Duffy-Gauss quadrature on the physical
// triangle. FEMM Henrotte vector formulation:
//
//   K[i, j] = (2 pi / mu) * r_i * r_j *
//             Integrate_T [ (d psi_i/dr)(d psi_j/dr)
//                         + (d psi_i/dz)(d psi_j/dz) ] / r dA
//
// Reference: w:/.../axifem/axifem_p2_triangle.py Phase B1c (sphere vs
// Stoll: -0.51% for 50 mm finite air box, Phase B3 Kelvin removes that).
// Axis DOFs (r_i = 0) auto-decouple (zero row/column via the r_i factor).
//
// Quadrature: 8x8 = 64-point Duffy-Gauss-Legendre. Exact for polynomial
// degree 15 in each variable; integrand has degree <=10, so 8-point is
// safe. Axis-touching cells have integrable 1/r near-axis singularity;
// the quadrature handles it numerically with degraded accuracy near axis.
// ---------------------------------------------------------------------------

// Module-static Duffy-Gauss 8x8 = 64 reference-triangle points (xi, eta, w).
// Computed at first call from numpy-equivalent 1D Gauss-Legendre on [0, 1].
// Weights sum to 1/2 = area of reference triangle.
namespace {
constexpr int N_GL = 8;
constexpr double _gl_nodes_01[N_GL] = {
    0.01985507175123188,  0.10166676129318664,  0.23723379504183551,
    0.40828267875217508,  0.59171732124782492,  0.76276620495816449,
    0.89833323870681336,  0.98014492824876812,
};
constexpr double _gl_weights_01[N_GL] = {
    0.05061426814518813,  0.11119051722668724,  0.15685332293894364,
    0.18134189168918100,  0.18134189168918100,  0.15685332293894364,
    0.11119051722668724,  0.05061426814518813,
};

// Iterate over Duffy-Gauss reference-triangle points. Inner body uses (xi,
// eta, weight). Weight already includes Duffy Jacobian (1 - u).
template <typename F>
inline void DuffyGauss8x8(F && body)
{
    for (int i = 0; i < N_GL; ++i) {
        double u = _gl_nodes_01[i];
        double u_w = _gl_weights_01[i];
        double one_minus_u = 1.0 - u;
        for (int k = 0; k < N_GL; ++k) {
            double v = _gl_nodes_01[k];
            double v_w = _gl_weights_01[k];
            double xi = u;
            double eta = v * one_minus_u;
            double w = u_w * v_w * one_minus_u;
            body(xi, eta, w);
        }
    }
}
}  // anonymous namespace

// Reference-triangle P2 Lagrange shape function and derivative helpers.
// Kept identical to the FE-side helpers (axi_henrotte_fe.cpp) so that the
// integrator uses the SAME geometric map as the FE itself (essential for
// curved meshes -- mesh.Curve(p>=2) makes the edge nodes lie off the
// straight-chord midpoint, and the FE's mid-edge nodes returned by the
// FESpace already include that curvature).
// NGSolve ET_TRIG convention (verified): V0 at ref (1,0), V1 at (0,1), V2 at (0,0).
// Local order [v0, v1, v2, m01, m12, m20]. Area coords lam_0=xi, lam_1=eta, lam_2=1-xi-eta.
static inline void IntP2RefShape(double xi, double eta, double N[6])
{
    double l0 = xi;
    double l1 = eta;
    double l2 = 1.0 - xi - eta;
    N[0] = l0 * (2.0 * l0 - 1.0);
    N[1] = l1 * (2.0 * l1 - 1.0);
    N[2] = l2 * (2.0 * l2 - 1.0);
    N[3] = 4.0 * l0 * l1;
    N[4] = 4.0 * l1 * l2;
    N[5] = 4.0 * l2 * l0;
}
static inline void IntP2RefDShape(double xi, double eta,
                                   double dN_dxi[6], double dN_deta[6])
{
    dN_dxi[0]  = 4.0 * xi - 1.0;
    dN_dxi[1]  = 0.0;
    dN_dxi[2]  = -3.0 + 4.0 * xi + 4.0 * eta;
    dN_dxi[3]  = 4.0 * eta;
    dN_dxi[4]  = -4.0 * eta;
    dN_dxi[5]  = 4.0 - 8.0 * xi - 4.0 * eta;
    dN_deta[0] = 0.0;
    dN_deta[1] = 4.0 * eta - 1.0;
    dN_deta[2] = -3.0 + 4.0 * xi + 4.0 * eta;
    dN_deta[3] = 4.0 * xi;
    dN_deta[4] = 4.0 - 4.0 * xi - 8.0 * eta;
    dN_deta[5] = -4.0 * xi;
}
static inline void IntP2GeomMap(const double rn[6], const double zn[6],
                                 double xi, double eta,
                                 double & rp, double & zp,
                                 double & dr_dxi, double & dr_deta,
                                 double & dz_dxi, double & dz_deta)
{
    double N[6], dNxi[6], dNeta[6];
    IntP2RefShape(xi, eta, N);
    IntP2RefDShape(xi, eta, dNxi, dNeta);
    rp = zp = 0.0;
    dr_dxi = dr_deta = dz_dxi = dz_deta = 0.0;
    for (int k = 0; k < 6; ++k) {
        rp      += N[k]    * rn[k];
        zp      += N[k]    * zn[k];
        dr_dxi  += dNxi[k] * rn[k];
        dr_deta += dNeta[k]* rn[k];
        dz_dxi  += dNxi[k] * zn[k];
        dz_deta += dNeta[k]* zn[k];
    }
}

void P2TriangleStiffness(const double rn[6], const double zn[6], double mu,
                          FlatMatrix<double> elmat)
{
    // Build 6x6 Vandermonde inverse so psi_i(r, z) = sum_j Vinv[j, i] m_j(r, z).
    Mat<6,6> Vand;
    for (int i = 0; i < 6; ++i) {
        double si = rn[i] * rn[i];
        Vand(i, 0) = 1.0;
        Vand(i, 1) = si;
        Vand(i, 2) = zn[i];
        Vand(i, 3) = si * si;
        Vand(i, 4) = si * zn[i];
        Vand(i, 5) = zn[i] * zn[i];
    }
    Mat<6,6> Vinv;
    CalcInverse(Vand, Vinv);

    Mat<6,6> K;
    K = 0.0;
    DuffyGauss8x8([&](double xi, double eta, double w) {
        // Per-quadrature-point geometric map + Jacobian (handles curved P2
        // meshes; reduces to the affine map for straight edges).
        double rp, zp, dr_dxi, dr_deta, dz_dxi, dz_deta;
        IntP2GeomMap(rn, zn, xi, eta, rp, zp,
                     dr_dxi, dr_deta, dz_dxi, dz_deta);
        if (rp <= EPS_AXIS) return;  // 1/r integrand near axis (integrable)
        double detJ = abs(dr_dxi * dz_deta - dr_deta * dz_dxi);
        double sp = rp * rp;
        double dm_dr[6] = { 0.0, 2.0 * rp, 0.0, 4.0 * rp * sp, 2.0 * rp * zp, 0.0 };
        double dm_dz[6] = { 0.0, 0.0, 1.0, 0.0, sp, 2.0 * zp };
        double dpsi_dr[6], dpsi_dz[6];
        for (int i = 0; i < 6; ++i) {
            double sr = 0.0, sz = 0.0;
            for (int j = 0; j < 6; ++j) {
                sr += Vinv(j, i) * dm_dr[j];
                sz += Vinv(j, i) * dm_dz[j];
            }
            dpsi_dr[i] = sr;
            dpsi_dz[i] = sz;
        }
        double factor = w * detJ / rp;
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                K(i, j) += factor * (dpsi_dr[i] * dpsi_dr[j]
                                    + dpsi_dz[i] * dpsi_dz[j]);
    });
    // Apply (2 pi / mu) * r_i * r_j pre-factor.
    double inv_mu = 1.0 / mu;
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            elmat(i, j) = 2.0 * PI * inv_mu * rn[i] * rn[j] * K(i, j);
}

// ---------------------------------------------------------------------------
// P2 triangle sigma-mass in V-DOF via Duffy-Gauss 8x8 quadrature.
// FEMM Henrotte: N_A_i = r_i * psi_i / r,
//   M[i, j] = 2 pi sigma r_i r_j * Integrate_T psi_i psi_j / r dA
// Axis DOFs auto-decouple via the r_i factor.
// ---------------------------------------------------------------------------
void P2TriangleSigmaMass(const double rn[6], const double zn[6], double sigma,
                          FlatMatrix<double> elmat)
{
    Mat<6,6> Vand;
    for (int i = 0; i < 6; ++i) {
        double si = rn[i] * rn[i];
        Vand(i, 0) = 1.0;
        Vand(i, 1) = si;
        Vand(i, 2) = zn[i];
        Vand(i, 3) = si * si;
        Vand(i, 4) = si * zn[i];
        Vand(i, 5) = zn[i] * zn[i];
    }
    Mat<6,6> Vinv;
    CalcInverse(Vand, Vinv);

    Mat<6,6> M;
    M = 0.0;
    DuffyGauss8x8([&](double xi, double eta, double w) {
        // Per-quadrature-point Lagrange geometric map (curved-P2-aware).
        double rp, zp, dr_dxi, dr_deta, dz_dxi, dz_deta;
        IntP2GeomMap(rn, zn, xi, eta, rp, zp,
                     dr_dxi, dr_deta, dz_dxi, dz_deta);
        if (rp <= EPS_AXIS) return;
        double detJ = abs(dr_dxi * dz_deta - dr_deta * dz_dxi);
        double sp = rp * rp;
        double m[6] = { 1.0, sp, zp, sp * sp, sp * zp, zp * zp };
        double psi[6];
        for (int i = 0; i < 6; ++i) {
            double s = 0.0;
            for (int j = 0; j < 6; ++j) s += Vinv(j, i) * m[j];
            psi[i] = s;
        }
        double factor = w * detJ / rp;
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                M(i, j) += factor * psi[i] * psi[j];
    });
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            elmat(i, j) = 2.0 * PI * sigma * rn[i] * rn[j]
                        * 0.5 * (M(i, j) + M(j, i));
}

// ---------------------------------------------------------------------------
// AxiHenrotteFE_Q2_Curved stiffness / sigma-mass via the element's own
// biquadratic geometric map + tensor Gauss-Legendre 8x8 on [0,1]^2.
// Mirrors P2TriangleStiffness/SigmaMass (V-DOF, 2 pi r_i r_j prefactor) but for
// the 9-node curved quad.  K[i,j] = (2 pi / mu) r_i r_j INT (grad psi_i . grad psi_j)/r.
// ---------------------------------------------------------------------------

// Biquadratic Lagrange geometric map (duplicate of the FE-side helper, kept here
// so the integrator uses the SAME map as the FE -- essential for curved meshes).
constexpr int Q2C_IX[9] = { 0, 2, 2, 0, 1, 2, 1, 0, 1 };
constexpr int Q2C_IY[9] = { 0, 0, 2, 2, 0, 1, 2, 1, 1 };
inline void IntQ2CurvedGeomMap(const double rn[9], const double zn[9],
                               double xi, double eta,
                               double & rp, double & zp,
                               double & dr_dxi, double & dr_deta,
                               double & dz_dxi, double & dz_deta)
{
    const double Lx[3]  = { 2*xi*xi - 3*xi + 1, -4*xi*xi + 4*xi, 2*xi*xi - xi };
    const double dLx[3] = { 4*xi - 3,           -8*xi + 4,        4*xi - 1 };
    const double Ly[3]  = { 2*eta*eta - 3*eta + 1, -4*eta*eta + 4*eta, 2*eta*eta - eta };
    const double dLy[3] = { 4*eta - 3,              -8*eta + 4,          4*eta - 1 };
    rp = zp = dr_dxi = dr_deta = dz_dxi = dz_deta = 0.0;
    for (int k = 0; k < 9; ++k) {
        double N   = Lx[Q2C_IX[k]]  * Ly[Q2C_IY[k]];
        double dNx = dLx[Q2C_IX[k]] * Ly[Q2C_IY[k]];
        double dNy = Lx[Q2C_IX[k]]  * dLy[Q2C_IY[k]];
        rp += N * rn[k];       zp += N * zn[k];
        dr_dxi += dNx * rn[k]; dr_deta += dNy * rn[k];
        dz_dxi += dNx * zn[k]; dz_deta += dNy * zn[k];
    }
}

// Physical-coord psi gradients at a quad point (uses the FE's scaling + coeffs).
inline void Q2CurvedPsiGrad(const AxiHenrotteFE_Q2_Curved & fe, double rp, double zp,
                            double dpsi_dr[9], double dpsi_dz[9])
{
    double rpr = rp / fe.L_r, zpr = (zp - fe.z_c) / fe.L_z, spr = rpr * rpr;
    double dmr[9] = { 0.0, 2*rpr, 4*rpr*spr, 0.0, 2*rpr*zpr,
                      4*rpr*spr*zpr, 0.0, 2*rpr*zpr*zpr, 4*rpr*spr*zpr*zpr };
    double dmz[9] = { 0.0, 0.0, 0.0, 1.0, spr, spr*spr, 2*zpr, 2*spr*zpr, 2*spr*spr*zpr };
    for (int i = 0; i < 9; ++i) {
        double gr = 0.0, gz = 0.0;
        for (int j = 0; j < 9; ++j) { gr += fe.coeffs[j][i]*dmr[j]; gz += fe.coeffs[j][i]*dmz[j]; }
        dpsi_dr[i] = gr / fe.L_r; dpsi_dz[i] = gz / fe.L_z;
    }
}
inline void Q2CurvedPsi(const AxiHenrotteFE_Q2_Curved & fe, double rp, double zp,
                        double psi[9])
{
    double rpr = rp / fe.L_r, zpr = (zp - fe.z_c) / fe.L_z, spr = rpr * rpr;
    double m[9] = { 1.0, spr, spr*spr, zpr, spr*zpr, spr*spr*zpr,
                    zpr*zpr, spr*zpr*zpr, spr*spr*zpr*zpr };
    for (int i = 0; i < 9; ++i) {
        double v = 0.0;
        for (int j = 0; j < 9; ++j) v += fe.coeffs[j][i]*m[j];
        psi[i] = v;
    }
}

void Q2CurvedStiffness(const AxiHenrotteFE_Q2_Curved & fe, double mu,
                       FlatMatrix<double> elmat)
{
    Mat<9,9> K; K = 0.0;
    for (int ix = 0; ix < N_GL; ++ix)
        for (int iy = 0; iy < N_GL; ++iy) {
            double xi = _gl_nodes_01[ix], eta = _gl_nodes_01[iy];
            double w = _gl_weights_01[ix] * _gl_weights_01[iy];
            double rp, zp, drx, dre, dzx, dze;
            IntQ2CurvedGeomMap(fe.rn, fe.zn, xi, eta, rp, zp, drx, dre, dzx, dze);
            if (rp <= EPS_AXIS) continue;
            double detJ = abs(drx * dze - dre * dzx);
            double gr[9], gz[9];
            Q2CurvedPsiGrad(fe, rp, zp, gr, gz);
            double factor = w * detJ / rp;
            for (int i = 0; i < 9; ++i)
                for (int j = 0; j < 9; ++j)
                    K(i, j) += factor * (gr[i]*gr[j] + gz[i]*gz[j]);
        }
    double inv_mu = 1.0 / mu;
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            elmat(i, j) = 2.0 * PI * inv_mu * fe.rn[i] * fe.rn[j] * K(i, j);
}

void Q2CurvedSigmaMass(const AxiHenrotteFE_Q2_Curved & fe, double sigma,
                       FlatMatrix<double> elmat)
{
    Mat<9,9> M; M = 0.0;
    for (int ix = 0; ix < N_GL; ++ix)
        for (int iy = 0; iy < N_GL; ++iy) {
            double xi = _gl_nodes_01[ix], eta = _gl_nodes_01[iy];
            double w = _gl_weights_01[ix] * _gl_weights_01[iy];
            double rp, zp, drx, dre, dzx, dze;
            IntQ2CurvedGeomMap(fe.rn, fe.zn, xi, eta, rp, zp, drx, dre, dzx, dze);
            if (rp <= EPS_AXIS) continue;
            double detJ = abs(drx * dze - dre * dzx);
            double psi[9]; Q2CurvedPsi(fe, rp, zp, psi);
            double factor = w * detJ / rp;
            for (int i = 0; i < 9; ++i)
                for (int j = 0; j < 9; ++j)
                    M(i, j) += factor * psi[i] * psi[j];
        }
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j)
            elmat(i, j) = 2.0 * PI * sigma * fe.rn[i] * fe.rn[j]
                        * 0.5 * (M(i, j) + M(j, i));
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
        Q2CurvedStiffness(*q2c, mu, elmat);
        return;
    }
    if (auto * p1 = dynamic_cast<const AxiHenrotteFE_P1_Triangle*>(&fel)) {
        P1TriangleStiffness(p1->r, p1->z, mu, elmat);
        return;
    }
    if (auto * p2 = dynamic_cast<const AxiHenrotteFE_P2_Triangle*>(&fel)) {
        P2TriangleStiffness(p2->r, p2->z, mu, elmat);
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
        Q2CurvedSigmaMass(*q2c, sigma, elmat);
        return;
    }
    if (auto * p1 = dynamic_cast<const AxiHenrotteFE_P1_Triangle*>(&fel)) {
        P1TriangleSigmaMass(p1->r, p1->z, sigma, elmat);
        return;
    }
    if (auto * p2 = dynamic_cast<const AxiHenrotteFE_P2_Triangle*>(&fel)) {
        P2TriangleSigmaMass(p2->r, p2->z, sigma, elmat);
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
// term, so axis-touching elements are integrable with the FULL Q2 9-monomial
// basis (unlike magnetic which drops 3 axis-incompatible monomials).
//
// We therefore bypass fe.Vinv (which assumes the magnetic axis-reduced basis)
// and compute the full 9x9 inverse Vandermonde in-place from element corner
// + s-midpoint edge midnode + face center coordinates.
//
// All matrices in q_heat_henrotte_generated.hpp already include the leading
// pi factor (from 2 pi r dr dz = pi ds dz).

namespace {

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
        "Nodal-T DOFs (no flux-function transformation).  Axis elements\n"
        "use the full 9-monomial Q2 basis (no 1/s factor in heat).")
        .def(py::init<shared_ptr<CoefficientFunction>>(),
             py::arg("k"));

    py::class_<AxiHenrotteHeatMassBFI, BilinearFormIntegrator,
               shared_ptr<AxiHenrotteHeatMassBFI>>(
        m, "AxiHenrotteHeatMassBFI",
        "Closed-form axisymmetric HEAT capacity (transient term):\n"
        "  m(T, v) = pi * Integrate( rho_c * T * v ) ds dz.\n"
        "Pass rho_c = rho * c_p [J/(m^3 K)] as a CoefficientFunction.")
        .def(py::init<shared_ptr<CoefficientFunction>>(),
             py::arg("rho_c"));
}

}  // namespace axifem
