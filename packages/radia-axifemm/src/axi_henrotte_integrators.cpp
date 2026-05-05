// axi_henrotte_integrators.cpp — Closed-form BilinearFormIntegrators.
//
// Element-matrix formulas ported verbatim from the validated Python prototype:
//   W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/axifemm/axifemm_quad.py
//     _element_matrices_quad_closed_form  (Q1 stiffness)
//     element_sigma_mass_quad             (Q1 sigma mass)
//   W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/axifemm/axifemm_core.py
//     element_matrices                    (P1 triangle stiffness, FEMM prob3big)
//   W:/30_CauerLadderNetwork/2026_04_01_長方形CLN/axifemm/sigma_mass.py
//     element_sigma_mass                  (P1 triangle sigma mass, Hammer 7-point)
//
// DOF convention (V-DOF): u_j = A_phi at vertex j. This matches the Python
// references (axifemm/disk_hiruma_quad.py and disk_hiruma.py).
// - Q1 closed form: M_V = T M_phi T with T_jj = 2 pi r_j; axis-vertex rows/
//   cols become zero (caller MUST Dirichlet axis edges).
// - P1 triangle stiffness: -pi factor applied to FEMM raw matrices, matches
//   axifemm_core.py post-multiplication.
// - P1 triangle sigma-mass: V-DOF Hammer 7-point integral of N_A_i N_A_j with
//   N_A_j = r_j psi_j / r (the A_phi basis at vertex j).
//
// The user MUST apply Dirichlet on every axis edge ("axis|...") for axis-
// touching elements; otherwise the V-DOF zero rows/cols make the system
// singular.

#include <comp.hpp>
#include <python_comp.hpp>
#include <pybind11/stl.h>
#include <cmath>
#include "axi_henrotte_integrators.hpp"
#include "axi_henrotte_fe.hpp"

namespace radia_axifemm {

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

// ---------------------------------------------------------------------------
// Build inverse Vandermonde for Q1 axis-aligned rectangle.
//   V[j, k] = (basis_k)(s_j, z_j) where basis_k in {1, s, z, s*z}
//   inv_V converts (phi_0, phi_1, phi_2, phi_3) to (a, b, c, d) coefficients.
// Vertex order: (sa, za), (sb, za), (sb, zb), (sa, zb)
// ---------------------------------------------------------------------------
void Q1InverseVandermonde(double ra, double rb, double za, double zb,
                          Mat<4,4> & inv_V)
{
    double sa = ra * ra;
    double sb = rb * rb;
    Mat<4,4> V;
    V(0,0) = 1; V(0,1) = sa; V(0,2) = za; V(0,3) = sa * za;
    V(1,0) = 1; V(1,1) = sb; V(1,2) = za; V(1,3) = sb * za;
    V(2,0) = 1; V(2,1) = sb; V(2,2) = zb; V(2,3) = sb * zb;
    V(3,0) = 1; V(3,1) = sa; V(3,2) = zb; V(3,3) = sa * zb;
    CalcInverse(V, inv_V);
}

// ---------------------------------------------------------------------------
// Q1 stiffness in (a, b, c, d) coefficient coords. Returns the 4x4 M_coef
// such that the energy quadratic form W = (1/2) coef^T M_coef coef.
// Inputs: sa = ra^2, sb = rb^2 (sa < sb), za < zb, mu = absolute permeability.
// ---------------------------------------------------------------------------
Mat<4,4> Q1StiffnessCoef(double sa, double sb, double za, double zb, double mu)
{
    Mat<4,4> M;
    M = 0.0;

    // W_z (only b, d coefficients contribute):
    //   M_coef[b,b] = (sb-sa)(zb-za) / (pi mu)
    //   M_coef[b,d] = (sb-sa)(zb^2-za^2) / (2 pi mu)
    //   M_coef[d,d] = (sb-sa)(zb^3-za^3) / (3 pi mu)
    double inv_pi_mu = 1.0 / (PI * mu);
    double Hz_bb = (sb - sa) * (zb - za) * inv_pi_mu;
    double Hz_bd = (sb - sa) * (zb*zb - za*za) * 0.5 * inv_pi_mu;
    double Hz_dd = (sb - sa) * (zb*zb*zb - za*za*za) / 3.0 * inv_pi_mu;
    M(1,1) += Hz_bb;
    M(1,3) += Hz_bd;
    M(3,1) += Hz_bd;
    M(3,3) += Hz_dd;

    // W_r (only c, d coefficients contribute):
    //   if sa > 0: log(sb/sa) term
    //   if sa = 0: shortcut — set Hr_cc, Hr_cd = 0 (axis Dirichlet handles it)
    double Hr_cc, Hr_cd, Hr_dd;
    if (sa < EPS_AXIS) {
        Hr_cc = 0.0;
        Hr_cd = 0.0;
        Hr_dd = (zb - za) * sb * sb / (8.0 * PI * mu);
    } else {
        double inv_4pi_mu = 1.0 / (4.0 * PI * mu);
        Hr_cc = (zb - za) * log(sb/sa) * inv_4pi_mu;
        Hr_cd = (zb - za) * (sb - sa) * inv_4pi_mu;
        Hr_dd = (zb - za) * (sb*sb - sa*sa) * 0.5 * inv_4pi_mu;
    }
    M(2,2) += Hr_cc;
    M(2,3) += Hr_cd;
    M(3,2) += Hr_cd;
    M(3,3) += Hr_dd;

    return M;
}

// ---------------------------------------------------------------------------
// Q1 sigma-mass in (a, b, c, d) coefficient coords.
//   M_coef[k, l] = (sigma / (4 pi)) * Integrate(basis_k basis_l / s, ds dz)
// Inputs: sa, sb, za, zb, sigma.
// ---------------------------------------------------------------------------
Mat<4,4> Q1SigmaMassCoef(double sa, double sb, double za, double zb, double sigma)
{
    double Sb_Sa_1 = sb - sa;
    double Sb_Sa_2 = 0.5 * (sb*sb - sa*sa);
    double Iz0 = zb - za;
    double Iz1 = 0.5 * (zb*zb - za*za);
    double Iz2 = (zb*zb*zb - za*za*za) / 3.0;

    Mat<4,4> I;
    I = 0.0;

    if (sa < EPS_AXIS) {
        // Axis-touching: only b, d couplings survive (axis vertices Dirichlet).
        I(1,1) = Iz0 * 0.5 * sb * sb;
        I(1,3) = Iz1 * 0.5 * sb * sb;
        I(3,1) = I(1,3);
        I(3,3) = Iz2 * 0.5 * sb * sb;
    } else {
        double log_sbsa = log(sb / sa);
        // basis indices: 0=(1), 1=(s), 2=(z), 3=(s*z)
        I(0,0) = Iz0 * log_sbsa;
        I(0,1) = Iz0 * Sb_Sa_1;             I(1,0) = I(0,1);
        I(0,2) = Iz1 * log_sbsa;             I(2,0) = I(0,2);
        I(0,3) = Iz1 * Sb_Sa_1;              I(3,0) = I(0,3);
        I(1,1) = Iz0 * Sb_Sa_2;
        I(1,2) = Iz1 * Sb_Sa_1;              I(2,1) = I(1,2);
        I(1,3) = Iz1 * Sb_Sa_2;              I(3,1) = I(1,3);
        I(2,2) = Iz2 * log_sbsa;
        I(2,3) = Iz2 * Sb_Sa_1;              I(3,2) = I(2,3);
        I(3,3) = Iz2 * Sb_Sa_2;
    }

    Mat<4,4> M = (sigma / (4.0 * PI)) * I;
    return M;
}

// ---------------------------------------------------------------------------
// Compute Q1 element matrix in V-DOF: elmat = T (inv_V^T M_coef inv_V) T,
// where T_jj = 2 pi r_j. This puts the matrix in the same DOF basis as the
// FEMM-style P1 triangle below, so the assembled global system is consistent.
// Axis-touching elements get zero rows/cols on axis vertices, eliminated by
// Dirichlet BC.
// ---------------------------------------------------------------------------
template <typename CoefBuilder>
void Q1ElementMatrix(double ra, double rb, double za, double zb,
                     CoefBuilder build_coef,
                     FlatMatrix<double> elmat)
{
    Mat<4,4> inv_V;
    Q1InverseVandermonde(ra, rb, za, zb, inv_V);

    Mat<4,4> M_coef = build_coef(ra*ra, rb*rb, za, zb);

    Mat<4,4> M_phi = Trans(inv_V) * M_coef * inv_V;

    // T = diag(2 pi r_0, 2 pi r_1, 2 pi r_2, 2 pi r_3) with r = (ra, rb, rb, ra).
    double r_nodes[4] = { ra, rb, rb, ra };
    double T[4];
    for (int j = 0; j < 4; ++j) T[j] = 2.0 * PI * r_nodes[j];

    for (int i = 0; i < 4; ++i)
        for (int j = 0; j < 4; ++j) {
            double v = T[i] * M_phi(i, j) * T[j];
            elmat(i, j) = 0.5 * (v + T[j] * M_phi(j, i) * T[i]);
        }
}

// ---------------------------------------------------------------------------
// P1 triangle stiffness (FEMM prob3big.cpp via axifemm_core.element_matrices).
// Inputs: rn[3], zn[3] vertex coords, mu (isotropic).
// Returns 3x3 in psi-DOF (= phi/V_FEMM = nodal value of FE function).
// ---------------------------------------------------------------------------
void P1TriangleStiffness(const double rn[3], const double zn[3], double mu,
                          FlatMatrix<double> elmat)
{
    // Geometric coefficients.
    double p[3] = { zn[1] - zn[2], zn[2] - zn[0], zn[0] - zn[1] };
    double q[3] = { rn[2] - rn[1], rn[0] - rn[2], rn[1] - rn[0] };
    double g[3] = { 0.5*(rn[1]+rn[2]), 0.5*(rn[0]+rn[2]), 0.5*(rn[0]+rn[1]) };

    double R = (rn[0] + rn[1] + rn[2]) / 3.0;
    double a_hat = 0.0;
    for (int j = 0; j < 3; ++j) a_hat += rn[j]*rn[j] * p[j] / (4.0 * R);

    int flag = 0;
    for (int j = 0; j < 3; ++j) if (abs(rn[j]) < EPS_AXIS) ++flag;

    double R_hat;
    if (flag == 2) {
        R_hat = R;
    } else if (flag == 1) {
        double r1, r2;
        if (abs(rn[0]) < EPS_AXIS)      { r1 = rn[1]; r2 = rn[2]; }
        else if (abs(rn[1]) < EPS_AXIS) { r1 = rn[2]; r2 = rn[0]; }
        else                            { r1 = rn[0]; r2 = rn[1]; }
        if (abs(r1 - r2) < EPS_AXIS)
            R_hat = 0.5 * r2;
        else
            R_hat = (r1 - r2) / (2.0 * (log(r1) - log(r2)));
    } else {
        // interior; q[i]=0 sub-cases first
        if (abs(q[0]) < EPS_AXIS)
            R_hat = q[1]*q[1] / (2.0 * (-q[1] + rn[0] * log(rn[0]/rn[2])));
        else if (abs(q[1]) < EPS_AXIS)
            R_hat = q[2]*q[2] / (2.0 * (-q[2] + rn[1] * log(rn[1]/rn[0])));
        else if (abs(q[2]) < EPS_AXIS)
            R_hat = q[0]*q[0] / (2.0 * (-q[0] + rn[2] * log(rn[2]/rn[1])));
        else
            R_hat = -(q[0]*q[1]*q[2]) /
                    (2.0 * (q[0]*rn[0]*log(rn[0])
                          + q[1]*rn[1]*log(rn[1])
                          + q[2]*rn[2]*log(rn[2])));
    }

    double K_r = -1.0 / (2.0 * a_hat * R);
    Mat<3,3> Mr;
    Mr = 0.0;
    for (int j = 0; j < 3; ++j)
        for (int k = j; k < 3; ++k) {
            Mr(j, k) = K_r * p[j] * rn[j] * p[k] * rn[k];
            Mr(k, j) = Mr(j, k);
        }
    // axis-singular diagonal boost (prob3big.cpp 250-251)
    double diag_sum = Mr(0,0) + Mr(1,1) + Mr(2,2);
    for (int j = 0; j < 3; ++j)
        if (abs(rn[j]) < EPS_AXIS)
            Mr(j, j) += diag_sum;

    double K_z = -1.0 / (2.0 * a_hat * R_hat);
    Mat<3,3> Mz;
    Mz = 0.0;
    for (int j = 0; j < 3; ++j)
        for (int k = j; k < 3; ++k) {
            Mz(j, k) = K_z * (q[j] * rn[j]) * (q[k] * rn[k])
                            * (g[j] / R) * (g[k] / R);
            Mz(k, j) = Mz(j, k);
        }

    // FEMM convention: Me = Mr/mu + Mz/mu (isotropic mu).
    // The Python prototype uses Mr / mu_z + Mz / mu_r; we collapse to one mu.
    double inv_mu = 1.0 / mu;
    Mat<3,3> Me = inv_mu * (Mr + Mz);

    // FEMM element matrix has FEMM-internal sign; flip via -pi prefactor so
    // that the assembled K is positive-definite and matches the Q1 closed form
    // convention.
    for (int i = 0; i < 3; ++i)
        for (int j = 0; j < 3; ++j)
            elmat(i, j) = -PI * Me(i, j);
}

// ---------------------------------------------------------------------------
// P1 triangle sigma-mass in V-DOF via Hammer 7-point quadrature.
// Verbatim port of W:/.../axifemm/sigma_mass.py:element_sigma_mass.
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
        Q1ElementMatrix(q1->r_a, q1->r_b, q1->z_a, q1->z_b,
                        [mu](double sa, double sb, double za, double zb) {
                            return Q1StiffnessCoef(sa, sb, za, zb, mu);
                        },
                        elmat);
        return;
    }
    if (auto * p1 = dynamic_cast<const AxiHenrotteFE_P1_Triangle*>(&fel)) {
        P1TriangleStiffness(p1->r, p1->z, mu, elmat);
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
        Q1ElementMatrix(q1->r_a, q1->r_b, q1->z_a, q1->z_b,
                        [sigma](double sa, double sb, double za, double zb) {
                            return Q1SigmaMassCoef(sa, sb, za, zb, sigma);
                        },
                        elmat);
        return;
    }
    if (auto * p1 = dynamic_cast<const AxiHenrotteFE_P1_Triangle*>(&fel)) {
        P1TriangleSigmaMass(p1->r, p1->z, sigma, elmat);
        return;
    }
    throw Exception(string("AxiHenrotteSigmaMassBFI: unsupported FE type ")
                    + typeid(fel).name());
}

// ---------------------------------------------------------------------------
// Python bindings
// ---------------------------------------------------------------------------

void ExportAxiHenrotteIntegrators(pybind11::module & m)
{
    namespace py = pybind11;

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
}

}  // namespace radia_axifemm
