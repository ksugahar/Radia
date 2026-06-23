// axi_henrotte_fe.cpp — Henrotte/Meeker axisymmetric FE implementation.
// Pybind11 bindings included via <python_comp.hpp> which sets up the
// pybind11 namespace alias and pulls in the dependent NGSolve headers.

#include <comp.hpp>
#include <python_comp.hpp>
#include <pybind11/stl.h>
#include "axi_henrotte_fe.hpp"

namespace axifem {

using namespace ngfem;

// ---------------------------------------------------------------------------
// AxiHenrotteFE_Q1_AxisAligned
// ---------------------------------------------------------------------------

AxiHenrotteFE_Q1_AxisAligned::AxiHenrotteFE_Q1_AxisAligned(
    double ra, double rb, double za, double zb)
  : AxiHenrotteBaseFE(4, 2),  // 4 DOFs, polynomial degree 2 (in r)
    r_a(ra), r_b(rb), z_a(za), z_b(zb) {}

void AxiHenrotteFE_Q1_AxisAligned::CalcShape(
    const IntegrationPoint & ip, BareSliceVector<> shape) const
{
    // Reference (xi, eta) in [0, 1]^2 -> physical (r, z) via linear mapping
    double xi = ip(0);
    double eta = ip(1);
    double r = r_a + (r_b - r_a) * xi;
    double z = z_a + (z_b - z_a) * eta;
    double s = r * r;

    // Solve Vandermonde-like 4x4 system: shape_i(r_j^2, z_j) = delta_ij
    // for the 4 vertices (s_a, z_a), (s_b, z_a), (s_b, z_b), (s_a, z_b).
    // shape_i = a_i + b_i*s + c_i*z + d_i*s*z
    // For axis-aligned rectangle the coefficients factor:
    //   shape_i = L_xi_i(s) * L_eta_i(z)
    //  where L_xi_i are linear interp in s between [s_a, s_b],
    //        L_eta_i are linear interp in z between [z_a, z_b].
    double s_a = r_a * r_a;
    double s_b = r_b * r_b;
    double L_s0 = (s_b - s) / (s_b - s_a);   // 1 at s=s_a
    double L_s1 = (s - s_a) / (s_b - s_a);   // 1 at s=s_b
    double L_z0 = (z_b - z) / (z_b - z_a);   // 1 at z=z_a
    double L_z1 = (z - z_a) / (z_b - z_a);   // 1 at z=z_b

    // Vertex order (matches NGSolve ET_QUAD):
    //   0: (s_a, z_a)
    //   1: (s_b, z_a)
    //   2: (s_b, z_b)
    //   3: (s_a, z_b)
    shape(0) = L_s0 * L_z0;
    shape(1) = L_s1 * L_z0;
    shape(2) = L_s1 * L_z1;
    shape(3) = L_s0 * L_z1;
}

void AxiHenrotteFE_Q1_AxisAligned::CalcDShape(
    const IntegrationPoint & ip, BareSliceMatrix<> dshape) const
{
    // Derivatives w.r.t. reference (xi, eta).
    // shape_i(r, z) with r = r_a + (r_b - r_a)*xi, z = z_a + (z_b - z_a)*eta.
    // dr/dxi  = r_b - r_a;  dz/deta = z_b - z_a;  cross terms zero.
    // d s/dxi = 2 r * (r_b - r_a).
    // dshape/dxi = ds/dxi * d(shape)/ds
    // dshape/deta = dz/deta * d(shape)/dz
    //
    // shape_i = L_s_i(s) * L_z_i(z) (factored), so:
    //   d(shape_i)/ds   = (dL_s/ds)_i * L_z_i
    //   d(shape_i)/dz   = L_s_i      * (dL_z/dz)_i

    double xi = ip(0);
    double eta = ip(1);
    double r = r_a + (r_b - r_a) * xi;
    double z = z_a + (z_b - z_a) * eta;
    double s = r * r;

    double s_a = r_a * r_a;
    double s_b = r_b * r_b;
    double inv_ds = 1.0 / (s_b - s_a);
    double inv_dz = 1.0 / (z_b - z_a);

    double L_s0 = (s_b - s) * inv_ds;
    double L_s1 = (s - s_a) * inv_ds;
    double dL_s0 = -inv_ds;
    double dL_s1 =  inv_ds;
    double L_z0 = (z_b - z) * inv_dz;
    double L_z1 = (z - z_a) * inv_dz;
    double dL_z0 = -inv_dz;
    double dL_z1 =  inv_dz;

    double ds_dxi = 2.0 * r * (r_b - r_a);
    double dz_deta = z_b - z_a;

    auto fill = [&](int idx, double dL_s, double L_s, double dL_z, double L_z) {
        dshape(idx, 0) = ds_dxi * dL_s * L_z;            // d/dxi
        dshape(idx, 1) = dz_deta * L_s * dL_z;           // d/deta
    };
    fill(0, dL_s0, L_s0, dL_z0, L_z0);
    fill(1, dL_s1, L_s1, dL_z0, L_z0);
    fill(2, dL_s1, L_s1, dL_z1, L_z1);
    fill(3, dL_s0, L_s0, dL_z1, L_z1);
}

// ---------------------------------------------------------------------------
// AxiHenrotteFE_Q2_AxisAligned
// ---------------------------------------------------------------------------

namespace {

constexpr double EPS_AXIS_FE = 1.0e-14;

// In-place Gauss-Jordan inverse of a small dense N x N matrix stored as
// double a[N][N]. Returns true on success, false on singular pivot. inv[][]
// must be allocated by caller (N x N).
template <int N>
bool InvertNxN(double a[N][N], double inv[N][N])
{
    double work[N][2*N];
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) work[i][j] = a[i][j];
        for (int j = 0; j < N; ++j) work[i][N + j] = (i == j) ? 1.0 : 0.0;
    }
    for (int col = 0; col < N; ++col) {
        // Partial pivot.
        int piv = col;
        double best = std::abs(work[col][col]);
        for (int i = col + 1; i < N; ++i)
            if (std::abs(work[i][col]) > best) { best = std::abs(work[i][col]); piv = i; }
        if (best < 1e-300) return false;
        if (piv != col)
            for (int j = 0; j < 2*N; ++j) std::swap(work[col][j], work[piv][j]);
        double diag = work[col][col];
        double inv_diag = 1.0 / diag;
        for (int j = 0; j < 2*N; ++j) work[col][j] *= inv_diag;
        for (int i = 0; i < N; ++i) {
            if (i == col) continue;
            double f = work[i][col];
            if (f == 0.0) continue;
            for (int j = 0; j < 2*N; ++j) work[i][j] -= f * work[col][j];
        }
    }
    for (int i = 0; i < N; ++i)
        for (int j = 0; j < N; ++j) inv[i][j] = work[i][N + j];
    return true;
}

// Q2 monomial basis (general, 9 entries) at (s, z).
inline void Q2MonomialsGeneral(double s, double z, double M[9])
{
    M[0] = 1.0;
    M[1] = s;
    M[2] = s * s;
    M[3] = z;
    M[4] = s * z;
    M[5] = s * s * z;
    M[6] = z * z;
    M[7] = s * z * z;
    M[8] = s * s * z * z;
}

inline void Q2MonomialsGeneralDs(double s, double z, double dM[9])
{
    dM[0] = 0.0;
    dM[1] = 1.0;
    dM[2] = 2.0 * s;
    dM[3] = 0.0;
    dM[4] = z;
    dM[5] = 2.0 * s * z;
    dM[6] = 0.0;
    dM[7] = z * z;
    dM[8] = 2.0 * s * z * z;
}

inline void Q2MonomialsGeneralDz(double s, double z, double dM[9])
{
    dM[0] = 0.0;
    dM[1] = 0.0;
    dM[2] = 0.0;
    dM[3] = 1.0;
    dM[4] = s;
    dM[5] = s * s;
    dM[6] = 2.0 * z;
    dM[7] = 2.0 * s * z;
    dM[8] = 2.0 * s * s * z;
}

// Q2 axis monomial basis (6 entries) at (s, z): {s, s^2, sz, s^2 z, sz^2, s^2 z^2}.
inline void Q2MonomialsAxis(double s, double z, double M[6])
{
    M[0] = s;
    M[1] = s * s;
    M[2] = s * z;
    M[3] = s * s * z;
    M[4] = s * z * z;
    M[5] = s * s * z * z;
}

inline void Q2MonomialsAxisDs(double s, double z, double dM[6])
{
    dM[0] = 1.0;
    dM[1] = 2.0 * s;
    dM[2] = z;
    dM[3] = 2.0 * s * z;
    dM[4] = z * z;
    dM[5] = 2.0 * s * z * z;
}

inline void Q2MonomialsAxisDz(double s, double z, double dM[6])
{
    dM[0] = 0.0;
    dM[1] = 0.0;
    dM[2] = s;
    dM[3] = s * s;
    dM[4] = 2.0 * s * z;
    dM[5] = 2.0 * s * s * z;
}

}  // anonymous namespace

AxiHenrotteFE_Q2_AxisAligned::AxiHenrotteFE_Q2_AxisAligned(
    double ra, double rb, double za, double zb)
  : AxiHenrotteBaseFE(9, 4),  // 9 DOFs, polynomial degree 4 in r (s^2 = r^4)
    r_a(ra), r_b(rb), z_a(za), z_b(zb),
    is_axis(ra < EPS_AXIS_FE)
{
    // Initialise Vinv to zero (axis nodes will stay zero in axis case).
    for (int i = 0; i < 9; ++i)
        for (int j = 0; j < 9; ++j) Vinv[i][j] = 0.0;

    double sa = ra * ra, sb = rb * rb;
    double sm = 0.5 * (sa + sb);
    double zm = 0.5 * (za + zb);

    // Node (s, z) coords matching JSON node_order_general.
    const double s_n[9] = { sa, sb, sb, sa, sm, sb, sm, sa, sm };
    const double z_n[9] = { za, za, zb, zb, za, zm, zb, zm, zm };

    if (!is_axis) {
        // Interior 9x9 Vandermonde V[i, j] = monomial_j evaluated at node i.
        double V[9][9];
        for (int i = 0; i < 9; ++i) {
            double M[9];
            Q2MonomialsGeneral(s_n[i], z_n[i], M);
            for (int j = 0; j < 9; ++j) V[i][j] = M[j];
        }
        double Vi[9][9];
        if (!InvertNxN<9>(V, Vi))
            throw Exception("AxiHenrotteFE_Q2_AxisAligned: singular interior Vandermonde");
        for (int i = 0; i < 9; ++i)
            for (int j = 0; j < 9; ++j) Vinv[i][j] = Vi[i][j];
        for (int k = 0; k < 9; ++k) nz_idx[k] = k;
        n_nz = 9;
    } else {
        // Axis-touching: 6 non-axis nodes are local indices {1, 2, 4, 5, 6, 8}.
        // (Axis nodes 0, 3, 7 have sa = 0 and contribute no DOF — A_phi(0) = 0.)
        const int idx[6] = { 1, 2, 4, 5, 6, 8 };
        double V[6][6];
        for (int i = 0; i < 6; ++i) {
            double M[6];
            Q2MonomialsAxis(s_n[idx[i]], z_n[idx[i]], M);
            for (int j = 0; j < 6; ++j) V[i][j] = M[j];
        }
        double Vi[6][6];
        if (!InvertNxN<6>(V, Vi))
            throw Exception("AxiHenrotteFE_Q2_AxisAligned: singular axis Vandermonde");
        // Embed 6x6 Vinv into 9x9 storage: row `idx[i]` <-> Vi[i, ...]; columns
        // 0..5 correspond to monomial_axis indices (different from monomial_general!).
        // Storage convention: Vinv[mono_general_index, local_node_index] would be
        // ambiguous in the axis case because the monomial index space differs.
        // Instead store Vinv[axis_mono_index, local_node_index] for axis case
        // (only first 6 rows used).
        for (int i = 0; i < 6; ++i)
            for (int j = 0; j < 6; ++j)
                Vinv[i][idx[j]] = Vi[i][j];
        for (int k = 0; k < 6; ++k) nz_idx[k] = idx[k];
        n_nz = 6;
    }
}

void AxiHenrotteFE_Q2_AxisAligned::CalcShape(
    const IntegrationPoint & ip, BareSliceVector<> shape) const
{
    // Reference (xi, eta) in [0, 1]^2 -> physical (r, z) -> (s, z).
    double xi = ip(0), eta = ip(1);
    double r = r_a + (r_b - r_a) * xi;
    double z = z_a + (z_b - z_a) * eta;
    double s = r * r;

    for (int j = 0; j < 9; ++j) shape(j) = 0.0;

    if (!is_axis) {
        double M[9];
        Q2MonomialsGeneral(s, z, M);
        // L_j(s, z) = sum_k Vinv[k, j] * M_k
        for (int j = 0; j < 9; ++j) {
            double v = 0.0;
            for (int k = 0; k < 9; ++k) v += Vinv[k][j] * M[k];
            shape(j) = v;
        }
    } else {
        double M[6];
        Q2MonomialsAxis(s, z, M);
        for (int q = 0; q < 6; ++q) {
            int j = nz_idx[q];          // global local-DOF index
            double v = 0.0;
            for (int k = 0; k < 6; ++k) v += Vinv[k][j] * M[k];
            shape(j) = v;
        }
        // Axis-side DOFs (0, 3, 7) stay at 0.
    }
}

void AxiHenrotteFE_Q2_AxisAligned::CalcDShape(
    const IntegrationPoint & ip, BareSliceMatrix<> dshape) const
{
    double xi = ip(0), eta = ip(1);
    double r = r_a + (r_b - r_a) * xi;
    double z = z_a + (z_b - z_a) * eta;
    double s = r * r;
    // ds/dxi  = 2 r * (r_b - r_a),   dz/deta = z_b - z_a
    double ds_dxi  = 2.0 * r * (r_b - r_a);
    double dz_deta = z_b - z_a;

    for (int j = 0; j < 9; ++j) { dshape(j, 0) = 0.0; dshape(j, 1) = 0.0; }

    if (!is_axis) {
        double dM_ds[9], dM_dz[9];
        Q2MonomialsGeneralDs(s, z, dM_ds);
        Q2MonomialsGeneralDz(s, z, dM_dz);
        for (int j = 0; j < 9; ++j) {
            double dL_ds = 0.0, dL_dz = 0.0;
            for (int k = 0; k < 9; ++k) {
                dL_ds += Vinv[k][j] * dM_ds[k];
                dL_dz += Vinv[k][j] * dM_dz[k];
            }
            dshape(j, 0) = ds_dxi  * dL_ds;
            dshape(j, 1) = dz_deta * dL_dz;
        }
    } else {
        double dM_ds[6], dM_dz[6];
        Q2MonomialsAxisDs(s, z, dM_ds);
        Q2MonomialsAxisDz(s, z, dM_dz);
        for (int q = 0; q < 6; ++q) {
            int j = nz_idx[q];
            double dL_ds = 0.0, dL_dz = 0.0;
            for (int k = 0; k < 6; ++k) {
                dL_ds += Vinv[k][j] * dM_ds[k];
                dL_dz += Vinv[k][j] * dM_dz[k];
            }
            dshape(j, 0) = ds_dxi  * dL_ds;
            dshape(j, 1) = dz_deta * dL_dz;
        }
    }
}

// ---------------------------------------------------------------------------
// AxiHenrotteFE_P1_Triangle
// ---------------------------------------------------------------------------

AxiHenrotteFE_P1_Triangle::AxiHenrotteFE_P1_Triangle(
    const double rs[3], const double zs[3])
  : AxiHenrotteBaseFE(3, 2)  // 3 DOFs, polynomial degree 2 (in r)
{
    for (int i = 0; i < 3; ++i) { r[i] = rs[i]; z[i] = zs[i]; }
    // Basis: shape_i(r, z) = alpha_i + beta_i*r^2 + gamma_i*z
    // such that shape_i(r_j, z_j) = delta_ij.
    // 3x3 system per i: [1, r_j^2, z_j] * (alpha, beta, gamma)^T_i = e_i
    double s[3] = { r[0]*r[0], r[1]*r[1], r[2]*r[2] };
    // Determinant
    double det = (s[1] - s[0]) * (z[2] - z[0]) - (z[1] - z[0]) * (s[2] - s[0]);
    if (det == 0.0) det = 1.0;  // degenerate — caller should avoid
    double inv_det = 1.0 / det;
    // Solve via Cramer for each i:
    //   shape_i = ... we use: (alpha, beta, gamma)_i = M^{-1} * e_i
    // M = [[1, s_0, z_0], [1, s_1, z_1], [1, s_2, z_2]]
    // Cofactors (transpose for inverse):
    auto col = [&](int i) {
        // Returns the i-th column of M^{-1}, i.e. (alpha_i, beta_i, gamma_i).
        // Compute as cofactors / det.
        // Indices: ip1, ip2 = (i+1)%3, (i+2)%3
        int ip1 = (i + 1) % 3;
        int ip2 = (i + 2) % 3;
        // Cofactor of (0, i): determinant of 2x2 minor with sign.
        double c0 = (s[ip1] * z[ip2] - s[ip2] * z[ip1]);  // alpha_i * det
        double c1 = -(z[ip2] - z[ip1]);                   // beta_i  * det
        double c2 = (s[ip2] - s[ip1]);                    // gamma_i * det
        // Sign correction for non-cyclic permutation:
        int sign = ((i % 2) == 0) ? 1 : -1;
        // Actually this is hand-derived; cleanest is to verify numerically below.
        (void)sign;
        alpha[i] = c0 * inv_det;
        beta[i]  = c1 * inv_det;
        gamma_[i] = c2 * inv_det;
    };
    col(0); col(1); col(2);
}

void AxiHenrotteFE_P1_Triangle::CalcShape(
    const IntegrationPoint & ip, BareSliceVector<> shape) const
{
    // Reference triangle: (xi, eta) with xi >= 0, eta >= 0, xi + eta <= 1.
    // Map to physical via standard P1 affine map:
    //   r = r_0 + (r_1 - r_0)*xi + (r_2 - r_0)*eta
    //   z = z_0 + (z_1 - z_0)*xi + (z_2 - z_0)*eta
    double xi = ip(0);
    double eta = ip(1);
    double rp = r[0] + (r[1] - r[0]) * xi + (r[2] - r[0]) * eta;
    double zp = z[0] + (z[1] - z[0]) * xi + (z[2] - z[0]) * eta;
    double sp = rp * rp;
    for (int i = 0; i < 3; ++i)
        shape(i) = alpha[i] + beta[i] * sp + gamma_[i] * zp;
}

void AxiHenrotteFE_P1_Triangle::CalcDShape(
    const IntegrationPoint & ip, BareSliceMatrix<> dshape) const
{
    double xi = ip(0);
    double eta = ip(1);
    double rp = r[0] + (r[1] - r[0]) * xi + (r[2] - r[0]) * eta;
    // d phi/dr = beta_i * 2 r ;  d phi/dz = gamma_i
    // dr/dxi = r_1 - r_0,    dr/deta = r_2 - r_0
    // dz/dxi = z_1 - z_0,    dz/deta = z_2 - z_0
    double dr_dxi = r[1] - r[0];
    double dr_deta = r[2] - r[0];
    double dz_dxi = z[1] - z[0];
    double dz_deta = z[2] - z[0];
    for (int i = 0; i < 3; ++i) {
        double dphi_dr = beta[i] * 2.0 * rp;
        double dphi_dz = gamma_[i];
        dshape(i, 0) = dphi_dr * dr_dxi + dphi_dz * dz_dxi;     // d/dxi
        dshape(i, 1) = dphi_dr * dr_deta + dphi_dz * dz_deta;   // d/deta
    }
}

// ---------------------------------------------------------------------------
// AxiHenrotteFE_P2_Triangle
// ---------------------------------------------------------------------------

AxiHenrotteFE_P2_Triangle::AxiHenrotteFE_P2_Triangle(
    const double rs[6], const double zs[6])
  : AxiHenrotteBaseFE(6, 4)  // 6 DOFs, polynomial degree 4 (in r, via r^4)
{
    for (int i = 0; i < 6; ++i) { r[i] = rs[i]; z[i] = zs[i]; }
    // Build 6x6 Vandermonde V[i, j] = m_j(r_i, z_i) where monomials are
    //   m = (1, s, z, s^2, s*z, z^2),  s = r^2.
    // Then Vinv = V^{-1}, and psi_i(r, z) = sum_j Vinv[j, i] * m_j(r, z).
    Mat<6,6> V;
    for (int i = 0; i < 6; ++i) {
        double s_i = r[i] * r[i];
        V(i, 0) = 1.0;
        V(i, 1) = s_i;
        V(i, 2) = z[i];
        V(i, 3) = s_i * s_i;
        V(i, 4) = s_i * z[i];
        V(i, 5) = z[i] * z[i];
    }
    Mat<6,6> Vinv_local;
    CalcInverse(V, Vinv_local);
    for (int i = 0; i < 6; ++i)
        for (int j = 0; j < 6; ++j)
            Vinv[i][j] = Vinv_local(i, j);
}

// Reference-triangle P2 Lagrange shape functions on NGSolve's ET_TRIG
// reference triangle convention (verified 2026-05-12 via
// examples/CLN/scripts/axifem/test_ngsolve_ref_tri_vertices.py):
//
//   Local DOF 0 (= mesh vertex V0) at ref (1, 0)
//   Local DOF 1 (= mesh vertex V1) at ref (0, 1)
//   Local DOF 2 (= mesh vertex V2) at ref (0, 0)
//   Local DOF 3 (= m01, midpoint V0-V1) at ref (0.5, 0.5)
//   Local DOF 4 (= m12, midpoint V1-V2) at ref (0, 0.5)
//   Local DOF 5 (= m20, midpoint V2-V0) at ref (0.5, 0)
//
// Area coords (NGSolve convention): lam_0 = xi, lam_1 = eta, lam_2 = 1-xi-eta.
//   N_0 = lam_0 (2 lam_0 - 1) = xi (2 xi - 1)            -- at V0
//   N_1 = lam_1 (2 lam_1 - 1) = eta (2 eta - 1)          -- at V1
//   N_2 = lam_2 (2 lam_2 - 1) = (1-xi-eta)(1-2xi-2eta)   -- at V2
//   N_3 = 4 lam_0 lam_1 = 4 xi eta                       -- at m01
//   N_4 = 4 lam_1 lam_2 = 4 eta (1-xi-eta)               -- at m12
//   N_5 = 4 lam_2 lam_0 = 4 (1-xi-eta) xi                -- at m20
//
// For a mesh that has been `.Curve(p>=2)`-d, the 6 nodes (r[k], z[k])
// returned by the FESpace are the curved-geometry positions of the
// corresponding NGSolve-reference nodes. The same shape functions then
// act as the geometric map (r(xi, eta), z(xi, eta)) = sum_k N_k (r[k], z[k]).
// For straight-edge meshes this reduces to the affine 3-vertex map.
static inline void P2RefShape(double xi, double eta, double N[6])
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

static inline void P2RefDShape(double xi, double eta,
                                double dN_dxi[6], double dN_deta[6])
{
    // dN/dxi
    dN_dxi[0]  = 4.0 * xi - 1.0;                       // d/dxi[xi(2xi-1)]
    dN_dxi[1]  = 0.0;                                   // d/dxi[eta(2eta-1)]
    dN_dxi[2]  = -3.0 + 4.0 * xi + 4.0 * eta;          // d/dxi[(1-xi-eta)(1-2xi-2eta)]
    dN_dxi[3]  = 4.0 * eta;                             // d/dxi[4 xi eta]
    dN_dxi[4]  = -4.0 * eta;                            // d/dxi[4 eta (1-xi-eta)]
    dN_dxi[5]  = 4.0 - 8.0 * xi - 4.0 * eta;            // d/dxi[4 (1-xi-eta) xi]
    // dN/deta
    dN_deta[0] = 0.0;
    dN_deta[1] = 4.0 * eta - 1.0;
    dN_deta[2] = -3.0 + 4.0 * xi + 4.0 * eta;
    dN_deta[3] = 4.0 * xi;
    dN_deta[4] = 4.0 - 4.0 * xi - 8.0 * eta;
    dN_deta[5] = -4.0 * xi;
}

// Map (xi, eta) -> (r, z) and its Jacobian using the 6-node Lagrange map.
static inline void P2GeomMap(const double r[6], const double z[6],
                              double xi, double eta,
                              double & rp, double & zp,
                              double & dr_dxi, double & dr_deta,
                              double & dz_dxi, double & dz_deta)
{
    double N[6], dNxi[6], dNeta[6];
    P2RefShape(xi, eta, N);
    P2RefDShape(xi, eta, dNxi, dNeta);
    rp = zp = 0.0;
    dr_dxi = dr_deta = dz_dxi = dz_deta = 0.0;
    for (int k = 0; k < 6; ++k) {
        rp      += N[k]    * r[k];
        zp      += N[k]    * z[k];
        dr_dxi  += dNxi[k] * r[k];
        dr_deta += dNeta[k]* r[k];
        dz_dxi  += dNxi[k] * z[k];
        dz_deta += dNeta[k]* z[k];
    }
}

void AxiHenrotteFE_P2_Triangle::CalcShape(
    const IntegrationPoint & ip, BareSliceVector<> shape) const
{
    // 6-node quadratic Lagrange geometric map: (xi, eta) -> (r, z).
    // For straight-edge meshes this reduces to the affine 3-vertex map.
    double xi = ip(0);
    double eta = ip(1);
    double rp, zp, dr_dxi, dr_deta, dz_dxi, dz_deta;
    P2GeomMap(r, z, xi, eta, rp, zp, dr_dxi, dr_deta, dz_dxi, dz_deta);
    double sp = rp * rp;
    // psi_i(rp, zp) = sum_j Vinv[j, i] * m_j with m = {1, s, z, s^2, sz, z^2}.
    double m[6] = { 1.0, sp, zp, sp * sp, sp * zp, zp * zp };
    for (int i = 0; i < 6; ++i) {
        double v = 0.0;
        for (int j = 0; j < 6; ++j) v += Vinv[j][i] * m[j];
        shape(i) = v;
    }
}

void AxiHenrotteFE_P2_Triangle::CalcDShape(
    const IntegrationPoint & ip, BareSliceMatrix<> dshape) const
{
    double xi = ip(0);
    double eta = ip(1);
    double rp, zp, dr_dxi, dr_deta, dz_dxi, dz_deta;
    P2GeomMap(r, z, xi, eta, rp, zp, dr_dxi, dr_deta, dz_dxi, dz_deta);
    double sp = rp * rp;
    // d m_j / d(r, z) at the physical point.
    double dm_dr[6] = { 0.0, 2.0 * rp, 0.0, 4.0 * rp * sp, 2.0 * rp * zp, 0.0 };
    double dm_dz[6] = { 0.0, 0.0, 1.0, 0.0, sp, 2.0 * zp };
    for (int i = 0; i < 6; ++i) {
        double dpsi_dr = 0.0, dpsi_dz = 0.0;
        for (int j = 0; j < 6; ++j) {
            dpsi_dr += Vinv[j][i] * dm_dr[j];
            dpsi_dz += Vinv[j][i] * dm_dz[j];
        }
        // Chain rule: d psi / d(xi, eta) = (d(r, z)/d(xi, eta))^T * d psi/d(r, z).
        dshape(i, 0) = dpsi_dr * dr_dxi + dpsi_dz * dz_dxi;
        dshape(i, 1) = dpsi_dr * dr_deta + dpsi_dz * dz_deta;
    }
}

// ---------------------------------------------------------------------------
// AxiHenrotteFE_Q2_Curved -- 9-node curved biquadratic Henrotte quad element.
// ---------------------------------------------------------------------------

namespace {
// FEMM/our Q2 node order -> 1D-Lagrange index (0=L@0, 1=L@0.5, 2=L@1) in xi, eta.
constexpr int Q2C_IX[9] = { 0, 2, 2, 0, 1, 2, 1, 0, 1 };
constexpr int Q2C_IY[9] = { 0, 0, 2, 2, 0, 1, 2, 1, 1 };

// Biquadratic Lagrange geometric map on the NGSolve quad ref [0,1]^2 through the
// 9 physical nodes; returns the physical point and the 2x2 geometric Jacobian.
inline void Q2CurvedGeomMap(const double rn[9], const double zn[9],
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
        rp += N * rn[k];      zp += N * zn[k];
        dr_dxi += dNx * rn[k]; dr_deta += dNy * rn[k];
        dz_dxi += dNx * zn[k]; dz_deta += dNy * zn[k];
    }
}

// Even-r-power monomials (scaled) m_j and their physical-coord derivatives.
inline void Q2CurvedMonomials(double sp, double zp, double m[9])
{
    m[0]=1.0; m[1]=sp; m[2]=sp*sp; m[3]=zp; m[4]=sp*zp;
    m[5]=sp*sp*zp; m[6]=zp*zp; m[7]=sp*zp*zp; m[8]=sp*sp*zp*zp;
}
}  // namespace

AxiHenrotteFE_Q2_Curved::AxiHenrotteFE_Q2_Curved(const double rs[9], const double zs[9])
  : AxiHenrotteBaseFE(9, 4)  // 9 DOFs, polynomial degree 4 (in r, via r^4)
{
    double rmax = 0.0, zmin = zs[0], zmax = zs[0];
    for (int i = 0; i < 9; ++i) {
        rn[i] = rs[i]; zn[i] = zs[i];
        double ar = rs[i] < 0 ? -rs[i] : rs[i];
        if (ar > rmax) rmax = ar;
        if (zs[i] < zmin) zmin = zs[i];
        if (zs[i] > zmax) zmax = zs[i];
    }
    L_r = rmax > 1e-30 ? rmax : 1e-30;
    z_c = 0.5 * (zmin + zmax);
    L_z = (0.5 * (zmax - zmin)) > 1e-30 ? 0.5 * (zmax - zmin) : 1e-30;

    Mat<9,9> V;
    for (int i = 0; i < 9; ++i) {
        double rpr = rn[i] / L_r, zpr = (zn[i] - z_c) / L_z, spr = rpr * rpr;
        double m[9]; Q2CurvedMonomials(spr, zpr, m);
        for (int j = 0; j < 9; ++j) V(i, j) = m[j];
    }
    Mat<9,9> Vi; CalcInverse(V, Vi);
    // psi_i = sum_j coeffs[j][i] m_j ;  V * C = I  =>  C = V^{-1}, coeffs[j][i]=Vi(j,i)
    for (int j = 0; j < 9; ++j)
        for (int i = 0; i < 9; ++i)
            coeffs[j][i] = Vi(j, i);
}

void AxiHenrotteFE_Q2_Curved::CalcShape(
    const IntegrationPoint & ip, BareSliceVector<> shape) const
{
    double rp, zp, a, b, c, d;
    Q2CurvedGeomMap(rn, zn, ip(0), ip(1), rp, zp, a, b, c, d);
    double rpr = rp / L_r, zpr = (zp - z_c) / L_z, spr = rpr * rpr;
    double m[9]; Q2CurvedMonomials(spr, zpr, m);
    for (int i = 0; i < 9; ++i) {
        double v = 0.0;
        for (int j = 0; j < 9; ++j) v += coeffs[j][i] * m[j];
        shape(i) = v;
    }
}

void AxiHenrotteFE_Q2_Curved::CalcDShape(
    const IntegrationPoint & ip, BareSliceMatrix<> dshape) const
{
    double rp, zp, dr_dxi, dr_deta, dz_dxi, dz_deta;
    Q2CurvedGeomMap(rn, zn, ip(0), ip(1), rp, zp,
                    dr_dxi, dr_deta, dz_dxi, dz_deta);
    double rpr = rp / L_r, zpr = (zp - z_c) / L_z, spr = rpr * rpr;
    // d m / d r', d m / d z'  (scaled), then chain to physical (/L_r, /L_z).
    double dmr[9] = { 0.0, 2*rpr, 4*rpr*spr, 0.0, 2*rpr*zpr,
                      4*rpr*spr*zpr, 0.0, 2*rpr*zpr*zpr, 4*rpr*spr*zpr*zpr };
    double dmz[9] = { 0.0, 0.0, 0.0, 1.0, spr, spr*spr, 2*zpr, 2*spr*zpr, 2*spr*spr*zpr };
    for (int i = 0; i < 9; ++i) {
        double dpsi_dr = 0.0, dpsi_dz = 0.0;
        for (int j = 0; j < 9; ++j) {
            dpsi_dr += coeffs[j][i] * dmr[j];
            dpsi_dz += coeffs[j][i] * dmz[j];
        }
        dpsi_dr /= L_r; dpsi_dz /= L_z;             // physical gradient
        // Chain rule to reference coords (NGSolve applies J^{-1} for the physical).
        dshape(i, 0) = dpsi_dr * dr_dxi + dpsi_dz * dz_dxi;   // d/dxi
        dshape(i, 1) = dpsi_dr * dr_deta + dpsi_dz * dz_deta; // d/deta
    }
}

// ---------------------------------------------------------------------------
// AxiHenrotteFE_Edge_Q1 / _Q2 -- 1D restriction of the parent quad's
// Lagrange basis to a single edge.  Used for Neumann RHS integration
// `LinearForm += q * v * ds(label)` in axisymmetric heat / magnetic
// solvers.  Basis values depend only on the edge's two endpoint
// coordinates (and the s-midpoint for Q2); no parent-quad info is
// required because the off-edge corners' Lagrange factors evaluate
// to 1 along the edge.
//
// CalcDShape is intentionally not implemented -- Neumann RHS does not
// need grad(v).  Throws if called, to surface unintended use.
// ---------------------------------------------------------------------------

namespace {

// Evaluate Lagrange basis at the 2 endpoints of a 1D segment in u-space,
// where u = s = r^2 (horizontal edge) or u = z (vertical edge).
//   L_0(u) = (u_1 - u) / (u_1 - u_0)
//   L_1(u) = (u - u_0) / (u_1 - u_0)
inline void LagrangeQ1Edge(double u0, double u1, double u, double L[2])
{
    double denom = u1 - u0;
    if (std::abs(denom) < 1e-300) { L[0] = L[1] = 0.0; return; }
    L[0] = (u1 - u) / denom;
    L[1] = (u  - u0) / denom;
}

// Q2 1D Lagrange: 3 nodes at u_0, u_1, u_m = (u_0 + u_1) / 2 (s-midpoint
// or z-midpoint convention).  L_0(u_0) = 1, L_1(u_1) = 1, L_m(u_m) = 1.
inline void LagrangeQ2Edge(double u0, double u1, double u, double L[3])
{
    double um = 0.5 * (u0 + u1);
    double d01 = u0 - u1, d0m = u0 - um;
    double d10 = u1 - u0, d1m = u1 - um;
    double dm0 = um - u0, dm1 = um - u1;
    if (std::abs(d01) < 1e-300 || std::abs(d0m) < 1e-300 ||
        std::abs(d1m) < 1e-300) {
        L[0] = L[1] = L[2] = 0.0;
        return;
    }
    L[0] = (u - u1) * (u - um) / (d01 * d0m);  // L at u0
    L[1] = (u - u0) * (u - um) / (d10 * d1m);  // L at u1
    L[2] = (u - u0) * (u - u1) / (dm0 * dm1);  // L at um
}

}  // namespace

void AxiHenrotteFE_Edge_Q1::CalcShape(
    const IntegrationPoint & ip, BareSliceVector<> shape) const
{
    // Segment reference parameter t in [0, 1].  Map to physical (r, z)
    // along the edge.
    double t = ip(0);
    double rp = (1.0 - t) * r0 + t * r1;
    double zp = (1.0 - t) * z0 + t * z1;
    bool horizontal = std::abs(z0 - z1) < 1e-12;
    double L[2];
    if (horizontal) {
        // Edge varies in r; basis varies in s = r^2.
        double s0 = r0 * r0;
        double s1 = r1 * r1;
        LagrangeQ1Edge(s0, s1, rp * rp, L);
    } else {
        // Edge varies in z (or oblique; we use z for axis-aligned grids).
        LagrangeQ1Edge(z0, z1, zp, L);
    }
    shape(0) = L[0];
    shape(1) = L[1];
}

void AxiHenrotteFE_Edge_Q1::CalcDShape(
    const IntegrationPoint & ip, BareSliceMatrix<> dshape) const
{
    throw Exception("AxiHenrotteFE_Edge_Q1::CalcDShape not implemented "
                    "(boundary gradient evaluation not supported; "
                    "Neumann RHS only needs the value trace)");
}

void AxiHenrotteFE_Edge_Q2::CalcShape(
    const IntegrationPoint & ip, BareSliceVector<> shape) const
{
    double t = ip(0);
    double rp = (1.0 - t) * r0 + t * r1;
    double zp = (1.0 - t) * z0 + t * z1;
    bool horizontal = std::abs(z0 - z1) < 1e-12;
    double L[3];
    if (horizontal) {
        double s0 = r0 * r0;
        double s1 = r1 * r1;
        LagrangeQ2Edge(s0, s1, rp * rp, L);
    } else {
        LagrangeQ2Edge(z0, z1, zp, L);
    }
    // GetDofNrs orders: vertex0, vertex1, edge midnode (NV + edge index).
    shape(0) = L[0];     // at vertex0
    shape(1) = L[1];     // at vertex1
    shape(2) = L[2];     // at edge midnode
}

void AxiHenrotteFE_Edge_Q2::CalcDShape(
    const IntegrationPoint & ip, BareSliceMatrix<> dshape) const
{
    throw Exception("AxiHenrotteFE_Edge_Q2::CalcDShape not implemented "
                    "(boundary gradient evaluation not supported; "
                    "Neumann RHS only needs the value trace)");
}

// ---------------------------------------------------------------------------
// Python bindings (Phase 2-A: just expose ctors for inspection)
// ---------------------------------------------------------------------------

void ExportAxiHenrotteFE(pybind11::module & m) {
    namespace py = pybind11;
    py::class_<AxiHenrotteFE_Q1_AxisAligned, ngfem::FiniteElement,
               std::shared_ptr<AxiHenrotteFE_Q1_AxisAligned>>(
        m, "AxiHenrotteFE_Q1_AxisAligned",
        "Q1 Henrotte FE on axis-aligned rectangle [r_a, r_b] x [z_a, z_b].\n"
        "Shape functions: a + b r^2 + c z + d r^2 z.")
        .def(py::init<double, double, double, double>(),
             py::arg("r_a"), py::arg("r_b"), py::arg("z_a"), py::arg("z_b"));

    py::class_<AxiHenrotteFE_Q2_AxisAligned, ngfem::FiniteElement,
               std::shared_ptr<AxiHenrotteFE_Q2_AxisAligned>>(
        m, "AxiHenrotteFE_Q2_AxisAligned",
        "Q2 Henrotte FE on axis-aligned rectangle. 9 DOFs (4 vertices + 4 edge\n"
        "midnodes + 1 face center) at s-midpoints (NOT physical r-midpoints).\n"
        "Axis-touching elements (r_a < EPS_AXIS) auto-switch to a 6-monomial\n"
        "restricted basis with zero shape functions on the axis-side nodes.")
        .def(py::init<double, double, double, double>(),
             py::arg("r_a"), py::arg("r_b"), py::arg("z_a"), py::arg("z_b"))
        .def_readonly("is_axis", &AxiHenrotteFE_Q2_AxisAligned::is_axis)
        .def_readonly("n_nz",    &AxiHenrotteFE_Q2_AxisAligned::n_nz);

    py::class_<AxiHenrotteFE_P1_Triangle, ngfem::FiniteElement,
               std::shared_ptr<AxiHenrotteFE_P1_Triangle>>(
        m, "AxiHenrotteFE_P1_Triangle",
        "P1 Henrotte FE on a triangle in (r, z).\n"
        "Shape functions: a + b r^2 + c z.")
        .def(py::init([](std::array<double, 3> rs, std::array<double, 3> zs) {
                return std::make_shared<AxiHenrotteFE_P1_Triangle>(rs.data(), zs.data());
            }),
             py::arg("rs"), py::arg("zs"));

    py::class_<AxiHenrotteFE_P2_Triangle, ngfem::FiniteElement,
               std::shared_ptr<AxiHenrotteFE_P2_Triangle>>(
        m, "AxiHenrotteFE_P2_Triangle",
        "P2 Henrotte FE on a triangle in (r, z) — 6 DOFs (3 vertex + 3 edge mid).\n"
        "Shape functions: {1, r^2, z, r^4, r^2 z, z^2} Lagrange interpolant.\n"
        "Node order: [v0, v1, v2, m01, m12, m20]. FEMM Henrotte vector\n"
        "convention: N_A_i = (r_i / r) psi_i. Axis vertices auto-decouple.\n"
        "Phase B1c reference: w:/.../axifem/axifem_p2_triangle.py")
        .def(py::init([](std::array<double, 6> rs, std::array<double, 6> zs) {
                return std::make_shared<AxiHenrotteFE_P2_Triangle>(rs.data(), zs.data());
            }),
             py::arg("rs"), py::arg("zs"));
}

}  // namespace axifem
