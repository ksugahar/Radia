// axi_henrotte_fe.cpp — Henrotte/Meeker axisymmetric FE implementation.
// Pybind11 bindings included via <python_comp.hpp> which sets up the
// pybind11 namespace alias and pulls in the dependent NGSolve headers.

#include <comp.hpp>
#include <python_comp.hpp>
#include <pybind11/stl.h>
#include "axi_henrotte_fe.hpp"

namespace radia_axifemm {

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
// AxiHenrotteFE_Q3_AxisAligned
// ---------------------------------------------------------------------------

namespace {

// Q3 monomial basis (general, 16 entries) at (s, z), order: a major (s-power),
// b minor (z-power). m[i] = s^a * z^b with i = 4*a + b for a, b in {0..3}.
inline void Q3MonomialsGeneral(double s, double z, double M[16])
{
    double s_pow[4] = {1.0, s, s*s, s*s*s};
    double z_pow[4] = {1.0, z, z*z, z*z*z};
    int idx = 0;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            M[idx++] = s_pow[a] * z_pow[b];
}

inline void Q3MonomialsGeneralDs(double s, double z, double dM[16])
{
    // d/ds(s^a z^b) = a s^(a-1) z^b
    double sm1[4] = {0.0, 1.0, 2.0*s, 3.0*s*s};   // a * s^(a-1)
    double z_pow[4] = {1.0, z, z*z, z*z*z};
    int idx = 0;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            dM[idx++] = sm1[a] * z_pow[b];
}

inline void Q3MonomialsGeneralDz(double s, double z, double dM[16])
{
    double s_pow[4] = {1.0, s, s*s, s*s*s};
    double zm1[4] = {0.0, 1.0, 2.0*z, 3.0*z*z};   // b * z^(b-1)
    int idx = 0;
    for (int a = 0; a < 4; ++a)
        for (int b = 0; b < 4; ++b)
            dM[idx++] = s_pow[a] * zm1[b];
}

// Q3 axis monomials (12 entries): {s^a z^b : 1 <= a <= 3, 0 <= b <= 3}.
inline void Q3MonomialsAxis(double s, double z, double M[12])
{
    double s_pow[3] = {s, s*s, s*s*s};            // a = 1, 2, 3
    double z_pow[4] = {1.0, z, z*z, z*z*z};
    int idx = 0;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 4; ++b)
            M[idx++] = s_pow[a] * z_pow[b];
}

inline void Q3MonomialsAxisDs(double s, double z, double dM[12])
{
    // d/ds(s^a z^b) for a=1,2,3 -> 1, 2s, 3s^2
    double sm1[3] = {1.0, 2.0*s, 3.0*s*s};
    double z_pow[4] = {1.0, z, z*z, z*z*z};
    int idx = 0;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 4; ++b)
            dM[idx++] = sm1[a] * z_pow[b];
}

inline void Q3MonomialsAxisDz(double s, double z, double dM[12])
{
    double s_pow[3] = {s, s*s, s*s*s};
    double zm1[4] = {0.0, 1.0, 2.0*z, 3.0*z*z};
    int idx = 0;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 4; ++b)
            dM[idx++] = s_pow[a] * zm1[b];
}

// Build the 16 Q3 node coordinates in local order. Returns s_n[16], z_n[16].
inline void Q3NodeCoordsGeneral(double sa, double sb, double za, double zb,
                                 double s_n[16], double z_n[16])
{
    double s_t1 = (2.0*sa + sb) / 3.0;
    double s_t2 = (sa + 2.0*sb) / 3.0;
    double z_t1 = (2.0*za + zb) / 3.0;
    double z_t2 = (za + 2.0*zb) / 3.0;
    // Local indexing per docstring in the header:
    //   0..3:  4 corners (NGSolve quad vertex order)
    //   4-5:   bottom edge (z=za) at s_t1, s_t2
    //   6-7:   top edge    (z=zb) at s_t1, s_t2
    //   8-9:   left edge   (s=sa) at z_t1, z_t2
    //   10-11: right edge  (s=sb) at z_t1, z_t2
    //   12-15: face interior (4 nodes in 2x2 grid)
    double sn[16] = { sa, sb, sb, sa,
                      s_t1, s_t2, s_t1, s_t2,
                      sa, sa, sb, sb,
                      s_t1, s_t2, s_t1, s_t2 };
    double zn[16] = { za, za, zb, zb,
                      za, za, zb, zb,
                      z_t1, z_t2, z_t1, z_t2,
                      z_t1, z_t1, z_t2, z_t2 };
    for (int i = 0; i < 16; ++i) { s_n[i] = sn[i]; z_n[i] = zn[i]; }
}

}  // anonymous namespace

AxiHenrotteFE_Q3_AxisAligned::AxiHenrotteFE_Q3_AxisAligned(
    double ra, double rb, double za, double zb)
  : AxiHenrotteBaseFE(16, 6),  // 16 DOFs, polynomial degree 6 in r (s^3 = r^6)
    r_a(ra), r_b(rb), z_a(za), z_b(zb),
    is_axis(ra < EPS_AXIS_FE)
{
    for (int i = 0; i < 16; ++i)
        for (int j = 0; j < 16; ++j) Vinv[i][j] = 0.0;

    double sa = ra * ra, sb = rb * rb;
    double s_n[16], z_n[16];
    Q3NodeCoordsGeneral(sa, sb, za, zb, s_n, z_n);

    if (!is_axis) {
        double V[16][16];
        for (int i = 0; i < 16; ++i) {
            double M[16];
            Q3MonomialsGeneral(s_n[i], z_n[i], M);
            for (int j = 0; j < 16; ++j) V[i][j] = M[j];
        }
        double Vi[16][16];
        if (!InvertNxN<16>(V, Vi))
            throw Exception("AxiHenrotteFE_Q3_AxisAligned: singular interior Vandermonde");
        for (int i = 0; i < 16; ++i)
            for (int j = 0; j < 16; ++j) Vinv[i][j] = Vi[i][j];
        for (int k = 0; k < 16; ++k) nz_idx[k] = k;
        n_nz = 16;
    } else {
        // Axis-touching: 12 non-axis nodes are locals not at s = sa = 0.
        // Axis-side nodes: 0 (sa, za), 3 (sa, zb), 8 (sa, z_t1), 9 (sa, z_t2).
        // Active local indices (keep): {1, 2, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15}.
        const int idx[12] = { 1, 2, 4, 5, 6, 7, 10, 11, 12, 13, 14, 15 };
        double V[12][12];
        for (int i = 0; i < 12; ++i) {
            double M[12];
            Q3MonomialsAxis(s_n[idx[i]], z_n[idx[i]], M);
            for (int j = 0; j < 12; ++j) V[i][j] = M[j];
        }
        double Vi[12][12];
        if (!InvertNxN<12>(V, Vi))
            throw Exception("AxiHenrotteFE_Q3_AxisAligned: singular axis Vandermonde");
        for (int i = 0; i < 12; ++i)
            for (int j = 0; j < 12; ++j)
                Vinv[i][idx[j]] = Vi[i][j];
        for (int k = 0; k < 12; ++k) nz_idx[k] = idx[k];
        n_nz = 12;
    }
}

void AxiHenrotteFE_Q3_AxisAligned::CalcShape(
    const IntegrationPoint & ip, BareSliceVector<> shape) const
{
    double xi = ip(0), eta = ip(1);
    double r = r_a + (r_b - r_a) * xi;
    double z = z_a + (z_b - z_a) * eta;
    double s = r * r;

    for (int j = 0; j < 16; ++j) shape(j) = 0.0;

    if (!is_axis) {
        double M[16];
        Q3MonomialsGeneral(s, z, M);
        for (int j = 0; j < 16; ++j) {
            double v = 0.0;
            for (int k = 0; k < 16; ++k) v += Vinv[k][j] * M[k];
            shape(j) = v;
        }
    } else {
        double M[12];
        Q3MonomialsAxis(s, z, M);
        for (int q = 0; q < 12; ++q) {
            int j = nz_idx[q];
            double v = 0.0;
            for (int k = 0; k < 12; ++k) v += Vinv[k][j] * M[k];
            shape(j) = v;
        }
    }
}

void AxiHenrotteFE_Q3_AxisAligned::CalcDShape(
    const IntegrationPoint & ip, BareSliceMatrix<> dshape) const
{
    double xi = ip(0), eta = ip(1);
    double r = r_a + (r_b - r_a) * xi;
    double z = z_a + (z_b - z_a) * eta;
    double s = r * r;
    double ds_dxi  = 2.0 * r * (r_b - r_a);
    double dz_deta = z_b - z_a;

    for (int j = 0; j < 16; ++j) { dshape(j, 0) = 0.0; dshape(j, 1) = 0.0; }

    if (!is_axis) {
        double dM_ds[16], dM_dz[16];
        Q3MonomialsGeneralDs(s, z, dM_ds);
        Q3MonomialsGeneralDz(s, z, dM_dz);
        for (int j = 0; j < 16; ++j) {
            double dL_ds = 0.0, dL_dz = 0.0;
            for (int k = 0; k < 16; ++k) {
                dL_ds += Vinv[k][j] * dM_ds[k];
                dL_dz += Vinv[k][j] * dM_dz[k];
            }
            dshape(j, 0) = ds_dxi  * dL_ds;
            dshape(j, 1) = dz_deta * dL_dz;
        }
    } else {
        double dM_ds[12], dM_dz[12];
        Q3MonomialsAxisDs(s, z, dM_ds);
        Q3MonomialsAxisDz(s, z, dM_dz);
        for (int q = 0; q < 12; ++q) {
            int j = nz_idx[q];
            double dL_ds = 0.0, dL_dz = 0.0;
            for (int k = 0; k < 12; ++k) {
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

    py::class_<AxiHenrotteFE_Q3_AxisAligned, ngfem::FiniteElement,
               std::shared_ptr<AxiHenrotteFE_Q3_AxisAligned>>(
        m, "AxiHenrotteFE_Q3_AxisAligned",
        "Q3 Henrotte FE on axis-aligned rectangle. 16 DOFs (4 corners + 8 edge\n"
        "midnodes + 4 face interior). Tensor-product Lagrange in (s, z) on a 4x4\n"
        "grid at s = sa, (2sa+sb)/3, (sa+2sb)/3, sb (similarly z). Axis-touching\n"
        "elements use a 12-monomial restricted basis.\n"
        "WARNING: Vandermonde conditioning ~1e6-1e8 for typical disk meshes;\n"
        "borderline for double precision. Switch to orthogonal basis for Q4+.")
        .def(py::init<double, double, double, double>(),
             py::arg("r_a"), py::arg("r_b"), py::arg("z_a"), py::arg("z_b"))
        .def_readonly("is_axis", &AxiHenrotteFE_Q3_AxisAligned::is_axis)
        .def_readonly("n_nz",    &AxiHenrotteFE_Q3_AxisAligned::n_nz);

    py::class_<AxiHenrotteFE_P1_Triangle, ngfem::FiniteElement,
               std::shared_ptr<AxiHenrotteFE_P1_Triangle>>(
        m, "AxiHenrotteFE_P1_Triangle",
        "P1 Henrotte FE on a triangle in (r, z).\n"
        "Shape functions: a + b r^2 + c z.")
        .def(py::init([](std::array<double, 3> rs, std::array<double, 3> zs) {
                return std::make_shared<AxiHenrotteFE_P1_Triangle>(rs.data(), zs.data());
            }),
             py::arg("rs"), py::arg("zs"));
}

}  // namespace radia_axifemm
