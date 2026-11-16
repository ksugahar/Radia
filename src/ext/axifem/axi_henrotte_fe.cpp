// axi_henrotte_fe.cpp — Henrotte/Meeker axisymmetric FE implementation.
// Pybind11 bindings included via <python_comp.hpp> which sets up the
// pybind11 namespace alias and pulls in the dependent NGSolve headers.

#include <comp.hpp>
#include <python_comp.hpp>
#include <pybind11/stl.h>
#include <algorithm>
#include <cmath>
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

void AxiHenrotteFE_Q1_AxisAligned::CalcVdof(const IntegrationPoint & ip, double * n,
                                            double (*dn)[2]) const
{
    // Reference (xi, eta) in [0, 1]^2 -> physical (r, z); psi bilinear in (s, z).
    const double r = r_a + (r_b - r_a) * ip(0);
    const double z = z_a + (z_b - z_a) * ip(1);
    const double s = r * r, sa = r_a * r_a, sb = r_b * r_b;
    const double is = 1.0 / (sb - sa), iz = 1.0 / (z_b - z_a);
    const double Ls[2] = { (sb - s) * is, (s - sa) * is }, dLs[2] = { -is, is };
    const double Lz[2] = { (z_b - z) * iz, (z - z_a) * iz }, dLz[2] = { -iz, iz };
    const int IS[4] = { 0, 1, 1, 0 }, IZ[4] = { 0, 0, 1, 1 };
    const double r_n[4] = { r_a, r_b, r_b, r_a };
    for (int j = 0; j < 4; ++j) {
        const double psi = Ls[IS[j]] * Lz[IZ[j]];
        const double ps = dLs[IS[j]] * Lz[IZ[j]], pz = Ls[IS[j]] * dLz[IZ[j]];
        double d[2];
        siso::VdofFromPsi(r_n[j], s, sb, psi, ps, pz, n[j], dn ? d : nullptr);
        if (dn) { dn[j][0] = d[0]; dn[j][1] = d[1]; }
    }
}

// ---------------------------------------------------------------------------
// AxiHenrotteFE_Q2_AxisAligned
// ---------------------------------------------------------------------------

namespace {

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
    is_axis(siso::OnAxisEdgeEnd(ra, ra, rb))
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

void AxiHenrotteFE_Q2_AxisAligned::CalcVdof(const IntegrationPoint & ip, double * n,
                                            double (*dn)[2]) const
{
    // Reference (xi, eta) in [0, 1]^2 -> physical (r, z) -> (s, z).
    const double r = r_a + (r_b - r_a) * ip(0);
    const double z = z_a + (z_b - z_a) * ip(1);
    const double s = r * r;
    const double sa = r_a * r_a, sb = r_b * r_b, sm = 0.5 * (sa + sb);
    const double r_n[9] = { r_a, r_b, r_b, r_a, std::sqrt(sm), r_b,
                            std::sqrt(sm), r_a, std::sqrt(sm) };
    double M[9], Ms[9], Mz[9];
    if (!is_axis) {
        Q2MonomialsGeneral(s, z, M);
        Q2MonomialsGeneralDs(s, z, Ms);
        Q2MonomialsGeneralDz(s, z, Mz);
    } else {
        Q2MonomialsAxis(s, z, M);
        Q2MonomialsAxisDs(s, z, Ms);
        Q2MonomialsAxisDz(s, z, Mz);
    }
    for (int j = 0; j < 9; ++j) { n[j] = 0.0; if (dn) dn[j][0] = dn[j][1] = 0.0; }
    for (int q = 0; q < n_nz; ++q) {
        const int j = nz_idx[q];   // local DOF index
        double psi = 0.0, ps = 0.0, pz = 0.0;
        for (int k = 0; k < n_nz; ++k) {
            psi += Vinv[k][j] * M[k];
            ps  += Vinv[k][j] * Ms[k];
            pz  += Vinv[k][j] * Mz[k];
        }
        double d[2];
        siso::VdofFromPsi(r_n[j], s, sb, psi, ps, pz, n[j], dn ? d : nullptr);
        if (dn) { dn[j][0] = d[0]; dn[j][1] = d[1]; }
    }
}
// ---------------------------------------------------------------------------
// s-isoparametric elements (axi_henrotte_siso.hpp)
// ---------------------------------------------------------------------------

AxiHenrotteFE_P1_Triangle::AxiHenrotteFE_P1_Triangle(const double rs[6], const double zs[6])
  : AxiHenrotteBaseFE(3, 2)
{
    geo.Init(rs, zs);
}

void AxiHenrotteFE_P1_Triangle::CalcVdof(const IntegrationPoint & ip, double * n,
                                         double (*dn)[2]) const
{
    geo.EvalVdof(ip(0), ip(1), n, dn);
}

AxiHenrotteFE_P2_Triangle::AxiHenrotteFE_P2_Triangle(const double rs[6], const double zs[6])
  : AxiHenrotteBaseFE(6, 4)
{
    geo.Init(rs, zs);
}

void AxiHenrotteFE_P2_Triangle::CalcVdof(const IntegrationPoint & ip, double * n,
                                         double (*dn)[2]) const
{
    geo.EvalVdof(ip(0), ip(1), n, dn);
}

AxiHenrotteFE_Q2_Curved::AxiHenrotteFE_Q2_Curved(const double rs[9], const double zs[9])
  : AxiHenrotteBaseFE(9, 4)
{
    geo.Init(rs, zs);
}

void AxiHenrotteFE_Q2_Curved::CalcVdof(const IntegrationPoint & ip, double * n,
                                       double (*dn)[2]) const
{
    geo.EvalVdof(ip(0), ip(1), n, dn);
}

// ---------------------------------------------------------------------------
// AxiHenrotteFE_Edge -- boundary trace of the V-DOF basis on one segment.
// The edge's (s, z) map uses siso::EdgeMidnode, the rule of the volume
// elements, and the edge parameter tau is matched to the segment point the way
// the volume evaluation matches it: through s on an s-midpoint edge (horizontal
// or ending on the axis: the evaluation radius is the mesh radius) and through
// z on a mesh-midpoint edge (tau = t there).  So the trace equals the volume
// value on both neighbours.  DOF order: vertex0, vertex1[, edge midnode].
// ---------------------------------------------------------------------------

void AxiHenrotteFE_Edge::CalcVdof(const IntegrationPoint & ip, double * n,
                                  double (*dn)[2]) const
{
    if (dn)
        throw Exception("AxiHenrotteFE_Edge: boundary gradient evaluation is not "
                        "supported (the boundary trace has values only)");
    // NGSolve's ET_SEGM reference has vertex 0 at x = 1 and vertex 1 at x = 0;
    // t runs from vertex 0 (t = 0) to vertex 1 (t = 1).
    const double t = 1.0 - ip(0);
    // Mesh radius along the segment (quadratic through the midpoint: exact for
    // straight and Curve(2) segments).
    const double rp = (1.0 - t) * (1.0 - 2.0 * t) * r0 + t * (2.0 * t - 1.0) * r1
                      + 4.0 * t * (1.0 - t) * rm;
    const double s0 = r0 * r0, s1 = r1 * r1;
    double sm, zmid;
    siso::MidKind kind;
    if (siso::CurvedEdge(r0, z0, r1, z1, rm, zm)
        && !siso::OnAxisEdgeEnd(r0, r0, r1) && !siso::OnAxisEdgeEnd(r1, r0, r1)) {
        // the rule of the volume elements: a curved edge off the axis takes the
        // mesh midpoint, and the volume target of its point at t is the (s, z)
        // edge point at t; one ending on the axis is s-matched (below)
        kind = siso::MidKind::Mesh;
        sm = rm * rm;
    } else {
        kind = siso::EdgeMidnode(r0, z0, r1, z1, rm, zm, sm, zmid);
    }
    // tau: s-matched on an s-midpoint edge (s is linear in tau there), else
    // the segment parameter itself (z is linear in it).
    const double tau = (kind != siso::MidKind::Mesh) ? (rp * rp - s0) / (s1 - s0) : t;
    // Edge (s, z) map and its tau derivative (quadratic Lagrange at 0, 1, 1/2).
    const double q0 = (1.0 - tau) * (1.0 - 2.0 * tau), q1 = tau * (2.0 * tau - 1.0);
    const double qm = 4.0 * tau * (1.0 - tau);
    const double s = q0 * s0 + q1 * s1 + qm * sm;
    const double ds = (4.0 * tau - 3.0) * s0 + (4.0 * tau - 1.0) * s1 + (4.0 - 8.0 * tau) * sm;
    double psi[3], dpsi[3], rn[3] = { r0, r1, std::sqrt(std::max(sm, 0.0)) };
    if (GetNDof() == 2) {
        psi[0] = 1.0 - tau; psi[1] = tau;
        dpsi[0] = -1.0;     dpsi[1] = 1.0;
    } else {
        psi[0] = q0; psi[1] = q1; psi[2] = qm;
        dpsi[0] = 4.0 * tau - 3.0; dpsi[1] = 4.0 * tau - 1.0; dpsi[2] = 4.0 - 8.0 * tau;
    }
    const double s_scale = std::max(std::max(s0, s1), 1e-300);
    for (int k = 0; k < GetNDof(); ++k) {
        // psi_s along the edge = (d psi / d tau) / (d s / d tau); only used at
        // the axis, where d s / d tau != 0 because the edge ends on the axis.
        const double ps = ds != 0.0 ? dpsi[k] / ds : 0.0;
        // evaluation radius = mesh radius, as in the volume (equal to s on
        // straight segments)
        siso::VdofFromPsi(rn[k], rp * rp, s_scale, psi[k], ps, 0.0, n[k], nullptr);
    }
}

// ---------------------------------------------------------------------------
// Python bindings (constructors for inspection)
// ---------------------------------------------------------------------------

void ExportAxiHenrotteFE(pybind11::module & m) {
    namespace py = pybind11;
    py::class_<AxiHenrotteFE_Q1_AxisAligned, ngfem::FiniteElement,
               std::shared_ptr<AxiHenrotteFE_Q1_AxisAligned>>(
        m, "AxiHenrotteFE_Q1_AxisAligned",
        "Q1 Henrotte FE on an axis-aligned rectangle [r_a, r_b] x [z_a, z_b]:\n"
        "psi bilinear in (s = r^2, z), basis N_i = r_i psi_i / r.")
        .def(py::init<double, double, double, double>(),
             py::arg("r_a"), py::arg("r_b"), py::arg("z_a"), py::arg("z_b"));

    py::class_<AxiHenrotteFE_Q2_AxisAligned, ngfem::FiniteElement,
               std::shared_ptr<AxiHenrotteFE_Q2_AxisAligned>>(
        m, "AxiHenrotteFE_Q2_AxisAligned",
        "Q2 Henrotte FE on an axis-aligned rectangle: psi biquadratic in\n"
        "(s = r^2, z) with s-midpoint nodes, basis N_i = r_i psi_i / r.\n"
        "Axis-touching elements use the 6 non-axis nodes.")
        .def(py::init<double, double, double, double>(),
             py::arg("r_a"), py::arg("r_b"), py::arg("z_a"), py::arg("z_b"))
        .def_readonly("is_axis", &AxiHenrotteFE_Q2_AxisAligned::is_axis)
        .def_readonly("n_nz",    &AxiHenrotteFE_Q2_AxisAligned::n_nz);

    py::class_<AxiHenrotteFE_P1_Triangle, ngfem::FiniteElement,
               std::shared_ptr<AxiHenrotteFE_P1_Triangle>>(
        m, "AxiHenrotteFE_P1_Triangle",
        "P1 Henrotte FE on a straight triangle: linear reference basis on the\n"
        "quadratic (s, z) map, basis N_i = r_i psi_i / r.")
        .def(py::init([](std::array<double, 3> rs, std::array<double, 3> zs) {
                double r6[6] = { rs[0], rs[1], rs[2], 0.5 * (rs[0] + rs[1]),
                                 0.5 * (rs[1] + rs[2]), 0.5 * (rs[2] + rs[0]) };
                double z6[6] = { zs[0], zs[1], zs[2], 0.5 * (zs[0] + zs[1]),
                                 0.5 * (zs[1] + zs[2]), 0.5 * (zs[2] + zs[0]) };
                return std::make_shared<AxiHenrotteFE_P1_Triangle>(r6, z6);
            }),
             py::arg("rs"), py::arg("zs"));

    py::class_<AxiHenrotteFE_P2_Triangle, ngfem::FiniteElement,
               std::shared_ptr<AxiHenrotteFE_P2_Triangle>>(
        m, "AxiHenrotteFE_P2_Triangle",
        "P2 Henrotte FE on a straight triangle, nodes [v0, v1, v2, m01, m12,\n"
        "m20]: quadratic reference basis on the (s, z) map, N_i = r_i psi_i / r.")
        .def(py::init([](std::array<double, 6> rs, std::array<double, 6> zs) {
                return std::make_shared<AxiHenrotteFE_P2_Triangle>(rs.data(), zs.data());
            }),
             py::arg("rs"), py::arg("zs"));
}

}  // namespace axifem
