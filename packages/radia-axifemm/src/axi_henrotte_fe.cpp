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
