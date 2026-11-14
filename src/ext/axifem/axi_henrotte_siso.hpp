// axi_henrotte_siso.hpp -- s-isoparametric geometry for the general Henrotte
// elements (P2 triangle, Q2 curved quad).
//
// The flux function psi = 2 pi r A_phi must behave like s = r^2 at the axis so
// that the 1/r stiffness stays finite.  The earlier general elements achieved
// this with physical monomials in (s, z) fitted at the nodes.  That basis is not
// conforming: on an edge that is oblique in (r, z) its trace is a quartic in the
// edge parameter, which the three shared edge DOFs cannot fix.
//
// Here the element geometry is a Lagrange map from the reference element to
// (s, z) through the nodes (r_k^2, z_k), and the basis is the reference Lagrange
// basis N_k.  A trace on an edge depends only on the edge's own nodes and the
// edge parameter, so neighbouring elements agree on every shared edge, and the
// regularity psi ~ s at the axis is kept.  A straight (r, z) edge whose midnode
// is the r-midpoint is reproduced exactly: s is quadratic and z linear in the
// edge parameter.
//
// The quadrature of the stiffness and sigma mass runs over the reference
// element through this (s, z) map.  Point evaluation from NGSolve (DiffOps)
// receives NGSolve reference coordinates: the point is mapped to (r, z) with the
// element's (r, z) Lagrange map through the same nodes, and the (s, z) map is
// inverted by Newton's method.  Both maps use the same nodes, so straight and
// Curve(2) elements are represented consistently; for higher curve orders the
// two maps differ by the geometric approximation error.

#ifndef RADIA_AXIFEM_AXI_HENROTTE_SISO_HPP
#define RADIA_AXIFEM_AXI_HENROTTE_SISO_HPP

#include <fem.hpp>
#include <cmath>
#include <string>

namespace axifem {
namespace siso {

using namespace ngfem;

// P2 triangle on NGSolve's ET_TRIG reference: V0 at (1,0), V1 at (0,1), V2 at
// (0,0); local order [v0, v1, v2, m01, m12, m20].
struct P2Tri {
    static constexpr int N = 6;
    static constexpr bool kTriangle = true;
    static constexpr int NE = 3;
    static constexpr int EDGE[3][3] = { {0, 1, 3}, {1, 2, 4}, {2, 0, 5} };  // a, b, mid
    static void Eval(double xi, double eta, double n[6], double dx[6], double dy[6])
    {
        double l0 = xi, l1 = eta, l2 = 1.0 - xi - eta;
        n[0] = l0 * (2.0 * l0 - 1.0);
        n[1] = l1 * (2.0 * l1 - 1.0);
        n[2] = l2 * (2.0 * l2 - 1.0);
        n[3] = 4.0 * l0 * l1;
        n[4] = 4.0 * l1 * l2;
        n[5] = 4.0 * l2 * l0;
        dx[0] = 4.0 * xi - 1.0;            dy[0] = 0.0;
        dx[1] = 0.0;                       dy[1] = 4.0 * eta - 1.0;
        dx[2] = -3.0 + 4.0 * xi + 4.0 * eta; dy[2] = -3.0 + 4.0 * xi + 4.0 * eta;
        dx[3] = 4.0 * eta;                 dy[3] = 4.0 * xi;
        dx[4] = -4.0 * eta;                dy[4] = 4.0 - 4.0 * xi - 8.0 * eta;
        dx[5] = 4.0 - 8.0 * xi - 4.0 * eta; dy[5] = -4.0 * xi;
    }
};

// Q2 quad on NGSolve's ET_QUAD reference [0,1]^2: corners 0-3, edge mids
// 4 (bottom), 5 (right), 6 (top), 7 (left), face centre 8.
struct Q2Quad {
    static constexpr int N = 9;
    static constexpr bool kTriangle = false;
    static constexpr int NE = 4;
    static constexpr int EDGE[4][3] = { {0, 1, 4}, {1, 2, 5}, {2, 3, 6}, {3, 0, 7} };
    static void Eval(double xi, double eta, double n[9], double dx[9], double dy[9])
    {
        static constexpr int IX[9] = { 0, 2, 2, 0, 1, 2, 1, 0, 1 };
        static constexpr int IY[9] = { 0, 0, 2, 2, 0, 1, 2, 1, 1 };
        const double Lx[3]  = { 2*xi*xi - 3*xi + 1, -4*xi*xi + 4*xi, 2*xi*xi - xi };
        const double dLx[3] = { 4*xi - 3,           -8*xi + 4,        4*xi - 1 };
        const double Ly[3]  = { 2*eta*eta - 3*eta + 1, -4*eta*eta + 4*eta, 2*eta*eta - eta };
        const double dLy[3] = { 4*eta - 3,              -8*eta + 4,          4*eta - 1 };
        for (int k = 0; k < 9; ++k) {
            n[k]  = Lx[IX[k]]  * Ly[IY[k]];
            dx[k] = dLx[IX[k]] * Ly[IY[k]];
            dy[k] = Lx[IX[k]]  * dLy[IY[k]];
        }
    }
};

// Point data of one element at one reference point.
template <class S>
struct PointData {
    double n[S::N];      // basis values
    double gs[S::N];     // d N / d s
    double gz[S::N];     // d N / d z
    double s, z;         // physical (s, z)
    double detJ;         // |det d(s, z) / d(ref)|
};

template <class S>
class Geometry {
public:
    // rg, zg: nodes sampled from the mesh transformation, used for the (r, z)
    //         map of NGSolve reference points.
    // s, z:   nodes of the (s, z) map; r = sqrt(s) is the DOF node radius.
    double rg[S::N], zg[S::N];
    double r[S::N], z[S::N], s[S::N];
    double s_scale, z_scale;   // max |s|; z extent of the element
    double s_tol, z_tol;       // Newton residual tolerances

    void Init(const double * rs, const double * zs)
    {
        double smax = 0.0, zmin = zs[0], zmax = zs[0];
        for (int k = 0; k < S::N; ++k) {
            rg[k] = rs[k]; zg[k] = zs[k];
            z[k] = zs[k]; s[k] = rs[k] * rs[k];
            smax = std::max(smax, s[k]);
            zmin = std::min(zmin, zs[k]); zmax = std::max(zmax, zs[k]);
        }
        s_scale = smax > 0.0 ? smax : 1.0;
        z_scale = zmax > zmin ? zmax - zmin : 1.0;
        // An edge with exactly one end on the axis: an r-midpoint midnode makes
        // s ~ t^2 along the edge, so the (s, z) map is singular at the axis and
        // psi ~ r there (infinite 1/r energy).  Put that midnode at the
        // s-midpoint instead; the rule depends only on the edge, so both
        // neighbours use the same node.
        const double axis_tol = 1e-24 * s_scale;
        bool moved = false;
        for (int e = 0; e < S::NE; ++e) {
            int a = S::EDGE[e][0], b = S::EDGE[e][1], m = S::EDGE[e][2];
            bool on_a = s[a] <= axis_tol, on_b = s[b] <= axis_tol;
            if (on_a != on_b) {
                s[m] = 0.5 * (s[a] + s[b]);
                z[m] = 0.5 * (z[a] + z[b]);
                moved = true;
            }
        }
        if (moved && S::N == 9)   // keep the face centre consistent (bilinear blend)
            s[8] = 0.5 * (s[4] + s[5] + s[6] + s[7]) - 0.25 * (s[0] + s[1] + s[2] + s[3]);
        for (int k = 0; k < S::N; ++k)
            r[k] = std::sqrt(std::max(s[k], 0.0));
        // Residual tolerance: a small fraction of the element extent, but
        // never below the round-off of the coordinates themselves (thin
        // boundary-layer elements far from the axis).
        double smin = s[0], smx = s[0], zabs = 0.0;
        for (int k = 0; k < S::N; ++k) {
            smin = std::min(smin, s[k]); smx = std::max(smx, s[k]);
            zabs = std::max(zabs, std::fabs(z[k]));
        }
        const double eps = 2.220446049250313e-16;
        s_tol = std::max(1e-12 * (smx - smin), 64.0 * eps * s_scale);
        z_tol = std::max(1e-12 * z_scale, 64.0 * eps * zabs);
    }

    // (r, z) Lagrange map through the sampled nodes and its Jacobian
    // J[a][b] = d x_a / d ref_b with x = (r, z).
    void MapRZ(double xi, double eta, double & rp, double & zp, double J[2][2]) const
    {
        double n[S::N], dx[S::N], dy[S::N];
        S::Eval(xi, eta, n, dx, dy);
        rp = zp = 0.0; J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
        for (int k = 0; k < S::N; ++k) {
            rp += n[k] * rg[k];       zp += n[k] * zg[k];
            J[0][0] += dx[k] * rg[k]; J[0][1] += dy[k] * rg[k];
            J[1][0] += dx[k] * zg[k]; J[1][1] += dy[k] * zg[k];
        }
    }

    // Basis, (s, z) gradients and geometry at a reference point of the (s, z) map.
    void AtRef(double xi, double eta, PointData<S> & pd) const
    {
        double dx[S::N], dy[S::N];
        S::Eval(xi, eta, pd.n, dx, dy);
        double sp = 0.0, zp = 0.0, a = 0.0, b = 0.0, c = 0.0, d = 0.0;
        for (int k = 0; k < S::N; ++k) {
            sp += pd.n[k] * s[k]; zp += pd.n[k] * z[k];
            a += dx[k] * s[k]; b += dy[k] * s[k];
            c += dx[k] * z[k]; d += dy[k] * z[k];
        }
        double det = a * d - b * c;
        if (!(std::fabs(det) > 0.0))
            throw Exception("H1Henrotte: degenerate (s, z) element map");
        // grad_(s,z) N = J^{-T} grad_ref N with J = [[a, b], [c, d]].
        double inv = 1.0 / det;
        for (int k = 0; k < S::N; ++k) {
            pd.gs[k] = ( d * dx[k] - c * dy[k]) * inv;
            pd.gz[k] = (-b * dx[k] + a * dy[k]) * inv;
        }
        pd.s = sp; pd.z = zp; pd.detJ = std::fabs(det);
    }

    // Reference point of the (s, z) map at physical (st, zt), starting from
    // (xi, eta).  Fails loudly if Newton's method does not converge.
    void Locate(double st, double zt, double & xi, double & eta) const
    {
        // Converged when the (s, z) residual is within s_tol, z_tol (see
        // Init).  The residual, not the step, is the criterion.
        const double xi0 = xi, eta0 = eta;
        for (int it = 0; it < 60; ++it) {
            double n[S::N], dx[S::N], dy[S::N];
            S::Eval(xi, eta, n, dx, dy);
            double sp = 0.0, zp = 0.0, a = 0.0, b = 0.0, c = 0.0, d = 0.0;
            for (int k = 0; k < S::N; ++k) {
                sp += n[k] * s[k]; zp += n[k] * z[k];
                a += dx[k] * s[k]; b += dy[k] * s[k];
                c += dx[k] * z[k]; d += dy[k] * z[k];
            }
            double fs = sp - st, fz = zp - zt;
            if (std::fabs(fs) <= s_tol && std::fabs(fz) <= z_tol)
                return;
            double det = a * d - b * c;
            if (!(std::fabs(det) > 0.0)) break;
            xi  -= ( d * fs - b * fz) / det;
            eta -= (-c * fs + a * fz) / det;
        }
        throw Exception("H1Henrotte: (s, z) element map inversion did not converge"
                        " at s=" + std::to_string(st) + " z=" + std::to_string(zt)
                        + " from ref (" + std::to_string(xi0) + ", "
                        + std::to_string(eta0) + "), r0=" + std::to_string(r[0])
                        + " z0=" + std::to_string(z[0]));
    }

    // Basis values and reference-coordinate derivatives for an NGSolve
    // reference point (DiffOp path).
    void AtNGSolveRef(double xi, double eta, double * shape, double (*dref)[2]) const
    {
        double rp, zp, J[2][2];
        MapRZ(xi, eta, rp, zp, J);
        double u = xi, v = eta;
        Locate(rp * rp, zp, u, v);
        PointData<S> pd;
        AtRef(u, v, pd);
        for (int k = 0; k < S::N; ++k) {
            if (shape) shape[k] = pd.n[k];
            if (dref) {
                double gr = 2.0 * rp * pd.gs[k];
                double gz = pd.gz[k];
                dref[k][0] = gr * J[0][0] + gz * J[1][0];
                dref[k][1] = gr * J[0][1] + gz * J[1][1];
            }
        }
    }

    // Reference quadrature: body(xi, eta, w).
    template <typename F>
    static void Quadrature(F && body)
    {
        static constexpr int NG = 8;
        static constexpr double x[NG] = {
            0.01985507175123188, 0.10166676129318664, 0.23723379504183551,
            0.40828267875217508, 0.59171732124782492, 0.76276620495816449,
            0.89833323870681336, 0.98014492824876812 };
        static constexpr double w[NG] = {
            0.05061426814518813, 0.11119051722668724, 0.15685332293894364,
            0.18134189168918100, 0.18134189168918100, 0.15685332293894364,
            0.11119051722668724, 0.05061426814518813 };
        for (int i = 0; i < NG; ++i)
            for (int j = 0; j < NG; ++j) {
                if (S::kTriangle)   // Duffy collapse, weights sum to 1/2
                    body(x[i], x[j] * (1.0 - x[i]), w[i] * w[j] * (1.0 - x[i]));
                else
                    body(x[i], x[j], w[i] * w[j]);
            }
    }

    // K_ij = (2 pi / mu) r_i r_j INT grad psi_i . grad psi_j / r dr dz
    //      = (2 pi / mu) r_i r_j INT (2 psi_i,s psi_j,s + psi_i,z psi_j,z / (2 s)) ds dz
    void Stiffness(double mu, FlatMatrix<double> elmat) const
    {
        double K[S::N][S::N] = {};
        Quadrature([&](double xi, double eta, double w) {
            PointData<S> pd;
            AtRef(xi, eta, pd);
            if (pd.s <= 1e-30 * s_scale) return;
            double f = w * pd.detJ;
            for (int i = 0; i < S::N; ++i)
                for (int j = 0; j < S::N; ++j)
                    K[i][j] += f * (2.0 * pd.gs[i] * pd.gs[j]
                                    + pd.gz[i] * pd.gz[j] / (2.0 * pd.s));
        });
        const double two_pi = 6.283185307179586476925;
        for (int i = 0; i < S::N; ++i)
            for (int j = 0; j < S::N; ++j)
                elmat(i, j) = two_pi / mu * r[i] * r[j] * 0.5 * (K[i][j] + K[j][i]);
    }

    // M_ij = 2 pi sigma r_i r_j INT psi_i psi_j / r dr dz
    //      = 2 pi sigma r_i r_j INT psi_i psi_j / (2 s) ds dz
    void SigmaMass(double sigma, FlatMatrix<double> elmat) const
    {
        double M[S::N][S::N] = {};
        Quadrature([&](double xi, double eta, double w) {
            PointData<S> pd;
            AtRef(xi, eta, pd);
            if (pd.s <= 1e-30 * s_scale) return;
            double f = w * pd.detJ / (2.0 * pd.s);
            for (int i = 0; i < S::N; ++i)
                for (int j = 0; j < S::N; ++j)
                    M[i][j] += f * pd.n[i] * pd.n[j];
        });
        const double two_pi = 6.283185307179586476925;
        for (int i = 0; i < S::N; ++i)
            for (int j = 0; j < S::N; ++j)
                elmat(i, j) = two_pi * sigma * r[i] * r[j] * 0.5 * (M[i][j] + M[j][i]);
    }
};

}  // namespace siso
}  // namespace axifem

#endif  // RADIA_AXIFEM_AXI_HENROTTE_SISO_HPP
