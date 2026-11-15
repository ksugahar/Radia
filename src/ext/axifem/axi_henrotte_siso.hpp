// axi_henrotte_siso.hpp -- s-isoparametric geometry of the general Henrotte
// elements (P1 and P2 triangles, Q2 general quads).
//
// Field convention (one for every evaluator and integrator of H1Henrotte):
// the DOFs are nodal values of A_phi and the basis is the Henrotte V-DOF basis
//     N_i = r_i psi_i / r,   psi = 2 pi r A_phi  (flux function),
// where psi_i is the element's scalar Lagrange basis.  psi ~ s = r^2 at the
// axis keeps the 1/r stiffness finite; axis nodes (r_i = 0) carry no basis.
//
// Geometry: a quadratic Lagrange map from the reference element to (s, z)
// through the element nodes.  The basis psi_i is the reference Lagrange basis
// (linear for P1, quadratic for P2 and Q2), so the trace on an edge depends
// only on that edge's nodes and the edge parameter: the space is conforming.
//
// Edge midnode rule (EdgeMidnode, shared by every family and by the boundary
// trace): the (s, z) midnode of an edge is the s-midpoint ((s_a+s_b)/2,
// (z_a+z_b)/2) when the edge is horizontal (z = const) or has exactly one end
// on the axis, and the mesh midpoint otherwise.
//   - mesh midpoint: s = r^2 is reproduced exactly along a straight edge, so
//     the element map equals the mesh map there.
//   - horizontal: the s-midpoint keeps the edge on its line z = const (the same
//     point set as the mesh edge) and makes the trace quadratic in s, the
//     trace of the closed-form axis-aligned quads.
//   - one end on the axis: a mesh midpoint would give s ~ t^2 along the edge
//     and psi ~ r at the axis (infinite 1/r energy).  The s-midpoint keeps
//     psi ~ s.  An oblique such edge becomes a parabola in (r, z): the element
//     occupies a slightly different region than the mesh element (an O(h)
//     shape change confined to elements with an oblique edge ending on the
//     axis).
// Only straight mesh elements are admitted (Curve(1)); curved, degenerate or
// folded (s, z) maps are refused.
//
// Point evaluation from NGSolve reference coordinates (EvalVdof): the target s
// is the mesh r^2 (the evaluation radius is always the mesh radius, so A_phi ~ r
// and B stay regular at the axis); the target z is the mesh z, corrected on the
// oblique axis edges.  On such an edge the mesh point at radius r is matched to
// the element edge point with the same s (edge parameter tau = u^2 for the
// fraction u from the axis end; the boundary trace uses the same rule), which
// lies (u - u^2) (z_off - z_axis) below in z; that shift times the edge's
// midnode function (zero on every other edge and on the axis) is added to the
// mesh z.  The (s, z) map is inverted at the target by Newton's method.  For
// every element without an oblique axis edge the target is the mesh point
// itself: values and gradients are those of the discrete field at the
// physical point.

#ifndef RADIA_AXIFEM_AXI_HENROTTE_SISO_HPP
#define RADIA_AXIFEM_AXI_HENROTTE_SISO_HPP

#include <fem.hpp>
#include <algorithm>
#include <cmath>
#include <string>

namespace axifem {
namespace siso {

using namespace ngfem;

// P1 basis on NGSolve's ET_TRIG reference: V0 at (1,0), V1 at (0,1), V2 at
// (0,0).  Used with the P2Tri geometry map (Geometry<P2Tri, P1Tri>): linear in
// the reference coordinates, nodes at the three vertices.
struct P1Tri {
    static constexpr int N = 3;
    static constexpr bool kTriangle = true;
    static void Eval(double xi, double eta, double n[3], double dx[3], double dy[3])
    {
        n[0] = xi; n[1] = eta; n[2] = 1.0 - xi - eta;
        dx[0] = 1.0;  dy[0] = 0.0;
        dx[1] = 0.0;  dy[1] = 1.0;
        dx[2] = -1.0; dy[2] = -1.0;
    }
};

// P2 triangle on NGSolve's ET_TRIG reference: V0 at (1,0), V1 at (0,1), V2 at
// (0,0); local order [v0, v1, v2, m01, m12, m20].
struct P2Tri {
    static constexpr int N = 6;
    static constexpr bool kTriangle = true;
    static constexpr int NE = 3;
    static constexpr int EDGE[3][3] = { {0, 1, 3}, {1, 2, 4}, {2, 0, 5} };  // a, b, mid
    static void NodeRef(int k, double & xi, double & eta)
    {
        static constexpr double X[6] = { 1, 0, 0, 0.5, 0, 0.5 };
        static constexpr double Y[6] = { 0, 1, 0, 0.5, 0.5, 0 };
        xi = X[k]; eta = Y[k];
    }
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
    static void NodeRef(int k, double & xi, double & eta)
    {
        static constexpr double X[9] = { 0, 1, 1, 0, 0.5, 1, 0.5, 0, 0.5 };
        static constexpr double Y[9] = { 0, 0, 1, 1, 0, 0.5, 1, 0.5, 0.5 };
        xi = X[k]; eta = Y[k];
    }
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

// Shared axis and horizontal-edge tests.  Both depend only on the edge, so the
// two neighbours of an edge and its boundary element classify it identically.
inline bool OnAxisEdgeEnd(double r, double ra, double rb)
{
    return std::fabs(r) <= 1e-12 * std::max(std::fabs(ra), std::fabs(rb));
}

inline bool HorizontalEdge(double ra, double za, double rb, double zb)
{
    return std::fabs(zb - za) <= 1e-12 * std::hypot(rb - ra, zb - za);
}

enum class MidKind { Mesh, SMidHorizontal, SMidAxis };

// Edge midnode rule (see the file comment).  rm, zm: mesh midpoint.
inline MidKind EdgeMidnode(double ra, double za, double rb, double zb,
                           double rm, double zm, double & sm, double & zmid)
{
    const bool on_a = OnAxisEdgeEnd(ra, ra, rb), on_b = OnAxisEdgeEnd(rb, ra, rb);
    if (HorizontalEdge(ra, za, rb, zb) || on_a != on_b) {
        sm = 0.5 * (ra * ra + rb * rb);
        zmid = 0.5 * (za + zb);
        return HorizontalEdge(ra, za, rb, zb) ? MidKind::SMidHorizontal : MidKind::SMidAxis;
    }
    sm = rm * rm;
    zmid = zm;
    return MidKind::Mesh;
}

// Henrotte V-DOF basis from the scalar basis: N = r_i psi / r with physical
// gradient dN/dr = r_i (2 psi_s - psi / s), dN/dz = r_i psi_z / r (psi_r =
// 2 r psi_s).  Near the axis (s <= 1e-8 s_scale) the leading term psi ~ psi_s s
// is used, N = r_i psi_s r and dN/dr = r_i psi_s, which avoids the
// cancellation in psi / s; on the axis N = 0, dN/dr = r_i psi_s, dN/dz = 0.
inline void VdofFromPsi(double ri, double s, double s_scale, double psi,
                        double ps, double pz, double & n, double * dn)
{
    if (ri == 0.0) {
        n = 0.0;
        if (dn) dn[0] = dn[1] = 0.0;
        return;
    }
    const double r = std::sqrt(std::max(s, 0.0));
    if (s <= 1e-8 * s_scale) {
        n = ri * ps * r;
        if (dn) { dn[0] = ri * ps; dn[1] = r > 0.0 ? ri * pz / r : 0.0; }
        return;
    }
    n = ri * psi / r;
    if (dn) { dn[0] = ri * (2.0 * ps - psi / s); dn[1] = ri * pz / r; }
}

// Point data of the basis B at one reference point.
template <class B>
struct PointData {
    double n[B::N];      // scalar basis psi_i
    double gs[B::N];     // d psi_i / d s
    double gz[B::N];     // d psi_i / d z
    double s, z;         // physical (s, z)
    double det;          // det d(s, z) / d(ref), signed
};

// G: geometry map (P2Tri or Q2Quad).  B: basis (P1Tri, P2Tri or Q2Quad), whose
// nodes are the first B::N geometry nodes.
template <class G, class B = G>
class Geometry {
public:
    static_assert(G::kTriangle == B::kTriangle, "basis and geometry element types differ");
    static_assert(B::N <= G::N, "basis nodes must be geometry nodes");
    double rg[G::N], zg[G::N];   // mesh nodes (straight element)
    double s[G::N], z[G::N];     // (s, z) map nodes after the edge midnode rule
    double target_dz[G::N];      // z correction of oblique axis edges (per midnode)
    bool has_target_dz = false;
    double r[B::N];              // DOF node radii
    double s_scale, z_scale;     // max s; z extent of the element
    double s_tol, z_tol;         // Newton residual tolerances

    void Init(const double * rs, const double * zs)
    {
        double rmax = 0.0, zmin = zs[0], zmax = zs[0];
        for (int k = 0; k < G::N; ++k) {
            rg[k] = rs[k]; zg[k] = zs[k];
            s[k] = rs[k] * rs[k]; z[k] = zs[k]; target_dz[k] = 0.0;
            rmax = std::max(rmax, std::fabs(rs[k]));
            zmin = std::min(zmin, zs[k]); zmax = std::max(zmax, zs[k]);
        }
        s_scale = rmax > 0.0 ? rmax * rmax : 1.0;
        z_scale = zmax > zmin ? zmax - zmin : 1.0;
        // Straight elements only: every mesh midnode is the chord midpoint and
        // the quad centre is the corner average.  A curved (Curve > 1) element
        // would make the (s, z) map and the mesh map disagree in the interior.
        auto require_straight = [&](double rm, double zm, double ra, double za,
                                    double rb, double zb) {
            const double len = std::hypot(rb - ra, zb - za);
            if (std::hypot(rm - 0.5 * (ra + rb), zm - 0.5 * (za + zb)) > 1e-10 * len)
                throw Exception("H1Henrotte: curved element at r0=" + std::to_string(rg[0])
                                + " z0=" + std::to_string(zg[0])
                                + "; only straight elements are supported (mesh.Curve(1))");
        };
        bool moved = false;
        for (int e = 0; e < G::NE; ++e) {
            const int a = G::EDGE[e][0], b = G::EDGE[e][1], m = G::EDGE[e][2];
            require_straight(rg[m], zg[m], rg[a], zg[a], rg[b], zg[b]);
            double sm, zm;
            const MidKind kind = EdgeMidnode(rg[a], zg[a], rg[b], zg[b], rg[m], zg[m], sm, zm);
            if (kind == MidKind::Mesh) continue;
            if (kind == MidKind::SMidAxis) {
                // N_m = 4 u (1 - u) on the edge; the matched point is
                // (u^2 - u) (z_off - z_axis) = -N_m (z_off - z_axis) / 4 away in z.
                const bool a_on_axis = OnAxisEdgeEnd(rg[a], rg[a], rg[b]);
                const double z_axis = a_on_axis ? zg[a] : zg[b];
                const double z_off = a_on_axis ? zg[b] : zg[a];
                target_dz[m] = -0.25 * (z_off - z_axis);
                has_target_dz = true;
            }
            s[m] = sm; z[m] = zm;
            moved = true;
        }
        if (G::N == 9) {
            require_straight(rg[8], zg[8], 0.5 * (rg[0] + rg[1]), 0.5 * (zg[0] + zg[1]),
                             0.5 * (rg[2] + rg[3]), 0.5 * (zg[2] + zg[3]));
            if (moved)   // keep the centre consistent with the edges (bilinear blend)
                s[8] = 0.5 * (s[4] + s[5] + s[6] + s[7]) - 0.25 * (s[0] + s[1] + s[2] + s[3]);
        }
        for (int k = 0; k < B::N; ++k)
            r[k] = std::sqrt(std::max(s[k], 0.0));
        // Residual tolerance: a small fraction of the element extent, but
        // never below the round-off of the coordinates themselves.
        double smin = s[0], smx = s[0], zabs = 0.0;
        for (int k = 0; k < G::N; ++k) {
            smin = std::min(smin, s[k]); smx = std::max(smx, s[k]);
            zabs = std::max(zabs, std::fabs(z[k]));
        }
        const double eps = 2.220446049250313e-16;
        s_tol = std::max(1e-12 * (smx - smin), 64.0 * eps * s_scale);
        z_tol = std::max(1e-12 * z_scale, 64.0 * eps * zabs);
        // Admission: the (s, z) map must be nondegenerate and keep the mesh
        // orientation (s = r^2 increases with r) at every quadrature point and
        // every node.
        double rp, zp, J[2][2];
        const double c = G::kTriangle ? 1.0 / 3.0 : 0.5;
        MapRZ(c, c, rp, zp, J);
        const double sgn = (J[0][0] * J[1][1] - J[0][1] * J[1][0]) > 0.0 ? 1.0 : -1.0;
        const double ds = smx > smin ? smx - smin : s_scale;
        const double det_min = 1e-12 * ds * z_scale;
        auto check = [&](double xi, double eta) {
            const double d = sgn * DetAt(xi, eta);
            if (!(d > det_min))
                throw Exception(std::string("H1Henrotte: ")
                                + (d < -det_min ? "folded" : "degenerate")
                                + " (s, z) element map at r0=" + std::to_string(rg[0])
                                + " z0=" + std::to_string(zg[0]));
        };
        Quadrature([&](double xi, double eta, double) { check(xi, eta); });
        for (int k = 0; k < G::N; ++k) {
            double xi, eta;
            G::NodeRef(k, xi, eta);
            check(xi, eta);
        }
    }

    double DetAt(double xi, double eta) const
    {
        double n[G::N], dx[G::N], dy[G::N];
        G::Eval(xi, eta, n, dx, dy);
        double a = 0.0, b = 0.0, c = 0.0, d = 0.0;
        for (int k = 0; k < G::N; ++k) {
            a += dx[k] * s[k]; b += dy[k] * s[k];
            c += dx[k] * z[k]; d += dy[k] * z[k];
        }
        return a * d - b * c;
    }

    // Mesh (r, z) map through the mesh nodes and its Jacobian.
    void MapRZ(double xi, double eta, double & rp, double & zp, double J[2][2]) const
    {
        double n[G::N], dx[G::N], dy[G::N];
        G::Eval(xi, eta, n, dx, dy);
        rp = zp = 0.0; J[0][0] = J[0][1] = J[1][0] = J[1][1] = 0.0;
        for (int k = 0; k < G::N; ++k) {
            rp += n[k] * rg[k];       zp += n[k] * zg[k];
            J[0][0] += dx[k] * rg[k]; J[0][1] += dy[k] * rg[k];
            J[1][0] += dx[k] * zg[k]; J[1][1] += dy[k] * zg[k];
        }
    }

    // Basis, (s, z) gradients and geometry at a reference point of the (s, z) map.
    void AtRef(double xi, double eta, PointData<B> & pd) const
    {
        double gn[G::N], gx[G::N], gy[G::N];
        G::Eval(xi, eta, gn, gx, gy);
        double sp = 0.0, zp = 0.0, a = 0.0, b = 0.0, c = 0.0, d = 0.0;
        for (int k = 0; k < G::N; ++k) {
            sp += gn[k] * s[k]; zp += gn[k] * z[k];
            a += gx[k] * s[k]; b += gy[k] * s[k];
            c += gx[k] * z[k]; d += gy[k] * z[k];
        }
        const double det = a * d - b * c;
        if (!(std::fabs(det) > 0.0))
            throw Exception("H1Henrotte: degenerate (s, z) element map");
        double dx[B::N], dy[B::N];
        B::Eval(xi, eta, pd.n, dx, dy);
        // grad_(s,z) psi = J^{-T} grad_ref psi with J = [[a, b], [c, d]].
        const double inv = 1.0 / det;
        for (int k = 0; k < B::N; ++k) {
            pd.gs[k] = ( d * dx[k] - c * dy[k]) * inv;
            pd.gz[k] = (-b * dx[k] + a * dy[k]) * inv;
        }
        pd.s = sp; pd.z = zp; pd.det = det;
    }

    // Reference point of the (s, z) map at (st, zt), starting from (xi, eta).
    // Fails loudly if Newton's method does not converge.
    void Locate(double st, double zt, double & xi, double & eta) const
    {
        const double xi0 = xi, eta0 = eta;
        for (int it = 0; it < 60; ++it) {
            double n[G::N], dx[G::N], dy[G::N];
            G::Eval(xi, eta, n, dx, dy);
            double sp = 0.0, zp = 0.0, a = 0.0, b = 0.0, c = 0.0, d = 0.0;
            for (int k = 0; k < G::N; ++k) {
                sp += n[k] * s[k]; zp += n[k] * z[k];
                a += dx[k] * s[k]; b += dy[k] * s[k];
                c += dx[k] * z[k]; d += dy[k] * z[k];
            }
            const double fs = sp - st, fz = zp - zt;
            if (std::fabs(fs) <= s_tol && std::fabs(fz) <= z_tol)
                return;
            const double det = a * d - b * c;
            if (!(std::fabs(det) > 0.0)) break;
            xi  -= ( d * fs - b * fz) / det;
            eta -= (-c * fs + a * fz) / det;
        }
        throw Exception("H1Henrotte: (s, z) element map inversion did not converge"
                        " at s=" + std::to_string(st) + " z=" + std::to_string(zt)
                        + " from ref (" + std::to_string(xi0) + ", "
                        + std::to_string(eta0) + "), r0=" + std::to_string(rg[0])
                        + " z0=" + std::to_string(zg[0]));
    }

    // V-DOF basis N_i and its physical (r, z) gradient at an NGSolve reference
    // point (see the file comment for the evaluation target).
    void EvalVdof(double xi, double eta, double * nv, double (*dnv)[2]) const
    {
        double rp, zp, J[2][2];
        MapRZ(xi, eta, rp, zp, J);
        double zt = zp;
        if (has_target_dz) {
            double n[G::N], dx[G::N], dy[G::N];
            G::Eval(xi, eta, n, dx, dy);
            for (int k = 0; k < G::N; ++k) zt += n[k] * target_dz[k];
        }
        double u = xi, v = eta;
        Locate(rp * rp, zt, u, v);
        PointData<B> pd;
        AtRef(u, v, pd);
        for (int k = 0; k < B::N; ++k) {
            double dn[2];
            // the target radius is exact; pd.s carries the Newton tolerance,
            // which dominates s itself near the axis
            VdofFromPsi(r[k], rp * rp, s_scale, pd.n[k], pd.gs[k], pd.gz[k], nv[k],
                        dnv ? dn : nullptr);
            if (dnv) { dnv[k][0] = dn[0]; dnv[k][1] = dn[1]; }
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
                if (G::kTriangle)   // Duffy collapse, weights sum to 1/2
                    body(x[i], x[j] * (1.0 - x[i]), w[i] * w[j] * (1.0 - x[i]));
                else
                    body(x[i], x[j], w[i] * w[j]);
            }
    }

    // K_ij = (2 pi / mu) INT grad N_i . grad N_j r dr dz with N_i = r_i psi_i / r
    //      = (2 pi / mu) r_i r_j INT (2 psi_i,s psi_j,s + psi_i,z psi_j,z / (2 s)) ds dz
    void Stiffness(double mu, FlatMatrix<double> elmat) const
    {
        double K[B::N][B::N] = {};
        Quadrature([&](double xi, double eta, double w) {
            PointData<B> pd;
            AtRef(xi, eta, pd);
            if (pd.s <= 1e-30 * s_scale) return;
            const double f = w * std::fabs(pd.det);
            for (int i = 0; i < B::N; ++i)
                for (int j = 0; j < B::N; ++j)
                    K[i][j] += f * (2.0 * pd.gs[i] * pd.gs[j]
                                    + pd.gz[i] * pd.gz[j] / (2.0 * pd.s));
        });
        const double two_pi = 6.283185307179586476925;
        for (int i = 0; i < B::N; ++i)
            for (int j = 0; j < B::N; ++j)
                elmat(i, j) = two_pi / mu * r[i] * r[j] * 0.5 * (K[i][j] + K[j][i]);
    }

    // M_ij = 2 pi sigma INT N_i N_j r dr dz
    //      = 2 pi sigma r_i r_j INT psi_i psi_j / (2 s) ds dz
    void SigmaMass(double sigma, FlatMatrix<double> elmat) const
    {
        double M[B::N][B::N] = {};
        Quadrature([&](double xi, double eta, double w) {
            PointData<B> pd;
            AtRef(xi, eta, pd);
            if (pd.s <= 1e-30 * s_scale) return;
            const double f = w * std::fabs(pd.det) / (2.0 * pd.s);
            for (int i = 0; i < B::N; ++i)
                for (int j = 0; j < B::N; ++j)
                    M[i][j] += f * pd.n[i] * pd.n[j];
        });
        const double two_pi = 6.283185307179586476925;
        for (int i = 0; i < B::N; ++i)
            for (int j = 0; j < B::N; ++j)
                elmat(i, j) = two_pi * sigma * r[i] * r[j] * 0.5 * (M[i][j] + M[j][i]);
    }
};

}  // namespace siso
}  // namespace axifem

#endif  // RADIA_AXIFEM_AXI_HENROTTE_SISO_HPP
