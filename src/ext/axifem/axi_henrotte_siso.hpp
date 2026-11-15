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
// Curved edges (Curve(2): midnode off the chord).  Off the axis they take the
// mesh midpoint; along such an edge (r(t)^2, z(t)) is quartic in s while the
// (s, z) edge is the quadratic s_q(t) through the same nodes, so the target s
// is corrected by the blended edge ratio rho_e(t) = 1 - s_q(t) / r(t)^2.  A
// curved edge with one end on the axis takes the s-midpoint like a straight
// one; its (s, z) edge is the straight segment between the end nodes, matched
// to the mesh edge through s (tau = (r(t)^2 - s_a) / (s_b - s_a)), and the
// target z is corrected by the blended shift dz_e(t) = z_a + tau (z_b - z_a) -
// z(t) (see Target below).  Degenerate or folded maps are refused.
//
// Point evaluation from NGSolve reference coordinates (EvalVdof): the target s
// is the mesh r^2 (the evaluation radius is always the mesh radius, so A_phi ~ r
// and B stay regular at the axis); the target z is the mesh z, corrected on the
// oblique axis edges.  On such an edge the mesh point at radius r is matched to
// the element edge point with the same s (edge parameter tau = u^2 for the
// fraction u from the axis end; the boundary trace uses the same rule), which
// lies (u - u^2) (z_off - z_axis) below in z; that shift times the edge's
// midnode function (zero on every other edge and on the axis) is added to the
// mesh z.  On curved edges the target s is r^2 (1 - sum_e B_e) with the
// blended edge ratios B_e (zero on every other edge and at the vertices), so
// the target point of a curved edge point at parameter t is the (s, z) edge
// point at the same t.  The (s, z) map is inverted at the target by Newton's method.  The
// basis on the mesh element is psi(Phi(r, z)) with this target map Phi; the
// map is the same from both sides of every edge, so the space stays conforming
// and the boundary trace equals the volume value.  The physical gradient is
// taken through Phi (chain rule), so it is the derivative of the value, and an
// element whose Phi is not (r^2, z) assembles its stiffness and mass by
// quadrature on the mesh element through the same evaluator.  Elsewhere Phi =
// (r^2, z) and the (s, z) quadrature is the same integral.

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

// A curved edge: its mesh midnode is off the chord.
inline bool CurvedEdge(double ra, double za, double rb, double zb, double rm, double zm)
{
    return std::hypot(rm - 0.5 * (ra + rb), zm - 0.5 * (za + zb))
           > 1e-10 * std::hypot(rb - ra, zb - za);
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
    if (s <= 1e-24 * s_scale) {
        // on the axis within the snap tolerance of the space (r <= 1e-12 of
        // the element): A = 0, dA/dz = 0; mesh maps of curved elements put
        // their axis points at |r| ~ 1e-16 rather than exactly 0
        n = 0.0;
        if (dn) { dn[0] = ri * ps; dn[1] = 0.0; }
        return;
    }
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
    bool curved[G::NE] = {};     // curved (Curve > 1) edges off the axis
    bool curved_axis[G::NE] = {};  // curved edges with one end on the axis
    bool has_curved = false;
    double r[B::N];              // DOF node radii
    double s_scale, z_scale;     // max s; z extent of the element
    double s_tol, z_tol;         // Newton residual tolerances

    void Init(const double * rs, const double * zs)
    {
        // (nodes within round-off of the axis arrive at r = 0: the space snaps
        // them with one mesh-wide tolerance, AxiHenrotteFESpace::SnapAxis)
        double rmax = 0.0, zmin = zs[0], zmax = zs[0];
        for (int k = 0; k < G::N; ++k) {
            rg[k] = rs[k]; zg[k] = zs[k];
            s[k] = rg[k] * rg[k]; z[k] = zs[k]; target_dz[k] = 0.0;
            rmax = std::max(rmax, std::fabs(rg[k]));
            zmin = std::min(zmin, zs[k]); zmax = std::max(zmax, zs[k]);
        }
        s_scale = rmax > 0.0 ? rmax * rmax : 1.0;
        z_scale = zmax > zmin ? zmax - zmin : 1.0;
        bool moved = false;
        for (int e = 0; e < G::NE; ++e) {
            const int a = G::EDGE[e][0], b = G::EDGE[e][1], m = G::EDGE[e][2];
            if (CurvedEdge(rg[a], zg[a], rg[b], zg[b], rg[m], zg[m])) {
                const bool on_a = OnAxisEdgeEnd(rg[a], rg[a], rg[b]);
                const bool on_b = OnAxisEdgeEnd(rg[b], rg[a], rg[b]);
                if (on_a && on_b)
                    throw Exception("H1Henrotte: curved edge on the axis at r0="
                                    + std::to_string(rg[0]) + " z0=" + std::to_string(zg[0]));
                // r > 0 inside the edge, r^2 monotone when one end is on the axis
                double prev = -1.0;
                for (int i = 0; i <= 64; ++i) {
                    const double ri = EdgeQuad(i / 64.0, rg[a], rg[b], rg[m]);
                    const double si = on_b ? EdgeQuad(1.0 - i / 64.0, rg[a], rg[b], rg[m]) : ri;
                    if ((i > 0 && i < 64 && !(ri > 0.0))
                        || ((on_a || on_b) && i > 0 && !(si * si > prev)))
                        throw Exception("H1Henrotte: curved edge reaches r <= 0 or turns back"
                                        " in r at r0=" + std::to_string(rg[0])
                                        + " z0=" + std::to_string(zg[0]));
                    prev = si * si;
                }
                has_curved = true;
                if (on_a || on_b) {
                    // s-midpoint, straight (s, z) edge; target z blended (Target)
                    curved_axis[e] = true;
                    s[m] = 0.5 * (s[a] + s[b]);
                    z[m] = 0.5 * (zg[a] + zg[b]);
                    moved = true;
                } else {
                    curved[e] = true;   // mesh midpoint; target s blended (Target)
                }
                continue;
            }
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
            if (!has_curved
                && CurvedEdge(0.5 * (rg[0] + rg[1]), 0.5 * (zg[0] + zg[1]),
                              0.5 * (rg[2] + rg[3]), 0.5 * (zg[2] + zg[3]), rg[8], zg[8]))
                throw Exception("H1Henrotte: quad centre off the corner average at r0="
                                + std::to_string(rg[0]) + " z0=" + std::to_string(zg[0])
                                + " with straight edges");
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
        // Curved elements: the target map (r, z) -> (s_t, z_t) must keep the
        // orientation; det d(s_t, z_t)/d(r, z) / (2 r) > 0 at every quadrature point.
        if (has_curved)
            Quadrature([&](double xi, double eta, double) {
                TargetPoint tp;
                Target(xi, eta, tp);
                const double d = tp.sr * tp.zz - tp.sz * tp.zr;
                if (!(d > 1e-6))
                    throw Exception("H1Henrotte: curved element too strongly curved for the"
                                    " (s, z) target map at r0=" + std::to_string(rg[0])
                                    + " z0=" + std::to_string(zg[0]));
            });
    }

    // Quadratic Lagrange along an edge at parameter t (nodes t = 0, 1, 1/2).
    static double EdgeQuad(double t, double fa, double fb, double fm)
    {
        return (1.0 - t) * (1.0 - 2.0 * t) * fa + t * (2.0 * t - 1.0) * fb
               + 4.0 * t * (1.0 - t) * fm;
    }
    static double EdgeQuadD(double t, double fa, double fb, double fm)
    {
        return (4.0 * t - 3.0) * fa + (4.0 * t - 1.0) * fb + (4.0 - 8.0 * t) * fm;
    }

    // Target point of the basis at a reference point.  s_t = r^2 (1 - Bsum),
    // Bsum = sum over curved edges e of w_e rho_e(t_e), rho_e(t) = 1 - s_q(t) /
    // r_e(t)^2 (zero at the edge's nodes), with the edge parameter t_e and a
    // blend weight w_e that is 1 on edge e and 0 on the other edges (triangle:
    // w = (l_a + l_b)^2, t = l_b / (l_a + l_b); quad: w = 1 - distance from the
    // edge, t along it).  z_t = z + sum_k n_k dz_k.  Stored derivatives stay
    // finite at the axis: sr = s_t,r / (2 r) = (1 - Bsum) - r Bsum_r / 2,
    // sz = s_t,z, zr = z_t,r / (2 r) (0 at r = 0), zz = z_t,z.
    struct TargetPoint {
        double r, z, st, zt;
        double sr, sz, zr, zz;
        double det;    // |det d(r, z) / d(ref)|
    };
    void Target(double xi, double eta, TargetPoint & tp) const
    {
        double rp, zp, J[2][2];
        MapRZ(xi, eta, rp, zp, J);
        const double dj = J[0][0] * J[1][1] - J[0][1] * J[1][0];
        auto phys = [&](double fx, double fy, double & fr, double & fz) {
            fr = ( fx * J[1][1] - fy * J[1][0]) / dj;
            fz = (-fx * J[0][1] + fy * J[0][0]) / dj;
        };
        double zt = zp, txi = J[1][0], teta = J[1][1];   // z_t and d z_t / d(ref)
        if (has_target_dz) {
            double n[G::N], dx[G::N], dy[G::N];
            G::Eval(xi, eta, n, dx, dy);
            for (int k = 0; k < G::N; ++k) {
                zt += n[k] * target_dz[k];
                txi += dx[k] * target_dz[k];
                teta += dy[k] * target_dz[k];
            }
        }
        double bs = 0.0, bx = 0.0, by = 0.0;
        if (has_curved) {
            for (int e = 0; e < G::NE; ++e) {
                if (!curved[e] && !curved_axis[e]) continue;
                const int a = G::EDGE[e][0], b = G::EDGE[e][1], m = G::EDGE[e][2];
                double t, tx, ty, w, wx, wy;
                if (G::kTriangle) {
                    double l[3], lx[3], ly[3];
                    P1Tri::Eval(xi, eta, l, lx, ly);
                    const double S = l[a] + l[b];
                    if (S <= 1e-14) continue;    // opposite vertex: w rho -> 0
                    const double Sx = lx[a] + lx[b], Sy = ly[a] + ly[b];
                    t = l[b] / S;
                    tx = (lx[b] * S - l[b] * Sx) / (S * S);
                    ty = (ly[b] * S - l[b] * Sy) / (S * S);
                    w = S * S; wx = 2.0 * S * Sx; wy = 2.0 * S * Sy;
                } else {
                    double ax, ay, ex, ey;
                    G::NodeRef(a, ax, ay); G::NodeRef(b, ex, ey);
                    const double dxe = ex - ax, dye = ey - ay;     // unit edge direction
                    const double nx = -dye, ny = dxe;
                    const double sg = ((0.5 - ax) * nx + (0.5 - ay) * ny) > 0.0 ? 1.0 : -1.0;
                    t = (xi - ax) * dxe + (eta - ay) * dye; tx = dxe; ty = dye;
                    w = 1.0 - sg * ((xi - ax) * nx + (eta - ay) * ny);
                    wx = -sg * nx; wy = -sg * ny;
                }
                const double re = EdgeQuad(t, rg[a], rg[b], rg[m]);
                const double red = EdgeQuadD(t, rg[a], rg[b], rg[m]);
                if (curved_axis[e]) {
                    // dz_e(t) = z_a + tau (z_b - z_a) - z(t), tau = (r^2 - s_a)/(s_b - s_a)
                    const double tau = (re * re - s[a]) / (s[b] - s[a]);
                    const double taud = 2.0 * re * red / (s[b] - s[a]);
                    const double dz = zg[a] + tau * (zg[b] - zg[a]) - EdgeQuad(t, zg[a], zg[b], zg[m]);
                    const double dzd = taud * (zg[b] - zg[a]) - EdgeQuadD(t, zg[a], zg[b], zg[m]);
                    zt += w * dz;
                    txi += wx * dz + w * dzd * tx;
                    teta += wy * dz + w * dzd * ty;
                    continue;
                }
                const double sq = EdgeQuad(t, s[a], s[b], s[m]);
                const double sqd = EdgeQuadD(t, s[a], s[b], s[m]);
                const double rho = 1.0 - sq / (re * re);
                const double rhod = -sqd / (re * re) + 2.0 * sq * red / (re * re * re);
                bs += w * rho;
                bx += wx * rho + w * rhod * tx;
                by += wy * rho + w * rhod * ty;
            }
        }
        double gr, gz, br = 0.0, bz = 0.0;
        phys(txi, teta, gr, gz);
        if (has_curved) phys(bx, by, br, bz);
        tp.r = rp; tp.z = zp; tp.det = std::fabs(dj);
        tp.st = rp * rp * (1.0 - bs);
        tp.zt = zt;
        tp.sr = (1.0 - bs) - 0.5 * rp * br;
        tp.sz = -rp * rp * bz;
        // on the axis (within the snap band of VdofFromPsi) psi_zt -> 0 like s,
        // so psi_zt dz_t/ds -> 0
        tp.zr = rp * rp > 1e-24 * s_scale ? gr / (2.0 * rp) : 0.0;
        tp.zz = gz;
    }

    // Elements whose target map is not (r^2, z) are integrated on the mesh
    // element through EvalVdof.
    bool PhysicalTarget() const { return has_target_dz || has_curved; }

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
        // At least one Newton step: near the axis the reference distance from
        // the axis edge is O(s) on the (s, z) map but O(r) at the mesh starting
        // point, whose residual O(r s_scale^0.5) falls below s_tol there; accepting it
        // would leave an O(1) error in psi_zt dz_t/ds = psi_zt dz_t/dr / (2r).
        // After one step the error is O(r^2), so that term errs by O(r).
        // Once within tolerance one more step is taken: the residual of a
        // quadratically converging step is far below s_tol, and psi_z / r near
        // the axis would otherwise carry the location error s_tol / r.
        bool polished = false;
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
            if (it > 0 && std::fabs(fs) <= s_tol && std::fabs(fz) <= z_tol) {
                if (polished) return;
                polished = true;
            }
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
    //
    // The basis on the mesh element is psi_i(Phi(r, z)) with the target map
    // Phi(r, z) = (r^2, z_t), z_t = z + sum_k n_k dz_k; its physical gradient is
    // taken through Phi (chain rule), so it is the derivative of the value:
    //   psi_r = 2 r psi_s + psi_zt dz_t/dr,   psi_z = psi_zt dz_t/dz.
    // Optionally returns the mesh point and |det d(r,z)/d(ref)|.
    void EvalVdof(double xi, double eta, double * nv, double (*dnv)[2],
                  double * rz = nullptr, double * mesh_det = nullptr) const
    {
        TargetPoint tp;
        Target(xi, eta, tp);
        const double rp = tp.r;
        if (rz) { rz[0] = tp.r; rz[1] = tp.z; }
        if (mesh_det) *mesh_det = tp.det;
        double u = xi, v = eta;
        Locate(tp.st, tp.zt, u, v);
        PointData<B> pd;
        AtRef(u, v, pd);
        for (int k = 0; k < B::N; ++k) {
            double dn[2];
            // psi_r / (2 r) and the physical psi_z through the target map; the
            // evaluation radius is the mesh radius (pd.s carries the Newton
            // tolerance, which dominates s near the axis)
            const double ps = pd.gs[k] * tp.sr + pd.gz[k] * tp.zr;
            const double pz = pd.gs[k] * tp.sz + pd.gz[k] * tp.zz;
            VdofFromPsi(r[k], rp * rp, s_scale, pd.n[k], ps, pz, nv[k],
                        dnv ? dn : nullptr);
            if (dnv) { dnv[k][0] = dn[0]; dnv[k][1] = dn[1]; }
        }
    }

    // Physical quadrature on the mesh element through EvalVdof, used where the
    // target map is not (r^2, z) (an oblique axis edge or curved edges):
    // body(n, dn, r, w) with w including |det| of the mesh map.
    template <typename F>
    void PhysicalQuadrature(F && body) const
    {
        // The target map makes the integrand non-polynomial in the reference
        // coordinates (the s-matched z shift, the blended edge ratio): the
        // 16-point rule (8 points leave ~1e-6 relative).
        auto at = [&](double xi, double eta, double w) {
            double n[B::N], dn[B::N][2], rz[2], det;
            EvalVdof(xi, eta, n, dn, rz, &det);
            body(n, dn, rz[0], w * det);
        };
        Quadrature<16>(at);
    }

    // Reference quadrature: body(xi, eta, w).
    template <int NG = 8, typename F>
    static void Quadrature(F && body)
    {
        static_assert(NG == 8 || NG == 16, "Gauss rule with 8 or 16 points");
        static constexpr double x8[8] = {
            0.01985507175123188, 0.10166676129318664, 0.23723379504183551,
            0.40828267875217508, 0.59171732124782492, 0.76276620495816449,
            0.89833323870681336, 0.98014492824876812 };
        static constexpr double w8[8] = {
            0.05061426814518813, 0.11119051722668724, 0.15685332293894364,
            0.18134189168918100, 0.18134189168918100, 0.15685332293894364,
            0.11119051722668724, 0.05061426814518813 };
        static constexpr double x16[16] = {
            0.005299532504175031, 0.0277124884633837, 0.06718439880608412,
            0.1222977958224985, 0.19106187779867811, 0.2709916111713863,
            0.35919822461037054, 0.4524937450811813, 0.5475062549188188,
            0.6408017753896295, 0.7290083888286136, 0.8089381222013219,
            0.8777022041775016, 0.9328156011939159, 0.9722875115366163,
            0.994700467495825 };
        static constexpr double w16[16] = {
            0.013576229705877088, 0.031126761969323728, 0.0475792558412463,
            0.062314485627767036, 0.07479799440828835, 0.08457825969750132,
            0.09130170752246182, 0.09472530522753432, 0.09472530522753432,
            0.09130170752246182, 0.08457825969750132, 0.07479799440828835,
            0.062314485627767036, 0.0475792558412463, 0.031126761969323728,
            0.013576229705877088 };
        const double * x = NG == 8 ? x8 : x16;
        const double * w = NG == 8 ? w8 : w16;
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
    //
    // Where the target map is (r^2, z) both forms are the same integral; an
    // element with an oblique axis edge is integrated in (r, z) through EvalVdof
    // instead, so the matrix is the form of the evaluated basis.
    void Stiffness(double mu, FlatMatrix<double> elmat) const
    {
        const double two_pi = 6.283185307179586476925;
        if (PhysicalTarget()) {
            double K[B::N][B::N] = {};
            PhysicalQuadrature([&](const double * n, const double (*dn)[2], double rp, double w) {
                double cr[B::N];
                for (int i = 0; i < B::N; ++i) cr[i] = dn[i][0] + n[i] / rp;
                for (int i = 0; i < B::N; ++i)
                    for (int j = 0; j < B::N; ++j)
                        K[i][j] += w * rp * (dn[i][1] * dn[j][1] + cr[i] * cr[j]);
            });
            for (int i = 0; i < B::N; ++i)
                for (int j = 0; j < B::N; ++j)
                    elmat(i, j) = two_pi / mu * 0.5 * (K[i][j] + K[j][i]);
            return;
        }
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
        for (int i = 0; i < B::N; ++i)
            for (int j = 0; j < B::N; ++j)
                elmat(i, j) = two_pi / mu * r[i] * r[j] * 0.5 * (K[i][j] + K[j][i]);
    }

    // M_ij = 2 pi sigma INT N_i N_j r dr dz
    //      = 2 pi sigma r_i r_j INT psi_i psi_j / (2 s) ds dz
    void SigmaMass(double sigma, FlatMatrix<double> elmat) const
    {
        if (PhysicalTarget()) {
            double M[B::N][B::N] = {};
            PhysicalQuadrature([&](const double * n, const double (*)[2], double rp, double w) {
                for (int i = 0; i < B::N; ++i)
                    for (int j = 0; j < B::N; ++j)
                        M[i][j] += w * rp * n[i] * n[j];
            });
            const double two_pi_sigma = 6.283185307179586476925 * sigma;
            for (int i = 0; i < B::N; ++i)
                for (int j = 0; j < B::N; ++j)
                    elmat(i, j) = two_pi_sigma * 0.5 * (M[i][j] + M[j][i]);
            return;
        }
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
