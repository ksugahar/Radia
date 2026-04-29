/*-------------------------------------------------------------------------
*
* File name:      rad_bem_galerkin.cpp
*
* Project:        RADIA
*
* Description:    Implementation of fast C++ Galerkin BEM assembler.
*                 See rad_bem_galerkin.h for API.
*
*                 Internal organization:
*                   - Triangle Gauss-Legendre rules (Stroud + Duffy)
*                   - Sauter-Schwab Duffy 4D rules (CV 2 / CE 5 / Id 6)
*                   - P2 Lagrange geometry evaluator
*                   - Per-pair assembly: regular / vertex / edge / identical
*                   - OpenMP parallel outer loop
*
-------------------------------------------------------------------------*/

#include "rad_bem_galerkin.h"

#include <vector>
#include <array>
#include <cmath>
#include <cstdint>
#include <stdexcept>
#include <algorithm>

#ifdef _OPENMP
#include <omp.h>
#endif

namespace radia { namespace bem {

namespace {  // file-local

constexpr double INV_4PI = 1.0 / (4.0 * 3.14159265358979323846);
constexpr double PI      = 3.14159265358979323846;

// ----------------------------------------------------------------------
// Triangle quadrature: Stroud 7-pt (deg 5) + 13-pt (deg 7) + Duffy GL n*n
// Each rule stores (xi, eta, w) with sum(w) = 1/2 (= area of T_ref).
// ----------------------------------------------------------------------
struct TriQuad {
    std::vector<double> xi;
    std::vector<double> eta;
    std::vector<double> w;
    int n() const { return static_cast<int>(xi.size()); }
};

TriQuad make_stroud_7() {
    TriQuad q;
    auto add = [&](double x, double e, double ww) {
        q.xi.push_back(x); q.eta.push_back(e); q.w.push_back(ww * 0.5);
    };
    add(1.0/3.0, 1.0/3.0, 0.225);
    add(0.0597158717, 0.4701420641, 0.1323941527);
    add(0.4701420641, 0.0597158717, 0.1323941527);
    add(0.4701420641, 0.4701420641, 0.1323941527);
    add(0.7974269853, 0.1012865073, 0.1259391805);
    add(0.1012865073, 0.7974269853, 0.1259391805);
    add(0.1012865073, 0.1012865073, 0.1259391805);
    return q;
}

TriQuad make_stroud_13() {
    TriQuad q;
    auto add = [&](double x, double e, double ww) {
        q.xi.push_back(x); q.eta.push_back(e); q.w.push_back(ww * 0.5);
    };
    add(0.333333333333333, 0.333333333333333, -0.149570044467670);
    add(0.260345966079038, 0.260345966079038,  0.175615257433204);
    add(0.260345966079038, 0.479308067841923,  0.175615257433204);
    add(0.479308067841923, 0.260345966079038,  0.175615257433204);
    add(0.065130102902216, 0.065130102902216,  0.053347235608839);
    add(0.065130102902216, 0.869739794195568,  0.053347235608839);
    add(0.869739794195568, 0.065130102902216,  0.053347235608839);
    add(0.312865496004875, 0.638444188569809,  0.077113760890257);
    add(0.638444188569809, 0.048690315425316,  0.077113760890257);
    add(0.048690315425316, 0.312865496004875,  0.077113760890257);
    add(0.312865496004875, 0.048690315425316,  0.077113760890257);
    add(0.638444188569809, 0.312865496004875,  0.077113760890257);
    add(0.048690315425316, 0.638444188569809,  0.077113760890257);
    return q;
}

// 1D Gauss-Legendre on [0, 1].  Nodes/weights for n = 1..16.
// Pre-tabulated up to n=12, computed via Newton-on-Legendre for any n.
// We embed values for n up to 16 directly.
void gl1d_on_unit_interval(int n, std::vector<double>& x,
                            std::vector<double>& w) {
    // Standard GL on [-1, 1], then map to [0, 1].
    // Here we hard-code up to n = 12 then bail out (for higher n we
    // could add a Golub-Welsch routine, but the BEM use case rarely
    // exceeds n_q = 10).
    static const std::vector<std::vector<double>> nodes_neg11 = {
        {0.0},
        {-0.5773502691896257, 0.5773502691896257},
        {-0.7745966692414834, 0.0, 0.7745966692414834},
        {-0.8611363115940526, -0.3399810435848563, 0.3399810435848563, 0.8611363115940526},
        {-0.9061798459386640, -0.5384693101056831, 0.0, 0.5384693101056831, 0.9061798459386640},
        {-0.9324695142031521, -0.6612093864662645, -0.2386191860831969,
          0.2386191860831969,  0.6612093864662645,  0.9324695142031521},
        {-0.9491079123427585, -0.7415311855993945, -0.4058451513773972, 0.0,
          0.4058451513773972,  0.7415311855993945,  0.9491079123427585},
        {-0.9602898564975363, -0.7966664774136267, -0.5255324099163290, -0.1834346424956498,
          0.1834346424956498,  0.5255324099163290,  0.7966664774136267,  0.9602898564975363},
        {-0.9681602395076261, -0.8360311073266358, -0.6133714327005904, -0.3242534234038089, 0.0,
          0.3242534234038089,  0.6133714327005904,  0.8360311073266358,  0.9681602395076261},
        {-0.9739065285171717, -0.8650633666889845, -0.6794095682990244, -0.4333953941292472,
         -0.1488743389816312,  0.1488743389816312,  0.4333953941292472,  0.6794095682990244,
          0.8650633666889845,  0.9739065285171717},
        {-0.9782286581460570, -0.8870625997680953, -0.7301520055740494, -0.5190961292068118,
         -0.2695431559523450,  0.0,                0.2695431559523450,  0.5190961292068118,
          0.7301520055740494,  0.8870625997680953,  0.9782286581460570},
        {-0.9815606342467192, -0.9041172563704749, -0.7699026741943047, -0.5873179542866175,
         -0.3678314989981802, -0.1252334085114689,  0.1252334085114689,  0.3678314989981802,
          0.5873179542866175,  0.7699026741943047,  0.9041172563704749,  0.9815606342467192},
    };
    static const std::vector<std::vector<double>> wts_neg11 = {
        {2.0},
        {1.0, 1.0},
        {0.5555555555555556, 0.8888888888888888, 0.5555555555555556},
        {0.3478548451374538, 0.6521451548625461, 0.6521451548625461, 0.3478548451374538},
        {0.2369268850561891, 0.4786286704993665, 0.5688888888888889, 0.4786286704993665,
         0.2369268850561891},
        {0.1713244923791704, 0.3607615730481386, 0.4679139345726910, 0.4679139345726910,
         0.3607615730481386, 0.1713244923791704},
        {0.1294849661688697, 0.2797053914892767, 0.3818300505051189, 0.4179591836734694,
         0.3818300505051189, 0.2797053914892767, 0.1294849661688697},
        {0.1012285362903763, 0.2223810344533745, 0.3137066458778873, 0.3626837833783620,
         0.3626837833783620, 0.3137066458778873, 0.2223810344533745, 0.1012285362903763},
        {0.0812743883615744, 0.1806481606948574, 0.2606106964029354, 0.3123470770400028,
         0.3302393550012598, 0.3123470770400028, 0.2606106964029354, 0.1806481606948574,
         0.0812743883615744},
        {0.0666713443086881, 0.1494513491505806, 0.2190863625159820, 0.2692667193099963,
         0.2955242247147529, 0.2955242247147529, 0.2692667193099963, 0.2190863625159820,
         0.1494513491505806, 0.0666713443086881},
        {0.0556685671161737, 0.1255803694649046, 0.1862902109277343, 0.2331937645919905,
         0.2628045445102467, 0.2729250867779006, 0.2628045445102467, 0.2331937645919905,
         0.1862902109277343, 0.1255803694649046, 0.0556685671161737},
        {0.0471753363865118, 0.1069393259953184, 0.1600783285433462, 0.2031674267230659,
         0.2334925365383548, 0.2491470458134028, 0.2491470458134028, 0.2334925365383548,
         0.2031674267230659, 0.1600783285433462, 0.1069393259953184, 0.0471753363865118},
    };
    if (n < 1 || n > (int)nodes_neg11.size()) {
        throw std::invalid_argument("gl1d_on_unit_interval: n out of range");
    }
    const auto& nn = nodes_neg11[n - 1];
    const auto& ww = wts_neg11[n - 1];
    x.resize(n); w.resize(n);
    for (int i = 0; i < n; ++i) {
        x[i] = 0.5 * (nn[i] + 1.0);
        w[i] = 0.5 * ww[i];
    }
}

// Duffy-collapsed n*n triangle Gauss-Legendre on T_ref.
// Mapping: (xi, eta) = (s*(1-t), s*t), Jacobian = s.
// Sum of weights = 1/2 (area of T_ref).  Exact for total degree 2n-1.
TriQuad make_duffy_n(int n) {
    std::vector<double> sx, sw, tx, tw;
    gl1d_on_unit_interval(n, sx, sw);
    gl1d_on_unit_interval(n, tx, tw);
    TriQuad q;
    q.xi.reserve(n * n); q.eta.reserve(n * n); q.w.reserve(n * n);
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            double s = sx[i], t = tx[j];
            q.xi.push_back(s * (1.0 - t));
            q.eta.push_back(s * t);
            q.w.push_back(s * sw[i] * tw[j]);
        }
    }
    return q;
}

TriQuad pick_tri_quad(int degree) {
    if (degree <= 5)  return make_stroud_7();
    if (degree <= 7)  return make_stroud_13();
    int n = (degree + 2) / 2;  // ceil((degree+1)/2)
    return make_duffy_n(n);
}

// ----------------------------------------------------------------------
// P2 Lagrange evaluator -- positions and surface element.
// p2_nodes ordering: [v0, v1, v2, mid01, mid12, mid20]
// Reference: corner 0 at (0,0), corner 1 at (1,0), corner 2 at (0,1).
// ----------------------------------------------------------------------
inline void p2_basis(double xi, double eta, double N[6]) {
    double L0 = 1.0 - xi - eta, L1 = xi, L2 = eta;
    N[0] = L0 * (2.0 * L0 - 1.0);
    N[1] = L1 * (2.0 * L1 - 1.0);
    N[2] = L2 * (2.0 * L2 - 1.0);
    N[3] = 4.0 * L0 * L1;
    N[4] = 4.0 * L1 * L2;
    N[5] = 4.0 * L0 * L2;
}

inline void p2_basis_grad(double xi, double eta,
                          double dN_dxi[6], double dN_deta[6]) {
    double L0 = 1.0 - xi - eta, L1 = xi, L2 = eta;
    dN_dxi[0]  = -(4.0 * L0 - 1.0);
    dN_deta[0] = -(4.0 * L0 - 1.0);
    dN_dxi[1]  =  (4.0 * L1 - 1.0);
    dN_deta[1] = 0.0;
    dN_dxi[2]  = 0.0;
    dN_deta[2] =  (4.0 * L2 - 1.0);
    dN_dxi[3]  = 4.0 * (L0 - L1);
    dN_deta[3] = -4.0 * L1;
    dN_dxi[4]  = 4.0 * L2;
    dN_deta[4] = 4.0 * L1;
    dN_dxi[5]  = -4.0 * L2;
    dN_deta[5] = 4.0 * (L0 - L2);
}

/**
 * Evaluate physical position r(xi, eta), surface area Jacobian
 * J = |dr/dxi x dr/deta|, and unit outward normal at one (xi, eta)
 * on a P2-curved triangle.
 *
 * @param p2  6-row, 3-col matrix of P2 nodes (flat ptr, row-major):
 *            p2[k*3 + d] = d-th coord of k-th P2 node.
 */
inline void p2_geom(const double* p2, double xi, double eta,
                    double r[3], double* J_out, double n_hat[3]) {
    double N[6], dNxi[6], dNeta[6];
    p2_basis(xi, eta, N);
    p2_basis_grad(xi, eta, dNxi, dNeta);
    double r_local[3] = {0, 0, 0};
    double dxi[3]  = {0, 0, 0};
    double deta[3] = {0, 0, 0};
    for (int k = 0; k < 6; ++k) {
        for (int d = 0; d < 3; ++d) {
            double pk = p2[k * 3 + d];
            r_local[d] += N[k]    * pk;
            dxi[d]     += dNxi[k] * pk;
            deta[d]    += dNeta[k] * pk;
        }
    }
    r[0] = r_local[0]; r[1] = r_local[1]; r[2] = r_local[2];
    // Cross product dxi x deta
    double cx = dxi[1] * deta[2] - dxi[2] * deta[1];
    double cy = dxi[2] * deta[0] - dxi[0] * deta[2];
    double cz = dxi[0] * deta[1] - dxi[1] * deta[0];
    double J = std::sqrt(cx*cx + cy*cy + cz*cz);
    *J_out = J;
    if (n_hat) {
        double inv = (J > 0) ? 1.0 / J : 0.0;
        n_hat[0] = cx * inv;
        n_hat[1] = cy * inv;
        n_hat[2] = cz * inv;
    }
}

// ----------------------------------------------------------------------
// Sauter-Schwab Duffy 4D quadrature -- direct port of NGSolve.bem
// intrules_SauterSchwab.cpp.  Each rule produces (xi_a, eta_a, xi_b,
// eta_b, w) per quadrature point, in the canonical T_ref^2 frame
// (NGSolve T_ref = [(0,0),(1,0),(0,1)], shared singular at (0,0) for
// vertex / [(0,0),(1,0)] edge for edge / whole tri for identical).
//
// Caller is responsible for permuting element corners to the canonical
// form and applying parity sign correction for normal-dependent kernels
// (see DL pair handler).
// ----------------------------------------------------------------------
struct SSQuad {
    std::vector<double> xi_a, eta_a, xi_b, eta_b, w;
    int n() const { return static_cast<int>(w.size()); }
};

void build_ss_common_vertex(int n_q, SSQuad& out) {
    std::vector<double> g, gw;
    gl1d_on_unit_interval(n_q, g, gw);
    out.xi_a.clear(); out.eta_a.clear();
    out.xi_b.clear(); out.eta_b.clear(); out.w.clear();
    int total = 2 * n_q * n_q * n_q * n_q;
    out.xi_a.reserve(total); out.eta_a.reserve(total);
    out.xi_b.reserve(total); out.eta_b.reserve(total);
    out.w.reserve(total);

    auto push = [&](double ip0, double ip1, double ip2, double ip3, double ww) {
        out.xi_a.push_back(ip0 - ip1);
        out.eta_a.push_back(ip1);
        out.xi_b.push_back(ip2 - ip3);
        out.eta_b.push_back(ip3);
        out.w.push_back(ww);
    };

    // 2 sub-cubes:
    //   1: ip = xi*(1, e1, e2, e2*e3),   J = xi^3 * e2
    //   2: ip = xi*(e2, e2*e3, 1, e1)
    for (int ie1 = 0; ie1 < n_q; ++ie1) {
        double e1 = g[ie1], we1 = gw[ie1];
        for (int ie2 = 0; ie2 < n_q; ++ie2) {
            double e2 = g[ie2], we2 = gw[ie2];
            for (int ie3 = 0; ie3 < n_q; ++ie3) {
                double e3 = g[ie3], we3 = gw[ie3];
                for (int ix = 0; ix < n_q; ++ix) {
                    double xi = g[ix], wxi = gw[ix];
                    double w_4d = we1 * we2 * we3 * wxi;
                    double J_per = xi * xi * xi * e2 * w_4d;
                    push(xi, xi*e1, xi*e2, xi*e2*e3, J_per);
                    push(xi*e2, xi*e2*e3, xi, xi*e1, J_per);
                }
            }
        }
    }
}

void build_ss_common_edge(int n_q, SSQuad& out) {
    std::vector<double> g, gw;
    gl1d_on_unit_interval(n_q, g, gw);
    out.xi_a.clear(); out.eta_a.clear();
    out.xi_b.clear(); out.eta_b.clear(); out.w.clear();
    int total = 5 * n_q * n_q * n_q * n_q;
    out.xi_a.reserve(total); out.eta_a.reserve(total);
    out.xi_b.reserve(total); out.eta_b.reserve(total);
    out.w.reserve(total);

    auto push = [&](double ip0, double ip1, double ip2, double ip3, double ww) {
        out.xi_a.push_back(ip0 - ip1);
        out.eta_a.push_back(ip1);
        out.xi_b.push_back(ip2 - ip3);
        out.eta_b.push_back(ip3);
        out.w.push_back(ww);
    };

    for (int ie1 = 0; ie1 < n_q; ++ie1) {
        double e1 = g[ie1], we1 = gw[ie1];
        for (int ie2 = 0; ie2 < n_q; ++ie2) {
            double e2 = g[ie2], we2 = gw[ie2];
            for (int ie3 = 0; ie3 < n_q; ++ie3) {
                double e3 = g[ie3], we3 = gw[ie3];
                for (int ix = 0; ix < n_q; ++ix) {
                    double xi = g[ix], wxi = gw[ix];
                    double w_4d = we1 * we2 * we3 * wxi;
                    double J_a = xi * xi * xi * e1 * e1 * w_4d;          // sub 0
                    double J_b = xi * xi * xi * e1 * e1 * e2 * w_4d;     // subs 1..4
                    push(xi*1.0,             xi*e1*e3,            xi*(1-e1*e2),     xi*e1*(1-e2),       J_a);
                    push(xi*1.0,             xi*e1,                xi*(1-e1*e2*e3),  xi*e1*e2*(1-e3),    J_b);
                    push(xi*(1-e1*e2),       xi*e1*(1-e2),         xi*1.0,           xi*e1*e2*e3,        J_b);
                    push(xi*(1-e1*e2*e3),    xi*e1*e2*(1-e3),      xi*1.0,           xi*e1,              J_b);
                    push(xi*(1-e1*e2*e3),    xi*e1*(1-e2*e3),      xi*1.0,           xi*e1*e2,           J_b);
                }
            }
        }
    }
}

void build_ss_identical(int n_q, SSQuad& out) {
    std::vector<double> g, gw;
    gl1d_on_unit_interval(n_q, g, gw);
    out.xi_a.clear(); out.eta_a.clear();
    out.xi_b.clear(); out.eta_b.clear(); out.w.clear();
    int total = 6 * n_q * n_q * n_q * n_q;
    out.xi_a.reserve(total); out.eta_a.reserve(total);
    out.xi_b.reserve(total); out.eta_b.reserve(total);
    out.w.reserve(total);

    auto push = [&](double ip0, double ip1, double ip2, double ip3, double ww) {
        out.xi_a.push_back(ip0 - ip1);
        out.eta_a.push_back(ip1);
        out.xi_b.push_back(ip2 - ip3);
        out.eta_b.push_back(ip3);
        out.w.push_back(ww);
    };

    for (int ie1 = 0; ie1 < n_q; ++ie1) {
        double e1 = g[ie1], we1 = gw[ie1];
        for (int ie2 = 0; ie2 < n_q; ++ie2) {
            double e2 = g[ie2], we2 = gw[ie2];
            for (int ie3 = 0; ie3 < n_q; ++ie3) {
                double e3 = g[ie3], we3 = gw[ie3];
                for (int ix = 0; ix < n_q; ++ix) {
                    double xi = g[ix], wxi = gw[ix];
                    double w_4d = we1 * we2 * we3 * wxi;
                    double J = xi * xi * xi * e1 * e1 * e2 * w_4d;
                    push(xi*1.0,            xi*(1-e1+e1*e2),     xi*(1-e1*e2*e3),  xi*(1-e1),            J);
                    push(xi*(1-e1*e2*e3),   xi*(1-e1),           xi*1.0,           xi*(1-e1+e1*e2),      J);
                    push(xi*1.0,            xi*e1*(1-e2+e2*e3),  xi*(1-e1*e2),     xi*e1*(1-e2),         J);
                    push(xi*(1-e1*e2),      xi*e1*(1-e2),        xi*1.0,           xi*e1*(1-e2+e2*e3),   J);
                    push(xi*(1-e1*e2*e3),   xi*e1*(1-e2*e3),     xi*1.0,           xi*e1*(1-e2),         J);
                    push(xi*1.0,            xi*e1*(1-e2),        xi*(1-e1*e2*e3),  xi*e1*(1-e2*e3),      J);
                }
            }
        }
    }
}

// ----------------------------------------------------------------------
// Permutation parity (cyclic vs non-cyclic perm of (0,1,2))
// ----------------------------------------------------------------------
inline int perm_parity_sign(int p0, int p1, int p2) {
    if ((p0==0 && p1==1 && p2==2) || (p0==1 && p1==2 && p2==0) ||
        (p0==2 && p1==0 && p2==1)) return +1;
    return -1;
}

// Re-order p2_nodes per a permutation of corners (and matching mid-edges).
// p2_in: 6 nodes (each 3 doubles).  perm is the new corner order.
// Mid-edge order in p2: [3=mid01, 4=mid12, 5=mid20] (between corners (0,1), (1,2), (2,0)).
// After permutation, new mid-edge k (between new corners k and (k+1)%3)
// = old mid-edge between perm[k] and perm[(k+1)%3].
inline void permute_p2(const double* p2_in, int p0, int p1, int p2,
                       double p2_out[18]) {
    int corners[3] = {p0, p1, p2};
    for (int k = 0; k < 3; ++k) {
        int src = corners[k];
        for (int d = 0; d < 3; ++d) p2_out[k*3 + d] = p2_in[src*3 + d];
    }
    // mid-edge between corners (a, b) (a < b) maps to:
    //   (0, 1) -> 3
    //   (1, 2) -> 4
    //   (0, 2) -> 5
    auto old_mid_idx = [](int a, int b) -> int {
        int lo = std::min(a, b), hi = std::max(a, b);
        if (lo == 0 && hi == 1) return 3;
        if (lo == 1 && hi == 2) return 4;
        if (lo == 0 && hi == 2) return 5;
        return -1;
    };
    for (int k = 0; k < 3; ++k) {
        int a = corners[k], b = corners[(k + 1) % 3];
        int src = old_mid_idx(a, b);
        for (int d = 0; d < 3; ++d) p2_out[(3 + k)*3 + d] = p2_in[src*3 + d];
    }
}

// Given two tris (by their global vertex IDs), return:
//   n_shared: number of shared vertices (0=regular, 1=vertex, 2=edge, 3=identical)
//   shared_pairs[i] = (la, lb) for i = 0..n_shared-1
inline int find_shared(const int64_t* tri_a, const int64_t* tri_b,
                        int la_arr[3], int lb_arr[3]) {
    int n = 0;
    for (int la = 0; la < 3; ++la) {
        for (int lb = 0; lb < 3; ++lb) {
            if (tri_a[la] == tri_b[lb]) {
                la_arr[n] = la; lb_arr[n] = lb;
                ++n;
            }
        }
    }
    return n;
}

// Per-pair Galerkin block accumulators.
// SL block (3x3) and DL block (3x3) for a pair (T_a, T_b).
// Each entry block[i_a, j_b] = int int L_i^a kernel L_j^b dS dS'.
struct PairBlock {
    double SL[9];
    double DL[9];
};

// Evaluate SL/DL at one quadrature point and accumulate 3x3 contributions.
// Inputs:
//   xi_a, eta_a: ref coords on T_a
//   xi_b, eta_b: ref coords on T_b
//   weight: SS or regular w_q (no kernel applied yet, no Jacobian -- those
//           are computed inside via p2_geom)
//   p2a, p2b: 6-node P2 data (18 doubles)
//   sl_out, dl_out: 3x3 accumulators (preallocated, length 9)
inline void accumulate_pair_quadpt(double xi_a, double eta_a,
                                   double xi_b, double eta_b,
                                   double weight,
                                   const double* p2a, const double* p2b,
                                   double sl_out[9], double dl_out[9]) {
    double r_a[3], r_b[3], n_hat_b[3];
    double J_a, J_b;
    p2_geom(p2a, xi_a, eta_a, r_a, &J_a, nullptr);
    p2_geom(p2b, xi_b, eta_b, r_b, &J_b, n_hat_b);

    double dx = r_a[0] - r_b[0];
    double dy = r_a[1] - r_b[1];
    double dz = r_a[2] - r_b[2];
    double r2 = dx*dx + dy*dy + dz*dz;
    if (r2 < 1e-300) r2 = 1e-300;
    double r = std::sqrt(r2);

    double L0_a = 1.0 - xi_a - eta_a;
    double L1_a = xi_a;
    double L2_a = eta_a;
    double L0_b = 1.0 - xi_b - eta_b;
    double L1_b = xi_b;
    double L2_b = eta_b;
    double La[3] = {L0_a, L1_a, L2_a};
    double Lb[3] = {L0_b, L1_b, L2_b};

    double JaJb = J_a * J_b;
    double w_sl = weight * JaJb / (4.0 * PI * r);

    double dot_n = dx * n_hat_b[0] + dy * n_hat_b[1] + dz * n_hat_b[2];
    double w_dl = weight * JaJb * dot_n / (4.0 * PI * r2 * r);

    for (int i = 0; i < 3; ++i) {
        double LiA_sl = w_sl * La[i];
        double LiA_dl = w_dl * La[i];
        for (int j = 0; j < 3; ++j) {
            sl_out[i * 3 + j] += LiA_sl * Lb[j];
            dl_out[i * 3 + j] += LiA_dl * Lb[j];
        }
    }
}

}  // namespace (file-local)


// =====================================================================
// Public entry point
// =====================================================================
void AssembleSLDL(
    const double* verts, int n_v,
    const int64_t* tris, int n_t,
    const double* p2_nodes,
    int regular_quad_degree,
    int singular_n_q,
    int n_threads,
    double* SL_out,
    double* DL_out)
{
    (void)verts;  // unused -- positions come from p2_nodes corners

    // Pre-build all quadrature rules
    TriQuad q_reg = pick_tri_quad(regular_quad_degree);
    SSQuad ss_cv, ss_ce, ss_id;
    build_ss_common_vertex(singular_n_q, ss_cv);
    build_ss_common_edge(singular_n_q,   ss_ce);
    build_ss_identical(singular_n_q,     ss_id);

    // Pre-cache regular-pair physical positions and weights per tri.
    // For each tri t: arrays of (n_qreg, 3) positions, (n_qreg) weights*Jac,
    //                 (n_qreg, 3) hat values.
    int nq_reg = q_reg.n();
    std::vector<double> reg_pts(static_cast<size_t>(n_t) * nq_reg * 3);
    std::vector<double> reg_w(static_cast<size_t>(n_t) * nq_reg);
    std::vector<double> reg_n(static_cast<size_t>(n_t) * nq_reg * 3);
    std::vector<double> reg_L(static_cast<size_t>(nq_reg) * 3);  // shared L_ref values
    for (int q = 0; q < nq_reg; ++q) {
        reg_L[q * 3 + 0] = 1.0 - q_reg.xi[q] - q_reg.eta[q];
        reg_L[q * 3 + 1] = q_reg.xi[q];
        reg_L[q * 3 + 2] = q_reg.eta[q];
    }

#ifdef _OPENMP
    if (n_threads > 0) omp_set_num_threads(n_threads);
#endif

#pragma omp parallel for schedule(dynamic, 8)
    for (int t = 0; t < n_t; ++t) {
        const double* p2 = &p2_nodes[t * 6 * 3];
        for (int q = 0; q < nq_reg; ++q) {
            double r[3], n_hat[3], J;
            p2_geom(p2, q_reg.xi[q], q_reg.eta[q], r, &J, n_hat);
            size_t base_pts = (static_cast<size_t>(t) * nq_reg + q) * 3;
            reg_pts[base_pts + 0] = r[0];
            reg_pts[base_pts + 1] = r[1];
            reg_pts[base_pts + 2] = r[2];
            reg_n[base_pts + 0] = n_hat[0];
            reg_n[base_pts + 1] = n_hat[1];
            reg_n[base_pts + 2] = n_hat[2];
            reg_w[static_cast<size_t>(t) * nq_reg + q] = q_reg.w[q] * J;
        }
    }

    // Zero output matrices
    std::fill(SL_out, SL_out + (size_t)n_v * n_v, 0.0);
    std::fill(DL_out, DL_out + (size_t)n_v * n_v, 0.0);

    // To allow lock-free parallel scatter, accumulate per-thread tiles.
    // Simpler: parallelize over a, accumulate into thread-local copies of
    // a SUBSET, then merge at end.  For low thread count and N up to ~5K,
    // the simple atomic-add path is fine.
    //
    // Choice: we use thread-local (n_v x n_v) buffers and reduce.  Memory:
    // n_threads * 8 * n_v * 2 = 16 * n_v * n_threads bytes.  For n_v=2477,
    // 8 threads, 2 mats: 2 * 2477^2 * 8 * 8 = 784 MB.  Too much.
    //
    // Alternative: serialize the scatter via #pragma omp atomic.  Bad
    // for cache.
    //
    // Best: accumulate into an "a-row" buffer of size 3*n_v per thread,
    // flush after each `a`.  Each row of (a) updates only 9 entries per
    // (a, b) pair = 9*n_t entries per a-row, all in rows {tris[a, 0..2]}.
    // We sequentially process each `a` (parallelize over `b` inside if
    // desired, or sequential for cache locality).
    //
    // For simplicity in this initial implementation, parallelize the
    // outer (a, b) loop and use atomic adds.  n_v ~ 2500 means n_v^2 ~
    // 6.25M doubles = 50 MB per matrix, fits easily; the atomic add
    // overhead is small relative to computation.

#pragma omp parallel for schedule(dynamic, 4)
    for (int a = 0; a < n_t; ++a) {
        const int64_t* Va = &tris[a * 3];
        const double* p2a = &p2_nodes[a * 6 * 3];
        const double* pts_a = &reg_pts[(size_t)a * nq_reg * 3];
        const double* w_a   = &reg_w[(size_t)a * nq_reg];

        for (int b = 0; b < n_t; ++b) {
            const int64_t* Vb = &tris[b * 3];
            const double* p2b = &p2_nodes[b * 6 * 3];

            int la_arr[3], lb_arr[3];
            int n_sh = find_shared(Va, Vb, la_arr, lb_arr);

            double SL_block[9] = {0,0,0,0,0,0,0,0,0};
            double DL_block[9] = {0,0,0,0,0,0,0,0,0};
            int parity_sign_b = +1;  // for DL only

            if (n_sh == 0) {
                // Regular: tensor-product Gauss
                const double* pts_b = &reg_pts[(size_t)b * nq_reg * 3];
                const double* w_b   = &reg_w[(size_t)b * nq_reg];
                const double* n_b   = &reg_n[(size_t)b * nq_reg * 3];
                for (int qa = 0; qa < nq_reg; ++qa) {
                    double xa[3] = {pts_a[qa*3], pts_a[qa*3+1], pts_a[qa*3+2]};
                    double La0 = reg_L[qa*3], La1 = reg_L[qa*3+1], La2 = reg_L[qa*3+2];
                    double wa = w_a[qa];
                    for (int qb = 0; qb < nq_reg; ++qb) {
                        double xb[3] = {pts_b[qb*3], pts_b[qb*3+1], pts_b[qb*3+2]};
                        double Lb0 = reg_L[qb*3], Lb1 = reg_L[qb*3+1], Lb2 = reg_L[qb*3+2];
                        double wb = w_b[qb];
                        double dx = xa[0] - xb[0];
                        double dy = xa[1] - xb[1];
                        double dz = xa[2] - xb[2];
                        double r2 = dx*dx + dy*dy + dz*dz;
                        double r = std::sqrt(r2);
                        double K_sl = INV_4PI / r;
                        double dot_n = dx * n_b[qb*3] + dy * n_b[qb*3+1] + dz * n_b[qb*3+2];
                        double K_dl = dot_n * INV_4PI / (r2 * r);
                        double w_pair = wa * wb;
                        double w_sl = K_sl * w_pair;
                        double w_dl = K_dl * w_pair;
                        // 3x3 outer product
                        double La[3] = {La0, La1, La2};
                        double Lb[3] = {Lb0, Lb1, Lb2};
                        for (int i = 0; i < 3; ++i) {
                            double Li_sl = w_sl * La[i];
                            double Li_dl = w_dl * La[i];
                            for (int j = 0; j < 3; ++j) {
                                SL_block[i*3+j] += Li_sl * Lb[j];
                                DL_block[i*3+j] += Li_dl * Lb[j];
                            }
                        }
                    }
                }
            } else {
                // Singular pair: pick SS rule + permutation
                int pa[3], pb[3];
                const SSQuad* qss = nullptr;
                int perm_la_to_orig[3], perm_lb_to_orig[3];

                if (n_sh == 1) {
                    int la = la_arr[0], lb = lb_arr[0];
                    pa[0] = la; pa[1] = (la + 1) % 3; pa[2] = (la + 2) % 3;
                    pb[0] = lb; pb[1] = (lb + 1) % 3; pb[2] = (lb + 2) % 3;
                    qss = &ss_cv;
                } else if (n_sh == 2) {
                    int la0 = la_arr[0], lb0 = lb_arr[0];
                    int la1 = la_arr[1], lb1 = lb_arr[1];
                    int apex_a = 0 + 1 + 2 - la0 - la1;
                    int apex_b = 0 + 1 + 2 - lb0 - lb1;
                    pa[0] = la0; pa[1] = la1; pa[2] = apex_a;
                    pb[0] = lb0; pb[1] = lb1; pb[2] = apex_b;
                    qss = &ss_ce;
                } else {  // n_sh == 3 identical
                    int perm_b_canon[3] = {-1, -1, -1};
                    for (int k = 0; k < 3; ++k) {
                        perm_b_canon[la_arr[k]] = lb_arr[k];
                    }
                    pa[0] = 0; pa[1] = 1; pa[2] = 2;
                    pb[0] = perm_b_canon[0]; pb[1] = perm_b_canon[1]; pb[2] = perm_b_canon[2];
                    qss = &ss_id;
                }

                for (int k = 0; k < 3; ++k) {
                    perm_la_to_orig[k] = pa[k];
                    perm_lb_to_orig[k] = pb[k];
                }
                parity_sign_b = perm_parity_sign(pb[0], pb[1], pb[2]);

                double p2a_perm[18], p2b_perm[18];
                permute_p2(p2a, pa[0], pa[1], pa[2], p2a_perm);
                permute_p2(p2b, pb[0], pb[1], pb[2], p2b_perm);

                int n_quad = qss->n();
                for (int q = 0; q < n_quad; ++q) {
                    accumulate_pair_quadpt(qss->xi_a[q], qss->eta_a[q],
                                            qss->xi_b[q], qss->eta_b[q],
                                            qss->w[q], p2a_perm, p2b_perm,
                                            SL_block, DL_block);
                }

                // DL gets parity sign correction (SL invariant under winding)
                if (parity_sign_b == -1) {
                    for (int k = 0; k < 9; ++k) DL_block[k] = -DL_block[k];
                }

                // Un-permute the 3x3 block:
                //  block_perm[i_perm, j_perm] -> SL_orig[pa[i_perm], pb[j_perm]]
                double SL_orig[9] = {0,0,0,0,0,0,0,0,0};
                double DL_orig[9] = {0,0,0,0,0,0,0,0,0};
                for (int i_p = 0; i_p < 3; ++i_p) {
                    for (int j_p = 0; j_p < 3; ++j_p) {
                        SL_orig[perm_la_to_orig[i_p] * 3 + perm_lb_to_orig[j_p]]
                            = SL_block[i_p * 3 + j_p];
                        DL_orig[perm_la_to_orig[i_p] * 3 + perm_lb_to_orig[j_p]]
                            = DL_block[i_p * 3 + j_p];
                    }
                }
                for (int k = 0; k < 9; ++k) {
                    SL_block[k] = SL_orig[k];
                    DL_block[k] = DL_orig[k];
                }
            }

            // Scatter into global SL_out, DL_out at rows Va[i], cols Vb[j]
            for (int i = 0; i < 3; ++i) {
                int gi = static_cast<int>(Va[i]);
                for (int j = 0; j < 3; ++j) {
                    int gj = static_cast<int>(Vb[j]);
                    double sv = SL_block[i*3 + j];
                    double dv = DL_block[i*3 + j];
                    if (sv != 0.0) {
#pragma omp atomic
                        SL_out[(size_t)gi * n_v + gj] += sv;
                    }
                    if (dv != 0.0) {
#pragma omp atomic
                        DL_out[(size_t)gi * n_v + gj] += dv;
                    }
                }
            }
        }  // b
    }  // a
}

}}  // namespace radia::bem
