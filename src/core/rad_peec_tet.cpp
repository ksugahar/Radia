// Tetrahedral volume PEEC (see rad_peec_tet.h).
#include "rad_peec_tet.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <cmath>
#include <cstdint>
#include <limits>
#include <stdexcept>
#include <string>
#include <unordered_map>

#include "rad_host_parallel.h"

#ifdef HAVE_LAPACK
#include "rad_blas.h"
#if defined(RADIA_BLAS_MKL)
#include <mkl_lapack.h>
#define RADIA_TET_GETRF zgetrf_
#define RADIA_TET_GETRS zgetrs_
#elif defined(RADIA_BLAS_NVPL)
#include <nvpl_lapack.h>
#define RADIA_TET_GETRF NVPL_LAPACK_zgetrf
#define RADIA_TET_GETRS NVPL_LAPACK_zgetrs
#endif
#endif

namespace radia {
namespace peec_tet {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kMu0Over4Pi = 1e-7;

using V3 = std::array<double, 3>;
inline V3 sub(const V3& a, const V3& b) { return {a[0] - b[0], a[1] - b[1], a[2] - b[2]}; }
inline V3 add(const V3& a, const V3& b) { return {a[0] + b[0], a[1] + b[1], a[2] + b[2]}; }
inline V3 mul(double s, const V3& a) { return {s * a[0], s * a[1], s * a[2]}; }
inline double dot(const V3& a, const V3& b) { return a[0] * b[0] + a[1] * b[1] + a[2] * b[2]; }
inline V3 cross(const V3& a, const V3& b) {
    return {a[1] * b[2] - a[2] * b[1], a[2] * b[0] - a[0] * b[2], a[0] * b[1] - a[1] * b[0]};
}
inline double norm(const V3& a) { return std::sqrt(dot(a, a)); }

// int_T 1/R dS and int_T R dS over a flat triangle (a, b, c), counterclockwise about the unit
// normal n, at the point x (Wilton et al. 1984; Graglia 1993 recursion for R).
void FaceIntegrals(const V3& x, const V3& a, const V3& b, const V3& c, const V3& n,
                   double& im1, double& ip1) {
    const double h = dot(sub(x, a), n);
    const double ah = std::fabs(h);
    const V3 rho = sub(x, mul(h, n));
    const V3 p[3] = {a, b, c};
    double sum_tl = 0.0, sum_beta = 0.0, sum_te = 0.0;
    for (int i = 0; i < 3; ++i) {
        const V3& pm = p[i];
        const V3& pp = p[(i + 1) % 3];
        const V3 e = sub(pp, pm);
        const double len = norm(e);
        const V3 lh = mul(1.0 / len, e);
        const V3 u = cross(lh, n);
        const double t = dot(sub(pm, rho), u);
        const double lp = dot(sub(pp, rho), lh), lm = dot(sub(pm, rho), lh);
        const double Rp = norm(sub(x, pp)), Rm = norm(sub(x, pm));
        const double R02 = t * t + h * h;
        if (std::fabs(t) <= 1e-12 * len)
            continue;                                   // x on the edge line: every term carries t
        // ln((R+ + l+)/(R- + l-)); where l < 0, R + l = R0^2/(R - l) avoids the cancellation
        // for points close to the edge line (R0^2 = t^2 + h^2 > 0 here).
        const double up = lp >= 0.0 ? Rp + lp : R02 / (Rp - lp);
        const double um = lm >= 0.0 ? Rm + lm : R02 / (Rm - lm);
        const double Li = std::log(up / um);
        sum_tl += t * Li;
        sum_beta += std::atan(t * lp / (R02 + ah * Rp)) - std::atan(t * lm / (R02 + ah * Rm));
        sum_te += t * 0.5 * (lp * Rp - lm * Rm + R02 * Li);
    }
    im1 = sum_tl - ah * sum_beta;
    ip1 = (h * h * im1 + sum_te) / 3.0;
}

struct Tet {
    V3 v[4];                      // vertices relative to the centroid (local frame)
    V3 centroid;                  // global centroid
    double volume = 0.0;
    double radius = 0.0;          // bounding radius about the centroid
    V3 n[4];                      // outward unit normal of the face opposite vertex k
    int face[4][3];               // that face's vertices, counterclockwise about n[k]
};

Tet MakeTet(const double* P, const int* idx, std::size_t e) {
    Tet T;
    V3 g[4];
    for (int k = 0; k < 4; ++k)
        g[k] = {P[3 * idx[k]], P[3 * idx[k] + 1], P[3 * idx[k] + 2]};
    T.centroid = mul(0.25, add(add(g[0], g[1]), add(g[2], g[3])));
    for (int k = 0; k < 4; ++k)
        T.v[k] = sub(g[k], T.centroid);
    T.volume = std::fabs(dot(sub(T.v[1], T.v[0]), cross(sub(T.v[2], T.v[0]), sub(T.v[3], T.v[0])))) / 6.0;
    double scale = 0.0;
    for (int k = 0; k < 4; ++k) {
        T.radius = std::max(T.radius, norm(T.v[k]));
        for (int j = 0; j < 4; ++j)
            scale = std::max(scale, norm(sub(T.v[k], T.v[j])));
    }
    if (!(T.volume > 1e-12 * scale * scale * scale) || !std::isfinite(T.volume))
        throw std::invalid_argument("peec_tet: tetrahedron " + std::to_string(e) + " is degenerate");
    for (int k = 0; k < 4; ++k) {
        int f[3], m = 0;
        for (int j = 0; j < 4; ++j)
            if (j != k)
                f[m++] = j;
        V3 nn = cross(sub(T.v[f[1]], T.v[f[0]]), sub(T.v[f[2]], T.v[f[0]]));
        if (dot(nn, sub(T.v[k], T.v[f[0]])) > 0.0) {
            std::swap(f[1], f[2]);
            nn = mul(-1.0, nn);
        }
        T.n[k] = mul(1.0 / norm(nn), nn);
        for (int j = 0; j < 3; ++j)
            T.face[k][j] = f[j];
    }
    return T;
}

// x and the returned p1 = int (y - c)/R dy in the local frame of T (relative to its centroid c).
void Potential(const V3& x, const Tet& T, double& p0, V3& p1) {
    p0 = 0.0;
    V3 flux = {0.0, 0.0, 0.0};
    for (int k = 0; k < 4; ++k) {
        const V3& a = T.v[T.face[k][0]];
        double im1, ip1;
        FaceIntegrals(x, a, T.v[T.face[k][1]], T.v[T.face[k][2]], T.n[k], im1, ip1);
        p0 += 0.5 * dot(sub(a, x), T.n[k]) * im1;     // div_y((y - x)/R) = 2/R
        flux = add(flux, mul(ip1, T.n[k]));           // int (y - x)/R dy = int_dT R n dS
    }
    p1 = add(mul(p0, x), flux);
}

// Collapsed Gauss rule (Duffy map) with n^3 points on a tetrahedron.
struct Rule { std::vector<V3> x; std::vector<double> w; };

void GaussLegendre01(int n, std::vector<double>& t, std::vector<double>& w) {
    static const double x2[] = {-0.5773502691896257645, 0.5773502691896257645};
    static const double w2[] = {1.0, 1.0};
    static const double x3[] = {-0.7745966692414833770, 0.0, 0.7745966692414833770};
    static const double w3[] = {5.0 / 9.0, 8.0 / 9.0, 5.0 / 9.0};
    static const double x4[] = {-0.8611363115940525752, -0.3399810435848562648,
                                0.3399810435848562648, 0.8611363115940525752};
    static const double w4[] = {0.3478548451374538574, 0.6521451548625461426,
                                0.6521451548625461426, 0.3478548451374538574};
    const double* xs = n == 2 ? x2 : n == 3 ? x3 : x4;
    const double* ws = n == 2 ? w2 : n == 3 ? w3 : w4;
    t.resize(n);
    w.resize(n);
    for (int i = 0; i < n; ++i) {
        t[i] = 0.5 * (xs[i] + 1.0);
        w[i] = 0.5 * ws[i];
    }
}

// Collapsed rule in the local frame; the apex and edge order follow the lexicographic order of
// the vertex coordinates, so a relabeled element gets the same points.
Rule MakeRule(const Tet& T, int n) {
    std::vector<double> t, w;
    GaussLegendre01(n, t, w);
    Rule r;
    int o[4] = {0, 1, 2, 3};
    std::sort(o, o + 4, [&](int a, int b) { return T.v[a] < T.v[b]; });
    const V3& v0 = T.v[o[0]];
    const V3 e1 = sub(T.v[o[1]], v0), e2 = sub(T.v[o[2]], v0), e3 = sub(T.v[o[3]], v0);
    for (int i = 0; i < n; ++i)
        for (int j = 0; j < n; ++j)
            for (int k = 0; k < n; ++k) {
                const double u = t[i], v = t[j], s = t[k];
                const double a = u, b = v * (1 - u), c = s * (1 - u) * (1 - v);
                r.x.push_back(add(v0, add(mul(a, e1), add(mul(b, e2), mul(c, e3)))));
                r.w.push_back(w[i] * w[j] * w[k] * (1 - u) * (1 - u) * (1 - v) * 6.0 * T.volume);
            }
    return r;
}

}  // namespace

void TetPotential(const double x[3], const double v[4][3], double& p0, double p1[3]) {
    double P[12];
    for (int k = 0; k < 4; ++k)
        for (int j = 0; j < 3; ++j)
            P[3 * k + j] = v[k][j];
    const int idx[4] = {0, 1, 2, 3};
    const Tet T = MakeTet(P, idx, 0);
    V3 q;
    Potential(sub({x[0], x[1], x[2]}, T.centroid), T, p0, q);
    for (int j = 0; j < 3; ++j)
        p1[j] = q[j] + p0 * T.centroid[j];
}

namespace {

// Directed moments of the pair (A outer, B inner), classified by centroid distance; x is taken
// relative to the centroid of A and y to that of B, so no moment carries the global position.
void DirectedMoments(const Tet& A, const Rule& nearA, const Rule& modA, const Tet& B, const Rule& modB,
                     double near_factor, double moderate_factor, PairMoments& M) {
    const V3 dab = sub(A.centroid, B.centroid);
    const double d = norm(dab);
    const double s = A.radius + B.radius;
    double S0 = 0.0, Sxp = 0.0;
    V3 Sx = {0, 0, 0}, Sp = {0, 0, 0};
    auto accumulate = [&](const V3& x, double w, double p0, const V3& p1) {
        S0 += w * p0;
        Sx = add(Sx, mul(w * p0, x));
        Sp = add(Sp, mul(w, p1));
        Sxp += w * dot(x, p1);
    };
    if (d < near_factor * s) {
        for (std::size_t q = 0; q < nearA.x.size(); ++q) {
            double p0;
            V3 p1;
            Potential(add(dab, nearA.x[q]), B, p0, p1);
            accumulate(nearA.x[q], nearA.w[q], p0, p1);
        }
    } else if (d < moderate_factor * s) {
        for (std::size_t q = 0; q < modA.x.size(); ++q) {
            double p0 = 0.0;
            V3 p1 = {0, 0, 0};
            for (std::size_t j = 0; j < modB.x.size(); ++j) {
                const double g = modB.w[j] / norm(sub(add(dab, modA.x[q]), modB.x[j]));
                p0 += g;
                p1 = add(p1, mul(g, modB.x[j]));
            }
            accumulate(modA.x[q], modA.w[q], p0, p1);
        }
    } else {
        const double g = B.volume / d;
        accumulate({0.0, 0.0, 0.0}, A.volume, g, {0.0, 0.0, 0.0});
    }
    M.S0 = S0;
    M.Sxp = Sxp;
    for (int c = 0; c < 3; ++c) {
        M.Sx[c] = Sx[c];
        M.Sp[c] = Sp[c];
    }
}

// Value of the directed pair for the bases (x - va) on A and (y - vb) on B.
inline double PairValue(const PairMoments& M, const V3& va, const V3& vb) {
    const V3 Sx = {M.Sx[0], M.Sx[1], M.Sx[2]}, Sp = {M.Sp[0], M.Sp[1], M.Sp[2]};
    return M.Sxp - dot(va, Sp) - dot(vb, Sx) + dot(va, vb) * M.S0;
}

// Validates the mesh and its DOF table and builds the tetrahedra.
std::vector<Tet> CheckedTets(const TetMesh& mesh, double near_factor, double moderate_factor, int near_order) {
    const std::size_t nt = mesh.n_tet, nd = mesh.n_dof;
    if (!mesh.points || !mesh.tets || !mesh.face_coef || !mesh.face_dof || nt == 0 || nd == 0)
        throw std::invalid_argument("peec_tet: empty mesh or DOF table");
    if (!(near_factor > 0.0) || !(moderate_factor >= near_factor) || near_order < 2 || near_order > 4)
        throw std::invalid_argument("peec_tet: invalid pair classes or near order");
    for (std::size_t i = 0; i < 3 * mesh.n_points; ++i)
        if (!std::isfinite(mesh.points[i]))
            throw std::invalid_argument("peec_tet: non-finite point coordinate");
    std::vector<Tet> T(nt);
    std::vector<int> seen(nd, 0);
    for (std::size_t e = 0; e < nt; ++e) {
        for (int k = 0; k < 4; ++k) {
            const int p = mesh.tets[4 * e + k];
            if (p < 0 || static_cast<std::size_t>(p) >= mesh.n_points)
                throw std::invalid_argument("peec_tet: tetrahedron " + std::to_string(e) + " references point "
                                            + std::to_string(p));
            const int d = mesh.face_dof[4 * e + k];
            if (d < -1 || d >= static_cast<int>(nd))
                throw std::invalid_argument("peec_tet: DOF " + std::to_string(d) + " out of range");
            if (d >= 0) {
                if (!std::isfinite(mesh.face_coef[4 * e + k]) || mesh.face_coef[4 * e + k] == 0.0)
                    throw std::invalid_argument("peec_tet: invalid basis coefficient");
                ++seen[d];
            }
        }
        T[e] = MakeTet(mesh.points, mesh.tets + 4 * e, e);
    }
    for (std::size_t d = 0; d < nd; ++d)
        if (seen[d] < 1 || seen[d] > 2)
            throw std::invalid_argument("peec_tet: DOF " + std::to_string(d) + " belongs to "
                                        + std::to_string(seen[d]) + " tetrahedra");
    return T;
}

std::atomic<std::uint64_t> g_kernel_id{0};

}  // namespace

std::vector<double> AssembleInductance(const TetMesh& mesh, double near_factor, double moderate_factor,
                                       int near_order) {
    const std::size_t nt = mesh.n_tet, nd = mesh.n_dof;
    const std::vector<Tet> T = CheckedTets(mesh, near_factor, moderate_factor, near_order);
    // Colors: tetrahedra of one color share no DOF, so their rows of L are written by one worker.
    std::vector<std::vector<std::size_t>> by_dof(nd);
    for (std::size_t e = 0; e < nt; ++e)
        for (int k = 0; k < 4; ++k)
            if (mesh.face_dof[4 * e + k] >= 0)
                by_dof[mesh.face_dof[4 * e + k]].push_back(e);
    std::vector<int> color(nt, -1);
    int n_color = 0;
    for (std::size_t e = 0; e < nt; ++e) {
        unsigned used = 0;
        for (int k = 0; k < 4; ++k) {
            const int d = mesh.face_dof[4 * e + k];
            if (d < 0)
                continue;
            for (std::size_t o : by_dof[d])
                if (color[o] >= 0)
                    used |= 1u << color[o];
        }
        int c = 0;
        while (used & (1u << c))
            ++c;
        color[e] = c;
        n_color = std::max(n_color, c + 1);
    }
    std::vector<Rule> near(nt), mod(nt);
    for (std::size_t e = 0; e < nt; ++e) {
        near[e] = MakeRule(T[e], near_order);
        mod[e] = MakeRule(T[e], 2);
    }
    std::vector<double> L(nd * nd, 0.0);
    for (int c = 0; c < n_color; ++c) {
        std::vector<std::size_t> members;
        for (std::size_t e = 0; e < nt; ++e)
            if (color[e] == c)
                members.push_back(e);
        radia::host::ParallelFor(radia::host::IntRange(members.size()), [&](std::size_t m) {
            const std::size_t e = members[m];
            const Tet& A = T[e];
            for (std::size_t f = 0; f < nt; ++f) {
                const Tet& B = T[f];
                PairMoments M;
                DirectedMoments(A, near[e], mod[e], B, mod[f], near_factor, moderate_factor, M);
                for (int k = 0; k < 4; ++k) {
                    const int di = mesh.face_dof[4 * e + k];
                    if (di < 0)
                        continue;
                    for (int kk = 0; kk < 4; ++kk) {
                        const int dj = mesh.face_dof[4 * f + kk];
                        if (dj < 0)
                            continue;
                        L[static_cast<std::size_t>(di) * nd + dj] += kMu0Over4Pi * mesh.face_coef[4 * e + k]
                                                                       * mesh.face_coef[4 * f + kk]
                                                                       * PairValue(M, A.v[k], B.v[kk]);
                    }
                }
            }
        });
    }
    for (std::size_t i = 0; i < nd; ++i)                 // average the two directed rules
        for (std::size_t j = i + 1; j < nd; ++j) {
            const double v = 0.5 * (L[i * nd + j] + L[j * nd + i]);
            L[i * nd + j] = L[j * nd + i] = v;
        }
    return L;
}

struct InductanceKernel::Impl {
    struct Owner { int tet = -1, k = -1; double c = 0.0; };
    std::vector<Tet> T;
    std::vector<Rule> near, mod;
    std::vector<std::array<Owner, 2>> owners;
    double near_factor = 3.0, moderate_factor = 10.0;
    std::uint64_t id = 0;

    bool IsQuadrature(int e, int f) const {
        return norm(sub(T[e].centroid, T[f].centroid)) < moderate_factor * (T[e].radius + T[f].radius);
    }
    // Near and moderate pairs dominate the cost and recur in up to 16 entries, so a bounded
    // per-thread cache keeps them; a cached value is bitwise the directly computed one.
    PairMoments Moments(int e, int f) const {
        PairMoments M;
        if (!IsQuadrature(e, f)) {
            DirectedMoments(T[e], near[e], mod[e], T[f], mod[f], near_factor, moderate_factor, M);
            return M;
        }
        struct Cache {
            std::uint64_t id = 0;
            std::unordered_map<std::uint64_t, PairMoments> map;
        };
        thread_local Cache cache;
        if (cache.id != id || cache.map.size() >= (std::size_t(1) << 20)) {
            cache.map.clear();
            cache.id = id;
        }
        const std::uint64_t key = static_cast<std::uint64_t>(e) * T.size() + static_cast<std::uint64_t>(f);
        const auto it = cache.map.find(key);
        if (it != cache.map.end())
            return it->second;
        DirectedMoments(T[e], near[e], mod[e], T[f], mod[f], near_factor, moderate_factor, M);
        cache.map.emplace(key, M);
        return M;
    }
};

InductanceKernel::InductanceKernel(const TetMesh& mesh, double near_factor, double moderate_factor,
                                   int near_order)
    : impl_(new Impl) {
    Impl& I = *impl_;
    I.T = CheckedTets(mesh, near_factor, moderate_factor, near_order);
    I.near_factor = near_factor;
    I.moderate_factor = moderate_factor;
    I.id = ++g_kernel_id;
    const std::size_t nt = mesh.n_tet;
    I.near.resize(nt);
    I.mod.resize(nt);
    for (std::size_t e = 0; e < nt; ++e) {
        I.near[e] = MakeRule(I.T[e], near_order);
        I.mod[e] = MakeRule(I.T[e], 2);
    }
    I.owners.assign(mesh.n_dof, {});
    for (std::size_t e = 0; e < nt; ++e)
        for (int k = 0; k < 4; ++k) {
            const int d = mesh.face_dof[4 * e + k];
            if (d < 0)
                continue;
            auto& o = I.owners[d];
            Impl::Owner& slot = o[0].tet < 0 ? o[0] : o[1];
            slot.tet = static_cast<int>(e);
            slot.k = k;
            slot.c = mesh.face_coef[4 * e + k];
        }
}

InductanceKernel::~InductanceKernel() = default;

std::size_t InductanceKernel::n_dof() const { return impl_->owners.size(); }

double InductanceKernel::Entry(int i, int j) const {
    const Impl& I = *impl_;
    if (i < 0 || j < 0 || static_cast<std::size_t>(i) >= I.owners.size()
        || static_cast<std::size_t>(j) >= I.owners.size())
        throw std::out_of_range("peec_tet::InductanceKernel::Entry: DOF out of range");
    double sum = 0.0;
    for (const auto& a : I.owners[i]) {
        if (a.tet < 0)
            continue;
        for (const auto& b : I.owners[j]) {
            if (b.tet < 0)
                continue;
            const V3& va = I.T[a.tet].v[a.k];
            const V3& vb = I.T[b.tet].v[b.k];
            const double forward = PairValue(I.Moments(a.tet, b.tet), va, vb);
            const double backward = PairValue(I.Moments(b.tet, a.tet), vb, va);
            sum += a.c * b.c * 0.5 * (forward + backward);
        }
    }
    return kMu0Over4Pi * sum;
}

void InductanceKernel::Support(int i, double center[3], double& radius) const {
    const Impl& I = *impl_;
    if (i < 0 || static_cast<std::size_t>(i) >= I.owners.size())
        throw std::out_of_range("peec_tet::InductanceKernel::Support: DOF out of range");
    const auto& o = I.owners[i];
    const Tet& A = I.T[o[0].tet];
    V3 c = {0, 0, 0};                                  // face centroid, local frame of A
    for (int j = 0; j < 4; ++j)
        if (j != o[0].k)
            c = add(c, mul(1.0 / 3.0, A.v[j]));
    radius = 0.0;
    for (const auto& a : o) {
        if (a.tet < 0)
            continue;
        const V3 shift = sub(I.T[a.tet].centroid, A.centroid);
        for (int j = 0; j < 4; ++j)
            radius = std::max(radius, norm(sub(add(shift, I.T[a.tet].v[j]), c)));
    }
    for (int j = 0; j < 3; ++j)
        center[j] = A.centroid[j] + c[j];
}

KKTSolution SolveKKT(std::size_t n_dof, const std::vector<double>& L, const Triplets& R, std::size_t n_con,
                     const Triplets& D, const std::vector<std::complex<double>>& b, double omega,
                     double max_bytes) {
#if !defined(HAVE_LAPACK) || !(defined(RADIA_BLAS_MKL) || defined(RADIA_BLAS_NVPL))
    (void)n_dof; (void)L; (void)R; (void)n_con; (void)D; (void)b; (void)omega; (void)max_bytes;
    throw std::runtime_error("peec_tet::SolveKKT: built without LAPACK");
#else
    const std::size_t N = n_dof + n_con;
    if (L.size() != n_dof * n_dof || b.size() != n_dof)
        throw std::invalid_argument("peec_tet::SolveKKT: L or b does not match n_dof");
    if (!(omega >= 0.0) || !std::isfinite(omega))
        throw std::invalid_argument("peec_tet::SolveKKT: omega must be finite and nonnegative");
    const double bytes = 2.0 * 16.0 * static_cast<double>(N) * static_cast<double>(N);   // A and its LU
    if (!(bytes <= max_bytes))
        throw std::invalid_argument("peec_tet::SolveKKT: dense system of " + std::to_string(N) + " unknowns needs "
                                    + std::to_string(bytes / 1073741824.0) + " GiB, above the "
                                    + std::to_string(max_bytes / 1073741824.0) + " GiB budget");
    std::vector<std::complex<double>> A(N * N, std::complex<double>(0.0, 0.0));   // column-major
    auto at = [&](std::size_t i, std::size_t j) -> std::complex<double>& { return A[j * N + i]; };
    for (std::size_t i = 0; i < n_dof; ++i)
        for (std::size_t j = 0; j < n_dof; ++j)
            at(i, j) = std::complex<double>(0.0, omega * L[i * n_dof + j]);
    if (R.row.size() != R.col.size() || R.row.size() != R.value.size()
        || D.row.size() != D.col.size() || D.row.size() != D.value.size())
        throw std::invalid_argument("peec_tet::SolveKKT: ragged triplets");
    for (std::size_t t = 0; t < R.row.size(); ++t) {
        if (R.row[t] < 0 || R.col[t] < 0 || static_cast<std::size_t>(R.row[t]) >= n_dof
            || static_cast<std::size_t>(R.col[t]) >= n_dof)
            throw std::invalid_argument("peec_tet::SolveKKT: R index out of range");
        at(R.row[t], R.col[t]) += R.value[t];
    }
    for (std::size_t t = 0; t < D.row.size(); ++t) {
        if (D.row[t] < 0 || D.col[t] < 0 || static_cast<std::size_t>(D.row[t]) >= n_con
            || static_cast<std::size_t>(D.col[t]) >= n_dof)
            throw std::invalid_argument("peec_tet::SolveKKT: D index out of range");
        at(n_dof + D.row[t], D.col[t]) += D.value[t];
        at(D.col[t], n_dof + D.row[t]) += D.value[t];
    }
    std::vector<std::complex<double>> LU = A, x(N, std::complex<double>(0.0, 0.0));
    std::copy(b.begin(), b.end(), x.begin());
    std::vector<int> ipiv(N);
    int ln = static_cast<int>(N), nrhs = 1, info = 0;
    char trans = 'N';
    RADIA_TET_GETRF(&ln, &ln, reinterpret_cast<radia_blas_complex_double*>(LU.data()), &ln, ipiv.data(), &info);
    if (info != 0)
        throw std::runtime_error("peec_tet::SolveKKT: zgetrf info " + std::to_string(info)
                                 + " (singular: unconstrained loop or ungrounded potential?)");
    RADIA_TET_GETRS(&trans, &ln, &nrhs, reinterpret_cast<radia_blas_complex_double*>(LU.data()), &ln,
                    ipiv.data(), reinterpret_cast<radia_blas_complex_double*>(x.data()), &ln, &info);
    if (info != 0)
        throw std::runtime_error("peec_tet::SolveKKT: zgetrs info " + std::to_string(info));
    double r2 = 0.0, b2 = 0.0;
    for (std::size_t i = 0; i < N; ++i) {
        std::complex<double> r = i < n_dof ? b[i] : std::complex<double>(0.0, 0.0);
        b2 += std::norm(r);
        for (std::size_t j = 0; j < N; ++j)
            r -= A[j * N + i] * x[j];
        r2 += std::norm(r);
    }
    KKTSolution sol;
    sol.true_relative_residual = std::sqrt(r2 / b2);
    if (!(sol.true_relative_residual <= 1e-6))
        throw std::runtime_error("peec_tet::SolveKKT: true relative residual "
                                 + std::to_string(sol.true_relative_residual) + " above 1e-6");
    sol.current.assign(x.begin(), x.begin() + n_dof);
    sol.potential.assign(x.begin() + n_dof, x.end());
    return sol;
#endif
}

}  // namespace peec_tet
}  // namespace radia
