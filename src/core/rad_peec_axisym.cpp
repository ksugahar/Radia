// Axisymmetric volume PEEC of a revolved conductor section (see rad_peec_axisym.h).
#include "rad_peec_axisym.h"

#include <algorithm>
#include <cmath>
#include <stdexcept>
#include <string>

#include "rad_host_parallel.h"

#ifdef HAVE_LAPACK
#include "rad_blas.h"
#if defined(RADIA_BLAS_MKL)
#include <mkl_lapack.h>
#define RADIA_AXI_GETRF zgetrf_
#define RADIA_AXI_GETRS zgetrs_
#elif defined(RADIA_BLAS_NVPL)
#include <nvpl_lapack.h>
#define RADIA_AXI_GETRF NVPL_LAPACK_zgetrf
#define RADIA_AXI_GETRS NVPL_LAPACK_zgetrs
#endif
#endif

namespace radia {
namespace peec_axisym {
namespace {

constexpr double kPi = 3.141592653589793238462643383279502884;
constexpr double kMu0 = 4e-7 * kPi;

// Carlson symmetric elliptic integrals (Carlson 1995, duplication theorem).
double CarlsonRF(double x, double y, double z)
{
    for (int it = 0; it < 200; ++it) {
        const double a = (x + y + z) / 3.0;
        const double dx = 1.0 - x / a, dy = 1.0 - y / a, dz = 1.0 - z / a;
        if (std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) < 1e-4) {
            const double e2 = dx * dy - dz * dz, e3 = dx * dy * dz;
            return (1.0 + (e2 / 24.0 - 0.1 - 3.0 * e3 / 44.0) * e2 + e3 / 14.0) / std::sqrt(a);
        }
        const double sx = std::sqrt(x), sy = std::sqrt(y), sz = std::sqrt(z);
        const double lam = sx * (sy + sz) + sy * sz;
        x = 0.25 * (x + lam); y = 0.25 * (y + lam); z = 0.25 * (z + lam);
    }
    throw std::runtime_error("CarlsonRF did not converge");
}

double CarlsonRD(double x, double y, double z)
{
    double sum = 0.0, fac = 1.0;
    for (int it = 0; it < 200; ++it) {
        const double a = (x + y + 3.0 * z) / 5.0;
        const double dx = 1.0 - x / a, dy = 1.0 - y / a, dz = 1.0 - z / a;
        if (std::max({std::abs(dx), std::abs(dy), std::abs(dz)}) < 1e-4) {
            const double ea = dx * dy, eb = dz * dz, ec = ea - eb, ed = ea - 6.0 * eb, ee = ed + ec + ec;
            constexpr double c1 = 3.0 / 14.0, c2 = 1.0 / 6.0, c3 = 9.0 / 22.0, c4 = 3.0 / 26.0;
            constexpr double c5 = 0.25 * c3, c6 = 1.5 * c4;
            return 3.0 * sum + fac * (1.0 + ed * (-c1 + c5 * ed - c6 * dz * ee)
                                      + dz * (c2 * ee + dz * (-c3 * ec + dz * c4 * ea)))
                                   / (a * std::sqrt(a));
        }
        const double sx = std::sqrt(x), sy = std::sqrt(y), sz = std::sqrt(z);
        const double lam = sx * (sy + sz) + sy * sz;
        sum += fac / (sz * (z + lam));
        fac *= 0.25;
        x = 0.25 * (x + lam); y = 0.25 * (y + lam); z = 0.25 * (z + lam);
    }
    throw std::runtime_error("CarlsonRD did not converge");
}

// Dunavant rules on barycentric coordinates (weights sum to 1).
struct Rule { int n; double b[7][3]; double w[7]; };

const Rule& R3()
{
    static const Rule r{3, {{2.0 / 3, 1.0 / 6, 1.0 / 6}, {1.0 / 6, 2.0 / 3, 1.0 / 6}, {1.0 / 6, 1.0 / 6, 2.0 / 3}},
                        {1.0 / 3, 1.0 / 3, 1.0 / 3}};
    return r;
}

const Rule& R6()
{
    constexpr double a = 0.445948490915965, b = 0.091576213509771;
    static const Rule r{6, {{1 - 2 * a, a, a}, {a, 1 - 2 * a, a}, {a, a, 1 - 2 * a},
                           {1 - 2 * b, b, b}, {b, 1 - 2 * b, b}, {b, b, 1 - 2 * b}},
                        {0.223381589678011, 0.223381589678011, 0.223381589678011,
                         0.109951743655322, 0.109951743655322, 0.109951743655322}};
    return r;
}

const Rule& R7()
{
    constexpr double c = 0.470142064105115, d = 0.101286507323456;
    static const Rule r{7, {{1.0 / 3, 1.0 / 3, 1.0 / 3}, {1 - 2 * c, c, c}, {c, 1 - 2 * c, c}, {c, c, 1 - 2 * c},
                           {1 - 2 * d, d, d}, {d, 1 - 2 * d, d}, {d, d, 1 - 2 * d}},
                        {0.225, 0.132394152788506, 0.132394152788506, 0.132394152788506,
                         0.125939180544827, 0.125939180544827, 0.125939180544827}};
    return r;
}

struct Tri { double v[3][2]; double area; double c[2]; double size; };

void RulePoints(const Tri& t, const Rule& rule, double (*out)[2])
{
    for (int q = 0; q < rule.n; ++q)
        for (int d = 0; d < 2; ++d)
            out[q][d] = rule.b[q][0] * t.v[0][d] + rule.b[q][1] * t.v[1][d] + rule.b[q][2] * t.v[2][d];
}

// int_T ln|x - y| dA_y in closed form (sum over the oriented edges).
double LogPotential(const double x[2], const Tri& t)
{
    double out = 0.0;
    for (int i = 0; i < 3; ++i) {
        const double* p = t.v[i];
        const double* q = t.v[(i + 1) % 3];
        const double ex = q[0] - p[0], ey = q[1] - p[1];
        const double length = std::hypot(ex, ey);
        const double ux = ex / length, uy = ey / length;
        const double cross = (p[0] - x[0]) * (q[1] - x[1]) - (p[1] - x[1]) * (q[0] - x[0]);
        const double h = std::abs(cross) / length;
        if (!(h > 0.0))
            continue;
        const double t1 = (p[0] - x[0]) * ux + (p[1] - x[1]) * uy;
        const double t2 = t1 + length;
        auto F = [h](double s) { return (h / 4) * (s * std::log(h * h + s * s) - 3 * s + 2 * h * std::atan(s / h)); };
        out += (cross > 0 ? 1.0 : -1.0) * (F(t2) - F(t1));
    }
    return out;
}

// Galerkin entry with target t and source s, singular part in closed form,
// normalized by both areas.
double Directed(const Tri& t, const Tri& s)
{
    const Rule& ro = R7();
    const Rule& ri = R6();
    double X[7][2], Y[6][2];
    RulePoints(t, ro, X);
    RulePoints(s, ri, Y);
    double logpart = 0.0, rem = 0.0;
    for (int q = 0; q < ro.n; ++q) {
        const double r = X[q][0];
        const double base = std::log(8 * r) - 2.0;
        logpart += ro.w[q] * kMu0 * r * (base * s.area - LogPotential(X[q], s));
        double inner = 0.0;
        for (int p = 0; p < ri.n; ++p) {
            const double dist = std::hypot(X[q][0] - Y[p][0], X[q][1] - Y[p][1]);
            const double subtr = kMu0 * r * (base - std::log(dist));
            inner += ri.w[p] * (FilamentMutual(X[q][0], X[q][1], Y[p][0], Y[p][1]) - subtr);
        }
        rem += ro.w[q] * inner;
    }
    const double total = t.area * logpart + t.area * s.area * rem;
    return total / (t.area * s.area);
}

double Moderate(const Tri& a, const Tri& b)
{
    const Rule& r = R3();
    double X[3][2], Y[3][2];
    RulePoints(a, r, X);
    RulePoints(b, r, Y);
    double v = 0.0;
    for (int i = 0; i < r.n; ++i)
        for (int j = 0; j < r.n; ++j)
            v += r.w[i] * r.w[j] * FilamentMutual(X[i][0], X[i][1], Y[j][0], Y[j][1]);
    return v;
}

std::vector<Tri> BuildTriangles(const RingSection& sec)
{
    if (!sec.points || !sec.triangles || sec.n_tri == 0)
        throw std::invalid_argument("AssembleRingInductance: empty section");
    std::vector<Tri> tris(sec.n_tri);
    for (std::size_t k = 0; k < sec.n_tri; ++k) {
        Tri& t = tris[k];
        for (int j = 0; j < 3; ++j) {
            const int idx = sec.triangles[3 * k + j];
            if (idx < 0 || static_cast<std::size_t>(idx) >= sec.n_points)
                throw std::invalid_argument("AssembleRingInductance: triangle " + std::to_string(k)
                                            + " references point " + std::to_string(idx));
            t.v[j][0] = sec.points[2 * idx];
            t.v[j][1] = sec.points[2 * idx + 1];
            if (!(t.v[j][0] > 0.0) || !std::isfinite(t.v[j][0]) || !std::isfinite(t.v[j][1]))
                throw std::invalid_argument("AssembleRingInductance: point " + std::to_string(idx)
                                            + " is not off the axis (r > 0, finite)");
        }
        const double signed_area = 0.5 * ((t.v[1][0] - t.v[0][0]) * (t.v[2][1] - t.v[0][1])
                                          - (t.v[2][0] - t.v[0][0]) * (t.v[1][1] - t.v[0][1]));
        if (signed_area < 0.0)                      // the log potential needs counterclockwise vertices
            for (int d = 0; d < 2; ++d)
                std::swap(t.v[1][d], t.v[2][d]);
        t.area = std::abs(signed_area);
        if (!(t.area > 0.0))
            throw std::invalid_argument("AssembleRingInductance: degenerate meridian triangle "
                                        + std::to_string(k));
        for (int d = 0; d < 2; ++d)
            t.c[d] = (t.v[0][d] + t.v[1][d] + t.v[2][d]) / 3.0;
        t.size = std::sqrt(2 * t.area);
    }
    return tris;
}

}  // namespace

double FilamentMutual(double r1, double z1, double r2, double z2)
{
    const double s = (r1 + r2) * (r1 + r2) + (z1 - z2) * (z1 - z2);
    const double m = 4.0 * r1 * r2 / s;
    const double m1 = ((r1 - r2) * (r1 - r2) + (z1 - z2) * (z1 - z2)) / s;
    const double k = std::sqrt(m);
    if (m < 0.5) {
        // (2/k - k) K - (2/k) E = (pi/16) k^3 2F1(3/2, 3/2; 3; m): a positive series, free of the
        // cancellation of the closed form at small m (far filaments).
        double term = 1.0, sum = 1.0;
        for (int n = 0; n < 200 && term > 1e-17 * sum; ++n) {
            term *= (1.5 + n) * (1.5 + n) / ((3.0 + n) * (1.0 + n)) * m;
            sum += term;
        }
        return kMu0 * std::sqrt(r1 * r2) * (kPi / 16.0) * m * k * sum;
    }
    const double K = CarlsonRF(0.0, m1, 1.0);
    const double E = K - (m / 3.0) * CarlsonRD(0.0, m1, 1.0);
    return kMu0 * std::sqrt(r1 * r2) * ((2 / k - k) * K - (2 / k) * E);
}

RingInductance AssembleRingInductance(const RingSection& section, double near_factor, double moderate_factor)
{
    if (!(near_factor > 0.0) || !(moderate_factor >= near_factor))
        throw std::invalid_argument("AssembleRingInductance: need 0 < near_factor <= moderate_factor");
    const std::vector<Tri> tris = BuildTriangles(section);
    const std::size_t n = tris.size();
    RingInductance out;
    out.n = n;
    out.M.assign(n * n, 0.0);
    out.area.resize(n);
    out.centroid.resize(2 * n);
    for (std::size_t k = 0; k < n; ++k) {
        out.area[k] = tris[k].area;
        out.centroid[2 * k] = tris[k].c[0];
        out.centroid[2 * k + 1] = tris[k].c[1];
    }
    double* M = out.M.data();
    // Row a fills the upper triangle (b >= a) and its mirror: no two rows write the same entry.
    radia::host::ParallelFor(radia::host::IntRange(n), [&](std::size_t a) {
        const Tri& ta = tris[a];
        for (std::size_t b = a; b < n; ++b) {
            const Tri& tb = tris[b];
            const double dist = std::hypot(ta.c[0] - tb.c[0], ta.c[1] - tb.c[1]);
            const double reach = ta.size + tb.size;
            double v;
            if (b == a || dist < near_factor * reach)
                v = 0.5 * (Directed(ta, tb) + Directed(tb, ta));
            else if (dist < moderate_factor * reach)
                v = Moderate(ta, tb);
            else
                v = FilamentMutual(ta.c[0], ta.c[1], tb.c[0], tb.c[1]);
            M[a * n + b] = v;
            M[b * n + a] = v;
        }
    });
    return out;
}

RingSolution SolveLoopVoltage(const RingInductance& L, double frequency_hz, double sigma)
{
#ifndef HAVE_LAPACK
    (void)L; (void)frequency_hz; (void)sigma;
    throw std::runtime_error("SolveLoopVoltage: built without LAPACK");
#else
    const std::size_t n = L.n;
    if (n == 0 || L.M.size() != n * n || L.area.size() != n || L.centroid.size() != 2 * n)
        throw std::invalid_argument("SolveLoopVoltage: inconsistent ring inductance");
    for (std::size_t i = 0; i < n; ++i)              // the column-major solve below relies on M = M^T
        for (std::size_t j = i + 1; j < n; ++j)
            if (!(L.M[i * n + j] == L.M[j * n + i]))
                throw std::invalid_argument("SolveLoopVoltage: M is not symmetric at ("
                                            + std::to_string(i) + ", " + std::to_string(j) + ")");
    if (!(frequency_hz > 0.0) || !(sigma > 0.0))
        throw std::invalid_argument("SolveLoopVoltage: frequency and sigma must be positive");
    const double omega = 2 * kPi * frequency_hz;
    RingSolution sol;
    sol.n = n;
    sol.resistance.resize(n);
    for (std::size_t k = 0; k < n; ++k)
        sol.resistance[k] = 2 * kPi * L.centroid[2 * k] / (sigma * L.area[k]);
    // Column-major Z for LAPACK; M is symmetric, so the layouts coincide.
    std::vector<std::complex<double>> Z(n * n);
    for (std::size_t i = 0; i < n * n; ++i)
        Z[i] = std::complex<double>(0.0, omega * L.M[i]);
    for (std::size_t k = 0; k < n; ++k)
        Z[k * n + k] += sol.resistance[k];
    std::vector<std::complex<double>> LU = Z;
    std::vector<std::complex<double>> x(n, std::complex<double>(1.0, 0.0));
    std::vector<int> ipiv(n);
    int ln = static_cast<int>(n), nrhs = 1, info = 0;
    char trans = 'N';
    RADIA_AXI_GETRF(&ln, &ln, reinterpret_cast<radia_blas_complex_double*>(LU.data()), &ln, ipiv.data(), &info);
    if (info != 0)
        throw std::runtime_error("SolveLoopVoltage: zgetrf info " + std::to_string(info));
    RADIA_AXI_GETRS(&trans, &ln, &nrhs, reinterpret_cast<radia_blas_complex_double*>(LU.data()), &ln,
                    ipiv.data(), reinterpret_cast<radia_blas_complex_double*>(x.data()), &ln, &info);
    if (info != 0)
        throw std::runtime_error("SolveLoopVoltage: zgetrs info " + std::to_string(info));
    double res2 = 0.0;
    for (std::size_t i = 0; i < n; ++i) {
        std::complex<double> r(-1.0, 0.0);
        for (std::size_t j = 0; j < n; ++j)
            r += Z[j * n + i] * x[j];
        res2 += std::norm(r);
    }
    sol.true_relative_residual = std::sqrt(res2 / static_cast<double>(n));
    if (!(sol.true_relative_residual <= 1e-6))
        throw std::runtime_error("SolveLoopVoltage: true relative residual "
                                 + std::to_string(sol.true_relative_residual) + " above 1e-6");
    std::complex<double> total(0.0, 0.0);
    for (const auto& c : x)
        total += c;
    sol.current = std::move(x);
    sol.impedance = 1.0 / total;
    sol.cell_loss.resize(n);
    const double scale = 1.0 / std::norm(total);   // per 1 A peak terminal current
    for (std::size_t k = 0; k < n; ++k)
        sol.cell_loss[k] = 0.5 * sol.resistance[k] * std::norm(sol.current[k]) * scale;
    return sol;
#endif
}

}  // namespace peec_axisym
}  // namespace radia
