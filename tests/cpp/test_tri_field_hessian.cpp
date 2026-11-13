// Guards: QuadTriFieldBasisHessian, the receiver Hessian d2H_k/dr_j dr_l of the
// degree <= 2 monomial surface-charge fields of a flat triangle, for receivers off
// the source plane (closed form) and far receivers (collapsed Gauss rule).
// Assumes: two references that share no code with the kernel under test:
// (1) a Richardson central difference of the exact production Jacobian
// QuadTriFieldDirectional, and (2) tensor-product Gauss quadrature of the analytic
// kernel derivative on a subdivided triangle, written here without Radia calls.
// DETECTOR_REL values are proposed regression detectors pending user approval,
// not adopted accuracy criteria.
#include "rad_hdiv_vim.h"
#include <chrono>
#include <cmath>
#include <cstdio>
#include <vector>

namespace {

constexpr double DETECTOR_REL = 1e-8;      // closed form / far rule vs references
constexpr double DETECTOR_SYMMETRY = 1e-10; // relative asymmetry and trace

using Tensor = double[10][3][3][3];

double norm_inf(const Tensor& t, int i)
{
    double m = 0.0;
    for (int k = 0; k < 3; ++k) for (int j = 0; j < 3; ++j) for (int l = 0; l < 3; ++l)
        m = std::fmax(m, std::fabs(t[i][k][j][l]));
    return m;
}

void monomial_coefficients(int index, double& sigma0, double s[3], double S[3][3])
{
    // PotentialMomentIndex order: 1, z, y, x, zz, yz, yy, xz, xy, xx.
    static const int alpha[10][3] = {{0,0,0},{0,0,1},{0,1,0},{1,0,0},{0,0,2},
                                     {0,1,1},{0,2,0},{1,0,1},{1,1,0},{2,0,0}};
    sigma0 = 0.0;
    for (int a = 0; a < 3; ++a) { s[a] = 0.0; for (int b = 0; b < 3; ++b) S[a][b] = 0.0; }
    const int* e = alpha[index];
    const int total = e[0] + e[1] + e[2];
    if (total == 0) sigma0 = 1.0;
    else if (total == 1) { for (int a = 0; a < 3; ++a) if (e[a]) s[a] = 1.0; }
    else {
        int axes[2], count = 0;
        for (int a = 0; a < 3; ++a) for (int p = 0; p < e[a]; ++p) axes[count++] = a;
        S[axes[0]][axes[1]] += 0.5; S[axes[1]][axes[0]] += 0.5;
    }
}

// d/dr_l of H_k for monomial `index`, from the exact production forward derivative.
void jacobian_column(const double V[3][3], const double r[3], int index, int l, double out[3])
{
    double sigma0, s[3], S[3][3];
    monomial_coefficients(index, sigma0, s, S);
    const double zero33[3][3] = {}, zero3[3] = {};
    double dr[3] = {0.0, 0.0, 0.0}; dr[l] = 1.0;
    double value[3];
    rad_hdiv::QuadTriFieldDirectional(V, zero33, r, dr, sigma0, 0.0, s, zero3, S, zero33,
                                      value, out);
}

// Richardson central difference of the exact Jacobian in r_j.
void finite_difference(const double V[3][3], const double r[3], double step, Tensor& out)
{
    for (int i = 0; i < 10; ++i)
        for (int j = 0; j < 3; ++j)
            for (int l = 0; l < 3; ++l) {
                double d[2][3];
                for (int level = 0; level < 2; ++level) {
                    const double delta = step / (1 << level);
                    double rp[3] = {r[0], r[1], r[2]}, rm[3] = {r[0], r[1], r[2]};
                    rp[j] += delta; rm[j] -= delta;
                    double jp[3], jm[3];
                    jacobian_column(V, rp, i, l, jp);
                    jacobian_column(V, rm, i, l, jm);
                    for (int k = 0; k < 3; ++k) d[level][k] = (jp[k] - jm[k]) / (2.0 * delta);
                }
                for (int k = 0; k < 3; ++k) out[i][k][j][l] = (4.0 * d[1][k] - d[0][k]) / 3.0;
            }
}

// Gauss-Legendre nodes on [0,1] by Newton iteration.
void gauss_legendre(int n, std::vector<double>& x, std::vector<double>& w)
{
    x.assign(n, 0.0); w.assign(n, 0.0);
    const double pi = 3.14159265358979323846;
    for (int i = 0; i < n; ++i) {
        double z = std::cos(pi * (i + 0.75) / (n + 0.5)), dp = 0.0;
        for (int it = 0; it < 100; ++it) {
            double p0 = 1.0, p1 = z;
            for (int k = 2; k <= n; ++k) { const double p2 = ((2*k-1)*z*p1 - (k-1)*p0)/k; p0 = p1; p1 = p2; }
            dp = n * (z * p1 - p0) / (z * z - 1.0);
            const double dz = p1 / dp; z -= dz;
            if (std::fabs(dz) < 1e-16) break;
        }
        x[i] = 0.5 * (1.0 - z); w[i] = 1.0 / ((1.0 - z * z) * dp * dp) * 2.0 * 0.5;
    }
}

// INT_T x^alpha d_j d_l [(r-x)_k/|r-x|^3] dS on 4^levels subtriangles, n x n collapsed Gauss each.
void quadrature(const double V[3][3], const double r[3], int levels, int n, Tensor& out)
{
    static const int alpha[10][3] = {{0,0,0},{0,0,1},{0,1,0},{1,0,0},{0,0,2},
                                     {0,1,1},{0,2,0},{1,0,1},{1,1,0},{2,0,0}};
    for (int i = 0; i < 10; ++i) for (int k = 0; k < 3; ++k) for (int j = 0; j < 3; ++j)
        for (int l = 0; l < 3; ++l) out[i][k][j][l] = 0.0;
    std::vector<double> gx, gw; gauss_legendre(n, gx, gw);
    struct Tri { double v[3][3]; };
    std::vector<Tri> tris(1);
    for (int a = 0; a < 3; ++a) for (int k = 0; k < 3; ++k) tris[0].v[a][k] = V[a][k];
    for (int level = 0; level < levels; ++level) {
        std::vector<Tri> next;
        for (const Tri& t : tris) {
            double mid[3][3];
            for (int a = 0; a < 3; ++a) for (int k = 0; k < 3; ++k)
                mid[a][k] = 0.5 * (t.v[a][k] + t.v[(a+1)%3][k]);
            Tri c[4];
            for (int k = 0; k < 3; ++k) {
                c[0].v[0][k]=t.v[0][k]; c[0].v[1][k]=mid[0][k]; c[0].v[2][k]=mid[2][k];
                c[1].v[0][k]=mid[0][k]; c[1].v[1][k]=t.v[1][k]; c[1].v[2][k]=mid[1][k];
                c[2].v[0][k]=mid[2][k]; c[2].v[1][k]=mid[1][k]; c[2].v[2][k]=t.v[2][k];
                c[3].v[0][k]=mid[0][k]; c[3].v[1][k]=mid[1][k]; c[3].v[2][k]=mid[2][k];
            }
            for (const Tri& child : c) next.push_back(child);
        }
        tris.swap(next);
    }
    for (const Tri& t : tris) {
        double e1[3], e2[3], cr[3];
        for (int k = 0; k < 3; ++k) { e1[k] = t.v[1][k]-t.v[0][k]; e2[k] = t.v[2][k]-t.v[0][k]; }
        cr[0]=e1[1]*e2[2]-e1[2]*e2[1]; cr[1]=e1[2]*e2[0]-e1[0]*e2[2]; cr[2]=e1[0]*e2[1]-e1[1]*e2[0];
        const double area2 = std::sqrt(cr[0]*cr[0]+cr[1]*cr[1]+cr[2]*cr[2]);
        for (int iu = 0; iu < n; ++iu) for (int iv = 0; iv < n; ++iv) {
            const double u = gx[iu], v = gx[iv];
            const double weight = gw[iu] * gw[iv] * area2 * (1.0 - u);
            double x[3], d[3];
            for (int k = 0; k < 3; ++k) { x[k] = t.v[0][k] + u*e1[k] + (1.0-u)*v*e2[k]; d[k] = r[k]-x[k]; }
            const double d2 = d[0]*d[0]+d[1]*d[1]+d[2]*d[2], R = std::sqrt(d2);
            const double R5 = d2*d2*R, R7 = R5*d2;
            for (int i = 0; i < 10; ++i) {
                double mono = 1.0;
                for (int a = 0; a < 3; ++a) for (int p = 0; p < alpha[i][a]; ++p) mono *= x[a];
                for (int k = 0; k < 3; ++k) for (int j = 0; j < 3; ++j) for (int l = 0; l < 3; ++l) {
                    const double kernel = 15.0*d[j]*d[k]*d[l]/R7
                        - 3.0*((k==l?d[j]:0.0)+(j==k?d[l]:0.0)+(j==l?d[k]:0.0))/R5;
                    out[i][k][j][l] += weight * mono * kernel;
                }
            }
        }
    }
}

double max_rel(const Tensor& a, const Tensor& b)
{
    double worst = 0.0;
    for (int i = 0; i < 10; ++i) {
        const double scale = std::fmax(norm_inf(b, i), 1e-300);
        for (int k = 0; k < 3; ++k) for (int j = 0; j < 3; ++j) for (int l = 0; l < 3; ++l)
            worst = std::fmax(worst, std::fabs(a[i][k][j][l] - b[i][k][j][l]) / scale);
    }
    return worst;
}

void symmetry(const Tensor& t, double& asym, double& trace)
{
    asym = trace = 0.0;
    for (int i = 0; i < 10; ++i) {
        const double scale = std::fmax(norm_inf(t, i), 1e-300);
        for (int k = 0; k < 3; ++k) {
            double tr = 0.0;
            for (int j = 0; j < 3; ++j) {
                tr += t[i][k][j][j];
                for (int l = 0; l < 3; ++l) {
                    const double v = t[i][k][j][l];
                    asym = std::fmax(asym, std::fabs(v - t[i][k][l][j]) / scale);
                    asym = std::fmax(asym, std::fabs(v - t[i][j][k][l]) / scale);
                    asym = std::fmax(asym, std::fabs(v - t[i][l][j][k]) / scale);
                }
            }
            trace = std::fmax(trace, std::fabs(tr) / scale);
        }
    }
}

}  // namespace

int main()
{
    // A generic triangle in 3D, away from the origin so every monomial is nontrivial.
    const double V[3][3] = {{0.31, -0.12, 0.47}, {1.08, 0.21, 0.29}, {0.44, 0.93, 0.71}};
    double e1[3], e2[3], n[3], c[3];
    for (int k = 0; k < 3; ++k) { e1[k] = V[1][k]-V[0][k]; e2[k] = V[2][k]-V[0][k]; c[k] = (V[0][k]+V[1][k]+V[2][k])/3.0; }
    n[0]=e1[1]*e2[2]-e1[2]*e2[1]; n[1]=e1[2]*e2[0]-e1[0]*e2[2]; n[2]=e1[0]*e2[1]-e1[1]*e2[0];
    const double nn = std::sqrt(n[0]*n[0]+n[1]*n[1]+n[2]*n[2]);
    for (double& x : n) x /= nn;
    const double size = std::sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);

    struct Case { const char* name; double point[3]; bool quadrature_ok; };
    std::vector<Case> cases;
    auto add = [&](const char* name, double along, double off1, double off2, bool q) {
        Case cs{name, {0, 0, 0}, q};
        for (int k = 0; k < 3; ++k) cs.point[k] = c[k] + along*size*n[k] + off1*size*e1[k]/size + off2*size*e2[k]/size;
        cases.push_back(cs);
    };
    add("centroid h=0.5", 0.5, 0.0, 0.0, true);
    add("centroid h=0.1", 0.1, 0.0, 0.0, true);
    add("centroid h=-0.02", -0.02, 0.0, 0.0, false);
    add("over vertex h=0.2", 0.2, -0.33, -0.33, true);
    add("outside edge h=0.3", 0.3, 0.6, 0.6, true);
    add("outside edge h=0.05", 0.05, 0.6, 0.6, true);
    add("closed form below switch", 1.5, 0.0, 0.0, true);
    add("far rule above switch", 1.56, 0.0, 0.0, true);
    add("far rule low height", 0.05, 0.9, 0.9, true);
    add("far rule distant", 6.0, 0.0, 0.0, true);

    double circumradius = 0.0;
    for (int a = 0; a < 3; ++a) {
        double d2 = 0.0;
        for (int k = 0; k < 3; ++k) d2 += (V[a][k]-c[k])*(V[a][k]-c[k]);
        circumradius = std::fmax(circumradius, std::sqrt(d2));
    }
    // The kernel switches to the Gauss rule at a ratio >= 2.
    auto far_ratio = [&](const double p[3]) {
        double d2 = 0.0;
        for (int k = 0; k < 3; ++k) d2 += (p[k]-c[k])*(p[k]-c[k]);
        return std::sqrt(d2) / circumradius;
    };
    int failures = 0;
    std::printf("{\n  \"cases\": [\n");
    for (std::size_t ci = 0; ci < cases.size(); ++ci) {
        const Case& cs = cases[ci];
        static Tensor hessian, fd, quad, quad_coarse;
        const bool ok = rad_hdiv::QuadTriFieldBasisHessian(V, cs.point, hessian);
        if (!ok) { std::fprintf(stderr, "%s rejected\n", cs.name); ++failures; continue; }
        double asym, trace; symmetry(hessian, asym, trace);
        finite_difference(V, cs.point, 2e-3 * size, fd);
        const double fd_rel = max_rel(hessian, fd);
        double q_rel = -1.0, q_conv = -1.0;
        if (cs.quadrature_ok) {
            quadrature(V, cs.point, 3, 24, quad_coarse);
            quadrature(V, cs.point, 4, 24, quad);
            q_rel = max_rel(hessian, quad);
            q_conv = max_rel(quad_coarse, quad);
        }
        const bool pass = asym <= DETECTOR_SYMMETRY && trace <= DETECTOR_SYMMETRY
            && fd_rel <= 1e-6 && (q_rel < 0.0 || (q_rel <= DETECTOR_REL && q_conv <= 1e-11));
        if (!pass) ++failures;
        std::printf("    {\"name\": \"%s\", \"far_ratio\": %.3f, \"asymmetry\": %.3e, \"trace\": %.3e, \"fd_rel\": %.3e, "
                    "\"quadrature_rel\": %.3e, \"quadrature_self_convergence\": %.3e, \"pass\": %s}%s\n",
                    cs.name, far_ratio(cs.point), asym, trace, fd_rel, q_rel, q_conv, pass ? "true" : "false",
                    ci + 1 < cases.size() ? "," : "");
    }
    std::printf("  ],\n  \"near_plane_scan\": [");
    // Informational: closed-form accuracy as an interior receiver approaches the
    // plane, against the finite-difference reference with a step below the height.
    {
        const double heights[] = {1e-2, 1e-3, 2e-4, 1e-5, 1e-6};
        for (int hi = 0; hi < 5; ++hi) {
            const double h = heights[hi] * size;
            double p[3];
            for (int k = 0; k < 3; ++k) p[k] = c[k] + 0.1*e1[k] + h*n[k];
            static Tensor t, f;
            const bool ok = rad_hdiv::QuadTriFieldBasisHessian(V, p, t);
            finite_difference(V, p, 0.2 * h, f);
            std::printf("%s{\"h_over_edge\": %.0e, \"accepted\": %s, \"fd_rel\": %.3e}",
                        hi ? ", " : "", heights[hi], ok ? "true" : "false", max_rel(t, f));
        }
    }
    // Near-plane receivers are outside this scope and must be rejected, not approximated.
    {
        double p[3]; for (int k = 0; k < 3; ++k) p[k] = c[k] + 1e-12 * n[k];
        static Tensor t;
        if (rad_hdiv::QuadTriFieldBasisHessian(V, p, t)) { std::fprintf(stderr, "on-plane accepted\n"); ++failures; }
    }
    // Setup cost: one Hessian call against one QuadTriFieldBasis call at the same receiver.
    {
        double p[3]; for (int k = 0; k < 3; ++k) p[k] = c[k] + 0.3 * size * n[k];
        static Tensor t; double basis[10][3]; const int reps = 20000; double sink = 0.0;
        auto t0 = std::chrono::steady_clock::now();
        for (int i = 0; i < reps; ++i) { p[0] += 1e-9; rad_hdiv::QuadTriFieldBasisHessian(V, p, t); sink += t[9][0][0][0]; }
        auto t1 = std::chrono::steady_clock::now();
        for (int i = 0; i < reps; ++i) { p[0] += 1e-9; rad_hdiv::QuadTriFieldBasis(V, p, basis); sink += basis[9][0]; }
        auto t2 = std::chrono::steady_clock::now();
        const double hess_us = std::chrono::duration<double, std::micro>(t1 - t0).count() / reps;
        const double field_us = std::chrono::duration<double, std::micro>(t2 - t1).count() / reps;
        std::printf("],\n  \"timing_us_per_call\": {\"hessian\": %.3f, \"field_basis\": %.3f, \"sink\": %.3e},\n",
                    hess_us, field_us, sink);
    }
    std::printf("  \"failures\": %d\n}\n", failures);
    return failures;
}
