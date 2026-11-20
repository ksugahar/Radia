// Oriented-box Galerkin partial inductance between two rectangular bars.
//
// One implementation for the CPU builder and the CUDA assembly: the outer
// box is integrated by tensor Gauss-Legendre quadrature, the inner box by
// the closed-form Newtonian prism potential; an outer box that crosses the
// inner box surface (where the potential is only C^1) is integrated over a
// composite of sub-boxes bisected near that surface. Plain arithmetic only,
// no allocation, so the same code compiles as host and device functions.
#pragma once

#include <cmath>

#if defined(__CUDACC__)
#define RADIA_PEEC_HD __host__ __device__
#else
#define RADIA_PEEC_HD
#endif

namespace radia_peec_box {

constexpr int kMaxOrder = 32;
constexpr int kSplitDepth = 3;
constexpr int kStack = 32;  // depth-first stack bound: 1 + kSplitDepth * 7 = 22

struct V3 {
    double x, y, z;
};

struct Frame {
    V3 c, e[3];
    double half[3];
};

// Trace counters; see PEECMatrixBuilder::OrientedBoxTrace.
struct Trace {
    bool composite = false;
    int splits = 0, leaves = 0, straddling_leaves = 0, contained_top = 0,
        contained_descendants = 0, cross_axes_skipped = 0, cross_axes_evaluated = 0;
};

// Gauss-Legendre rules n = 1..order on [-1, 1], packed: rule n starts at n(n-1)/2.
struct Rules {
    const double* x;
    const double* w;
    int order;
};

RADIA_PEEC_HD inline int RuleOffset(int n) { return n * (n - 1) / 2; }

// Gauss-Legendre nodes and weights on [-1, 1] for n = 1..order, packed by
// RuleOffset (host only; the device receives the packed tables).
inline void GaussRules(int order, double* x, double* w) {
    const double pi = 3.141592653589793238462643383279502884;
    for (int n = 1; n <= order; ++n) {
        double* xn = x + RuleOffset(n);
        double* wn = w + RuleOffset(n);
        if (n == 1) { xn[0] = 0.0; wn[0] = 2.0; continue; }
        for (int i = 0; i < n; ++i) {
            double z = std::cos(pi * (i + 0.75) / (n + 0.5)), dp = 0.0;
            for (int it = 0; it < 100; ++it) {
                double p0 = 1.0, p1 = 0.0;
                for (int k = 1; k <= n; ++k) {
                    const double p2 = p1;
                    p1 = p0;
                    p0 = ((2 * k - 1) * z * p1 - (k - 1) * p2) / k;
                }
                dp = n * (z * p0 - p1) / (z * z - 1.0);
                const double dz = p0 / dp;
                z -= dz;
                if (std::abs(dz) < 1e-15) break;
            }
            xn[i] = -z;
            wn[i] = 2.0 / ((1.0 - z * z) * dp * dp);
        }
    }
}

RADIA_PEEC_HD inline double Dot(const V3& a, const V3& b) { return a.x * b.x + a.y * b.y + a.z * b.z; }

// ln(u + r) without cancellation when u < 0 (q = r^2 - u^2 > 0).
RADIA_PEEC_HD inline double LogUPlusR(double u, double r, double q) {
    return u >= 0.0 ? log(u + r) : log(q / (r - u));
}

// Corner primitive F of int dV / |r| over a rectangular prism.
RADIA_PEEC_HD inline double PrismCorner(double x, double y, double z) {
    const double x2 = x * x, y2 = y * y, z2 = z * z;
    const double r = sqrt(x2 + y2 + z2);
    if (r == 0.0) return 0.0;
    double t = 0.0;
    if (x != 0.0 && y != 0.0) t += x * y * LogUPlusR(z, r, x2 + y2);
    if (y != 0.0 && z != 0.0) t += y * z * LogUPlusR(x, r, y2 + z2);
    if (z != 0.0 && x != 0.0) t += z * x * LogUPlusR(y, r, z2 + x2);
    if (x != 0.0) t -= 0.5 * x2 * atan(y * z / (x * r));
    if (y != 0.0) t -= 0.5 * y2 * atan(z * x / (y * r));
    if (z != 0.0) t -= 0.5 * z2 * atan(x * y / (z * r));
    return t;
}

// Newtonian potential int_box dV / |p - r|, box bounds given relative to p.
RADIA_PEEC_HD inline double PrismPotential(const double lo[3], const double hi[3]) {
    double v = 0.0;
    for (int i = 0; i < 2; ++i)
        for (int j = 0; j < 2; ++j)
            for (int k = 0; k < 2; ++k) {
                const int lows = 3 - (i + j + k);
                const double sgn = (lows % 2 == 0) ? 1.0 : -1.0;
                v += sgn * PrismCorner(i ? hi[0] : lo[0], j ? hi[1] : lo[1], k ? hi[2] : lo[2]);
            }
    return v;
}

// Frame of a bar: e0 along the current, e1 the width axis (given, or the
// default normal of e0), e2 = e0 x e1; half extents (length, width, height)/2.
RADIA_PEEC_HD inline Frame FrameFromBar(const V3& c, const V3& d, const V3* width_dir, double length,
                                        double width, double height) {
    V3 ew;
    if (width_dir) {
        ew = *width_dir;
    } else {
        const V3 ref = fabs(d.x) < 0.9 ? V3{1, 0, 0} : V3{0, 1, 0};
        const double x = ref.y * d.z - ref.z * d.y, y = ref.z * d.x - ref.x * d.z,
                     z = ref.x * d.y - ref.y * d.x;
        const double n = sqrt(x * x + y * y + z * z);
        ew = V3{x / n, y / n, z / n};
    }
    Frame f;
    f.c = c;
    f.e[0] = d;
    f.e[1] = ew;
    f.e[2] = V3{d.y * ew.z - d.z * ew.y, d.z * ew.x - d.x * ew.z, d.x * ew.y - d.y * ew.x};
    f.half[0] = 0.5 * length;
    f.half[1] = 0.5 * width;
    f.half[2] = 0.5 * height;
    return f;
}

struct Cell {
    double u[3], h[3];
    int depth;
};

RADIA_PEEC_HD inline V3 CellCentre(const Frame& f, const double u[3]) {
    return V3{f.c.x + u[0] * f.e[0].x + u[1] * f.e[1].x + u[2] * f.e[2].x,
              f.c.y + u[0] * f.e[0].y + u[1] * f.e[1].y + u[2] * f.e[2].y,
              f.c.z + u[0] * f.e[0].z + u[1] * f.e[1].z + u[2] * f.e[2].z};
}

// True when the cell crosses the inner box surface: it overlaps the inner box
// with positive measure without lying inside it. Touching counts as separated.
RADIA_PEEC_HD inline bool Straddles(const Frame& fo, const Cell& cell, const Frame& fn, Trace* trace) {
    const V3 cw = CellCentre(fo, cell.u);
    bool contained = true;
    for (int s = 0; s < 8 && contained; ++s) {
        const double a = (s & 1) ? cell.h[0] : -cell.h[0], b = (s & 2) ? cell.h[1] : -cell.h[1],
                     c = (s & 4) ? cell.h[2] : -cell.h[2];
        const V3 p{cw.x + a * fo.e[0].x + b * fo.e[1].x + c * fo.e[2].x - fn.c.x,
                   cw.y + a * fo.e[0].y + b * fo.e[1].y + c * fo.e[2].y - fn.c.y,
                   cw.z + a * fo.e[0].z + b * fo.e[1].z + c * fo.e[2].z - fn.c.z};
        for (int k = 0; k < 3; ++k)
            if (fabs(Dot(p, fn.e[k])) > fn.half[k]) contained = false;
    }
    if (contained) {
        if (trace) ++(cell.depth == 0 ? trace->contained_top : trace->contained_descendants);
        return false;
    }
    const double tol =
        1e-12 * (sqrt(cell.h[0] * cell.h[0] + cell.h[1] * cell.h[1] + cell.h[2] * cell.h[2]) +
                 sqrt(fn.half[0] * fn.half[0] + fn.half[1] * fn.half[1] + fn.half[2] * fn.half[2]));
    const V3 d{fn.c.x - cw.x, fn.c.y - cw.y, fn.c.z - cw.z};
    auto separated = [&](const V3& ax) {
        double ra = 0.0, rb = 0.0;
        for (int k = 0; k < 3; ++k) {
            ra += cell.h[k] * fabs(Dot(fo.e[k], ax));
            rb += fn.half[k] * fabs(Dot(fn.e[k], ax));
        }
        return fabs(Dot(d, ax)) >= ra + rb - tol;
    };
    for (int k = 0; k < 3; ++k)
        if (separated(fo.e[k]) || separated(fn.e[k])) return false;
    for (int a = 0; a < 3; ++a)
        for (int b = 0; b < 3; ++b) {
            const V3& x = fo.e[a];
            const V3& y = fn.e[b];
            const V3 z{x.y * y.z - x.z * y.y, x.z * y.x - x.x * y.z, x.x * y.y - x.y * y.x};
            const double nz = sqrt(Dot(z, z));
            if (!(nz > 1e-12)) {
                if (trace) ++trace->cross_axes_skipped;
                continue;
            }
            if (trace) ++trace->cross_axes_evaluated;
            if (separated(V3{z.x / nz, z.y / nz, z.z / nz})) return false;
        }
    return true;
}

// Potential of the inner box at p (relative to its centre).
RADIA_PEEC_HD inline double InnerPotential(const Frame& fn, double px, double py, double pz) {
    double lo[3], hi[3];
    for (int k = 0; k < 3; ++k) {
        const double pk = px * fn.e[k].x + py * fn.e[k].y + pz * fn.e[k].z;
        lo[k] = -fn.half[k] - pk;
        hi[k] = fn.half[k] - pk;
    }
    return PrismPotential(lo, hi);
}

// Distance from p to the inner box surface (from outside or inside).
RADIA_PEEC_HD inline double SurfaceDistance(const Frame& fn, const V3& p) {
    const V3 r{p.x - fn.c.x, p.y - fn.c.y, p.z - fn.c.z};
    double out = 0.0, in = HUGE_VAL;
    for (int k = 0; k < 3; ++k) {
        const double q = fabs(Dot(r, fn.e[k]));
        const double o = fmax(q - fn.half[k], 0.0);
        out += o * o;
        in = fmin(in, fn.half[k] - q);
    }
    return out > 0.0 ? sqrt(out) : fmax(in, 0.0);
}

// Canonical orientation: the 1/r integral depends on the box as a point set,
// so reversing a bar flips the entry's sign exactly.
RADIA_PEEC_HD inline void Canonical(Frame& f) {
    const V3& d = f.e[0];
    const double lead = d.x != 0.0 ? d.x : (d.y != 0.0 ? d.y : d.z);
    if (lead < 0.0) {
        f.e[0] = V3{-d.x, -d.y, -d.z};
        f.e[2] = V3{-f.e[2].x, -f.e[2].y, -f.e[2].z};
    }
}

// Lexicographic order on (volume, centre, half extents, frame): the outer
// box is chosen geometrically, never by index.
RADIA_PEEC_HD inline bool KeyLess(const Frame& a, const Frame& b) {
    const double ka[16] = {a.half[0] * a.half[1] * a.half[2], a.c.x, a.c.y, a.c.z, a.half[0], a.half[1],
                           a.half[2], a.e[0].x, a.e[0].y, a.e[0].z, a.e[1].x, a.e[1].y, a.e[1].z,
                           a.e[2].x, a.e[2].y, a.e[2].z};
    const double kb[16] = {b.half[0] * b.half[1] * b.half[2], b.c.x, b.c.y, b.c.z, b.half[0], b.half[1],
                           b.half[2], b.e[0].x, b.e[0].y, b.e[0].z, b.e[1].x, b.e[1].y, b.e[1].z,
                           b.e[2].x, b.e[2].y, b.e[2].z};
    for (int k = 0; k < 16; ++k) {
        if (ka[k] < kb[k]) return true;
        if (kb[k] < ka[k]) return false;
    }
    return false;
}

RADIA_PEEC_HD inline double CellIntegral(const Frame& fo, const Frame& fn, const Cell& cell,
                                         const int n[3], const Rules& q) {
    const double *x0 = q.x + RuleOffset(n[0]), *w0 = q.w + RuleOffset(n[0]);
    const double *x1 = q.x + RuleOffset(n[1]), *w1 = q.w + RuleOffset(n[1]);
    const double *x2 = q.x + RuleOffset(n[2]), *w2 = q.w + RuleOffset(n[2]);
    double acc = 0.0;
    for (int a = 0; a < n[0]; ++a)
        for (int b = 0; b < n[1]; ++b)
            for (int c = 0; c < n[2]; ++c) {
                const double u = cell.u[0] + cell.h[0] * x0[a], v = cell.u[1] + cell.h[1] * x1[b],
                             w = cell.u[2] + cell.h[2] * x2[c];
                acc += w0[a] * w1[b] * w2[c] *
                       InnerPotential(fn, fo.c.x + u * fo.e[0].x + v * fo.e[1].x + w * fo.e[2].x - fn.c.x,
                                      fo.c.y + u * fo.e[0].y + v * fo.e[1].y + w * fo.e[2].y - fn.c.y,
                                      fo.c.z + u * fo.e[0].z + v * fo.e[1].z + w * fo.e[2].z - fn.c.z);
            }
    return acc;
}

// Partial inductance [H] between bars fi and fj (frames as built, before
// canonical orientation). Returns 0 for exactly orthogonal currents.
RADIA_PEEC_HD inline double Inductance(Frame fi, Frame fj, bool self, const Rules& q, Trace* trace) {
    const double dot = fi.e[0].x * fj.e[0].x + fi.e[0].y * fj.e[0].y + fi.e[0].z * fj.e[0].z;
    if (!self && dot == 0.0) return 0.0;
    Canonical(fi);
    Canonical(fj);
    const bool swap = KeyLess(fj, fi);
    const Frame& fo = swap ? fj : fi;  // outer (quadrature)
    const Frame& fn = swap ? fi : fj;  // inner (closed form)
    const double mu0_4pi = 1.0e-7;
    const double area_n = 4.0 * fn.half[1] * fn.half[2];
    const Cell top{{0.0, 0.0, 0.0}, {fo.half[0], fo.half[1], fo.half[2]}, 0};
    if (!self && Straddles(fo, top, fn, trace)) {
        if (trace) trace->composite = true;
        const int leaf = (q.order + 1) / 2 > 1 ? (q.order + 1) / 2 : 1;
        const double alpha = 0.5 * q.order;
        Cell stack[kStack];
        int sp = 0;
        stack[sp++] = top;
        double integral = 0.0;  // int_o phi_n dV
        while (sp > 0) {
            const Cell cell = stack[--sp];
            const bool cross = Straddles(fo, cell, fn, trace);
            if (cross && cell.depth < kSplitDepth) {
                if (trace) ++trace->splits;
                const double hm = fmax(cell.h[0], fmax(cell.h[1], cell.h[2]));
                int axes[3], na = 0;
                for (int k = 0; k < 3; ++k)
                    if (cell.h[k] >= 0.5 * hm) axes[na++] = k;
                for (int m = 0; m < (1 << na); ++m) {
                    Cell child = cell;
                    child.depth = cell.depth + 1;
                    for (int t = 0; t < na; ++t) {
                        const int k = axes[t];
                        child.h[k] = 0.5 * cell.h[k];
                        child.u[k] += ((m >> t) & 1) ? child.h[k] : -child.h[k];
                    }
                    stack[sp++] = child;
                }
                continue;
            }
            if (trace) {
                ++trace->leaves;
                if (cross) ++trace->straddling_leaves;
            }
            int n[3] = {leaf, leaf, leaf};
            if (!cross) {
                const double ds = SurfaceDistance(fn, CellCentre(fo, cell.u));
                if (ds > 0.0)
                    for (int k = 0; k < 3; ++k) {
                        const double m = fmin(static_cast<double>(leaf), ceil(alpha * 2.0 * cell.h[k] / ds));
                        n[k] = m > 1.0 ? static_cast<int>(m) : 1;
                    }
            }
            integral += CellIntegral(fo, fn, cell, n, q) * cell.h[0] * cell.h[1] * cell.h[2];
        }
        const double area_o = 4.0 * fo.half[1] * fo.half[2];
        return mu0_4pi * dot * integral / (area_o * area_n);
    }
    const double dx = fo.c.x - fn.c.x, dy = fo.c.y - fn.c.y, dz = fo.c.z - fn.c.z;
    const double dist = sqrt(dx * dx + dy * dy + dz * dz);
    int n[3];
    for (int k = 0; k < 3; ++k) {
        const double ext = 2.0 * fo.half[k];
        double m = static_cast<double>(q.order);
        if (!self && dist > 0.0) m = fmin(m, ceil(0.5 * q.order * ext / dist));
        n[k] = m > 1.0 ? static_cast<int>(m) : 1;
    }
    const double mean = CellIntegral(fo, fn, top, n, q) / 8.0;
    return mu0_4pi * (self ? 1.0 : dot) * (2.0 * fo.half[0]) * mean / area_n;
}

}  // namespace radia_peec_box
