// PRIMA reduction of a dense branch PEEC port model (see rad_peec_prima.h).
#include "rad_peec_prima.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <string>

#ifdef HAVE_LAPACK
#include "rad_blas.h"
#if defined(RADIA_BLAS_MKL)
#include <mkl_lapack.h>
#define RADIA_PRIMA_POTRF dpotrf_
#define RADIA_PRIMA_POTRS dpotrs_
#define RADIA_PRIMA_SYEV dsyev_
#elif defined(RADIA_BLAS_NVPL)
#include <nvpl_lapack.h>
#define RADIA_PRIMA_POTRF NVPL_LAPACK_dpotrf
#define RADIA_PRIMA_POTRS NVPL_LAPACK_dpotrs
#define RADIA_PRIMA_SYEV NVPL_LAPACK_dsyev
#endif
#endif

namespace radia {
namespace peec_prima {
namespace {

using Vec = std::vector<double>;

double Dot(const Vec& a, const Vec& b) {
    double s = 0.0;
    for (std::size_t i = 0; i < a.size(); ++i) s += a[i] * b[i];
    return s;
}
double Norm(const Vec& a) { return std::sqrt(Dot(a, a)); }

void Validate(const PortModelInput& in) {
    const std::size_t N = in.n_branch;
    if (N == 0) throw std::invalid_argument("peec_prima: no branches");
    if (in.L.size() != N * N) throw std::invalid_argument("peec_prima: L must be n_branch x n_branch");
    if (in.R.size() != N) throw std::invalid_argument("peec_prima: R must have n_branch entries");
    if (in.branch_nodes.size() != 2 * N)
        throw std::invalid_argument("peec_prima: branch_nodes must be n_branch x 2");
    for (std::size_t i = 0; i < N; ++i) {
        if (!(std::isfinite(in.R[i]) && in.R[i] > 0.0))
            throw std::invalid_argument("peec_prima: R[" + std::to_string(i) + "] must be finite and > 0");
        for (std::size_t j = 0; j < N; ++j) {
            const double v = in.L[i * N + j];
            if (!std::isfinite(v)) throw std::invalid_argument("peec_prima: L has a non-finite entry");
            if (v != in.L[j * N + i])
                throw std::invalid_argument("peec_prima: L is not exactly symmetric at ("
                                            + std::to_string(i) + ", " + std::to_string(j) + ")");
        }
    }
    const auto nn = static_cast<long long>(in.n_nodes);
    for (std::size_t b = 0; b < N; ++b) {
        const int u = in.branch_nodes[2 * b], v = in.branch_nodes[2 * b + 1];
        if (u < 0 || v < 0 || u >= nn || v >= nn)
            throw std::invalid_argument("peec_prima: branch " + std::to_string(b) + " node out of range");
        if (u == v) throw std::invalid_argument("peec_prima: branch " + std::to_string(b) + " is a self-loop");
    }
    if (in.port_plus.empty() || in.port_plus.size() != in.port_minus.size())
        throw std::invalid_argument("peec_prima: need >= 1 port with matching plus/minus lists");
    for (std::size_t p = 0; p < in.port_plus.size(); ++p) {
        const int a = in.port_plus[p], c = in.port_minus[p];
        if (a < 0 || c < 0 || a >= nn || c >= nn || a == c)
            throw std::invalid_argument("peec_prima: port " + std::to_string(p) + " has invalid nodes");
    }
}

}  // namespace

#ifndef HAVE_LAPACK
ReducedModel Reduce(const PortModelInput&, std::size_t, double) {
    throw std::runtime_error("peec_prima::Reduce requires LAPACK (build without HAVE_LAPACK)");
}
#else

ReducedModel Reduce(const PortModelInput& in, std::size_t order, double s0) {
    Validate(in);
    if (order == 0) throw std::invalid_argument("peec_prima: order must be >= 1");
    if (!(std::isfinite(s0) && s0 >= 0.0))
        throw std::invalid_argument("peec_prima: expansion point s0 must be real, finite and >= 0");

    const std::size_t N = in.n_branch;
    const int ground = in.port_minus[0];
    std::vector<int> idx(in.n_nodes, -1);
    std::size_t nr = 0;
    for (std::size_t k = 0; k < in.n_nodes; ++k)
        if (static_cast<int>(k) != ground) idx[k] = static_cast<int>(nr++);
    const std::size_t dim = N + nr;
    const std::size_t P = in.port_plus.size();
    if (order > dim) order = dim;

    // Branch ends in reduced node numbering (-1 = ground): A[u, b] = -1, A[v, b] = +1.
    std::vector<int> bu(N), bv(N);
    for (std::size_t b = 0; b < N; ++b) {
        bu[b] = idx[in.branch_nodes[2 * b]];
        bv[b] = idx[in.branch_nodes[2 * b + 1]];
    }
    auto apply_At = [&](const double* V, double* out) {   // out(N) = A^T V
        for (std::size_t b = 0; b < N; ++b)
            out[b] = (bv[b] >= 0 ? V[bv[b]] : 0.0) - (bu[b] >= 0 ? V[bu[b]] : 0.0);
    };
    auto apply_A = [&](const double* I, double* out) {    // out(nr) = A I
        std::fill(out, out + nr, 0.0);
        for (std::size_t b = 0; b < N; ++b) {
            if (bv[b] >= 0) out[bv[b]] += I[b];
            if (bu[b] >= 0) out[bu[b]] -= I[b];
        }
    };
    auto apply_L = [&](const double* x, double* out) {
        for (std::size_t i = 0; i < N; ++i) {
            const double* row = in.L.data() + i * N;
            double s = 0.0;
            for (std::size_t j = 0; j < N; ++j) s += row[j] * x[j];
            out[i] = s;
        }
    };

    // Full magnetic energy admission: L itself must be positive definite.
    // The reduced-E check below cannot see an energy hole of L that the
    // Krylov basis misses (e.g. a zero-incidence loop with negative energy).
    int n = static_cast<int>(N), info = 0;
    char uplo = 'L';
    {
        std::vector<double> Lc(in.L);
        RADIA_PRIMA_POTRF(&uplo, &n, Lc.data(), &n, &info);
        if (info != 0)
            throw std::runtime_error("peec_prima: L is not positive definite (full magnetic energy "
                                     "not admissible, dpotrf info " + std::to_string(info) + ")");
    }

    // Cholesky of K = R + s0 L (symmetric, so the row-major buffer is its own transpose).
    std::vector<double> K(N * N);
    for (std::size_t i = 0; i < N * N; ++i) K[i] = s0 * in.L[i];
    for (std::size_t i = 0; i < N; ++i) K[i * N + i] += in.R[i];
    RADIA_PRIMA_POTRF(&uplo, &n, K.data(), &n, &info);
    if (info != 0)
        throw std::runtime_error("peec_prima: R + s0 L is not positive definite (dpotrf info "
                                 + std::to_string(info) + ")");

    // Schur complement S = A K^{-1} A^T through Y = K^{-1} A^T (column-major N x nr).
    std::vector<double> Y(N * nr, 0.0);
    for (std::size_t b = 0; b < N; ++b) {
        if (bv[b] >= 0) Y[b + N * bv[b]] += 1.0;
        if (bu[b] >= 0) Y[b + N * bu[b]] -= 1.0;
    }
    int nrhs = static_cast<int>(nr);
    if (nr == 0) throw std::invalid_argument("peec_prima: the model has no non-ground node");
    RADIA_PRIMA_POTRS(&uplo, &n, &nrhs, K.data(), &n, Y.data(), &n, &info);
    if (info != 0) throw std::runtime_error("peec_prima: dpotrs info " + std::to_string(info));
    std::vector<double> S(nr * nr, 0.0);
    for (std::size_t c = 0; c < nr; ++c) {
        const double* y = Y.data() + N * c;
        double* s = S.data() + nr * c;
        for (std::size_t b = 0; b < N; ++b) {
            if (bv[b] >= 0) s[bv[b]] += y[b];
            if (bu[b] >= 0) s[bu[b]] -= y[b];
        }
    }
    Y.clear(); Y.shrink_to_fit();
    int m = static_cast<int>(nr);
    RADIA_PRIMA_POTRF(&uplo, &m, S.data(), &m, &info);
    if (info != 0)
        throw std::runtime_error("peec_prima: Schur complement A K^-1 A^T is not positive definite "
                                 "(dpotrf info " + std::to_string(info)
                                 + "); is the branch graph connected?");

    double worst_residual = 0.0;
    // Solve M w = [f; g], M = [[K, A^T], [A, 0]], with one refinement step when needed.
    auto saddle_once = [&](const Vec& rhs, Vec& w) {
        Vec t(rhs.begin(), rhs.begin() + N);
        int one = 1;
        RADIA_PRIMA_POTRS(&uplo, &n, &one, K.data(), &n, t.data(), &n, &info);   // K^-1 f
        if (info != 0) throw std::runtime_error("peec_prima: dpotrs info " + std::to_string(info));
        Vec Vn(nr);
        apply_A(t.data(), Vn.data());
        for (std::size_t r = 0; r < nr; ++r) Vn[r] -= rhs[N + r];
        RADIA_PRIMA_POTRS(&uplo, &m, &one, S.data(), &m, Vn.data(), &m, &info);
        if (info != 0) throw std::runtime_error("peec_prima: dpotrs info " + std::to_string(info));
        Vec AtV(N);
        apply_At(Vn.data(), AtV.data());
        for (std::size_t b = 0; b < N; ++b) AtV[b] = rhs[b] - AtV[b];
        RADIA_PRIMA_POTRS(&uplo, &n, &one, K.data(), &n, AtV.data(), &n, &info);
        if (info != 0) throw std::runtime_error("peec_prima: dpotrs info " + std::to_string(info));
        w.assign(dim, 0.0);
        std::copy(AtV.begin(), AtV.end(), w.begin());
        std::copy(Vn.begin(), Vn.end(), w.begin() + N);
    };
    auto residual = [&](const Vec& w, const Vec& rhs, Vec& r) {   // r = rhs - M w
        r.assign(dim, 0.0);
        Vec Lw(N), AtV(N), Aw(nr);
        apply_L(w.data(), Lw.data());
        apply_At(w.data() + N, AtV.data());
        apply_A(w.data(), Aw.data());
        for (std::size_t b = 0; b < N; ++b)
            r[b] = rhs[b] - (in.R[b] * w[b] + s0 * Lw[b] + AtV[b]);
        for (std::size_t k = 0; k < nr; ++k) r[N + k] = rhs[N + k] - Aw[k];
    };
    auto solve_saddle = [&](const Vec& rhs) {
        Vec w, r, dw;
        saddle_once(rhs, w);
        residual(w, rhs, r);
        const double rn = Norm(rhs);
        double rel = Norm(r) / std::max(rn, std::numeric_limits<double>::min());
        if (rel > 1e-12) {
            saddle_once(r, dw);
            for (std::size_t i = 0; i < dim; ++i) w[i] += dw[i];
            residual(w, rhs, r);
            rel = Norm(r) / std::max(rn, std::numeric_limits<double>::min());
        }
        if (!(rel <= RELATIVE_LIMIT))
            throw std::runtime_error("peec_prima: shifted saddle solve true residual "
                                     + std::to_string(rel) + " exceeds 1e-6");
        worst_residual = std::max(worst_residual, rel);
        return w;
    };
    // (G + s0 E)^{-1} r with G + s0 E = diag(I, -I) M.
    auto solve_pencil = [&](Vec r) {
        for (std::size_t k = N; k < dim; ++k) r[k] = -r[k];
        return solve_saddle(r);
    };

    // Port columns of B = [0; P] in reduced numbering.
    std::vector<Vec> Bcol(P, Vec(dim, 0.0));
    for (std::size_t p = 0; p < P; ++p) {
        if (idx[in.port_plus[p]] >= 0) Bcol[p][N + idx[in.port_plus[p]]] += 1.0;
        if (idx[in.port_minus[p]] >= 0) Bcol[p][N + idx[in.port_minus[p]]] -= 1.0;
    }

    // Block Arnoldi, one column at a time (band form): candidates are the
    // start block and then the operator applied to each accepted column in
    // order; Gram-Schmidt twice; a candidate whose remainder falls below
    // DEFLATION of its own norm is deflated.
    constexpr double DEFLATION = 1e-10;
    ReducedModel out;
    std::vector<Vec> basis;
    auto accept = [&](Vec w, int source) {
        const double before = Norm(w);
        if (!(before > 0.0)) { out.deflated.push_back(source); return; }
        for (int pass = 0; pass < 2; ++pass)
            for (const Vec& v : basis) {
                const double h = Dot(w, v);
                for (std::size_t i = 0; i < dim; ++i) w[i] -= h * v[i];
            }
        const double after = Norm(w);
        if (after <= DEFLATION * before) { out.deflated.push_back(source); return; }
        for (double& x : w) x /= after;
        basis.push_back(std::move(w));
    };
    for (std::size_t p = 0; p < P; ++p)
        accept(solve_pencil(Bcol[p]), -1 - static_cast<int>(p));
    if (basis.size() > order)
        throw std::invalid_argument("peec_prima: order " + std::to_string(order)
                                    + " is below the independent port start-block rank "
                                    + std::to_string(basis.size()));
    for (std::size_t j = 0; j < basis.size() && basis.size() < order; ++j) {
        Vec Ev(dim, 0.0);
        apply_L(basis[j].data(), Ev.data());
        accept(solve_pencil(Ev), static_cast<int>(j));
    }
    const std::size_t q = basis.size();
    if (q == 0) throw std::runtime_error("peec_prima: empty Krylov basis (zero port excitation)");

    // Congruence projection.
    out.q = q; out.n_ports = P; out.n_full = dim;
    out.V.assign(dim * q, 0.0);
    for (std::size_t c = 0; c < q; ++c)
        for (std::size_t i = 0; i < dim; ++i) out.V[i * q + c] = basis[c][i];
    std::vector<Vec> EV(q, Vec(dim, 0.0)), GV(q, Vec(dim, 0.0));
    for (std::size_t c = 0; c < q; ++c) {
        const Vec& v = basis[c];
        apply_L(v.data(), EV[c].data());
        Vec AtV(N), Av(nr);
        apply_At(v.data() + N, AtV.data());
        apply_A(v.data(), Av.data());
        for (std::size_t b = 0; b < N; ++b) GV[c][b] = in.R[b] * v[b] + AtV[b];
        for (std::size_t k = 0; k < nr; ++k) GV[c][N + k] = -Av[k];
    }
    out.E.assign(q * q, 0.0); out.G.assign(q * q, 0.0); out.B.assign(q * P, 0.0);
    for (std::size_t r = 0; r < q; ++r) {
        for (std::size_t c = 0; c < q; ++c) {
            out.E[r * q + c] = Dot(basis[r], EV[c]);
            out.G[r * q + c] = Dot(basis[r], GV[c]);
        }
        for (std::size_t p = 0; p < P; ++p) out.B[r * P + p] = Dot(basis[r], Bcol[p]);
    }
    for (std::size_t r = 0; r < q; ++r)
        for (std::size_t c = r + 1; c < q; ++c) {
            const double a = 0.5 * (out.E[r * q + c] + out.E[c * q + r]);
            out.E[r * q + c] = out.E[c * q + r] = a;
        }

    // Reduced magnetic energy must be >= 0 up to roundoff.
    std::vector<double> Ew(out.E), w(q);
    int iq = static_cast<int>(q), lwork = -1;
    char jobz = 'N';
    double wq = 0.0;
    RADIA_PRIMA_SYEV(&jobz, &uplo, &iq, Ew.data(), &iq, w.data(), &wq, &lwork, &info);
    lwork = std::max(1, static_cast<int>(wq));
    std::vector<double> work(lwork);
    RADIA_PRIMA_SYEV(&jobz, &uplo, &iq, Ew.data(), &iq, w.data(), work.data(), &lwork, &info);
    if (info != 0) throw std::runtime_error("peec_prima: dsyev info " + std::to_string(info));
    const double scale = std::max(std::fabs(w.front()), std::fabs(w.back()));
    const double tol = static_cast<double>(dim) * std::numeric_limits<double>::epsilon() * scale;
    if (w.front() < -tol)
        throw std::runtime_error("peec_prima: reduced inductance has negative energy (min eigenvalue "
                                 + std::to_string(w.front()) + " H, roundoff scale "
                                 + std::to_string(tol) + " H); is L positive semidefinite?");
    out.min_E_eig_over_max = scale > 0.0 ? w.front() / scale : 0.0;
    out.max_solve_residual = worst_residual;
    return out;
}

#endif  // HAVE_LAPACK

}  // namespace peec_prima
}  // namespace radia
