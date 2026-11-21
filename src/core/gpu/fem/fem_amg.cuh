/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
// Classical AMG hierarchy built on the device: strength, parallel maximal
// independent set coarsening, direct interpolation, Galerkin products and the
// coarse inverse. The host sees control scalars only.
#include "fem_api.hpp"
#include "fem_device.cuh"

namespace radia_fem {

__global__ void k_strength(int n, const int* rp, const int* ci, const double* av, double theta, int* strong) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double max_neg = 0, max_abs = 0;
    for (int k = rp[i]; k < rp[i + 1]; ++k) {
        if (ci[k] == i) continue;
        max_neg = fmax(max_neg, -av[k]);
        max_abs = fmax(max_abs, fabs(av[k]));
    }
    for (int k = rp[i]; k < rp[i + 1]; ++k) {
        int s = 0;
        if (ci[k] != i) {
            if (max_neg > 0) s = -av[k] >= theta * max_neg;
            else if (max_abs > 0) s = fabs(av[k]) >= theta * max_abs;
        }
        strong[k] = s;
    }
}

// Union of the strength graph and its transpose. Columns are sorted per row.
__global__ void k_adjacent(int n, const int* rp, const int* ci, const int* strong, int* adjacent) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    for (int k = rp[i]; k < rp[i + 1]; ++k) {
        const int j = ci[k];
        int low = rp[j], high = rp[j + 1] - 1, reverse = 0;
        while (low <= high) {
            const int mid = (low + high) / 2;
            if (ci[mid] == i) { reverse = strong[mid]; break; }
            if (ci[mid] < i) low = mid + 1;
            else high = mid - 1;
        }
        adjacent[k] = strong[k] | reverse;
    }
}

__global__ void k_strong_degree(int n, const int* rp, const int* ci, const int* strong, int* degree) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    for (int k = rp[i]; k < rp[i + 1]; ++k)
        if (strong[k]) atomicAdd(degree + ci[k], 1);
}

// SplitMix64 gives the same tie-break on every device and host.
__host__ __device__ inline double hash01(uint64_t i) {
    uint64_t z = i + 0x9E3779B97F4A7C15ULL;
    z = (z ^ (z >> 30)) * 0xBF58476D1CE4E5B9ULL;
    z = (z ^ (z >> 27)) * 0x94D049BB133111EBULL;
    z ^= z >> 31;
    return double(z >> 11) * (1. / 9007199254740992.);
}

enum CoarseState { UNDECIDED = 0, COARSE = 1, FINE = 2, ISOLATED = 3 };

__global__ void k_pmis_init(int n, const int* rp, const int* adjacent, const int* degree, double* measure, int* state) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    measure[i] = degree[i] + hash01(uint64_t(i));
    int connected = 0;
    for (int k = rp[i]; k < rp[i + 1]; ++k) connected |= adjacent[k];
    state[i] = connected ? UNDECIDED : ISOLATED;
}
__global__ void k_pmis_select(int n, const int* rp, const int* ci, const int* adjacent, const double* measure,
                              const int* state, int* next) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    next[i] = state[i];
    if (state[i] != UNDECIDED) return;
    for (int k = rp[i]; k < rp[i + 1]; ++k) {
        const int j = ci[k];
        if (adjacent[k] && state[j] == UNDECIDED && (measure[j] > measure[i] || (measure[j] == measure[i] && j > i)))
            return;
    }
    next[i] = COARSE;
}
__global__ void k_pmis_mark(int n, const int* rp, const int* ci, const int* strong, const int* next, int* state, int* undecided) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int value = next[i];
    if (value == UNDECIDED) {
        for (int k = rp[i]; k < rp[i + 1]; ++k)
            if (strong[k] && next[ci[k]] == COARSE) { value = FINE; break; }
        if (value == UNDECIDED) atomicAdd(undecided, 1);
    }
    state[i] = value;
}
__global__ void k_pmis_finish(int n, const int* state, int* counts) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) counts[i] = state[i] == COARSE || state[i] == UNDECIDED;
}
__global__ void k_interpolation_counts(int n, const int* rp, const int* ci, const int* strong, const int* flags, const int* state, int* counts) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    if (flags[i]) { counts[i] = 1; return; }
    if (state[i] == ISOLATED) { counts[i] = 0; return; }
    int total = 0, any_coarse = 0;
    for (int k = rp[i]; k < rp[i + 1]; ++k) {
        if (ci[k] == i || !flags[ci[k]]) continue;
        any_coarse = 1;
        total += strong[k];
    }
    counts[i] = total ? total : any_coarse;
}
__global__ void k_interpolation_fill(int n, const int* rp, const int* ci, const double* av, const int* strong, const int* flags,
                                     const int* coarse_index, const int* prp, int* pci, double* pv) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n || prp[i + 1] == prp[i]) return;
    int position = prp[i];
    if (flags[i]) { pci[position] = coarse_index[i]; pv[position] = 1.; return; }
    double denominator = 0;
    int strong_coarse = 0;
    for (int k = rp[i]; k < rp[i + 1]; ++k) {
        const int j = ci[k];
        if (j == i) denominator += av[k];
        else if (strong[k] && flags[j]) ++strong_coarse;
        else denominator += av[k];
    }
    if (fabs(denominator) < 1e-30) denominator = 1.;
    if (strong_coarse) {
        for (int k = rp[i]; k < rp[i + 1]; ++k) {
            const int j = ci[k];
            if (j != i && strong[k] && flags[j]) {
                pci[position] = coarse_index[j];
                pv[position] = -av[k] / denominator;
                ++position;
            }
        }
        return;
    }
    int best = -1;
    double largest = -1;
    for (int k = rp[i]; k < rp[i + 1]; ++k) {
        const int j = ci[k];
        if (j != i && flags[j] && fabs(av[k]) > largest) { largest = fabs(av[k]); best = j; }
    }
    pci[position] = coarse_index[best];
    pv[position] = 1.;
}
__global__ void k_dense_clear(long long count, double* a) {
    const long long i = (long long)(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < count) a[i] = 0;
}
__global__ void k_dense_fill(int n, const int* rp, const int* ci, const double* av, double* a, double* diag) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double d = 0;
    for (int k = rp[i]; k < rp[i + 1]; ++k) {
        a[size_t(i) * n + ci[k]] = av[k];
        if (ci[k] == i) d = av[k];
    }
    diag[i] = d;
}
__global__ void k_gj_snapshot(int n, int pivot, const double* a, double* row, double* column) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    row[i] = a[size_t(pivot) * n + i];
    column[i] = a[size_t(i) * n + pivot];
}
// One Gauss-Jordan sweep. A pivot at or below the relative limit marks a kernel
// direction of a positive semidefinite coarse operator: that unknown stays zero.
__global__ void k_gj_sweep(int n, int pivot, double limit, const double* row, const double* column, const double* diag, double* a, int* pinned) {
    const long long index = (long long)(blockIdx.x) * blockDim.x + threadIdx.x;
    if (index >= (long long)(n) * n) return;
    const int i = int(index / n), j = int(index % n);
    const double p = row[pivot];
    if (!(p > limit * diag[pivot])) {
        if (i == pivot || j == pivot) a[index] = 0;
        if (index == 0) atomicAdd(pinned, 1);
        return;
    }
    const double d = 1. / p;
    if (i == pivot && j == pivot) a[index] = d;
    else if (i == pivot) a[index] = row[j] * d;
    else if (j == pivot) a[index] = -column[i] * d;
    else a[index] -= column[i] * row[j] * d;
}
__global__ void k_verify_dense(long long count, const double* a, int* bad) {
    const long long i = (long long)(blockIdx.x) * blockDim.x + threadIdx.x;
    if (i < count && !isfinite(a[i])) atomicExch(bad, 1);
}

struct AmgOptions {
    double theta = 0.25;
    int max_levels = 25;
    int min_coarse = 50;
    double stall = 0.;
    int dense_max = 2048;      // PROPOSED: the largest device dense coarse inverse.
    double pin_limit = 1e-12;  // PROPOSED: relative pivot limit of that inverse.
    int coarse_sweeps = 10;    // Above dense_max the coarsest level uses l1-Jacobi sweeps.
};

struct AmgLevel {
    Csr a, p, r;
    Dev<double> l1, x, b, residual, tmp;
};

class Amg {
public:
    AmgReport build(Ledger& ledger, Csr fine, const AmgOptions& options, const char* name) {
        options_ = options;
        AmgReport report;
        Csr current = std::move(fine);
        for (int level = 0;; ++level) {
            levels_.emplace_back();
            AmgLevel& l = levels_.back();
            l.a = std::move(current);
            allocate_level(l, name);
            report.rows.push_back(l.a.n);
            report.nonzeros.push_back(l.a.nnz);
            if (level + 1 >= options.max_levels || l.a.n <= options.min_coarse) break;
            int coarse_rows = 0;
            Csr p = coarsen(ledger, l.a, options.theta, coarse_rows, name);
            if (coarse_rows == 0 || coarse_rows >= l.a.n) break;
            if (options.stall > 0 && l.a.n <= options.min_coarse * 10 && coarse_rows > options.stall * l.a.n) break;
            Csr r = transpose(ledger, p.view(), p.nnz, name);
            Csr ap = multiply_csr(ledger, l.a.view(), p.view(), name);
            current = multiply_csr(ledger, r.view(), ap.view(), name);
            l.p = std::move(p);
            l.r = std::move(r);
        }
        AmgLevel& last = levels_.back();
        if (last.a.n <= options.dense_max) {
            factor_dense(ledger, last, name);
            report.coarse = "device_dense_inverse";
            report.pinned = pinned_;
        } else {
            report.coarse = "l1_jacobi_sweeps";
        }
        return report;
    }
    int rows() const { return levels_.front().a.n; }
    double* input() { return levels_.front().b.p; }
    double* output() { return levels_.front().x.p; }
    void cycle(size_t k = 0) {
        AmgLevel& l = levels_[k];
        const int n = l.a.n;
        const int blocks = blocks_of(n);
        if (k + 1 == levels_.size()) {
            if (dense_.p) {
                dense_coarse<<<blocks, TH>>>(n, dense_.p, l.b.p, l.x.p);
            } else {
                initial_relax<<<blocks, TH>>>(n, l.b.p, l.l1.p, l.x.p);
                for (int s = 1; s < options_.coarse_sweeps; ++s) relax(l);
            }
            return;
        }
        initial_relax<<<blocks, TH>>>(n, l.b.p, l.l1.p, l.x.p);
        multiply_residual<<<blocks, TH>>>(n, l.a.rowptr.p, l.a.columns.p, l.a.values.p, l.x.p, l.b.p, l.residual.p);
        AmgLevel& next = levels_[k + 1];
        multiply<<<blocks_of(next.a.n), TH>>>(next.a.n, l.r.rowptr.p, l.r.columns.p, l.r.values.p, l.residual.p, next.b.p);
        cycle(k + 1);
        multiply_add<<<blocks, TH>>>(n, l.p.rowptr.p, l.p.columns.p, l.p.values.p, next.x.p, 1., l.x.p);
        relax(l);
    }

private:
    static void allocate_level(AmgLevel& l, const char* name) {
        const size_t n = size_t(l.a.n);
        l.l1.allocate(n, name);
        l.x.allocate(n, name);
        l.b.allocate(n, name);
        l.residual.allocate(n, name);
        l.tmp.allocate(n, name);
        k_l1_norms<<<blocks_of(l.a.n), TH>>>(l.a.n, l.a.rowptr.p, l.a.values.p, l.l1.p);
        ck_kernels(name);
    }
    static void relax(AmgLevel& l) {
        multiply_relax<<<blocks_of(l.a.n), TH>>>(l.a.n, l.a.rowptr.p, l.a.columns.p, l.a.values.p, l.x.p, l.b.p, l.l1.p, l.tmp.p);
        std::swap(l.x.p, l.tmp.p);
    }
    static Csr coarsen(Ledger& ledger, const Csr& a, double theta, int& coarse_rows, const char* name) {
        const int n = a.n, nnz = a.nnz, blocks = blocks_of(n);
        const int* rp = a.rowptr.p;
        const int* ci = a.columns.p;
        Dev<int> strong(size_t(nnz), name), adjacent(size_t(nnz), name), degree(size_t(n), name), state(size_t(n), name), next(size_t(n), name);
        Dev<double> measure(size_t(n), name);
        Dev<int> undecided(1, name), flags(size_t(n), name), coarse_index;
        k_strength<<<blocks, TH>>>(n, rp, ci, a.values.p, theta, strong.p);
        k_adjacent<<<blocks, TH>>>(n, rp, ci, strong.p, adjacent.p);
        degree.zero();
        k_strong_degree<<<blocks, TH>>>(n, rp, ci, strong.p, degree.p);
        k_pmis_init<<<blocks, TH>>>(n, rp, adjacent.p, degree.p, measure.p, state.p);
        for (int iteration = 0; iteration < 100; ++iteration) {
            undecided.zero();
            k_pmis_select<<<blocks, TH>>>(n, rp, ci, adjacent.p, measure.p, state.p, next.p);
            k_pmis_mark<<<blocks, TH>>>(n, rp, ci, strong.p, next.p, state.p, undecided.p);
            if (read_scalar(ledger, undecided.p) == 0) break;
        }
        // Any row left undecided after the iteration cap becomes coarse.
        k_pmis_finish<<<blocks, TH>>>(n, state.p, flags.p);
        coarse_rows = rowptr_from_counts(ledger, flags.p, n, coarse_index, name);
        Csr p;
        if (coarse_rows == 0 || coarse_rows >= n) return p;
        Dev<int> counts(size_t(n), name);
        k_interpolation_counts<<<blocks, TH>>>(n, rp, ci, strong.p, flags.p, state.p, counts.p);
        p.n = n;
        p.m = coarse_rows;
        p.nnz = rowptr_from_counts(ledger, counts.p, n, p.rowptr, name);
        p.columns.allocate(size_t(p.nnz), name);
        p.values.allocate(size_t(p.nnz), name);
        k_interpolation_fill<<<blocks, TH>>>(n, rp, ci, a.values.p, strong.p, flags.p, coarse_index.p, p.rowptr.p, p.columns.p, p.values.p);
        ck_kernels(name);
        return p;
    }
    void factor_dense(Ledger& ledger, const AmgLevel& l, const char* name) {
        const int n = l.a.n;
        const long long count = (long long)(n) * n;
        dense_.allocate(size_t(count), name);
        Dev<double> diag(size_t(n), name), row(size_t(n), name), column(size_t(n), name);
        Dev<int> pinned(1, name), bad(1, name);
        pinned.zero();
        bad.zero();
        k_dense_clear<<<blocks_of(count), TH>>>(count, dense_.p);
        k_dense_fill<<<blocks_of(n), TH>>>(n, l.a.rowptr.p, l.a.columns.p, l.a.values.p, dense_.p, diag.p);
        for (int k = 0; k < n; ++k) {
            k_gj_snapshot<<<blocks_of(n), TH>>>(n, k, dense_.p, row.p, column.p);
            k_gj_sweep<<<blocks_of(count), TH>>>(n, k, options_.pin_limit, row.p, column.p, diag.p, dense_.p, pinned.p);
        }
        k_verify_dense<<<blocks_of(count), TH>>>(count, dense_.p, bad.p);
        ck_kernels(name);
        check_flag(ledger, bad, "device AMG: the dense coarse inverse is not finite");
        pinned_ = read_scalar(ledger, pinned.p);
    }
    AmgOptions options_;
    std::vector<AmgLevel> levels_;
    Dev<double> dense_;
    int pinned_ = 0;
};

}  // namespace radia_fem
