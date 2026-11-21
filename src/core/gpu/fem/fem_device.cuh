/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
// Device primitives of the FEM route: owned arrays, the byte ledger, sorted
// triple plans and CSR algebra. No MATLAB, Python or NGSolve types, no x86
// intrinsics or MKL: this compiles for sm_75/sm_80 with MSVC and sm_100 with g++.
// Everything that crosses to or from the host goes through the ledger.
#include <cub/device/device_radix_sort.cuh>
#include <cub/device/device_run_length_encode.cuh>
#include <cub/device/device_scan.cuh>
#include <cuda_runtime.h>
#include <algorithm>
#include <climits>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

namespace radia_fem {

constexpr int TH = 256;
inline int blocks_of(long long n) { return int((std::max<long long>(n, 1) + TH - 1) / TH); }

inline void ck(cudaError_t status, const char* where) {
    if (status != cudaSuccess)
        throw std::runtime_error(std::string(where) + ": " + cudaGetErrorString(status));
}
inline void ck_kernels(const char* where) { ck(cudaGetLastError(), where); }

struct DeviceMemory {
    size_t live = 0, peak = 0;
};
inline DeviceMemory& device_memory() { static DeviceMemory memory; return memory; }

// Every byte that crosses is recorded by payload and direction.
class Ledger {
public:
    void add(const std::string& payload, const char* direction, uint64_t bytes) {
        rows_[{payload, direction}] += bytes;
    }
    const std::map<std::pair<std::string, std::string>, uint64_t>& rows() const { return rows_; }
private:
    std::map<std::pair<std::string, std::string>, uint64_t> rows_;
};

// The same reserve as the other resident device routes.
inline void check_reserve(size_t bytes, const char* name) {
    size_t available = 0, total = 0;
    ck(cudaMemGetInfo(&available, &total), "cudaMemGetInfo");
    const size_t reserve = std::max<size_t>(2ULL << 30, total / 5);
    if (available <= reserve || bytes > available - reserve)
        throw std::runtime_error(std::string("device memory budget exceeded by ") + name + ": requested "
            + std::to_string(bytes) + " bytes, free " + std::to_string(available) + ", reserve " + std::to_string(reserve));
}

template <class T> class Dev {
public:
    Dev() = default;
    Dev(size_t count, const char* name) { allocate(count, name); }
    Dev(const Dev&) = delete;
    Dev& operator=(const Dev&) = delete;
    Dev(Dev&& other) noexcept : p(other.p), n(other.n) { other.p = nullptr; other.n = 0; }
    Dev& operator=(Dev&& other) noexcept {
        if (this != &other) { release(); p = other.p; n = other.n; other.p = nullptr; other.n = 0; }
        return *this;
    }
    ~Dev() { release(); }
    void allocate(size_t count, const char* name) {
        release();
        const size_t bytes = std::max<size_t>(count, 1) * sizeof(T);
        check_reserve(bytes, name);
        ck(cudaMalloc(&p, bytes), name);
        n = count;
        auto& memory = device_memory();
        memory.live += bytes;
        memory.peak = std::max(memory.peak, memory.live);
    }
    void release() noexcept {
        if (p) {
            cudaFree(p);
            device_memory().live -= std::max<size_t>(n, 1) * sizeof(T);
        }
        p = nullptr;
        n = 0;
    }
    void zero() { ck(cudaMemset(p, 0, std::max<size_t>(n, 1) * sizeof(T)), "cudaMemset"); }
    T* p = nullptr;
    size_t n = 0;
};

template <class T> Dev<T> upload(Ledger& ledger, const char* payload, const T* host, size_t count) {
    Dev<T> result(count, payload);
    if (count) ck(cudaMemcpy(result.p, host, count * sizeof(T), cudaMemcpyHostToDevice), payload);
    ledger.add(payload, "H2D", count * sizeof(T));
    return result;
}
template <class T> void download(Ledger& ledger, const char* payload, const T* device, T* host, size_t count) {
    if (count) ck(cudaMemcpy(host, device, count * sizeof(T), cudaMemcpyDeviceToHost), payload);
    ledger.add(payload, "D2H", count * sizeof(T));
}
// One control scalar for host-side flow, never a matrix entry.
template <class T> T read_scalar(Ledger& ledger, const T* device) {
    T value{};
    ck_kernels("control scalar");
    ck(cudaMemcpy(&value, device, sizeof(T), cudaMemcpyDeviceToHost), "control scalar");
    ledger.add("scalar.control", "D2H", sizeof(T));
    return value;
}
inline void check_flag(Ledger& ledger, const Dev<int>& flag, const char* message) {
    if (read_scalar(ledger, flag.p)) throw std::runtime_error(message);
}

inline int bits_for(long long values) {
    int bits = 1;
    while ((1LL << bits) < values) ++bits;
    return bits;
}

__global__ void k_widen(int n, const int* counts, long long* out) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) out[i] = counts[i];
    else if (i == n) out[i] = 0;
}
__global__ void k_narrow(int n, const long long* in, int* out) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i <= n) out[i] = int(in[i]);
}
// counts -> row pointers on the device; only the total crosses as a control scalar.
inline int rowptr_from_counts(Ledger& ledger, const int* counts, int n, Dev<int>& rowptr, const char* name) {
    Dev<long long> wide(size_t(n) + 1, name), sums(size_t(n) + 1, name);
    k_widen<<<blocks_of(n + 1), TH>>>(n, counts, wide.p);
    size_t bytes = 0;
    ck(cub::DeviceScan::ExclusiveSum(nullptr, bytes, wide.p, sums.p, n + 1), name);
    Dev<unsigned char> scratch(bytes, name);
    ck(cub::DeviceScan::ExclusiveSum(scratch.p, bytes, wide.p, sums.p, n + 1), name);
    const long long total = read_scalar(ledger, sums.p + n);
    if (total < 0 || total > INT_MAX - TH)
        throw std::runtime_error(std::string(name) + ": row pointer exceeds int32 capacity, total " + std::to_string(total));
    rowptr.allocate(size_t(n) + 1, name);
    k_narrow<<<blocks_of(n + 1), TH>>>(n, sums.p, rowptr.p);
    ck_kernels(name);
    return int(total);
}

__global__ void k_iota(int n, int* out) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) out[i] = i;
}
__global__ void k_decode_keys(int nu, int colbits, const uint64_t* keys, int* counts, int* columns) {
    int u = blockIdx.x * blockDim.x + threadIdx.x;
    if (u >= nu) return;
    atomicAdd(counts + int(keys[u] >> colbits), 1);
    columns[u] = int(keys[u] & ((uint64_t(1) << colbits) - 1));
}

// Triples sorted by (row, column). order[k] names the source slot of sorted
// entry k; slots starts[u]..starts[u+1] of order sum into nonzero u. The row
// and column come from the key, so the pattern is built on the device only.
struct Plan {
    int rows = 0, cols = 0, triples = 0, nnz = 0;
    Dev<int> order, starts, rowptr, columns;
};
inline Plan make_plan(Ledger& ledger, Dev<uint64_t>& keys, int triples, int rows, int cols, const char* name) {
    Plan plan;
    plan.rows = rows; plan.cols = cols; plan.triples = triples;
    const int colbits = bits_for(cols), totalbits = colbits + bits_for(rows);
    if (totalbits > 62) throw std::runtime_error(std::string(name) + ": key width exceeded");
    Dev<int> index(triples, name);
    plan.order.allocate(triples, name);
    Dev<uint64_t> sorted(triples, name);
    k_iota<<<blocks_of(triples), TH>>>(triples, index.p);
    size_t bytes = 0;
    ck(cub::DeviceRadixSort::SortPairs(nullptr, bytes, keys.p, sorted.p, index.p, plan.order.p, triples, 0, totalbits), name);
    {
        Dev<unsigned char> scratch(bytes, name);
        ck(cub::DeviceRadixSort::SortPairs(scratch.p, bytes, keys.p, sorted.p, index.p, plan.order.p, triples, 0, totalbits), name);
    }
    index.release();
    Dev<uint64_t> unique(triples, name);
    Dev<int> runs(triples, name), count(1, name);
    bytes = 0;
    ck(cub::DeviceRunLengthEncode::Encode(nullptr, bytes, sorted.p, unique.p, runs.p, count.p, triples), name);
    {
        Dev<unsigned char> scratch(bytes, name);
        ck(cub::DeviceRunLengthEncode::Encode(scratch.p, bytes, sorted.p, unique.p, runs.p, count.p, triples), name);
    }
    plan.nnz = read_scalar(ledger, count.p);
    int check = rowptr_from_counts(ledger, runs.p, plan.nnz, plan.starts, name);
    if (check != triples) throw std::runtime_error(std::string(name) + ": run lengths do not cover the triples");
    Dev<int> rowcount(rows, name);
    rowcount.zero();
    plan.columns.allocate(plan.nnz, name);
    k_decode_keys<<<blocks_of(plan.nnz), TH>>>(plan.nnz, colbits, unique.p, rowcount.p, plan.columns.p);
    check = rowptr_from_counts(ledger, rowcount.p, rows, plan.rowptr, name);
    if (check != plan.nnz) throw std::runtime_error(std::string(name) + ": row counts do not cover the pattern");
    ck_kernels(name);
    return plan;
}

template <class S>
__global__ void k_reduce(int nnz, const int* starts, const int* order, const S* source, S* out) {
    int u = blockIdx.x * blockDim.x + threadIdx.x;
    if (u >= nnz) return;
    S sum = source[order[starts[u]]];
    for (int k = starts[u] + 1; k < starts[u + 1]; ++k) sum += source[order[k]];
    out[u] = sum;
}
template <class S> void reduce(const Plan& plan, const S* source, S* out) {
    k_reduce<S><<<blocks_of(plan.nnz), TH>>>(plan.nnz, plan.starts.p, plan.order.p, source, out);
}

struct CsrView {
    int n = 0, m = 0;
    const int* rowptr = nullptr;
    const int* columns = nullptr;
    const double* values = nullptr;
};
struct Csr {
    int n = 0, m = 0, nnz = 0;
    Dev<int> rowptr, columns;
    Dev<double> values;
    CsrView view() const { return CsrView{n, m, rowptr.p, columns.p, values.p}; }
};

// Shared AMS/AMG primitives (multiply, relax, residual, mask, dense action).
#include "../ams_amg_kernels.cuh"

__global__ void k_row_of_entry(int n, const int* rowptr, int* rowof) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    for (int k = rowptr[i]; k < rowptr[i + 1]; ++k) rowof[k] = i;
}
__global__ void k_count_columns(int nnz, const int* columns, int* counts) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < nnz) atomicAdd(counts + columns[k], 1);
}
__global__ void k_permute_int(int count, const int* order, const int* source, int* out) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < count) out[k] = source[order[k]];
}
__global__ void k_permute_double(int count, const int* order, const double* source, double* out) {
    int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < count) out[k] = source[order[k]];
}

// Stable sort by column keeps each transposed row ordered.
struct TransposeMap {
    Dev<int> order, rowptr, columns;
    int nnz = 0;
};
inline TransposeMap transpose_map(Ledger& ledger, int n, int m, int nnz, const int* rowptr, const int* columns, const char* name) {
    TransposeMap map;
    map.nnz = nnz;
    Dev<int> keys(nnz, name), sortedkeys(nnz, name), index(nnz, name), rowof(nnz, name), counts(m, name);
    map.order.allocate(nnz, name);
    k_iota<<<blocks_of(nnz), TH>>>(nnz, index.p);
    k_row_of_entry<<<blocks_of(n), TH>>>(n, rowptr, rowof.p);
    ck(cudaMemcpy(keys.p, columns, size_t(nnz) * sizeof(int), cudaMemcpyDeviceToDevice), name);
    size_t bytes = 0;
    const int bits = bits_for(m);
    ck(cub::DeviceRadixSort::SortPairs(nullptr, bytes, keys.p, sortedkeys.p, index.p, map.order.p, nnz, 0, bits), name);
    {
        Dev<unsigned char> scratch(bytes, name);
        ck(cub::DeviceRadixSort::SortPairs(scratch.p, bytes, keys.p, sortedkeys.p, index.p, map.order.p, nnz, 0, bits), name);
    }
    counts.zero();
    k_count_columns<<<blocks_of(nnz), TH>>>(nnz, columns, counts.p);
    int total = rowptr_from_counts(ledger, counts.p, m, map.rowptr, name);
    if (total != nnz) throw std::runtime_error(std::string(name) + ": transpose lost entries");
    map.columns.allocate(nnz, name);
    k_permute_int<<<blocks_of(nnz), TH>>>(nnz, map.order.p, rowof.p, map.columns.p);
    ck_kernels(name);
    return map;
}
inline Csr transpose(Ledger& ledger, const CsrView& a, int nnz, const char* name) {
    TransposeMap map = transpose_map(ledger, a.n, a.m, nnz, a.rowptr, a.columns, name);
    Csr t;
    t.n = a.m; t.m = a.n; t.nnz = nnz;
    t.values.allocate(nnz, name);
    k_permute_double<<<blocks_of(nnz), TH>>>(nnz, map.order.p, a.values, t.values.p);
    t.rowptr = std::move(map.rowptr);
    t.columns = std::move(map.columns);
    ck_kernels(name);
    return t;
}

// Expand-sort-compress product on the device. An exact zero of A contributes
// nothing, so identity rows of a constrained operator stay decoupled.
__global__ void k_product_counts(int n, const int* arp, const int* aci, const double* av, const int* brp, int* counts) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int total = 0;
    for (int k = arp[i]; k < arp[i + 1]; ++k)
        if (av[k] != 0.) total += brp[aci[k] + 1] - brp[aci[k]];
    counts[i] = total;
}
__global__ void k_product_expand(int n, int colbits, const int* offsets, const int* arp, const int* aci, const double* av,
                                 const int* brp, const int* bci, const double* bv, uint64_t* keys, double* values) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    int position = offsets[i];
    for (int k = arp[i]; k < arp[i + 1]; ++k) {
        const double a = av[k];
        if (a == 0.) continue;
        for (int l = brp[aci[k]]; l < brp[aci[k] + 1]; ++l, ++position) {
            keys[position] = (uint64_t(i) << colbits) | uint64_t(bci[l]);
            values[position] = a * bv[l];
        }
    }
}
__global__ void k_segment_sum(int nnz, const int* starts, const double* values, double* out) {
    int u = blockIdx.x * blockDim.x + threadIdx.x;
    if (u >= nnz) return;
    double sum = 0;
    for (int k = starts[u]; k < starts[u + 1]; ++k) sum += values[k];
    out[u] = sum;
}
inline Csr multiply_csr(Ledger& ledger, const CsrView& a, const CsrView& b, const char* name) {
    if (a.m != b.n) throw std::runtime_error(std::string(name) + ": product dimensions disagree");
    Dev<int> counts(a.n, name), offsets;
    k_product_counts<<<blocks_of(a.n), TH>>>(a.n, a.rowptr, a.columns, a.values, b.rowptr, counts.p);
    const int triples = rowptr_from_counts(ledger, counts.p, a.n, offsets, name);
    counts.release();
    const int colbits = bits_for(b.m), totalbits = colbits + bits_for(a.n);
    Dev<uint64_t> keys(triples, name), sortedkeys(triples, name);
    Dev<double> values(triples, name), sortedvalues(triples, name);
    k_product_expand<<<blocks_of(a.n), TH>>>(a.n, colbits, offsets.p, a.rowptr, a.columns, a.values,
                                             b.rowptr, b.columns, b.values, keys.p, values.p);
    offsets.release();
    size_t bytes = 0;
    ck(cub::DeviceRadixSort::SortPairs(nullptr, bytes, keys.p, sortedkeys.p, values.p, sortedvalues.p, triples, 0, totalbits), name);
    {
        Dev<unsigned char> scratch(bytes, name);
        ck(cub::DeviceRadixSort::SortPairs(scratch.p, bytes, keys.p, sortedkeys.p, values.p, sortedvalues.p, triples, 0, totalbits), name);
    }
    keys.release(); values.release();
    Dev<uint64_t> unique(triples, name);
    Dev<int> runs(triples, name), count(1, name), starts;
    bytes = 0;
    ck(cub::DeviceRunLengthEncode::Encode(nullptr, bytes, sortedkeys.p, unique.p, runs.p, count.p, triples), name);
    {
        Dev<unsigned char> scratch(bytes, name);
        ck(cub::DeviceRunLengthEncode::Encode(scratch.p, bytes, sortedkeys.p, unique.p, runs.p, count.p, triples), name);
    }
    sortedkeys.release();
    const int nnz = read_scalar(ledger, count.p);
    int check = rowptr_from_counts(ledger, runs.p, nnz, starts, name);
    if (check != triples) throw std::runtime_error(std::string(name) + ": product runs lost entries");
    Csr c;
    c.n = a.n; c.m = b.m; c.nnz = nnz;
    c.values.allocate(nnz, name);
    c.columns.allocate(nnz, name);
    k_segment_sum<<<blocks_of(nnz), TH>>>(nnz, starts.p, sortedvalues.p, c.values.p);
    Dev<int> rowcount(a.n, name);
    rowcount.zero();
    k_decode_keys<<<blocks_of(nnz), TH>>>(nnz, colbits, unique.p, rowcount.p, c.columns.p);
    check = rowptr_from_counts(ledger, rowcount.p, a.n, c.rowptr, name);
    if (check != nnz) throw std::runtime_error(std::string(name) + ": product rows lost entries");
    ck_kernels(name);
    return c;
}

// Row-wise operations used by every route.
__global__ void k_axpy(int n, double alpha, const double* x, double* y) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) y[i] += alpha * x[i];
}
__global__ void k_xpay(int n, double alpha, const double* x, double* y) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) y[i] = x[i] + alpha * y[i];
}
__global__ void k_verify_finite(int n, const double* a, int* flag) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n && !isfinite(a[i])) atomicExch(flag, 1);
}
__global__ void k_l1_norms(int n, const int* rowptr, const double* values, double* l1) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double sum = 0;
    for (int k = rowptr[i]; k < rowptr[i + 1]; ++k) sum += fabs(values[k]);
    l1[i] = sum > 0 ? sum : 1.;
}
// Identity rows on constrained DOFs, with the constrained columns removed.
__global__ void k_constrain_values(int n, const int* rowptr, const int* columns, const int* free_dofs, const double* in, double* out) {
    int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    for (int k = rowptr[i]; k < rowptr[i + 1]; ++k) {
        const int j = columns[k];
        out[k] = free_dofs[i] ? (free_dofs[j] ? in[k] : 0.) : (i == j ? 1. : 0.);
    }
}

}  // namespace radia_fem
