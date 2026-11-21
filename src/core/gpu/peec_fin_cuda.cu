// Device-resident fin PEEC solve; see peec_fin_host.hpp.
#include "peec_fin_host.hpp"
#include "../peec_oriented_box.hpp"
#include <cublas_v2.h>
#include <cuComplex.h>
#include <cuda_runtime.h>
#include <cusolverDn.h>
#include <algorithm>
#include <chrono>
#include <cmath>
#include <limits>
#include <sstream>
#include <stdexcept>
#include <string>

namespace radia_peec_fin {
namespace {
namespace pb = radia_peec_box;
constexpr double kRelativeLimit = 1e-6;

void check(cudaError_t c, const char* what) {
    if (c != cudaSuccess) throw std::runtime_error(std::string(what) + ": " + cudaGetErrorString(c));
}
void check(cublasStatus_t c, const char* what) {
    if (c != CUBLAS_STATUS_SUCCESS) throw std::runtime_error(std::string(what) + ": cuBLAS status " + std::to_string(int(c)));
}
void check(cusolverStatus_t c, const char* what) {
    if (c != CUSOLVER_STATUS_SUCCESS) throw std::runtime_error(std::string(what) + ": cuSOLVER status " + std::to_string(int(c)));
}
template <class T> struct Buf {
    T* p = nullptr;
    std::size_t n = 0;
    explicit Buf(std::size_t count) : n(count) {
        check(cudaMalloc(reinterpret_cast<void**>(&p), (count ? count : 1) * sizeof(T)), "cudaMalloc");
    }
    Buf(const Buf&) = delete;
    Buf& operator=(const Buf&) = delete;
    ~Buf() { if (p) cudaFree(p); }
    void release() { if (p) cudaFree(p); p = nullptr; }
};
struct Handles {
    cublasHandle_t blas = nullptr;
    cusolverDnHandle_t solver = nullptr;
    Handles() { check(cublasCreate(&blas), "cublasCreate"); check(cusolverDnCreate(&solver), "cusolverDnCreate"); }
    ~Handles() { if (solver) cusolverDnDestroy(solver); if (blas) cublasDestroy(blas); }
};
using Clock = std::chrono::steady_clock;
double seconds(Clock::time_point a) { return std::chrono::duration<double>(Clock::now() - a).count(); }

struct Bars {
    const double *c, *d, *wd, *len, *wid, *hei;
    const std::uint8_t* has_wd;
};

__device__ pb::Frame bar_frame(const Bars& b, int i) {
    const pb::V3 c{b.c[3 * i], b.c[3 * i + 1], b.c[3 * i + 2]};
    const pb::V3 d{b.d[3 * i], b.d[3 * i + 1], b.d[3 * i + 2]};
    const pb::V3 w{b.wd[3 * i], b.wd[3 * i + 1], b.wd[3 * i + 2]};
    return pb::FrameFromBar(c, d, b.has_wd[i] ? &w : nullptr, b.len[i], b.wid[i], b.hei[i]);
}

// One thread per upper-triangle entry (i <= j), column-major symmetric L.
__global__ void assemble_l(int m, Bars b, pb::Rules q, double* L, long long n_pairs) {
    for (long long t = blockIdx.x * (long long)blockDim.x + threadIdx.x; t < n_pairs;
         t += (long long)gridDim.x * blockDim.x) {
        // t -> (i, j), j >= i, row-major upper triangle
        const double mm = 2.0 * m + 1.0;
        long long i = (long long)floor((mm - sqrt(mm * mm - 8.0 * (double)t)) / 2.0);
        auto start = [m](long long r) { return r * (2LL * m - r + 1) / 2; };
        while (i > 0 && start(i) > t) --i;
        while (start(i + 1) <= t) ++i;
        const long long j = i + (t - start(i));
        const double v = pb::Inductance(bar_frame(b, (int)i), bar_frame(b, (int)j), i == j, q, nullptr);
        L[i + j * (long long)m] = v;
        L[j + i * (long long)m] = v;
    }
}

__global__ void scan_l(long long n, const double* L, double* out /* [0]=max|L|, [1]=nonfinite */) {
    for (long long t = blockIdx.x * (long long)blockDim.x + threadIdx.x; t < n; t += (long long)gridDim.x * blockDim.x) {
        const double v = L[t];
        if (!isfinite(v)) atomicAdd(out + 1, 1.0);
        else {
            unsigned long long* a = reinterpret_cast<unsigned long long*>(out);
            unsigned long long old = *a, assumed;
            const double av = fabs(v);
            do {
                assumed = old;
                if (__longlong_as_double(assumed) >= av) break;
                old = atomicCAS(a, assumed, __double_as_longlong(av));
            } while (assumed != old);
        }
    }
}

__global__ void form_z(long long n2, int m, double omega, const double* L, const double* R, cuDoubleComplex* Z) {
    for (long long t = blockIdx.x * (long long)blockDim.x + threadIdx.x; t < n2; t += (long long)gridDim.x * blockDim.x) {
        const long long i = t % m, j = t / m;
        Z[t] = make_cuDoubleComplex(i == j ? R[i] : 0.0, omega * L[t]);
    }
}

// B = A_r^T (m x n_red, column-major): +1 at the from node, -1 at the to node.
__global__ void incidence_t(int m, int n_red, const int* nodes, const int* map, cuDoubleComplex* B) {
    for (int f = blockIdx.x * blockDim.x + threadIdx.x; f < m; f += gridDim.x * blockDim.x) {
        const int a = map[nodes[2 * f]], c = map[nodes[2 * f + 1]];
        if (a >= 0) B[f + (long long)a * m] = cuCadd(B[f + (long long)a * m], make_cuDoubleComplex(1.0, 0.0));
        if (c >= 0) B[f + (long long)c * m] = cuCsub(B[f + (long long)c * m], make_cuDoubleComplex(1.0, 0.0));
    }
}

double device_norm(cublasHandle_t h, int n, const cuDoubleComplex* x) {
    double r = 0.0;
    check(cublasDznrm2(h, n, x, 1, &r), "cublasDznrm2");
    return r;
}

// Spectrum summary on the device: [0] non-finite count, [1] |lambda| <= tol count.
__global__ void scan_spectrum(int m, const double* ev, double tol, unsigned long long* out) {
    for (int t = blockIdx.x * blockDim.x + threadIdx.x; t < m; t += gridDim.x * blockDim.x) {
        const double v = ev[t];
        if (!isfinite(v)) atomicAdd(out, 1ULL);
        else if (fabs(v) <= tol) atomicAdd(out + 1, 1ULL);
    }
}

int download_info(const int* info, std::size_t& bytes) {
    int h = 0;
    check(cudaMemcpy(&h, info, sizeof h, cudaMemcpyDeviceToHost), "download info");
    bytes += sizeof h;
    return h;
}

// Energy gate (same rule as fin_volume.require_energy_admissible) on a
// device-resident column-major symmetric L. Only status words and scalar
// spectrum summaries come back to the host.
void energy_gate(const Handles& h, const double* L, int m, GateResult& g, std::size_t& bytes) {
    const std::size_t mm = std::size_t(m) * std::size_t(m);
    {
        Buf<double> stat(2);
        check(cudaMemset(stat.p, 0, 2 * sizeof(double)), "memset");
        scan_l<<<2048, 256>>>((long long)mm, L, stat.p);
        check(cudaGetLastError(), "scan_l launch");
        check(cudaDeviceSynchronize(), "scan_l");
        double s[2];
        check(cudaMemcpy(s, stat.p, sizeof s, cudaMemcpyDeviceToHost), "download");
        bytes += sizeof s;
        if (s[1] > 0) throw std::runtime_error("fin PEEC GPU: inductance matrix has non-finite entries");
        g.max_abs_L = s[0];
    }
    Buf<double> Lc(mm);
    check(cudaMemcpy(Lc.p, L, mm * sizeof(double), cudaMemcpyDeviceToDevice), "copy L");
    Buf<int> info(1);
    int lwork = 0;
    check(cusolverDnDpotrf_bufferSize(h.solver, CUBLAS_FILL_MODE_LOWER, m, Lc.p, m, &lwork), "potrf size");
    {
        Buf<double> work(lwork);
        check(cusolverDnDpotrf(h.solver, CUBLAS_FILL_MODE_LOWER, m, Lc.p, m, work.p, lwork, info.p), "potrf");
    }
    const int pinfo = download_info(info.p, bytes);
    if (pinfo < 0) throw std::runtime_error("fin PEEC GPU: potrf argument error " + std::to_string(pinfo));
    g.cholesky = pinfo == 0;
    if (g.cholesky) return;
    // Not positive definite to working precision: decide on the spectrum.
    check(cudaMemcpy(Lc.p, L, mm * sizeof(double), cudaMemcpyDeviceToDevice), "copy L");
    Buf<double> ev(m);
    check(cusolverDnDsyevd_bufferSize(h.solver, CUSOLVER_EIG_MODE_NOVECTOR, CUBLAS_FILL_MODE_LOWER, m,
                                      Lc.p, m, ev.p, &lwork), "syevd size");
    {
        Buf<double> work(lwork);
        check(cusolverDnDsyevd(h.solver, CUSOLVER_EIG_MODE_NOVECTOR, CUBLAS_FILL_MODE_LOWER, m, Lc.p, m,
                               ev.p, work.p, lwork, info.p), "syevd");
    }
    const int einfo = download_info(info.p, bytes);
    if (einfo != 0)
        throw std::runtime_error("fin PEEC GPU: eigensolver did not converge (syevd info " +
                                 std::to_string(einfo) + ")");
    double lo = 0.0, hi = 0.0;  // ascending eigenvalues: extremes only
    check(cudaMemcpy(&lo, ev.p, sizeof lo, cudaMemcpyDeviceToHost), "download");
    check(cudaMemcpy(&hi, ev.p + (m - 1), sizeof hi, cudaMemcpyDeviceToHost), "download");
    bytes += 2 * sizeof(double);
    if (!std::isfinite(lo) || !std::isfinite(hi))
        throw std::runtime_error("fin PEEC GPU: eigensolver returned non-finite extremes");
    const double scale = std::max(std::abs(lo), std::abs(hi));
    const double tol = m * std::numeric_limits<double>::epsilon() * scale;
    Buf<unsigned long long> cnt(2);
    check(cudaMemset(cnt.p, 0, 2 * sizeof(unsigned long long)), "memset");
    scan_spectrum<<<256, 256>>>(m, ev.p, tol, cnt.p);
    check(cudaGetLastError(), "scan_spectrum launch");
    unsigned long long c2[2];
    check(cudaMemcpy(c2, cnt.p, sizeof c2, cudaMemcpyDeviceToHost), "download");
    bytes += sizeof c2;
    if (c2[0] > 0) throw std::runtime_error("fin PEEC GPU: eigensolver returned non-finite eigenvalues");
    if (lo < -tol) {
        std::ostringstream msg;
        msg << "fin PEEC GPU: inductance matrix has negative energy: min eigenvalue " << lo
            << " H (" << lo / scale << " of the spectral scale; roundoff scale " << tol << " H)";
        throw std::runtime_error(msg.str());
    }
    g.min_eig_over_max = scale > 0 ? lo / scale : 0.0;
    g.n_null = int(c2[1]);
}

}  // namespace

GateResult energy_gate_host(const double* L_host, int m) {
    if (m <= 0) throw std::invalid_argument("fin PEEC GPU: empty matrix");
    Handles h;
    Buf<double> L(std::size_t(m) * m);
    check(cudaMemcpy(L.p, L_host, std::size_t(m) * m * sizeof(double), cudaMemcpyHostToDevice), "upload");
    GateResult g;
    std::size_t bytes = 0;
    energy_gate(h, L.p, m, g, bytes);
    return g;
}

Output solve_host(const Input& in) {
    const int m = in.m;
    if (m <= 0 || in.order < 1 || in.order > pb::kMaxOrder || in.n_nodes <= 0 || in.n_reduced <= 0 ||
        !(in.frequency > 0.0))
        throw std::invalid_argument("fin PEEC GPU: invalid sizes, order or frequency");
    Output out;
    Handles h;
    const std::size_t mm = std::size_t(m) * std::size_t(m);
    auto up = [&](auto* dst, const auto* src, std::size_t count) {
        const std::size_t bytes = count * sizeof(*src);
        check(cudaMemcpy(dst, src, bytes, cudaMemcpyHostToDevice), "upload");
        out.bytes_to_device += bytes;
    };
    // Non-matrix inputs.
    Buf<double> c(3 * m), d(3 * m), wd(3 * m), len(m), wid(m), hei(m), R(m);
    Buf<std::uint8_t> has(m);
    const int nr = in.order * (in.order + 1) / 2;
    Buf<double> rx(nr), rw(nr);
    Buf<int> nodes(2 * m), map(in.n_nodes);
    up(c.p, in.center, 3 * m); up(d.p, in.direction, 3 * m); up(wd.p, in.width_dir, 3 * m);
    up(len.p, in.length, m); up(wid.p, in.width, m); up(hei.p, in.height, m); up(R.p, in.resistance, m);
    up(has.p, in.has_width_dir, m); up(rx.p, in.rule_x, nr); up(rw.p, in.rule_w, nr);
    up(nodes.p, in.branch_nodes, 2 * m); up(map.p, in.node_map, in.n_nodes);

    // 1. L on the device.
    auto t0 = Clock::now();
    Buf<double> L(mm);
    const long long n_pairs = (long long)m * (m + 1) / 2;
    const Bars bars{c.p, d.p, wd.p, len.p, wid.p, hei.p, has.p};
    assemble_l<<<4096, 128>>>(m, bars, pb::Rules{rx.p, rw.p, in.order}, L.p, n_pairs);
    check(cudaGetLastError(), "assemble_l launch");
    check(cudaDeviceSynchronize(), "assemble_l");
    out.t_assemble_s = seconds(t0);

    // 2. Energy gate on the device.
    t0 = Clock::now();
    {
        GateResult g;
        energy_gate(h, L.p, m, g, out.bytes_to_host);
        out.cholesky = g.cholesky; out.min_eig_over_max = g.min_eig_over_max;
        out.n_null = g.n_null; out.max_abs_L = g.max_abs_L;
    }
    out.t_energy_s = seconds(t0);

    // 3. Z = R + j omega L, LU on the device.
    t0 = Clock::now();
    const double omega = 2.0 * 3.141592653589793238462643383279502884 * in.frequency;
    Buf<cuDoubleComplex> Z(mm);
    form_z<<<4096, 256>>>((long long)mm, m, omega, L.p, R.p, Z.p);
    check(cudaGetLastError(), "form_z");
    Buf<cuDoubleComplex> Zf(mm);
    check(cudaMemcpy(Zf.p, Z.p, mm * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice), "copy Z");
    int lwork = 0;
    check(cusolverDnZgetrf_bufferSize(h.solver, m, m, Zf.p, m, &lwork), "getrf size");
    Buf<cuDoubleComplex> work(lwork);
    Buf<int> piv(m), info(1);
    check(cusolverDnZgetrf(h.solver, m, m, Zf.p, m, work.p, piv.p, info.p), "getrf");
    int hinfo = 0;
    check(cudaMemcpy(&hinfo, info.p, sizeof hinfo, cudaMemcpyDeviceToHost), "download");
    out.bytes_to_host += sizeof hinfo;
    if (hinfo != 0) throw std::runtime_error("fin PEEC GPU: Z is singular (getrf info " + std::to_string(hinfo) + ")");
    work.release();
    out.t_factor_s = seconds(t0);

    // 4. X = Z^-1 A_r^T, Y = A_r X, V = Y^-1 J, I = X V.
    t0 = Clock::now();
    const int n = in.n_reduced;
    Buf<cuDoubleComplex> B(std::size_t(m) * n), X(std::size_t(m) * n), Y(std::size_t(n) * n);
    check(cudaMemset(B.p, 0, std::size_t(m) * n * sizeof(cuDoubleComplex)), "memset");
    incidence_t<<<256, 256>>>(m, n, nodes.p, map.p, B.p);
    check(cudaGetLastError(), "incidence");
    check(cudaMemcpy(X.p, B.p, std::size_t(m) * n * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice), "copy");
    check(cusolverDnZgetrs(h.solver, CUBLAS_OP_N, m, n, Zf.p, m, piv.p, X.p, m, info.p), "getrs");
    if (download_info(info.p, out.bytes_to_host) != 0) throw std::runtime_error("fin PEEC GPU: getrs argument error");
    const cuDoubleComplex one = make_cuDoubleComplex(1, 0), zero = make_cuDoubleComplex(0, 0),
                          mone = make_cuDoubleComplex(-1, 0);
    check(cublasZgemm(h.blas, CUBLAS_OP_T, CUBLAS_OP_N, n, n, m, &one, B.p, m, X.p, m, &zero, Y.p, n), "A X");
    Zf.release();
    Buf<cuDoubleComplex> V(n), I(m), r(m);
    Buf<cuDoubleComplex> J(n);
    up(reinterpret_cast<std::complex<double>*>(J.p), in.injection, n);
    check(cudaMemcpy(V.p, J.p, n * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice), "copy J");
    check(cusolverDnZgetrf_bufferSize(h.solver, n, n, Y.p, n, &lwork), "getrf size");
    Buf<cuDoubleComplex> w2(lwork);
    Buf<int> piv2(n);
    check(cusolverDnZgetrf(h.solver, n, n, Y.p, n, w2.p, piv2.p, info.p), "getrf Y");
    hinfo = download_info(info.p, out.bytes_to_host);
    if (hinfo != 0)
        throw std::runtime_error("fin PEEC GPU: nodal admittance is singular (getrf info " + std::to_string(hinfo) + ")");
    check(cusolverDnZgetrs(h.solver, CUBLAS_OP_N, n, 1, Y.p, n, piv2.p, V.p, n, info.p), "getrs Y");
    if (download_info(info.p, out.bytes_to_host) != 0) throw std::runtime_error("fin PEEC GPU: getrs Y argument error");
    check(cublasZgemv(h.blas, CUBLAS_OP_N, m, n, &one, X.p, m, V.p, 1, &zero, I.p, 1), "X V");
    // True residual of the branch law: Z I - A^T V.
    check(cublasZgemv(h.blas, CUBLAS_OP_N, m, n, &one, B.p, m, V.p, 1, &zero, r.p, 1), "A^T V");
    const double rhs_norm = device_norm(h.blas, m, r.p);
    check(cublasZgemv(h.blas, CUBLAS_OP_N, m, m, &one, Z.p, m, I.p, 1, &mone, r.p, 1), "Z I - A^T V");
    out.z_true_residual = device_norm(h.blas, m, r.p) / (rhs_norm > 0 ? rhs_norm : 1.0);
    if (!(out.z_true_residual <= kRelativeLimit))
        throw std::runtime_error("fin PEEC GPU: true residual " + std::to_string(out.z_true_residual) + " > 1e-6");
    // True KCL residual: A_r I - J (A_r I = B^T I).
    {
        Buf<cuDoubleComplex> k(n);
        check(cudaMemcpy(k.p, J.p, n * sizeof(cuDoubleComplex), cudaMemcpyDeviceToDevice), "copy J");
        const double jn = device_norm(h.blas, n, J.p);
        check(cublasZgemv(h.blas, CUBLAS_OP_T, m, n, &one, B.p, m, I.p, 1, &mone, k.p, 1), "A I - J");
        out.kcl_true_residual = device_norm(h.blas, n, k.p) / (jn > 0 ? jn : 1.0);
        if (!(out.kcl_true_residual <= kRelativeLimit))
            throw std::runtime_error("fin PEEC GPU: KCL residual " + std::to_string(out.kcl_true_residual) + " > 1e-6");
    }
    // I^H L I and I^H R I.
    {
        Buf<double> Ir(m), Ii(m), tmp(m);
        // Real and imaginary parts via strided copies.
        check(cublasDcopy(h.blas, m, reinterpret_cast<const double*>(I.p), 2, Ir.p, 1), "copy re");
        check(cublasDcopy(h.blas, m, reinterpret_cast<const double*>(I.p) + 1, 2, Ii.p, 1), "copy im");
        const double done = 1.0, dzero = 0.0;
        double a = 0.0, b2 = 0.0;
        check(cublasDsymv(h.blas, CUBLAS_FILL_MODE_LOWER, m, &done, L.p, m, Ir.p, 1, &dzero, tmp.p, 1), "L re");
        check(cublasDdot(h.blas, m, Ir.p, 1, tmp.p, 1, &a), "dot");
        check(cublasDsymv(h.blas, CUBLAS_FILL_MODE_LOWER, m, &done, L.p, m, Ii.p, 1, &dzero, tmp.p, 1), "L im");
        check(cublasDdot(h.blas, m, Ii.p, 1, tmp.p, 1, &b2), "dot");
        out.current_L_current = a + b2;
    }
    out.current.resize(m);
    out.voltage.resize(n);
    check(cudaMemcpy(out.current.data(), I.p, m * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost), "download I");
    check(cudaMemcpy(out.voltage.data(), V.p, n * sizeof(cuDoubleComplex), cudaMemcpyDeviceToHost), "download V");
    out.bytes_to_host += (std::size_t(m) + n) * sizeof(cuDoubleComplex);
    for (int f = 0; f < m; ++f) out.current_R_current += in.resistance[f] * std::norm(out.current[f]);
    out.t_solve_s = seconds(t0);
    return out;
}

}  // namespace radia_peec_fin
