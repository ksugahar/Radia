/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
// Krylov iterations on the device. Vectors stay on the device; dot products cross
// as control scalars and are recorded. Real systems use conjugate gradients and
// complex-symmetric systems use COCR with unconjugated products, as the CPU route.
#include <cublas_v2.h>
#include <cmath>
#include <complex>
#include <cstdio>
#include <functional>
#include <string>
#include "fem_api.hpp"
#include "fem_device.cuh"

namespace radia_fem {

__global__ void k_copy(int n, const double* x, double* y) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n) y[i] = x[i];
}
__global__ void k_caxpy(int n, double re, double im, const double* xr, const double* xi, double* yr, double* yi) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    yr[i] += re * xr[i] - im * xi[i];
    yi[i] += re * xi[i] + im * xr[i];
}
__global__ void k_cxpay(int n, double re, double im, const double* xr, const double* xi, double* yr, double* yi) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const double ar = yr[i], ai = yi[i];
    yr[i] = xr[i] + re * ar - im * ai;
    yi[i] = xi[i] + re * ai + im * ar;
}

class Blas {
public:
    explicit Blas(Ledger& ledger) : ledger_(ledger) {
        const cublasStatus_t status = cublasCreate(&handle_);
        if (status != CUBLAS_STATUS_SUCCESS) throw std::runtime_error("cuBLAS handle creation failed");
    }
    ~Blas() { cublasDestroy(handle_); }
    Blas(const Blas&) = delete;
    Blas& operator=(const Blas&) = delete;
    double dot(int n, const double* x, const double* y) {
        double value = 0;
        ck_kernels("dot");
        if (cublasDdot(handle_, n, x, 1, y, 1, &value) != CUBLAS_STATUS_SUCCESS) throw std::runtime_error("cuBLAS dot failed");
        ledger_.add("scalar.krylov", "D2H", sizeof(double));
        return value;
    }
    std::complex<double> dot(int n, const double* xr, const double* xi, const double* yr, const double* yi) {
        return {dot(n, xr, yr) - dot(n, xi, yi), dot(n, xr, yi) + dot(n, xi, yr)};
    }

private:
    Ledger& ledger_;
    cublasHandle_t handle_ = nullptr;
};

using RealAction = std::function<void(const double*, double*)>;
using ComplexAction = std::function<void(const double*, const double*, double*, double*)>;

inline std::string not_converged(const char* method, int max_iterations, double relative, double tolerance) {
    char text[128];
    std::snprintf(text, sizeof text, "device %s did not converge in %d iterations: relative residual %.3e, requested %.3e", method, max_iterations, relative,
                  tolerance);
    return text;
}

inline SolveReport pcg(Ledger& ledger, int n, const RealAction& A, const RealAction& M, const double* b, double* x,
                       double tolerance, int max_iterations) {
    Blas blas(ledger);
    const int blocks = blocks_of(n);
    Dev<double> r(size_t(n), "cg"), z(size_t(n), "cg"), p(size_t(n), "cg"), q(size_t(n), "cg");
    SolveReport report;
    k_copy<<<blocks, TH>>>(n, b, r.p);
    const double b_norm = std::sqrt(blas.dot(n, b, b));
    report.initial_norm = b_norm;
    if (!(b_norm > 0) || !std::isfinite(b_norm))
        throw std::runtime_error("device CG: the right-hand side is zero or nonfinite on the free DOFs");
    M(r.p, z.p);
    k_copy<<<blocks, TH>>>(n, z.p, p.p);
    double rho = blas.dot(n, r.p, z.p);
    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        A(p.p, q.p);
        const double curvature = blas.dot(n, p.p, q.p);
        if (!(curvature > 0) || !(rho > 0) || !std::isfinite(rho))
            throw std::runtime_error("device CG breakdown at iteration " + std::to_string(iteration) + ": rho " + std::to_string(rho)
                                     + " curvature " + std::to_string(curvature));
        const double alpha = rho / curvature;
        k_axpy<<<blocks, TH>>>(n, alpha, p.p, x);
        k_axpy<<<blocks, TH>>>(n, -alpha, q.p, r.p);
        const double r_norm = std::sqrt(blas.dot(n, r.p, r.p));
        report.history.push_back(r_norm / b_norm);
        report.iterations = iteration;
        report.final_norm = r_norm;
        if (!std::isfinite(r_norm)) throw std::runtime_error("device CG: nonfinite residual");
        if (r_norm <= tolerance * b_norm) return report;
        M(r.p, z.p);
        const double rho_next = blas.dot(n, r.p, z.p);
        k_xpay<<<blocks, TH>>>(n, rho_next / rho, z.p, p.p);
        rho = rho_next;
    }
    throw std::runtime_error(not_converged("CG", max_iterations, report.final_norm / b_norm, tolerance));
}

inline SolveReport cocr(Ledger& ledger, int n, const ComplexAction& A, const ComplexAction& M, const double* br, const double* bi,
                        double* xr, double* xi, double tolerance, int max_iterations) {
    Blas blas(ledger);
    const int blocks = blocks_of(n);
    const char* name = "cocr";
    Dev<double> rr(size_t(n), name), ri(size_t(n), name), tr(size_t(n), name), ti(size_t(n), name), pr(size_t(n), name), pi(size_t(n), name),
        qr(size_t(n), name), qi(size_t(n), name), ur(size_t(n), name), ui(size_t(n), name), sr(size_t(n), name), si(size_t(n), name);
    SolveReport report;
    k_copy<<<blocks, TH>>>(n, br, rr.p);
    k_copy<<<blocks, TH>>>(n, bi, ri.p);
    const double b_norm = std::sqrt(blas.dot(n, br, br) + blas.dot(n, bi, bi));
    report.initial_norm = b_norm;
    if (!(b_norm > 0) || !std::isfinite(b_norm))
        throw std::runtime_error("device COCR: the right-hand side is zero or nonfinite on the free DOFs");
    M(rr.p, ri.p, tr.p, ti.p);
    k_copy<<<blocks, TH>>>(n, tr.p, pr.p);
    k_copy<<<blocks, TH>>>(n, ti.p, pi.p);
    A(pr.p, pi.p, qr.p, qi.p);
    std::complex<double> rho = blas.dot(n, tr.p, ti.p, qr.p, qi.p);
    for (int iteration = 1; iteration <= max_iterations; ++iteration) {
        M(qr.p, qi.p, ur.p, ui.p);
        const std::complex<double> mu = blas.dot(n, ur.p, ui.p, qr.p, qi.p);
        if (!(std::abs(mu) > 0) || !std::isfinite(std::abs(mu)))
            throw std::runtime_error("device COCR breakdown at iteration " + std::to_string(iteration));
        const std::complex<double> alpha = rho / mu;
        k_caxpy<<<blocks, TH>>>(n, alpha.real(), alpha.imag(), pr.p, pi.p, xr, xi);
        k_caxpy<<<blocks, TH>>>(n, -alpha.real(), -alpha.imag(), qr.p, qi.p, rr.p, ri.p);
        k_caxpy<<<blocks, TH>>>(n, -alpha.real(), -alpha.imag(), ur.p, ui.p, tr.p, ti.p);
        const double r_norm = std::sqrt(blas.dot(n, rr.p, rr.p) + blas.dot(n, ri.p, ri.p));
        report.history.push_back(r_norm / b_norm);
        report.iterations = iteration;
        report.final_norm = r_norm;
        if (!std::isfinite(r_norm)) throw std::runtime_error("device COCR: nonfinite residual");
        if (r_norm <= tolerance * b_norm) return report;
        A(tr.p, ti.p, sr.p, si.p);
        const std::complex<double> rho_next = blas.dot(n, tr.p, ti.p, sr.p, si.p);
        const std::complex<double> beta = rho_next / rho;
        if (!std::isfinite(std::abs(beta)))
            throw std::runtime_error("device COCR breakdown at iteration " + std::to_string(iteration));
        k_cxpay<<<blocks, TH>>>(n, beta.real(), beta.imag(), tr.p, ti.p, pr.p, pi.p);
        k_cxpay<<<blocks, TH>>>(n, beta.real(), beta.imag(), sr.p, si.p, qr.p, qi.p);
        rho = rho_next;
    }
    throw std::runtime_error(not_converged("COCR", max_iterations, report.final_norm / b_norm, tolerance));
}

}  // namespace radia_fem
