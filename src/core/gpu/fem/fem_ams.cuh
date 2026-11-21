/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
// Auxiliary-space Maxwell preconditioner built on the device. The discrete
// gradient and the Nedelec interpolation come from the mesh edge list and the
// vertex coordinates; each Galerkin operator is a device product and each
// auxiliary hierarchy is a device AMG.
#include "fem_amg.cuh"

namespace radia_fem {

// Edge row: first the low vertex with -1, then the high vertex with +1. Columns
// stay sorted because the low vertex number comes first.
__global__ void k_edge_pattern(int nedge, const int* edge_vertices, int* rp, int* ci) {
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e > nedge) return;
    rp[e] = 2 * e;
    if (e < nedge) {
        ci[2 * e] = edge_vertices[2 * e];
        ci[2 * e + 1] = edge_vertices[2 * e + 1];
    }
}
__global__ void k_gradient_values(int nedge, double* values) {
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= nedge) return;
    values[2 * e] = -1.;
    values[2 * e + 1] = 1.;
}
__global__ void k_interpolation_values(int nedge, int direction, const int* edge_vertices, const double* coordinates, double* values) {
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= nedge) return;
    const double delta = coordinates[size_t(edge_vertices[2 * e + 1]) * 3 + direction] - coordinates[size_t(edge_vertices[2 * e]) * 3 + direction];
    values[2 * e] = 0.5 * delta;
    values[2 * e + 1] = 0.5 * delta;
}
// A row with no nonzero value keeps unit diagonal when it has a diagonal entry.
__global__ void k_fix_zero_rows(int n, const int* rp, const int* ci, double* av) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    double sum = 0;
    for (int k = rp[i]; k < rp[i + 1]; ++k) sum += fabs(av[k]);
    if (sum > 0) return;
    for (int k = rp[i]; k < rp[i + 1]; ++k)
        if (ci[k] == i) av[k] = 1.;
}

struct AmsSpace {
    Csr p, pt;
    Amg amg;
    void add(const double* residual, double alpha, double* output) {
        multiply<<<blocks_of(pt.n), TH>>>(pt.n, pt.rowptr.p, pt.columns.p, pt.values.p, residual, amg.input());
        amg.cycle();
        multiply_add<<<blocks_of(p.n), TH>>>(p.n, p.rowptr.p, p.columns.p, p.values.p, amg.output(), alpha, output);
    }
};

class Ams {
public:
    // Takes the real operator, sets identity rows on the constrained DOFs and
    // builds the four auxiliary spaces from it.
    std::vector<AmgReport> build(Ledger& ledger, const CsrView& source, int nnz, const int* free_dofs, const double* values,
                                 int nedge, int nv, const int* edge_vertices, const double* coordinates, int min_coarse, double theta) {
        if (source.n != nedge)
            throw std::runtime_error("device AMS: the operator must have one row per mesh edge, rows " + std::to_string(source.n)
                                     + " edges " + std::to_string(nedge));
        const char* name = "AMS setup";
        n_ = source.n;
        free_ = free_dofs;
        a_.n = a_.m = n_;
        a_.nnz = nnz;
        a_.rowptr.allocate(size_t(n_) + 1, name);
        a_.columns.allocate(size_t(nnz), name);
        a_.values.allocate(size_t(nnz), name);
        ck(cudaMemcpy(a_.rowptr.p, source.rowptr, (size_t(n_) + 1) * sizeof(int), cudaMemcpyDeviceToDevice), name);
        ck(cudaMemcpy(a_.columns.p, source.columns, size_t(nnz) * sizeof(int), cudaMemcpyDeviceToDevice), name);
        k_constrain_values<<<blocks_of(n_), TH>>>(n_, source.rowptr, source.columns, free_dofs, values, a_.values.p);
        l1_.allocate(size_t(n_), name);
        k_l1_norms<<<blocks_of(n_), TH>>>(n_, a_.rowptr.p, a_.values.p, l1_.p);
        x_.allocate(size_t(n_), name);
        b_.allocate(size_t(n_), name);
        r_.allocate(size_t(n_), name);
        tmp_.allocate(size_t(n_), name);
        ck_kernels(name);
        AmgOptions options;
        options.theta = theta;
        options.min_coarse = min_coarse;
        options.stall = 0.8;
        std::vector<AmgReport> reports;
        for (int s = 0; s < 4; ++s) {
            AmsSpace& space = spaces_[s];
            Csr& p = space.p;
            p.n = nedge;
            p.m = nv;
            p.nnz = 2 * nedge;
            p.rowptr.allocate(size_t(nedge) + 1, name);
            p.columns.allocate(size_t(p.nnz), name);
            p.values.allocate(size_t(p.nnz), name);
            k_edge_pattern<<<blocks_of(nedge + 1), TH>>>(nedge, edge_vertices, p.rowptr.p, p.columns.p);
            if (s == 0) k_gradient_values<<<blocks_of(nedge), TH>>>(nedge, p.values.p);
            else k_interpolation_values<<<blocks_of(nedge), TH>>>(nedge, s - 1, edge_vertices, coordinates, p.values.p);
            ck_kernels(name);
            space.pt = transpose(ledger, p.view(), p.nnz, name);
            Csr ap = multiply_csr(ledger, a_.view(), p.view(), name);
            Csr galerkin = multiply_csr(ledger, space.pt.view(), ap.view(), name);
            k_fix_zero_rows<<<blocks_of(galerkin.n), TH>>>(galerkin.n, galerkin.rowptr.p, galerkin.columns.p, galerkin.values.p);
            ck_kernels(name);
            options.constant_kernel = s == 0;
            reports.push_back(space.amg.build(ledger, std::move(galerkin), options, name));
        }
        return reports;
    }
    // The 01210 cycle of the CPU route: fine smoothing, gradient correction, the
    // three nodal corrections from one residual, gradient correction, smoothing.
    void apply(const double* rhs, double* output) {
        const int blocks = blocks_of(n_);
        ck(cudaMemcpy(b_.p, rhs, size_t(n_) * sizeof(double), cudaMemcpyDeviceToDevice), "AMS input");
        mask<<<blocks, TH>>>(n_, free_, b_.p);
        initial_relax<<<blocks, TH>>>(n_, b_.p, l1_.p, x_.p);
        gradient();
        residual();
        for (int s = 1; s < 4; ++s) spaces_[s].add(r_.p, 1., x_.p);
        gradient();
        multiply_relax<<<blocks, TH>>>(n_, a_.rowptr.p, a_.columns.p, a_.values.p, x_.p, b_.p, l1_.p, tmp_.p);
        std::swap(x_.p, tmp_.p);
        mask<<<blocks, TH>>>(n_, free_, x_.p);
        ck(cudaMemcpy(output, x_.p, size_t(n_) * sizeof(double), cudaMemcpyDeviceToDevice), "AMS output");
        ck_kernels("AMS apply");
    }

    int rows() const { return n_; }
    int columns() const { return spaces_[0].p.m; }
    // The gradient (space 0) or an interpolation (space 1 to 3) or its transpose on a device vector.
    void transfer(int space, bool transpose, const double* x, double* y) const {
        const Csr& p = transpose ? spaces_[space].pt : spaces_[space].p;
        multiply<<<blocks_of(p.n), TH>>>(p.n, p.rowptr.p, p.columns.p, p.values.p, x, y);
    }

private:
    void residual() {
        multiply_residual<<<blocks_of(n_), TH>>>(n_, a_.rowptr.p, a_.columns.p, a_.values.p, x_.p, b_.p, r_.p);
    }
    void gradient() {
        residual();
        spaces_[0].add(r_.p, 1., x_.p);
    }
    int n_ = 0;
    const int* free_ = nullptr;
    Csr a_;
    Dev<double> l1_, x_, b_, r_, tmp_;
    AmsSpace spaces_[4];
};

}  // namespace radia_fem
