/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
// Balancing domain decomposition by constraints for p>1 spaces, built on the
// device. Element matrices come from the device forms; the interior inverse, the
// harmonic extension and the wirebasket Schur complement are batched per element
// and reduced through device plans. The wirebasket is solved by device AMG or AMS.
// This is a Radia-owned preconditioner, compared with validated routes on identical
// inputs; it is not a copy of another implementation's numbers.
#include "fem_amg.cuh"
#include "fem_ams.cuh"
#include "fem_elements.cuh"
#include "fem_krylov.cuh"

namespace radia_fem {

struct cplx {
    double re, im;
    __host__ __device__ cplx& operator+=(const cplx& o) { re += o.re; im += o.im; return *this; }
};
__host__ __device__ inline cplx operator*(const cplx& a, const cplx& b) { return cplx{a.re * b.re - a.im * b.im, a.re * b.im + a.im * b.re}; }
__host__ __device__ inline cplx operator-(const cplx& a, const cplx& b) { return cplx{a.re - b.re, a.im - b.im}; }
template <class S> __host__ __device__ inline S make_scalar(double re, double im);
template <> __host__ __device__ inline double make_scalar<double>(double re, double) { return re; }
template <> __host__ __device__ inline cplx make_scalar<cplx>(double re, double im) { return cplx{re, im}; }
__host__ __device__ inline double magnitude(double a) { return fabs(a); }
__host__ __device__ inline double magnitude(const cplx& a) { return sqrt(a.re * a.re + a.im * a.im); }
__host__ __device__ inline double invert_scalar(double a) { return 1. / a; }
__host__ __device__ inline cplx invert_scalar(const cplx& a) {
    const double d = a.re * a.re + a.im * a.im;
    return cplx{a.re / d, -a.im / d};
}
__host__ __device__ inline double scale_by(double a, double w) { return a * w; }
__host__ __device__ inline cplx scale_by(const cplx& a, double w) { return cplx{a.re * w, a.im * w}; }

template <class S>
__device__ inline S element_entry(const GroupView& g, const int* free_dofs, const int* local, int c, int e, int a, int b) {
    if (!free_dofs[local[a]] || !free_dofs[local[b]]) return make_scalar<S>(a == b ? 1. : 0., 0.);
    const double* coef = g.coef + size_t(e) * COLUMNS;
    S value = make_scalar<S>(0., 0.);
    if (coef[STIFF_RE] != 0. || coef[STIFF_IM] != 0.)
        value += make_scalar<S>(coef[STIFF_RE], coef[STIFF_IM]) * make_scalar<S>(form_entry(g.stiff, g.nl, c, e, a, b), 0.);
    if (coef[MASS_RE] != 0. || coef[MASS_IM] != 0.)
        value += make_scalar<S>(coef[MASS_RE], coef[MASS_IM]) * make_scalar<S>(form_entry(g.mass, g.nl, c, e, a, b), 0.);
    return value;
}

// One thread per element of the batch. Interior inverse in place, harmonic
// extension, Schur complement and the interior weights.
template <class S>
__global__ void k_bddc_element(GroupView g, int first, int count, int ni, int nw, const int* free_dofs, const int* int_local, const int* wb_local,
                               S* scratch, S* inner, S* harmonic, S* schur, double* weights, int* bad) {
    const int index = blockIdx.x * blockDim.x + threadIdx.x;
    if (index >= count) return;
    const int e = first + index;
    const int* local = g.dofs + size_t(e) * g.nl;
    const int c = g.cls[e];
    S* d = inner + size_t(e) * ni * ni;
    S* he = harmonic + size_t(e) * ni * nw;
    S* se = schur + size_t(e) * nw * nw;
    double* w = weights + size_t(e) * ni;
    S* aiw = scratch + size_t(index) * ni * nw;
    for (int r = 0; r < ni; ++r) {
        for (int s = 0; s < ni; ++s) d[size_t(r) * ni + s] = element_entry<S>(g, free_dofs, local, c, e, int_local[r], int_local[s]);
        w[r] = magnitude(d[size_t(r) * ni + r]);
        for (int k = 0; k < nw; ++k) aiw[size_t(r) * nw + k] = element_entry<S>(g, free_dofs, local, c, e, int_local[r], wb_local[k]);
    }
    for (int k = 0; k < ni; ++k) {
        const S pivot = d[size_t(k) * ni + k];
        const double size = magnitude(pivot);
        if (!(size > 0) || !isfinite(size)) { atomicExch(bad, 1); return; }
        const S inverse = invert_scalar(pivot);
        d[size_t(k) * ni + k] = make_scalar<S>(1., 0.);
        for (int j = 0; j < ni; ++j) d[size_t(k) * ni + j] = d[size_t(k) * ni + j] * inverse;
        for (int i = 0; i < ni; ++i) {
            if (i == k) continue;
            const S factor = d[size_t(i) * ni + k];
            d[size_t(i) * ni + k] = make_scalar<S>(0., 0.);
            for (int j = 0; j < ni; ++j) d[size_t(i) * ni + j] = d[size_t(i) * ni + j] - factor * d[size_t(k) * ni + j];
        }
    }
    for (int r = 0; r < ni; ++r)
        for (int k = 0; k < nw; ++k) {
            S sum = make_scalar<S>(0., 0.);
            for (int s = 0; s < ni; ++s) sum += d[size_t(r) * ni + s] * aiw[size_t(s) * nw + k];
            he[size_t(r) * nw + k] = make_scalar<S>(0., 0.) - sum;
        }
    for (int a = 0; a < nw; ++a)
        for (int b = 0; b < nw; ++b) {
            S value = element_entry<S>(g, free_dofs, local, c, e, wb_local[a], wb_local[b]);
            for (int r = 0; r < ni; ++r) value += aiw[size_t(r) * nw + a] * he[size_t(r) * nw + b];
            se[size_t(a) * nw + b] = value;
        }
    for (int r = 0; r < ni; ++r) {
        for (int k = 0; k < nw; ++k) he[size_t(r) * nw + k] = scale_by(he[size_t(r) * nw + k], w[r]);
        for (int s = 0; s < ni; ++s) d[size_t(r) * ni + s] = scale_by(d[size_t(r) * ni + s], w[r] * w[s]);
    }
}

__global__ void k_bddc_keys(int ne, int nl, int nr, int nc, int colbits, const int* dofs, const int* rl, const int* cl, const int* rowmap,
                            const int* colmap, uint64_t* keys) {
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= ne) return;
    const int* local = dofs + size_t(e) * nl;
    uint64_t* out = keys + size_t(e) * nr * nc;
    for (int r = 0; r < nr; ++r) {
        const int row = rowmap ? rowmap[local[rl[r]]] : local[rl[r]];
        for (int c = 0; c < nc; ++c) {
            const int column = colmap ? colmap[local[cl[c]]] : local[cl[c]];
            out[size_t(r) * nc + c] = (uint64_t(row) << colbits) | uint64_t(column);
        }
    }
}
__global__ void k_bddc_weight_keys(int ne, int nl, int ni, const int* dofs, const int* int_local, uint64_t* keys) {
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= ne) return;
    const int* local = dofs + size_t(e) * nl;
    for (int r = 0; r < ni; ++r) keys[size_t(e) * ni + r] = uint64_t(local[int_local[r]]) << 1;
}
// The split must agree with the DOF classification on every element.
__global__ void k_bddc_check_split(int ne, int nl, int nw, int ni, const int* dofs, const int* wb_of_dof, const int* wb_local, const int* int_local, int* bad) {
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= ne) return;
    const int* local = dofs + size_t(e) * nl;
    for (int k = 0; k < nw; ++k)
        if (wb_of_dof[local[wb_local[k]]] < 0) atomicExch(bad, 1);
    for (int k = 0; k < ni; ++k)
        if (wb_of_dof[local[int_local[k]]] >= 0) atomicExch(bad, 1);
}
__global__ void k_bddc_wb_numbering(int n, const int* wb_of_dof, const int* free_dofs, int* wb_to_dof, int* free_wb) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n || wb_of_dof[i] < 0) return;
    wb_to_dof[wb_of_dof[i]] = i;
    free_wb[wb_of_dof[i]] = free_dofs[i];
}
__global__ void k_gather(int nw, const int* wb_to_dof, const double* x, double* out) {
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < nw) out[k] = x[wb_to_dof[k]];
}
__global__ void k_scatter_add(int nw, const int* wb_to_dof, const double* x, double* out) {
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < nw) out[wb_to_dof[k]] += x[k];
}
__global__ void k_split_cplx(int nnz, const cplx* in, double* re, double* im) {
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < nnz) { re[k] = in[k].re; im[k] = in[k].im; }
}
__global__ void k_scale_rows(int n, const int* rp, const double* weights, double* values) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const double s = weights[i] > 0 ? 1. / weights[i] : 0.;
    for (int k = rp[i]; k < rp[i + 1]; ++k) values[k] *= s;
}
__global__ void k_scale_both(int n, const int* rp, const int* ci, const double* weights, double* values) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i >= n) return;
    const double s = weights[i] > 0 ? 1. / weights[i] : 0.;
    for (int k = rp[i]; k < rp[i + 1]; ++k) {
        const double t = weights[ci[k]] > 0 ? 1. / weights[ci[k]] : 0.;
        values[k] *= s * t;
    }
}
// The weight plan has one column: a row holds its summed weight or nothing.
__global__ void k_expand_weights(int n, const int* rp, const double* compact, double* weights) {
    const int i = blockIdx.x * blockDim.x + threadIdx.x;
    if (i < n && rp[i + 1] > rp[i]) weights[i] = compact[rp[i]];
}
__global__ void k_add_values(int nnz, const double* a, const double* b, double* out) {
    const int k = blockIdx.x * blockDim.x + threadIdx.x;
    if (k < nnz) out[k] = a[k] + b[k];
}

// One sparsity pattern with real and optional imaginary values.
struct CMat {
    int n = 0, m = 0, nnz = 0;
    Dev<int> rowptr, columns;
    Dev<double> re, im;
};
inline void mult(const CMat& a, const double* xr, const double* xi, double* yr, double* yi) {
    if (a.nnz == 0) {
        ck(cudaMemset(yr, 0, size_t(a.n) * sizeof(double)), "BDDC product");
        if (yi) ck(cudaMemset(yi, 0, size_t(a.n) * sizeof(double)), "BDDC product");
        return;
    }
    const int blocks = blocks_of(a.n);
    multiply<<<blocks, TH>>>(a.n, a.rowptr.p, a.columns.p, a.re.p, xr, yr);
    if (a.im.p) {
        multiply_add<<<blocks, TH>>>(a.n, a.rowptr.p, a.columns.p, a.im.p, xi, -1., yr);
        multiply<<<blocks, TH>>>(a.n, a.rowptr.p, a.columns.p, a.re.p, xi, yi);
        multiply_add<<<blocks, TH>>>(a.n, a.rowptr.p, a.columns.p, a.im.p, xr, 1., yi);
    }
}
inline void mult_add(const CMat& a, const double* xr, const double* xi, double alpha, double* yr, double* yi) {
    if (a.nnz == 0) return;
    const int blocks = blocks_of(a.n);
    multiply_add<<<blocks, TH>>>(a.n, a.rowptr.p, a.columns.p, a.re.p, xr, alpha, yr);
    if (a.im.p) {
        multiply_add<<<blocks, TH>>>(a.n, a.rowptr.p, a.columns.p, a.im.p, xi, -alpha, yr);
        multiply_add<<<blocks, TH>>>(a.n, a.rowptr.p, a.columns.p, a.re.p, xi, alpha, yi);
        multiply_add<<<blocks, TH>>>(a.n, a.rowptr.p, a.columns.p, a.im.p, xr, alpha, yi);
    }
}

inline void fill_matrix(Plan& plan, const double* all, CMat& out, const char* name) {
    out.nnz = plan.nnz;
    if (!plan.nnz) return;
    out.re.allocate(size_t(plan.nnz), name);
    reduce(plan, all, out.re.p);
    out.rowptr = std::move(plan.rowptr);
    out.columns = std::move(plan.columns);
    ck_kernels(name);
}
inline void fill_matrix(Plan& plan, const cplx* all, CMat& out, const char* name) {
    out.nnz = plan.nnz;
    if (!plan.nnz) return;
    Dev<cplx> compact(size_t(plan.nnz), name);
    reduce(plan, all, compact.p);
    out.re.allocate(size_t(plan.nnz), name);
    out.im.allocate(size_t(plan.nnz), name);
    k_split_cplx<<<blocks_of(plan.nnz), TH>>>(plan.nnz, compact.p, out.re.p, out.im.p);
    out.rowptr = std::move(plan.rowptr);
    out.columns = std::move(plan.columns);
    ck_kernels(name);
}

class Bddc {
public:
    std::vector<AmgReport> build(Ledger& ledger, Assembled& op, bool complex_values, const std::vector<BddcGroupInput>& split, const int* wb_of_dof,
                                 int nw, int nedge, int nv, const int* edge_vertices, const double* coordinates, int cycles, int min_coarse, double theta) {
        if (int(split.size()) != op.groups.count)
            throw std::runtime_error("device BDDC: the split must name every element group");
        if (nw <= 0 || cycles < 1) throw std::runtime_error("device BDDC: invalid wirebasket size or cycle count");
        if (nedge == 0 && complex_values) throw std::runtime_error("device BDDC: a complex system needs the AMS wirebasket");
        n_ = op.n;
        nw_ = nw;
        cycles_ = cycles;
        free_ = op.free_dofs.p;
        complex_ = complex_values;
        if (nedge > 0) {
            if (nw != nedge) throw std::runtime_error("device BDDC: the AMS wirebasket must be the mesh edges, wirebasket " + std::to_string(nw)
                                                      + " edges " + std::to_string(nedge));
            for (int e = 0; e < nedge; ++e)
                if (wb_of_dof[e] != e) throw std::runtime_error("device BDDC: wirebasket unknown " + std::to_string(e) + " is not edge " + std::to_string(e));
        }
        wb_of_dof_ = upload(ledger, "mesh.wb_of_dof", wb_of_dof, size_t(n_));
        wb_to_dof_.allocate(size_t(nw), "BDDC");
        free_wb_.allocate(size_t(nw), "BDDC");
        ck(cudaMemset(wb_to_dof_.p, 0xff, size_t(nw) * sizeof(int)), "BDDC");
        k_bddc_wb_numbering<<<blocks_of(n_), TH>>>(n_, wb_of_dof_.p, free_, wb_to_dof_.p, free_wb_.p);
        std::vector<Dev<int>> wb_local, int_local;
        Dev<int> bad(1, "BDDC");
        bad.zero();
        for (size_t g = 0; g < split.size(); ++g) {
            wb_local.push_back(upload(ledger, "mesh.bddc_split", split[g].wb_local, size_t(split[g].nw)));
            int_local.push_back(upload(ledger, "mesh.bddc_split", split[g].int_local, size_t(split[g].ni)));
            const GroupView& view = op.groups.g[g];
            if (split[g].nw + split[g].ni != view.nl)
                throw std::runtime_error("device BDDC: the split does not cover the element DOFs");
            k_bddc_check_split<<<blocks_of(view.ne), TH>>>(view.ne, view.nl, split[g].nw, split[g].ni, view.dofs, wb_of_dof_.p, wb_local[g].p, int_local[g].p, bad.p);
        }
        ck_kernels("BDDC split");
        check_flag(ledger, bad, "device BDDC: the element split disagrees with the DOF classification");
        if (complex_values) build_typed<cplx>(ledger, op, split, wb_local, int_local);
        else build_typed<double>(ledger, op, split, wb_local, int_local);
        allocate_work();
        return build_coarse(ledger, nedge, nv, edge_vertices, coordinates, min_coarse, theta);
    }
    void apply(const double* xr, const double* xi, double* yr, double* yi) {
        const int blocks = blocks_of(n_), wblocks = blocks_of(nw_);
        k_gather<<<wblocks, TH>>>(nw_, wb_to_dof_.p, xr, twr_.p);
        if (complex_) k_gather<<<wblocks, TH>>>(nw_, wb_to_dof_.p, xi, twi_.p);
        mult_add(ht_, xr, xi, 1., twr_.p, twi_.p);
        mask<<<wblocks, TH>>>(nw_, free_wb_.p, twr_.p);
        if (complex_) mask<<<wblocks, TH>>>(nw_, free_wb_.p, twi_.p);
        ywr_.zero();
        if (complex_) ywi_.zero();
        k_copy<<<wblocks, TH>>>(nw_, twr_.p, rwr_.p);
        if (complex_) k_copy<<<wblocks, TH>>>(nw_, twi_.p, rwi_.p);
        for (int s = 0; s < cycles_; ++s) {
            coarse(rwr_.p, zwr_.p);
            k_axpy<<<wblocks, TH>>>(nw_, 1., zwr_.p, ywr_.p);
            if (complex_) {
                coarse(rwi_.p, zwi_.p);
                k_axpy<<<wblocks, TH>>>(nw_, 1., zwi_.p, ywi_.p);
            }
            if (s + 1 == cycles_) break;
            mult(s_, ywr_.p, ywi_.p, swr_.p, swi_.p);
            k_copy<<<wblocks, TH>>>(nw_, twr_.p, rwr_.p);
            k_axpy<<<wblocks, TH>>>(nw_, -1., swr_.p, rwr_.p);
            mask<<<wblocks, TH>>>(nw_, free_wb_.p, rwr_.p);
            if (complex_) {
                k_copy<<<wblocks, TH>>>(nw_, twi_.p, rwi_.p);
                k_axpy<<<wblocks, TH>>>(nw_, -1., swi_.p, rwi_.p);
                mask<<<wblocks, TH>>>(nw_, free_wb_.p, rwi_.p);
            }
        }
        mult(inner_, xr, xi, yr, yi);
        mult_add(h_, ywr_.p, ywi_.p, 1., yr, yi);
        k_scatter_add<<<wblocks, TH>>>(nw_, wb_to_dof_.p, ywr_.p, yr);
        mask<<<blocks, TH>>>(n_, free_, yr);
        if (complex_) {
            k_scatter_add<<<wblocks, TH>>>(nw_, wb_to_dof_.p, ywi_.p, yi);
            mask<<<blocks, TH>>>(n_, free_, yi);
        }
        ck_kernels("BDDC apply");
    }

private:
    static long long checked_total(long long total, const char* name) {
        if (total >= INT_MAX - TH) throw std::runtime_error(std::string("device BDDC: ") + name + " triple count exceeds int32 capacity");
        return total;
    }
    template <class S>
    void build_typed(Ledger& ledger, const Assembled& op, const std::vector<BddcGroupInput>& split, const std::vector<Dev<int>>& wb_local,
                     const std::vector<Dev<int>>& int_local) {
        const char* name = "BDDC setup";
        long long total_h = 0, total_d = 0, total_s = 0, total_w = 0;
        std::vector<long long> off_h, off_d, off_s, off_w;
        for (size_t g = 0; g < split.size(); ++g) {
            const long long ne = op.groups.g[g].ne, ni = split[g].ni, nw = split[g].nw;
            off_h.push_back(total_h); off_d.push_back(total_d); off_s.push_back(total_s); off_w.push_back(total_w);
            total_h += ne * ni * nw;
            total_d += ne * ni * ni;
            total_s += ne * nw * nw;
            total_w += ne * ni;
        }
        checked_total(total_h, "harmonic");
        checked_total(total_d, "interior");
        checked_total(total_s, "Schur");
        const int bits_n = bits_for(n_), bits_w = bits_for(nw_);
        Plan plan_h, plan_d, plan_s, plan_w;
        const auto build_plan = [&](long long total, int rows, int cols, int colbits, int which, Plan& plan) {
            if (total == 0) return;
            Dev<uint64_t> keys(size_t(total), name);
            for (size_t g = 0; g < split.size(); ++g) {
                const GroupView& view = op.groups.g[g];
                if (which == 0)
                    k_bddc_keys<<<blocks_of(view.ne), TH>>>(view.ne, view.nl, split[g].ni, split[g].nw, colbits, view.dofs, int_local[g].p, wb_local[g].p,
                                                            nullptr, wb_of_dof_.p, keys.p + off_h[g]);
                else if (which == 1)
                    k_bddc_keys<<<blocks_of(view.ne), TH>>>(view.ne, view.nl, split[g].ni, split[g].ni, colbits, view.dofs, int_local[g].p, int_local[g].p,
                                                            nullptr, nullptr, keys.p + off_d[g]);
                else if (which == 2)
                    k_bddc_keys<<<blocks_of(view.ne), TH>>>(view.ne, view.nl, split[g].nw, split[g].nw, colbits, view.dofs, wb_local[g].p, wb_local[g].p,
                                                            wb_of_dof_.p, wb_of_dof_.p, keys.p + off_s[g]);
                else
                    k_bddc_weight_keys<<<blocks_of(view.ne), TH>>>(view.ne, view.nl, split[g].ni, view.dofs, int_local[g].p, keys.p + off_w[g]);
            }
            ck_kernels(name);
            plan = make_plan(ledger, keys, int(total), rows, cols, name);
        };
        build_plan(total_h, n_, nw_, bits_w, 0, plan_h);
        build_plan(total_d, n_, n_, bits_n, 1, plan_d);
        build_plan(total_s, nw_, nw_, bits_w, 2, plan_s);
        build_plan(total_w, n_, 1, 1, 3, plan_w);
        Dev<S> all_h(size_t(total_h), name), all_d(size_t(total_d), name), all_s(size_t(total_s), name);
        Dev<double> all_w(size_t(total_w), name);
        Dev<int> bad(1, name);
        bad.zero();
        for (size_t g = 0; g < split.size(); ++g) {
            const GroupView& view = op.groups.g[g];
            const int ni = split[g].ni, nw = split[g].nw;
            // Short launches keep every kernel well inside a device watchdog. Elements of a group
            // cost the same, so one block per processor leaves none idle for a launch.
            const int batch = int(std::min<long long>(view.ne, 64LL * device_processors()));
            Dev<S> scratch(size_t(batch) * ni * nw, name);
            for (int first = 0; first < view.ne; first += batch) {
                const int count = std::min(batch, view.ne - first);
                k_bddc_element<S><<<(count + 63) / 64, 64>>>(view, first, count, ni, nw, free_, int_local[g].p, wb_local[g].p, scratch.p,
                                                            all_d.p + off_d[g], all_h.p + off_h[g], all_s.p + off_s[g], all_w.p + off_w[g], bad.p);
                ck_kernels(name);
            }
        }
        check_flag(ledger, bad, "device BDDC: an element interior block has a zero or nonfinite pivot");
        Dev<double> weights(size_t(n_), name);
        weights.zero();
        if (total_w > 0) {
            Dev<double> compact(size_t(plan_w.nnz), name);
            reduce(plan_w, all_w.p, compact.p);
            k_expand_weights<<<blocks_of(n_), TH>>>(n_, plan_w.rowptr.p, compact.p, weights.p);
            ck_kernels(name);
        }
        all_w.release();
        fill_matrix(plan_h, all_h.p, h_, name);
        all_h.release();
        h_.n = n_; h_.m = nw_;
        fill_matrix(plan_d, all_d.p, inner_, name);
        all_d.release();
        inner_.n = n_; inner_.m = n_;
        fill_matrix(plan_s, all_s.p, s_, name);
        all_s.release();
        s_.n = nw_; s_.m = nw_;
        if (h_.nnz) {
            k_scale_rows<<<blocks_of(n_), TH>>>(n_, h_.rowptr.p, weights.p, h_.re.p);
            if (h_.im.p) k_scale_rows<<<blocks_of(n_), TH>>>(n_, h_.rowptr.p, weights.p, h_.im.p);
            TransposeMap map = transpose_map(ledger, n_, nw_, h_.nnz, h_.rowptr.p, h_.columns.p, name);
            ht_.n = nw_; ht_.m = n_; ht_.nnz = h_.nnz;
            ht_.re.allocate(size_t(h_.nnz), name);
            k_permute_double<<<blocks_of(h_.nnz), TH>>>(h_.nnz, map.order.p, h_.re.p, ht_.re.p);
            if (h_.im.p) {
                ht_.im.allocate(size_t(h_.nnz), name);
                k_permute_double<<<blocks_of(h_.nnz), TH>>>(h_.nnz, map.order.p, h_.im.p, ht_.im.p);
            }
            ht_.rowptr = std::move(map.rowptr);
            ht_.columns = std::move(map.columns);
        } else {
            ht_.n = nw_; ht_.m = n_;
        }
        if (inner_.nnz) {
            k_scale_both<<<blocks_of(n_), TH>>>(n_, inner_.rowptr.p, inner_.columns.p, weights.p, inner_.re.p);
            if (inner_.im.p) k_scale_both<<<blocks_of(n_), TH>>>(n_, inner_.rowptr.p, inner_.columns.p, weights.p, inner_.im.p);
        }
        ck_kernels(name);
        if (s_.nnz == 0) throw std::runtime_error("device BDDC: the wirebasket Schur complement is empty");
    }
    std::vector<AmgReport> build_coarse(Ledger& ledger, int nedge, int nv, const int* edge_vertices, const double* coordinates, int min_coarse, double theta) {
        const char* name = "BDDC coarse";
        Dev<double> surrogate(size_t(s_.nnz), name);
        if (s_.im.p) k_add_values<<<blocks_of(s_.nnz), TH>>>(s_.nnz, s_.re.p, s_.im.p, surrogate.p);
        else ck(cudaMemcpy(surrogate.p, s_.re.p, size_t(s_.nnz) * sizeof(double), cudaMemcpyDeviceToDevice), name);
        ck_kernels(name);
        if (nedge > 0) {
            CsrView view{nw_, nw_, s_.rowptr.p, s_.columns.p, surrogate.p};
            ams_ = std::make_unique<Ams>();
            return ams_->build(ledger, view, s_.nnz, free_wb_.p, surrogate.p, nedge, nv, edge_vertices, coordinates, min_coarse, theta);
        }
        Csr a;
        a.n = a.m = nw_;
        a.nnz = s_.nnz;
        a.rowptr.allocate(size_t(nw_) + 1, name);
        a.columns.allocate(size_t(s_.nnz), name);
        a.values.allocate(size_t(s_.nnz), name);
        ck(cudaMemcpy(a.rowptr.p, s_.rowptr.p, (size_t(nw_) + 1) * sizeof(int), cudaMemcpyDeviceToDevice), name);
        ck(cudaMemcpy(a.columns.p, s_.columns.p, size_t(s_.nnz) * sizeof(int), cudaMemcpyDeviceToDevice), name);
        k_constrain_values<<<blocks_of(nw_), TH>>>(nw_, s_.rowptr.p, s_.columns.p, free_wb_.p, surrogate.p, a.values.p);
        ck_kernels(name);
        AmgOptions options;
        options.theta = theta;
        options.min_coarse = min_coarse;
        options.stall = 0.8;
        amg_ = std::make_unique<Amg>();
        return {amg_->build(ledger, std::move(a), options, name)};
    }
    void allocate_work() {
        const size_t w = size_t(nw_);
        const char* name = "BDDC work";
        twr_.allocate(w, name); ywr_.allocate(w, name); rwr_.allocate(w, name); zwr_.allocate(w, name); swr_.allocate(w, name);
        if (complex_) { twi_.allocate(w, name); ywi_.allocate(w, name); rwi_.allocate(w, name); zwi_.allocate(w, name); swi_.allocate(w, name); }
    }
    void coarse(const double* r, double* z) {
        if (ams_) { ams_->apply(r, z); return; }
        const int blocks = blocks_of(nw_);
        k_copy<<<blocks, TH>>>(nw_, r, amg_->input());
        amg_->cycle();
        k_copy<<<blocks, TH>>>(nw_, amg_->output(), z);
        mask<<<blocks, TH>>>(nw_, free_wb_.p, z);
    }
    int n_ = 0, nw_ = 0, cycles_ = 1;
    bool complex_ = false;
    const int* free_ = nullptr;
    Dev<int> wb_of_dof_, wb_to_dof_, free_wb_;
    CMat h_, ht_, inner_, s_;
    std::unique_ptr<Ams> ams_;
    std::unique_ptr<Amg> amg_;
    Dev<double> twr_, twi_, ywr_, ywi_, rwr_, rwi_, zwr_, zwi_, swr_, swi_;
};

}  // namespace radia_fem
