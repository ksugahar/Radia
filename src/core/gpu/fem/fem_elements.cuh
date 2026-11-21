/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
// Element forms on the device. The host supplies mesh arrays and reference
// tabulations; metrics, the sparsity pattern and every matrix value are computed
// here. Values are summed per nonzero through the triple plan, so the result does
// not depend on thread order and no floating-point atomics are needed.
#include "fem_api.hpp"
#include "fem_device.cuh"

namespace radia_fem {

constexpr int MAX_GROUPS = 8;

// Six symmetric components (00,01,02,11,12,22) with the off-diagonal reference
// terms summed into one coefficient.
struct FormView {
    int nq = 0;
    int affine = 0;
    const double* metric = nullptr;  // ne * (affine ? 1 : nq) * 6
    const double* table = nullptr;   // affine: ncls * nl * nl * 6; else ncls * nq * nl * 3
};
struct GroupView {
    int ne = 0, nl = 0, first_triple = 0, reserved = 0;
    const int* dofs = nullptr;
    const int* cls = nullptr;
    const double* coef = nullptr;
    FormView stiff, mass;
};
struct Groups {
    int count = 0;
    GroupView g[MAX_GROUPS];
};

__device__ inline bool metric_at(int law, int dim, const double* J, double* g) {
    double c00 = 0, c01 = 0, c02 = 0, c11 = 0, c12 = 0, c22 = 0;
    for (int i = 0; i < 3; ++i) {
        const double a = J[i * dim], b = J[i * dim + 1], d = dim == 3 ? J[i * dim + 2] : 0.;
        c00 += a * a; c01 += a * b; c02 += a * d; c11 += b * b; c12 += b * d; c22 += d * d;
    }
    for (int k = 0; k < 6; ++k) g[k] = 0;
    double det, i00, i01, i02, i11, i12, i22;
    if (dim == 3) {
        i00 = c11 * c22 - c12 * c12;
        i01 = c02 * c12 - c01 * c22;
        i02 = c01 * c12 - c02 * c11;
        det = c00 * i00 + c01 * i01 + c02 * i02;
        i11 = c00 * c22 - c02 * c02;
        i12 = c01 * c02 - c00 * c12;
        i22 = c00 * c11 - c01 * c01;
    } else {
        det = c00 * c11 - c01 * c01;
        i00 = c11; i01 = -c01; i02 = 0; i11 = c00; i12 = 0; i22 = 0;
    }
    if (!(det > 0) || !isfinite(det)) return false;
    const double measure = sqrt(det), scale = measure / det;
    if (law == LAW_VALUE) {
        g[0] = measure;
    } else if (law == LAW_CURL) {
        g[0] = c00 / measure; g[1] = c01 / measure; g[2] = c02 / measure;
        g[3] = c11 / measure; g[4] = c12 / measure; g[5] = c22 / measure;
    } else {  // LAW_GRAD and LAW_COVARIANT share measure * inverse(J^T J).
        g[0] = i00 * scale; g[1] = i01 * scale; g[2] = i02 * scale;
        g[3] = i11 * scale; g[4] = i12 * scale; g[5] = i22 * scale;
    }
    return true;
}

// A straight simplex has one metric and no weight in it: the weight is in the
// reference tensor. Every other element has one weighted metric per point.
__global__ void k_metric(int ne, int points, int ngeo, int dim, int law, const double* nodes, const double* dgeo,
                         const double* jacobian, const double* weights, double* metric, int* bad) {
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= ne) return;
    for (int q = 0; q < points; ++q) {
        double J[9];
        for (int i = 0; i < 3 * dim; ++i) J[i] = 0;
        if (jacobian) {
            const double* row = jacobian + (size_t(e) * points + q) * 3 * dim;
            for (int i = 0; i < 3 * dim; ++i) J[i] = row[i];
        } else {
            for (int n = 0; n < ngeo; ++n) {
                const double* X = nodes + (size_t(e) * ngeo + n) * 3;
                const double* d = dgeo + (size_t(q) * ngeo + n) * dim;
                for (int i = 0; i < 3; ++i)
                    for (int j = 0; j < dim; ++j) J[i * dim + j] += X[i] * d[j];
            }
        }
        double* g = metric + (size_t(e) * points + q) * 6;
        if (!metric_at(law, dim, J, g)) atomicExch(bad, 1);
        const double w = weights ? weights[q] : 1.;
        for (int k = 0; k < 6; ++k) g[k] *= w;
    }
}

__global__ void k_reference_tensor(int ncls, int nl, int nq, const double* table, const double* weights, double* tensor) {
    const long long index = (long long)(blockIdx.x) * blockDim.x + threadIdx.x;
    const long long total = (long long)(ncls) * nl * nl;
    if (index >= total) return;
    const int c = int(index / (nl * nl)), ab = int(index % (nl * nl)), a = ab / nl, b = ab % nl;
    double t[6] = {0, 0, 0, 0, 0, 0};
    for (int q = 0; q < nq; ++q) {
        const double* da = table + ((size_t(c) * nq + q) * nl + a) * 3;
        const double* db = table + ((size_t(c) * nq + q) * nl + b) * 3;
        const double w = weights[q];
        t[0] += w * da[0] * db[0];
        t[1] += w * (da[0] * db[1] + da[1] * db[0]);
        t[2] += w * (da[0] * db[2] + da[2] * db[0]);
        t[3] += w * da[1] * db[1];
        t[4] += w * (da[1] * db[2] + da[2] * db[1]);
        t[5] += w * da[2] * db[2];
    }
    for (int k = 0; k < 6; ++k) tensor[index * 6 + k] = t[k];
}

__device__ inline double form_entry(const FormView& f, int nl, int c, int e, int a, int b) {
    if (f.affine) {
        const double* t = f.table + (size_t(c) * nl * nl + size_t(a) * nl + b) * 6;
        const double* g = f.metric + size_t(e) * 6;
        return g[0] * t[0] + g[1] * t[1] + g[2] * t[2] + g[3] * t[3] + g[4] * t[4] + g[5] * t[5];
    }
    double sum = 0;
    for (int q = 0; q < f.nq; ++q) {
        const double* da = f.table + ((size_t(c) * f.nq + q) * nl + a) * 3;
        const double* db = f.table + ((size_t(c) * f.nq + q) * nl + b) * 3;
        const double* g = f.metric + (size_t(e) * f.nq + q) * 6;
        sum += g[0] * da[0] * db[0] + g[1] * (da[0] * db[1] + da[1] * db[0]) + g[2] * (da[0] * db[2] + da[2] * db[0])
             + g[3] * da[1] * db[1] + g[4] * (da[1] * db[2] + da[2] * db[1]) + g[5] * da[2] * db[2];
    }
    return sum;
}

__global__ void k_triple_keys(int ne, int nl, int colbits, const int* dofs, uint64_t* keys) {
    const int e = blockIdx.x * blockDim.x + threadIdx.x;
    if (e >= ne) return;
    const int* local = dofs + size_t(e) * nl;
    uint64_t* out = keys + size_t(e) * nl * nl;
    for (int a = 0; a < nl; ++a)
        for (int b = 0; b < nl; ++b)
            out[size_t(a) * nl + b] = (uint64_t(local[a]) << colbits) | uint64_t(local[b]);
}

__global__ void k_assemble_form(Groups groups, int nnz, const int* starts, const int* order, int stiff_column,
                                int mass_column, double* out) {
    const int u = blockIdx.x * blockDim.x + threadIdx.x;
    if (u >= nnz) return;
    double sum = 0;
    for (int k = starts[u]; k < starts[u + 1]; ++k) {
        const int slot = order[k];
        int index = 0;
        while (index + 1 < groups.count && slot >= groups.g[index + 1].first_triple) ++index;
        const GroupView& group = groups.g[index];
        const int local = slot - group.first_triple, nl = group.nl, nl2 = nl * nl;
        const int e = local / nl2, ab = local - e * nl2, a = ab / nl, b = ab - a * nl;
        const double cs = group.coef[size_t(e) * COLUMNS + stiff_column];
        const double cm = group.coef[size_t(e) * COLUMNS + mass_column];
        const int c = group.cls[e];
        if (cs != 0.) sum += cs * form_entry(group.stiff, nl, c, e, a, b);
        if (cm != 0.) sum += cm * form_entry(group.mass, nl, c, e, a, b);
    }
    out[u] = sum;
}

struct FormStorage {
    Dev<double> metric, table;
    FormView view() const {
        FormView v;
        v.metric = metric.p;
        v.table = table.p;
        return v;
    }
};
struct GroupStorage {
    Dev<int> dofs, cls;
    Dev<double> coef;
    FormStorage stiff, mass;
    int ne = 0, nl = 0, ncls = 0;
};

inline void build_form(Ledger& ledger, const GroupInput& input, const Dev<double>& nodes, const FormInput& form,
                       const char* prefix, FormStorage& out, FormView& view) {
    view = FormView{};
    if (form.law == LAW_NONE) return;
    if (input.dim == 2 && form.law == LAW_CURL)
        throw std::runtime_error("FEM device form: a boundary group has no curl form");
    const std::string tag(prefix);
    if (form.nq <= 0 || !form.weights || !form.table)
        throw std::runtime_error("FEM device form: quadrature and tables are required for " + tag);
    if (form.affine && (!nodes.p || form.jacobian))
        throw std::runtime_error("FEM device form: an affine form needs mesh nodes and no Jacobian table");
    if (!form.jacobian && (!nodes.p || !form.dgeo))
        throw std::runtime_error("FEM device form: mesh nodes with geometry derivatives or a Jacobian table are required");
    const int points = form.affine ? 1 : form.nq;
    const size_t table_values = size_t(input.ncls) * form.nq * input.nl * 3;
    Dev<double> weights = upload(ledger, ("table.weights." + tag).c_str(), form.weights, size_t(form.nq));
    Dev<double> table = upload(ledger, ("table.shape." + tag).c_str(), form.table, table_values);
    Dev<double> dgeo, jacobian;
    if (form.jacobian) {
        jacobian = upload(ledger, ("mesh.jacobian." + tag).c_str(), form.jacobian,
                          size_t(input.ne) * points * 3 * input.dim);
    } else {
        dgeo = upload(ledger, ("table.geometry." + tag).c_str(), form.dgeo, size_t(points) * input.ngeo * input.dim);
    }
    out.metric.allocate(size_t(input.ne) * points * 6, "metric");
    Dev<int> bad(1, "metric flag");
    bad.zero();
    k_metric<<<blocks_of(input.ne), TH>>>(input.ne, points, input.ngeo, input.dim, form.law, nodes.p, dgeo.p,
                                          jacobian.p, form.affine ? nullptr : weights.p, out.metric.p, bad.p);
    ck_kernels("metric");
    check_flag(ledger, bad, "FEM device form: an element has a non-positive or nonfinite measure");
    view.nq = points;
    view.affine = form.affine ? 1 : 0;
    if (form.affine) {
        const long long total = (long long)(input.ncls) * input.nl * input.nl;
        out.table.allocate(size_t(total) * 6, "reference tensor");
        k_reference_tensor<<<blocks_of(total), TH>>>(input.ncls, input.nl, form.nq, table.p, weights.p, out.table.p);
        ck_kernels("reference tensor");
    } else {
        out.table = std::move(table);
    }
    view.metric = out.metric.p;
    view.table = out.table.p;
}

struct Assembled {
    int n = 0;
    Plan plan;
    Groups groups;
    std::vector<GroupStorage> storage;
    Dev<int> free_dofs;
    Dev<double> re, im, aux;
    CsrView view(const double* values) const { return CsrView{n, n, plan.rowptr.p, plan.columns.p, values}; }
};

inline void assemble(Ledger& ledger, Assembled& out, int ndof, bool complex_values,
                     const std::vector<GroupInput>& inputs, const int* free_dofs) {
    if (inputs.empty() || int(inputs.size()) > MAX_GROUPS)
        throw std::runtime_error("FEM device assembly: unsupported number of element groups " + std::to_string(inputs.size()));
    if (ndof <= 0) throw std::runtime_error("FEM device assembly: empty space");
    out.n = ndof;
    out.free_dofs = upload(ledger, "mesh.free_dofs", free_dofs, size_t(ndof));
    long long total_triples = 0;
    out.storage.resize(inputs.size());
    out.groups.count = int(inputs.size());
    for (size_t g = 0; g < inputs.size(); ++g) {
        const GroupInput& input = inputs[g];
        if (input.ne <= 0 || input.nl <= 0 || input.ncls <= 0 || !input.dofs || !input.cls || !input.coef)
            throw std::runtime_error("FEM device assembly: incomplete element group");
        if (input.dim != 2 && input.dim != 3)
            throw std::runtime_error("FEM device assembly: unsupported element dimension " + std::to_string(input.dim));
        const long long triples = (long long)(input.ne) * input.nl * input.nl;
        if (total_triples + triples >= INT_MAX - TH)
            throw std::runtime_error("FEM device assembly: triple count exceeds int32 capacity");
        for (long long i = 0; i < (long long)(input.ne) * input.nl; ++i)
            if (input.dofs[i] < 0 || input.dofs[i] >= ndof)
                throw std::runtime_error("FEM device assembly: DOF number outside the space " + std::to_string(input.dofs[i]));
        for (int e = 0; e < input.ne; ++e) {
            if (input.cls[e] < 0 || input.cls[e] >= input.ncls)
                throw std::runtime_error("FEM device assembly: orientation class outside the table");
            const double* coef = input.coef + size_t(e) * COLUMNS;
            for (int k = 0; k < COLUMNS; ++k)
                if (!std::isfinite(coef[k])) throw std::runtime_error("FEM device assembly: nonfinite coefficient");
            if (input.stiff.law == LAW_NONE && (coef[STIFF_RE] != 0. || coef[STIFF_IM] != 0. || coef[STIFF_AUX] != 0.))
                throw std::runtime_error("FEM device assembly: stiffness coefficient without a stiffness form");
            if (input.mass.law == LAW_NONE && (coef[MASS_RE] != 0. || coef[MASS_IM] != 0. || coef[MASS_AUX] != 0.))
                throw std::runtime_error("FEM device assembly: mass coefficient without a mass form");
        }
        GroupStorage& store = out.storage[g];
        store.ne = input.ne; store.nl = input.nl; store.ncls = input.ncls;
        store.dofs = upload(ledger, "mesh.dofs", input.dofs, size_t(input.ne) * input.nl);
        store.cls = upload(ledger, "mesh.classes", input.cls, size_t(input.ne));
        store.coef = upload(ledger, "coef.elements", input.coef, size_t(input.ne) * COLUMNS);
        GroupView& view = out.groups.g[g];
        view.ne = input.ne; view.nl = input.nl; view.first_triple = int(total_triples);
        view.dofs = store.dofs.p; view.cls = store.cls.p; view.coef = store.coef.p;
        Dev<double> nodes;
        if (input.nodes) nodes = upload(ledger, "mesh.nodes", input.nodes, size_t(input.ne) * input.ngeo * 3);
        build_form(ledger, input, nodes, input.stiff, "stiff", store.stiff, view.stiff);
        build_form(ledger, input, nodes, input.mass, "mass", store.mass, view.mass);
        total_triples += triples;
    }
    const int triples = int(total_triples);
    const int colbits = bits_for(ndof);
    Dev<uint64_t> keys(size_t(triples), "triple keys");
    for (size_t g = 0; g < inputs.size(); ++g) {
        const GroupView& view = out.groups.g[g];
        k_triple_keys<<<blocks_of(view.ne), TH>>>(view.ne, view.nl, colbits, view.dofs, keys.p + view.first_triple);
    }
    ck_kernels("triple keys");
    out.plan = make_plan(ledger, keys, triples, ndof, ndof, "assembly plan");
    const int nnz = out.plan.nnz;
    out.re.allocate(size_t(nnz), "operator re");
    k_assemble_form<<<blocks_of(nnz), TH>>>(out.groups, nnz, out.plan.starts.p, out.plan.order.p, STIFF_RE, MASS_RE, out.re.p);
    if (complex_values) {
        out.im.allocate(size_t(nnz), "operator im");
        k_assemble_form<<<blocks_of(nnz), TH>>>(out.groups, nnz, out.plan.starts.p, out.plan.order.p, STIFF_IM, MASS_IM, out.im.p);
    }
    out.aux.allocate(size_t(nnz), "operator aux");
    k_assemble_form<<<blocks_of(nnz), TH>>>(out.groups, nnz, out.plan.starts.p, out.plan.order.p, STIFF_AUX, MASS_AUX, out.aux.p);
    ck_kernels("assembly");
}

}  // namespace radia_fem
