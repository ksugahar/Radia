/* SPDX-License-Identifier: MPL-2.0 */
// Device implementation behind the Python module of the FEM route. All kernels
// live in the fem_*.cuh headers; this file owns lifetime, stage timing and the
// mapping of the host API to device operations.
#include <cmath>
#include <chrono>
#include "fem_api.hpp"
#include "fem_device.cuh"
#include "fem_elements.cuh"
#include "fem_amg.cuh"
#include "fem_ams.cuh"
#include "fem_krylov.cuh"
#include "fem_bddc.cuh"

namespace radia_fem {

namespace {
class Stage {
public:
    Stage(std::vector<std::pair<std::string, double>>& table, const char* name)
        : table_(table), name_(name), start_(std::chrono::steady_clock::now()) {}
    ~Stage() {
        cudaDeviceSynchronize();
        const std::chrono::duration<double> elapsed = std::chrono::steady_clock::now() - start_;
        table_.emplace_back(name_, elapsed.count());
    }
private:
    std::vector<std::pair<std::string, double>>& table_;
    const char* name_;
    std::chrono::steady_clock::time_point start_;
};
}  // namespace

class System::Impl {
public:
    Ledger ledger;
    Assembled op;
    bool complex_values = false;
    std::unique_ptr<Amg> amg;
    std::unique_ptr<Ams> ams;
    std::unique_ptr<Bddc> bddc;
    Dev<int> edge_vertices;
    Dev<double> coordinates;
    std::vector<std::pair<std::string, double>> timing;

    const double* values(int which) const {
        if (which == 0) return op.re.p;
        if (which == 1 && op.im.p) return op.im.p;
        if (which == 2) return op.aux.p;
        throw std::runtime_error("FEM device route: no operator with index " + std::to_string(which));
    }
    void action(const double* x, double* y) {
        const int n = op.n, blocks = blocks_of(n);
        multiply<<<blocks, TH>>>(n, op.plan.rowptr.p, op.plan.columns.p, op.re.p, x, y);
        mask<<<blocks, TH>>>(n, op.free_dofs.p, y);
    }
    void complex_action(const double* xr, const double* xi, double* yr, double* yi) {
        const int n = op.n, blocks = blocks_of(n);
        const int* rp = op.plan.rowptr.p;
        const int* ci = op.plan.columns.p;
        multiply<<<blocks, TH>>>(n, rp, ci, op.re.p, xr, yr);
        multiply_add<<<blocks, TH>>>(n, rp, ci, op.im.p, xi, -1., yr);
        multiply<<<blocks, TH>>>(n, rp, ci, op.re.p, xi, yi);
        multiply_add<<<blocks, TH>>>(n, rp, ci, op.im.p, xr, 1., yi);
        mask<<<blocks, TH>>>(n, op.free_dofs.p, yr);
        mask<<<blocks, TH>>>(n, op.free_dofs.p, yi);
    }
    void precondition(const double* r, double* z) {
        const int n = op.n, blocks = blocks_of(n);
        if (amg) {
            k_copy<<<blocks, TH>>>(n, r, amg->input());
            amg->cycle();
            k_copy<<<blocks, TH>>>(n, amg->output(), z);
            mask<<<blocks, TH>>>(n, op.free_dofs.p, z);
        } else if (ams) {
            ams->apply(r, z);
        } else if (bddc) {
            bddc->apply(r, nullptr, z, nullptr);
        } else {
            throw std::runtime_error("FEM device route: no preconditioner was built");
        }
    }
    void complex_precondition(const double* rr, const double* ri, double* zr, double* zi) {
        if (ams) {
            ams->apply(rr, zr);
            ams->apply(ri, zi);
        } else if (bddc) {
            bddc->apply(rr, ri, zr, zi);
        } else {
            throw std::runtime_error("FEM device route: a complex system needs AMS or BDDC");
        }
    }
    void require_one_preconditioner() const {
        if (amg || ams || bddc) throw std::runtime_error("FEM device route: a preconditioner was already built");
    }
    void upload_edges(int nedge, int nv, const int* edges, const double* coords) {
        for (int e = 0; e < nedge; ++e)
            if (edges[2 * e] < 0 || edges[2 * e + 1] >= nv || edges[2 * e] >= edges[2 * e + 1])
                throw std::runtime_error("FEM device route: edge " + std::to_string(e) + " is not ordered low to high within the vertices");
        for (int i = 0; i < 3 * nv; ++i)
            if (!std::isfinite(coords[i])) throw std::runtime_error("FEM device route: nonfinite vertex coordinate");
        edge_vertices = upload(ledger, "mesh.edge_vertices", edges, size_t(nedge) * 2);
        coordinates = upload(ledger, "mesh.coordinates", coords, size_t(nv) * 3);
    }
};

System::System(int ndof, bool complex_values, const std::vector<GroupInput>& groups, const int* free_dofs)
    : impl_(new Impl) {
    for (int i = 0; i < ndof; ++i)
        if (free_dofs[i] != 0 && free_dofs[i] != 1) throw std::runtime_error("FEM device route: the free DOF mask must hold 0 or 1");
    impl_->complex_values = complex_values;
    device_memory().peak = device_memory().live;
    Stage stage(impl_->timing, "assembly_seconds");
    assemble(impl_->ledger, impl_->op, ndof, complex_values, groups, free_dofs);
}
System::~System() = default;

AmgReport System::setup_amg(int min_coarse, double theta) {
    Impl& m = *impl_;
    m.require_one_preconditioner();
    if (m.complex_values) throw std::runtime_error("FEM device route: AMG applies to a real system");
    Stage stage(m.timing, "preconditioner_seconds");
    Csr a;
    const int n = m.op.n, nnz = m.op.plan.nnz;
    a.n = a.m = n;
    a.nnz = nnz;
    a.rowptr.allocate(size_t(n) + 1, "AMG operator");
    a.columns.allocate(size_t(nnz), "AMG operator");
    a.values.allocate(size_t(nnz), "AMG operator");
    ck(cudaMemcpy(a.rowptr.p, m.op.plan.rowptr.p, (size_t(n) + 1) * sizeof(int), cudaMemcpyDeviceToDevice), "AMG operator");
    ck(cudaMemcpy(a.columns.p, m.op.plan.columns.p, size_t(nnz) * sizeof(int), cudaMemcpyDeviceToDevice), "AMG operator");
    k_constrain_values<<<blocks_of(n), TH>>>(n, m.op.plan.rowptr.p, m.op.plan.columns.p, m.op.free_dofs.p, m.op.aux.p, a.values.p);
    ck_kernels("AMG operator");
    AmgOptions options;
    options.theta = theta;
    options.min_coarse = min_coarse;
    auto hierarchy = std::make_unique<Amg>();
    AmgReport report = hierarchy->build(m.ledger, std::move(a), options, "AMG setup");
    m.amg = std::move(hierarchy);
    return report;
}

std::vector<AmgReport> System::setup_ams(int nedge, int nv, const int* edges, const double* coords, int min_coarse, double theta,
                                         const AmsOptions&) {
    Impl& m = *impl_;
    m.require_one_preconditioner();
    Stage stage(m.timing, "preconditioner_seconds");
    m.upload_edges(nedge, nv, edges, coords);
    auto route = std::make_unique<Ams>();
    std::vector<AmgReport> reports = route->build(m.ledger, m.op.view(m.op.aux.p), m.op.plan.nnz, m.op.free_dofs.p, m.op.aux.p, nedge, nv,
                                                  m.edge_vertices.p, m.coordinates.p, min_coarse, theta);
    m.ams = std::move(route);
    return reports;
}

std::vector<AmgReport> System::setup_bddc(const std::vector<BddcGroupInput>& split, const int* wb_of_dof, int nw, int nedge, int nv,
                                          const int* edges, const double* coords, int cycles, int min_coarse, double theta) {
    Impl& m = *impl_;
    m.require_one_preconditioner();
    Stage stage(m.timing, "preconditioner_seconds");
    if (nedge > 0) m.upload_edges(nedge, nv, edges, coords);
    auto route = std::make_unique<Bddc>();
    std::vector<AmgReport> reports = route->build(m.ledger, m.op, m.complex_values, split, wb_of_dof, nw, nedge, nv,
                                                  m.edge_vertices.p, m.coordinates.p, cycles, min_coarse, theta);
    m.bddc = std::move(route);
    return reports;
}

void System::update_tensors(const std::vector<const double*>& tensors) {
    Impl& m = *impl_;
    if (tensors.size() != m.op.storage.size())
        throw std::runtime_error("FEM device route: one tensor entry per element group is required");
    m.amg.reset();
    m.ams.reset();
    m.bddc.reset();
    Stage stage(m.timing, "tensor_update_seconds");
    for (size_t g = 0; g < tensors.size(); ++g) {
        if (!tensors[g]) continue;
        FormStorage& form = m.op.storage[g].stiff;
        if (!form.geometry.p) throw std::runtime_error("FEM device route: this group was not built with a point tensor");
        const int ne = m.op.storage[g].ne, points = m.op.groups.g[g].stiff.nq;
        Dev<double> tensor = upload(m.ledger, "coef.tensor.stiff", tensors[g], size_t(ne) * points * 6);
        k_tensor_metric<<<blocks_of(ne), TH>>>(ne, points, form.geometry.p, tensor.p, form.metric.p);
        ck_kernels("tensor metric");
    }
    assemble_values(m.op, m.complex_values);
}

int System::nonzeros() const { return impl_->op.plan.nnz; }

void System::probe_action(int which, const double* x, double* y) {
    Impl& m = *impl_;
    const int n = m.op.n;
    const double* values = m.values(which);
    Dev<double> in = upload(m.ledger, "probe.input", x, size_t(n));
    Dev<double> out(size_t(n), "probe");
    multiply<<<blocks_of(n), TH>>>(n, m.op.plan.rowptr.p, m.op.plan.columns.p, values, in.p, out.p);
    ck_kernels("probe action");
    download(m.ledger, "probe.output", out.p, y, size_t(n));
}

void System::probe_preconditioner(const double* x_re, const double* x_im, double* y_re, double* y_im) {
    Impl& m = *impl_;
    const int n = m.op.n;
    Dev<double> in_re = upload(m.ledger, "probe.input", x_re, size_t(n));
    Dev<double> out_re(size_t(n), "probe");
    if (m.complex_values) {
        Dev<double> in_im = upload(m.ledger, "probe.input", x_im, size_t(n));
        Dev<double> out_im(size_t(n), "probe");
        m.complex_precondition(in_re.p, in_im.p, out_re.p, out_im.p);
        ck_kernels("probe preconditioner");
        download(m.ledger, "probe.output", out_re.p, y_re, size_t(n));
        download(m.ledger, "probe.output", out_im.p, y_im, size_t(n));
    } else {
        m.precondition(in_re.p, out_re.p);
        ck_kernels("probe preconditioner");
        download(m.ledger, "probe.output", out_re.p, y_re, size_t(n));
    }
}

SolveReport System::solve(const double* rhs_re, const double* rhs_im, const double* lift_re, const double* lift_im,
                          double tolerance, int max_iterations, double relative_limit, double* x_re, double* x_im) {
    Impl& m = *impl_;
    const int n = m.op.n, blocks = blocks_of(n);
    const int* free_dofs = m.op.free_dofs.p;
    if (!(tolerance > 0) || max_iterations <= 0) throw std::runtime_error("FEM device solve: invalid tolerance or iteration limit");
    Stage stage(m.timing, "solve_seconds");
    Dev<double> br = upload(m.ledger, "rhs.re", rhs_re, size_t(n));
    Dev<double> xr(size_t(n), "solve"), yr(size_t(n), "solve"), gr(size_t(n), "solve");
    xr.zero();
    gr.zero();
    if (lift_re) ck(cudaMemcpy(gr.p, lift_re, size_t(n) * sizeof(double), cudaMemcpyHostToDevice), "lift.re");
    m.ledger.add("lift.re", "H2D", lift_re ? size_t(n) * sizeof(double) : 0);
    SolveReport report;
    if (!m.complex_values) {
        // b = f - A g on the free DOFs, x restricted to the free DOFs.
        multiply<<<blocks, TH>>>(n, m.op.plan.rowptr.p, m.op.plan.columns.p, m.op.re.p, gr.p, yr.p);
        k_axpy<<<blocks, TH>>>(n, -1., yr.p, br.p);
        mask<<<blocks, TH>>>(n, free_dofs, br.p);
        report = pcg(m.ledger, n, [&](const double* a, double* b) { m.action(a, b); },
                     [&](const double* a, double* b) { m.precondition(a, b); }, br.p, xr.p, tolerance, max_iterations);
        m.action(xr.p, yr.p);
        k_axpy<<<blocks, TH>>>(n, -1., br.p, yr.p);
        Blas blas(m.ledger);
        report.true_residual = std::sqrt(blas.dot(n, yr.p, yr.p)) / report.initial_norm;
        if (!(report.true_residual <= relative_limit))
            throw std::runtime_error("FEM device solve: true relative residual " + std::to_string(report.true_residual)
                                     + " exceeds the limit " + std::to_string(relative_limit));
        k_axpy<<<blocks, TH>>>(n, 1., gr.p, xr.p);
        ck_kernels("solve");
        download(m.ledger, "result.re", xr.p, x_re, size_t(n));
        return report;
    }
    if (!rhs_im || !x_im) throw std::runtime_error("FEM device solve: a complex system needs complex data");
    Dev<double> bi = upload(m.ledger, "rhs.im", rhs_im, size_t(n));
    Dev<double> xi(size_t(n), "solve"), yi(size_t(n), "solve"), gi(size_t(n), "solve");
    xi.zero();
    gi.zero();
    if (lift_im) ck(cudaMemcpy(gi.p, lift_im, size_t(n) * sizeof(double), cudaMemcpyHostToDevice), "lift.im");
    m.ledger.add("lift.im", "H2D", lift_im ? size_t(n) * sizeof(double) : 0);
    const int* rp = m.op.plan.rowptr.p;
    const int* ci = m.op.plan.columns.p;
    multiply<<<blocks, TH>>>(n, rp, ci, m.op.re.p, gr.p, yr.p);
    multiply_add<<<blocks, TH>>>(n, rp, ci, m.op.im.p, gi.p, -1., yr.p);
    multiply<<<blocks, TH>>>(n, rp, ci, m.op.re.p, gi.p, yi.p);
    multiply_add<<<blocks, TH>>>(n, rp, ci, m.op.im.p, gr.p, 1., yi.p);
    k_axpy<<<blocks, TH>>>(n, -1., yr.p, br.p);
    k_axpy<<<blocks, TH>>>(n, -1., yi.p, bi.p);
    mask<<<blocks, TH>>>(n, free_dofs, br.p);
    mask<<<blocks, TH>>>(n, free_dofs, bi.p);
    report = cocr(m.ledger, n, [&](const double* a, const double* b, double* c, double* d) { m.complex_action(a, b, c, d); },
                  [&](const double* a, const double* b, double* c, double* d) { m.complex_precondition(a, b, c, d); },
                  br.p, bi.p, xr.p, xi.p, tolerance, max_iterations);
    m.complex_action(xr.p, xi.p, yr.p, yi.p);
    k_axpy<<<blocks, TH>>>(n, -1., br.p, yr.p);
    k_axpy<<<blocks, TH>>>(n, -1., bi.p, yi.p);
    Blas blas(m.ledger);
    report.true_residual = std::sqrt(blas.dot(n, yr.p, yr.p) + blas.dot(n, yi.p, yi.p)) / report.initial_norm;
    if (!(report.true_residual <= relative_limit))
        throw std::runtime_error("FEM device solve: true relative residual " + std::to_string(report.true_residual)
                                 + " exceeds the limit " + std::to_string(relative_limit));
    k_axpy<<<blocks, TH>>>(n, 1., gr.p, xr.p);
    k_axpy<<<blocks, TH>>>(n, 1., gi.p, xi.p);
    ck_kernels("solve");
    download(m.ledger, "result.re", xr.p, x_re, size_t(n));
    download(m.ledger, "result.im", xi.p, x_im, size_t(n));
    return report;
}

std::vector<LedgerRow> System::ledger() const {
    std::vector<LedgerRow> rows;
    for (const auto& item : impl_->ledger.rows()) {
        LedgerRow row;
        row.payload = item.first.first;
        row.direction = item.first.second;
        row.bytes = item.second;
        rows.push_back(row);
    }
    return rows;
}
uint64_t System::peak_device_bytes() const { return device_memory().peak; }
std::vector<std::pair<std::string, double>> System::timings() const { return impl_->timing; }

std::string device_name() {
    int device = 0;
    ck(cudaGetDevice(&device), "cudaGetDevice");
    cudaDeviceProp properties;
    ck(cudaGetDeviceProperties(&properties, device), "cudaGetDeviceProperties");
    return properties.name;
}

}  // namespace radia_fem
