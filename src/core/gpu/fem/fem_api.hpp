/* SPDX-License-Identifier: MPL-2.0 */
#pragma once
// Host-side boundary of the FEM device route. It carries mesh arrays, reference
// tabulations, right-hand sides and results only; no matrix crosses this API.
#include <cstdint>
#include <memory>
#include <string>
#include <utility>
#include <vector>

namespace radia_fem {

// How a reference quantity maps to physical space.
enum Law { LAW_NONE = 0, LAW_GRAD = 1, LAW_VALUE = 2, LAW_CURL = 3, LAW_COVARIANT = 4 };

// Per-element coefficient columns, six per element.
enum Column { STIFF_RE = 0, STIFF_IM = 1, STIFF_AUX = 2, MASS_RE = 3, MASS_IM = 4, MASS_AUX = 5, COLUMNS = 6 };

struct FormInput {
    int law = LAW_NONE;
    int nq = 0;
    bool affine = false;
    const double* weights = nullptr;   // nq
    const double* table = nullptr;     // ncls * nq * nl * 3
    const double* dgeo = nullptr;      // nq * ngeo * dim, with nodes
    const double* jacobian = nullptr;  // ne * nq * 3 * dim, on curved elements
};

struct GroupInput {
    int dim = 3, ne = 0, nl = 0, ncls = 0, ngeo = 0;
    const int* dofs = nullptr;      // ne * nl
    const int* cls = nullptr;       // ne
    const double* nodes = nullptr;  // ne * ngeo * 3
    const double* coef = nullptr;   // ne * COLUMNS
    FormInput stiff, mass;
};

struct LedgerRow {
    std::string payload, direction;
    uint64_t bytes = 0;
};

struct SolveReport {
    int iterations = 0;
    double initial_norm = 0, final_norm = 0, true_residual = 0;
    std::vector<double> history;
};

struct AmgReport {
    std::vector<int> rows;
    std::vector<int> nonzeros;
    int pinned = 0;
    std::string coarse;
};

struct AmsOptions {
    int cycles = 1;
    double omega = 1.;
};

// Element split of the BDDC preconditioner. Local positions are uniform across a
// group and are checked on the host against the element DOF numbers.
struct BddcGroupInput {
    int nw = 0, ni = 0;
    const int* wb_local = nullptr;
    const int* int_local = nullptr;
};

class System {
public:
    System(int ndof, bool complex_values, const std::vector<GroupInput>& groups, const int* free_dofs);
    ~System();
    System(const System&) = delete;
    System& operator=(const System&) = delete;

    // Each preconditioner is built on the device from the assembled operators.
    AmgReport setup_amg(int min_coarse, double theta);
    std::vector<AmgReport> setup_ams(int nedge, int nv, const int* edge_vertices, const double* coordinates,
                                     int min_coarse, double theta, const AmsOptions& options);
    std::vector<AmgReport> setup_bddc(const std::vector<BddcGroupInput>& split, const int* wb_of_dof, int nw,
                                      int nedge, int nv, const int* edge_vertices, const double* coordinates,
                                      int cycles, int min_coarse, double theta);

    // Action of an assembled operator on a host probe vector. Which is 0 for the
    // real part, 1 for the imaginary part and 2 for the auxiliary real operator.
    void probe_action(int which, const double* x, double* y);
    void probe_preconditioner(const double* x_re, const double* x_im, double* y_re, double* y_im);
    int nonzeros() const;

    SolveReport solve(const double* rhs_re, const double* rhs_im, const double* lift_re, const double* lift_im,
                      double tolerance, int max_iterations, double relative_limit, double* x_re, double* x_im);

    std::vector<LedgerRow> ledger() const;
    uint64_t peak_device_bytes() const;
    std::vector<std::pair<std::string, double>> timings() const;

private:
    class Impl;
    std::unique_ptr<Impl> impl_;
};

// Host-readable name of the device this process uses.
std::string device_name();

}  // namespace radia_fem
