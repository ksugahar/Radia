#pragma once
// Device-resident fin PEEC port solve: oriented-box L assembled on the GPU,
// energy gate, Z = R + j omega L factorised, nodal reduction and port solve,
// all on the device. Only bar geometry, R diagonal and incidence go up; only
// branch currents, node voltages and scalar diagnostics come back.
#include <complex>
#include <cstddef>
#include <cstdint>
#include <vector>

namespace radia_peec_fin {

struct Input {
    int m = 0;                       // bars
    const double* center = nullptr;  // (m, 3) m
    const double* direction = nullptr;   // (m, 3) unit
    const double* width_dir = nullptr;   // (m, 3) unit, used where has_width_dir
    const std::uint8_t* has_width_dir = nullptr;  // (m,)
    const double* length = nullptr;  // (m,)
    const double* width = nullptr;
    const double* height = nullptr;
    const double* resistance = nullptr;  // (m,) ohm
    int order = 0;                   // oriented-box order, 1..32
    const double* rule_x = nullptr;  // packed Gauss rules n = 1..order
    const double* rule_w = nullptr;
    double frequency = 0.0;
    const int* branch_nodes = nullptr;  // (m, 2) from, to
    int n_nodes = 0;
    const int* node_map = nullptr;   // (n_nodes,) reduced index, -1 = ground
    int n_reduced = 0;
    const std::complex<double>* injection = nullptr;  // (n_reduced,) A
};

struct Output {
    std::vector<std::complex<double>> current;   // (m,) branch currents, A
    std::vector<std::complex<double>> voltage;   // (n_reduced,) node voltages, V
    double current_R_current = 0.0;   // Re(I^H R I)
    double current_L_current = 0.0;   // Re(I^H L I)
    bool cholesky = false;
    double min_eig_over_max = 0.0;    // set when Cholesky fails and the spectrum admits L
    int n_null = 0;
    double max_abs_L = 0.0;
    double z_true_residual = 0.0;     // ||Z I - A^T V|| / ||A^T V|| on the device
    double kcl_true_residual = 0.0;   // ||A I - J|| / ||J|| on the device
    double t_assemble_s = 0.0, t_energy_s = 0.0, t_factor_s = 0.0, t_solve_s = 0.0;
    std::size_t bytes_to_device = 0, bytes_to_host = 0;
};

struct GateResult {
    bool cholesky = false;
    double min_eig_over_max = 0.0;
    int n_null = 0;
    double max_abs_L = 0.0;
};

Output solve_host(const Input& in);
// The device energy gate on a host-supplied symmetric matrix (test entry; uploads one matrix).
GateResult energy_gate_host(const double* L_column_major, int m);

}  // namespace radia_peec_fin
