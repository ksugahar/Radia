#pragma once

#include <vector>

namespace radia { namespace ih {

struct CSRMatrix {
    int n = 0;
    std::vector<int> row_ptr;
    std::vector<int> col;
    std::vector<double> value;
};

struct ThermalState {
    std::vector<double> temperature_K;
    double time_s = 0.0;
    double previous_angle_rad = 0.0;
};

struct ThermalStepOptions {
    double dt_s = 0.0;
    double conductivity_scale = 1.0;
    double convection_W_per_m2K = 0.0;
    double ambient_temperature_K = 293.15;
    double tolerance = 1.0e-10;
    int max_iterations = 500;
};

// Advance M*T' + K*T = f(q, Tamb) by backward Euler.  M and K are assembled
// by the checked NGSolve mesh contract; this function owns only the native
// state update and linear solve used by the Simulink block.
void advance_thermal(const CSRMatrix& mass, const CSRMatrix& stiffness,
                     const CSRMatrix* convection,
                     const std::vector<double>& source_W,
                     const std::vector<double>& cell_weights,
                     double angle_now_rad, const ThermalStepOptions& options,
                     ThermalState& state);

}}  // namespace radia::ih
