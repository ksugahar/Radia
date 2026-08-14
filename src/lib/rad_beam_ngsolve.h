#pragma once

#include "rad_beam_dynamics.h"
#include "rad_beam_transfer.h"

#include <gridfunction.hpp>

#include <array>
#include <memory>
#include <string>
#include <vector>

namespace radia::beam {

using Vector3 = std::array<double, 3>;

struct GridFunctionLinearizationOptions {
    double magnetic_rigidity_t_m = 0.0;
    double sample_radius_m = 1.0e-3;
    Vector3 initial_horizontal{1.0, 0.0, 0.0};
    double curvature_sign = 1.0;
    double gradient_sign = 1.0;
    unsigned multipole_order = 1;
    unsigned maximum_map_order = 1;
    double maximum_step_m = 1.0e-3;
    std::size_t maximum_steps = 1000000;
};

struct GridFunctionSegmentLinearization {
    Vector3 reference_position_m{};
    Vector3 horizontal_axis{};
    Vector3 vertical_axis{};
    Vector3 tangent_axis{};
    Vector3 center_field_local_t{};
    Vector3 fitted_center_field_local_t{};
    // Row-major [field component Bx, By, Bs][derivative x, y].
    std::array<double, 6> field_gradient_local_t_per_m{};
    // Row-major [sample][global xyz] and [sample][local Bx, By, Bs].
    std::array<double, 27> sample_positions_m{};
    std::array<double, 27> sample_fields_local_t{};
    double curvature_per_m = 0.0;
    double normal_gradient_per_m2 = 0.0;
    double skew_gradient_per_m2 = 0.0;
    double transverse_divergence_t_per_m = 0.0;
    double transverse_curl_mismatch_t_per_m = 0.0;
    double center_fit_bias_t = 0.0;
    double rms_fit_residual_t = 0.0;
    double maximum_fit_residual_t = 0.0;
    std::size_t fit_rank = 3;
    double scaled_design_condition = 1.5;
    TransverseMagneticMultipoleExpansion multipoles;
    double multipole_rms_fit_residual_t = 0.0;
    double multipole_maximum_fit_residual_t = 0.0;
    std::size_t multipole_fit_rank = 2;
    double multipole_scaled_design_condition = 1.0;
    DynamicsJet6 dynamics_jet;
    Matrix6 a_per_m;
};

struct GridFunctionTransferReport6 {
    double magnetic_rigidity_t_m = 0.0;
    double sample_radius_m = 0.0;
    unsigned multipole_order = 1;
    VariationalReport6 transfer;
    std::vector<GridFunctionSegmentLinearization> linearizations;
};

// Evaluate a real three-component NGSolve GridFunction directly at nine
// transverse points per reference station. The local least-squares field jet
// is converted to the same combined-function linear generator used by the
// Radia accelerator-design APIs. NGSolve remains responsible for point search,
// element transformations, orientation, and GridFunction evaluation.
GridFunctionTransferReport6 PropagateGridFunctionLinearMap(
    const std::shared_ptr<ngcomp::GridFunction>& field,
    const std::vector<double>& segment_lengths_m,
    const std::vector<Vector3>& reference_positions_m,
    const std::vector<Vector3>& reference_tangents,
    const std::vector<std::string>& names,
    const GridFunctionLinearizationOptions& options);

GridFunctionTransferReport6 PropagateGridFunctionMultipoleMap(
    const std::shared_ptr<ngcomp::GridFunction>& field,
    const std::vector<double>& segment_lengths_m,
    const std::vector<Vector3>& reference_positions_m,
    const std::vector<Vector3>& reference_tangents,
    const std::vector<std::string>& names,
    const GridFunctionLinearizationOptions& options);

// Direct point-evaluation field for validating an expansion-based map against
// the solved NGSolve field. NGSolve retains point search and mapped evaluation.
class NGSolveGridFunctionField final : public Field {
public:
    explicit NGSolveGridFunctionField(
        std::shared_ptr<ngcomp::GridFunction> field);

    FieldSample Evaluate(const Vec3& position_m, double time_s,
                         const FieldRequest& request = {}) const override;
    std::string TypeName() const override;
    const std::shared_ptr<ngcomp::GridFunction>& GridFunction() const;

private:
    std::shared_ptr<ngcomp::GridFunction> field_;
};

}  // namespace radia::beam
