#pragma once

#include <array>
#include <cstddef>
#include <string>
#include <vector>

namespace radia::beam {

inline constexpr std::size_t kPhaseSpaceDimension = 6;

struct Matrix6 {
    std::array<double, 36> values{};

    double& operator()(std::size_t row, std::size_t column) {
        return values[row * kPhaseSpaceDimension + column];
    }
    double operator()(std::size_t row, std::size_t column) const {
        return values[row * kPhaseSpaceDimension + column];
    }
};

// One output index followed by two symmetric input indices.
struct Tensor3Map6 {
    std::array<double, 216> values{};

    double& operator()(std::size_t output, std::size_t first,
                       std::size_t second) {
        return values[(output * kPhaseSpaceDimension + first) *
                      kPhaseSpaceDimension + second];
    }
    double operator()(std::size_t output, std::size_t first,
                      std::size_t second) const {
        return values[(output * kPhaseSpaceDimension + first) *
                      kPhaseSpaceDimension + second];
    }
};

// One output index followed by three symmetric input indices.
struct Tensor4Map6 {
    std::array<double, 1296> values{};

    double& operator()(std::size_t output, std::size_t first,
                       std::size_t second, std::size_t third) {
        return values[((output * kPhaseSpaceDimension + first) *
                       kPhaseSpaceDimension + second) *
                      kPhaseSpaceDimension + third];
    }
    double operator()(std::size_t output, std::size_t first,
                      std::size_t second, std::size_t third) const {
        return values[((output * kPhaseSpaceDimension + first) *
                       kPhaseSpaceDimension + second) *
                      kPhaseSpaceDimension + third];
    }
};

// Factorial convention:
//   u_out = R u + 1/2 T[u,u] + 1/6 U[u,u,u] + O(u^4).
struct TaylorMap6 {
    Matrix6 r;
    Tensor3Map6 t;
    Tensor4Map6 u;
};

// Local equation jet:
//   u' = A u + 1/2 F2[u,u] + 1/6 F3[u,u,u] + O(u^4).
struct DynamicsJet6 {
    Matrix6 a_per_m;
    Tensor3Map6 f2_per_m;
    Tensor4Map6 f3_per_m;
};

struct DynamicsSegment6 {
    double length_m = 0.0;
    DynamicsJet6 jet;
    std::string name;
};

struct VariationalOptions {
    unsigned maximum_order = 3;
    double maximum_step_m = 1.0e-3;
    std::size_t maximum_steps = 1000000;
    std::size_t maximum_region_pairs = 100000;
    double input_symmetry_tolerance = 1.0e-12;
};

struct TransferStation6 {
    double path_length_m = 0.0;
    std::size_t boundary_index = 0;
    TaylorMap6 map_from_start;
    Matrix6 r_to_end;
};

struct RegionNonlinearContribution6 {
    std::size_t region_index = 0;
    std::string name;
    double s_begin_m = 0.0;
    double s_end_m = 0.0;
    Tensor3Map6 t_at_end;
    Tensor4Map6 u_direct_at_end;
    Tensor4Map6 u_local_cascade_at_end;
};

struct RegionPairNonlinearContribution6 {
    std::size_t upstream_region = 0;
    std::size_t downstream_region = 0;
    Tensor4Map6 u_cascade_at_end;
    double maximum_absolute_entry = 0.0;
};

struct VariationalDiagnostics {
    std::size_t integration_steps = 0;
    double r_composition_error = 0.0;
    double t_reconstruction_error = 0.0;
    double u_reconstruction_error = 0.0;
    double t_input_symmetry_defect = 0.0;
    double u_input_symmetry_defect = 0.0;
};

struct VariationalReport6 {
    unsigned maximum_order = 0;
    TaylorMap6 endpoint_map;
    std::vector<TransferStation6> stations;
    std::vector<RegionNonlinearContribution6> regions;
    std::vector<RegionPairNonlinearContribution6> region_pairs;
    VariationalDiagnostics diagnostics;
};

Matrix6 IdentityMatrix6();
TaylorMap6 IdentityTaylorMap6();

Matrix6 Multiply(const Matrix6& left, const Matrix6& right);
TaylorMap6 ComposeTaylorMaps(const TaylorMap6& outer,
                             const TaylorMap6& inner,
                             unsigned maximum_order = 3);

TaylorMap6 IntegrateConstantJet(const DynamicsJet6& jet, double length_m,
                                unsigned maximum_order = 3);

VariationalReport6 PropagateVariationalMap(
    const std::vector<DynamicsSegment6>& segments,
    const VariationalOptions& options = {});

double MaximumAbsoluteEntry(const Matrix6& value);
double MaximumAbsoluteEntry(const Tensor3Map6& value);
double MaximumAbsoluteEntry(const Tensor4Map6& value);
double MaximumAbsoluteDifference(const Matrix6& left, const Matrix6& right);
double MaximumAbsoluteDifference(const Tensor3Map6& left,
                                 const Tensor3Map6& right);
double MaximumAbsoluteDifference(const Tensor4Map6& left,
                                 const Tensor4Map6& right);
double InputSymmetryDefect(const Tensor3Map6& value);
double InputSymmetryDefect(const Tensor4Map6& value);

}  // namespace radia::beam
