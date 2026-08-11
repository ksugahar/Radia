#include "rad_beam_transfer.h"

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

namespace {

using radia::beam::DynamicsJet6;
using radia::beam::DynamicsSegment6;
using radia::beam::Matrix6;
using radia::beam::PropagateVariationalMap;
using radia::beam::VariationalOptions;

void RequireClose(double actual, double expected, double tolerance,
                  const std::string& name) {
    if (std::abs(actual - expected) > tolerance) {
        std::ostringstream message;
        message << std::setprecision(17) << name << ": actual=" << actual
                << " expected=" << expected
                << " error=" << std::abs(actual - expected)
                << " tolerance=" << tolerance;
        throw std::runtime_error(message.str());
    }
}

void TestDrift() {
    DynamicsJet6 jet;
    jet.a_per_m(0, 1) = 1.0;
    jet.a_per_m(2, 3) = 1.0;
    jet.a_per_m(4, 5) = 1.0;
    VariationalOptions options;
    options.maximum_order = 3;
    options.maximum_step_m = 0.13;
    const auto report = PropagateVariationalMap(
        {{2.0, jet, "drift"}}, options);
    const Matrix6& r = report.endpoint_map.r;
    RequireClose(r(0, 0), 1.0, 1.0e-15, "drift R11");
    RequireClose(r(0, 1), 2.0, 1.0e-14, "drift R12");
    RequireClose(r(2, 3), 2.0, 1.0e-14, "drift R34");
    RequireClose(r(4, 5), 2.0, 1.0e-14, "drift R56");
    RequireClose(radia::beam::MaximumAbsoluteEntry(report.endpoint_map.t),
                 0.0, 0.0, "drift T");
    RequireClose(radia::beam::MaximumAbsoluteEntry(report.endpoint_map.u),
                 0.0, 0.0, "drift U");
    RequireClose(report.diagnostics.r_composition_error, 0.0, 2.0e-14,
                 "drift composition");
}

void TestNormalQuadrupole() {
    const double strength = 1.7;
    const double length = 0.8;
    const double root = std::sqrt(strength);
    DynamicsJet6 jet;
    jet.a_per_m(0, 1) = 1.0;
    jet.a_per_m(1, 0) = -strength;
    jet.a_per_m(2, 3) = 1.0;
    jet.a_per_m(3, 2) = strength;
    VariationalOptions options;
    options.maximum_order = 1;
    options.maximum_step_m = 0.001;
    const auto report = PropagateVariationalMap(
        {{length, jet, "normal_quadrupole"}}, options);
    const Matrix6& r = report.endpoint_map.r;
    const double phase = root * length;
    RequireClose(r(0, 0), std::cos(phase), 3.0e-13, "quad R11");
    RequireClose(r(0, 1), std::sin(phase) / root, 3.0e-13, "quad R12");
    RequireClose(r(1, 0), -root * std::sin(phase), 5.0e-13, "quad R21");
    RequireClose(r(1, 1), std::cos(phase), 3.0e-13, "quad R22");
    RequireClose(r(2, 2), std::cosh(phase), 5.0e-13, "quad R33");
    RequireClose(r(2, 3), std::sinh(phase) / root, 5.0e-13, "quad R34");
    RequireClose(r(3, 2), root * std::sinh(phase), 8.0e-13, "quad R43");
    RequireClose(r(3, 3), std::cosh(phase), 5.0e-13, "quad R44");
}

void TestNonlinearAttribution() {
    const double f2_upstream = 2.0;
    const double f2_downstream = 1.5;
    const double f3_direct = -0.7;
    const double first_length = 0.3;
    const double second_length = 0.4;
    const double third_length = 0.2;

    DynamicsJet6 first;
    first.f2_per_m(1, 0, 0) = f2_upstream;
    DynamicsJet6 second;
    second.f2_per_m(2, 0, 1) = f2_downstream;
    second.f2_per_m(2, 1, 0) = f2_downstream;
    DynamicsJet6 third;
    third.f3_per_m(3, 0, 0, 0) = f3_direct;

    VariationalOptions options;
    options.maximum_order = 3;
    options.maximum_step_m = 1.0;
    const auto report = PropagateVariationalMap(
        {{first_length, first, "upstream_sextupole"},
         {second_length, second, "downstream_sextupole"},
         {third_length, third, "direct_octupole"}}, options);

    const double expected_t_upstream = f2_upstream * first_length;
    const double expected_t_downstream = f2_downstream * second_length;
    const double expected_pair_u = 3.0 * f2_upstream * first_length *
                                   f2_downstream * second_length;
    const double expected_direct_u = f3_direct * third_length;
    RequireClose(report.endpoint_map.t(1, 0, 0), expected_t_upstream,
                 2.0e-14, "upstream T100");
    RequireClose(report.endpoint_map.t(2, 0, 1), expected_t_downstream,
                 2.0e-14, "downstream T201");
    RequireClose(report.endpoint_map.t(2, 1, 0), expected_t_downstream,
                 2.0e-14, "downstream T210");
    RequireClose(report.endpoint_map.u(2, 0, 0, 0), expected_pair_u,
                 3.0e-14, "ordered sextupole cascade");
    RequireClose(report.endpoint_map.u(3, 0, 0, 0), expected_direct_u,
                 2.0e-14, "direct cubic term");

    if (report.region_pairs.size() != 1)
        throw std::runtime_error("expected one ordered nonlinear region pair");
    const auto& pair = report.region_pairs.front();
    if (pair.upstream_region != 0 || pair.downstream_region != 1)
        throw std::runtime_error("nonlinear region pair has wrong ordering");
    RequireClose(pair.u_cascade_at_end(2, 0, 0, 0), expected_pair_u,
                 3.0e-14, "pair-attributed U2000");
    RequireClose(report.regions[2].u_direct_at_end(3, 0, 0, 0),
                 expected_direct_u, 2.0e-14, "region direct U3000");
    RequireClose(report.diagnostics.t_reconstruction_error, 0.0, 3.0e-14,
                 "T reconstruction");
    RequireClose(report.diagnostics.u_reconstruction_error, 0.0, 4.0e-14,
                 "U reconstruction");
    RequireClose(report.diagnostics.t_input_symmetry_defect, 0.0, 1.0e-14,
                 "T symmetry");
    RequireClose(report.diagnostics.u_input_symmetry_defect, 0.0, 1.0e-14,
                 "U symmetry");
}

void TestSingleRegionCascadeAcrossSubsteps() {
    DynamicsJet6 jet;
    jet.f2_per_m(1, 0, 0) = 2.0;
    jet.f2_per_m(2, 0, 1) = 1.5;
    jet.f2_per_m(2, 1, 0) = 1.5;
    VariationalOptions options;
    options.maximum_order = 3;
    options.maximum_step_m = 0.1;
    const auto report = PropagateVariationalMap(
        {{1.0, jet, "combined_nonlinear_region"}}, options);

    RequireClose(report.endpoint_map.u(2, 0, 0, 0), 4.5, 2.0e-13,
                 "single-region cubic cascade");
    RequireClose(report.regions[0].u_local_cascade_at_end(2, 0, 0, 0),
                 4.5, 2.0e-13, "single-region attributed cascade");
    if (!report.region_pairs.empty())
        throw std::runtime_error(
            "one physical region must not produce a self-pair attribution");
    RequireClose(report.diagnostics.u_reconstruction_error, 0.0, 3.0e-13,
                 "single-region U reconstruction");
}

void TestRejectsNonsymmetricJet() {
    DynamicsJet6 invalid;
    invalid.f2_per_m(0, 0, 1) = 1.0;
    VariationalOptions options;
    options.maximum_order = 2;
    bool rejected = false;
    try {
        (void)PropagateVariationalMap({{1.0, invalid, "invalid"}}, options);
    } catch (const std::invalid_argument&) {
        rejected = true;
    }
    if (!rejected)
        throw std::runtime_error("nonsymmetric F2 was not rejected");
}

}  // namespace

int main() {
    try {
        TestDrift();
        TestNormalQuadrupole();
        TestNonlinearAttribution();
        TestSingleRegionCascadeAcrossSubsteps();
        TestRejectsNonsymmetricJet();
        std::cout << "rad_beam_transfer: all tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "rad_beam_transfer: FAILED: " << error.what() << '\n';
        return 1;
    }
}
