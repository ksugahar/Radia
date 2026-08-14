#include "rad_beam_transfer.h"

#include <algorithm>
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
using radia::beam::HamiltonianJet6;
using radia::beam::Matrix6;
using radia::beam::PropagateVariationalMap;
using radia::beam::TransverseMagneticMultipoleExpansion;
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

void TestMultipoleExpansionBuildsChromaticDynamicsJet() {
    TransverseMagneticMultipoleExpansion expansion;
    expansion.order = 3;
    expansion.normal_t_per_m_power = {0.0, 2.4, 5.0, -7.0};
    expansion.skew_t_per_m_power = {0.0, -0.6, 1.5, 2.0};
    const double rigidity = 3.0;
    const auto jet = radia::beam::BuildParaxialMagneticDynamicsJet(
        expansion, rigidity, 1.0, 1.0, 3);

    RequireClose(jet.a_per_m(1, 0), -0.8, 1.0e-15,
                 "multipole normal quadrupole");
    RequireClose(jet.a_per_m(1, 2), -0.2, 1.0e-15,
                 "multipole skew quadrupole");
    RequireClose(jet.a_per_m(3, 0), -0.2, 1.0e-15,
                 "multipole skew vertical force");
    RequireClose(jet.a_per_m(3, 2), 0.8, 1.0e-15,
                 "multipole normal vertical force");
    RequireClose(jet.f2_per_m(1, 0, 0), -10.0 / 3.0, 1.0e-15,
                 "normal sextupole xx");
    RequireClose(jet.f2_per_m(1, 0, 2), 1.0, 1.0e-15,
                 "skew sextupole xy");
    RequireClose(jet.f2_per_m(1, 0, 5), 0.8, 1.0e-15,
                 "quadrupole chromatic term");
    RequireClose(jet.f2_per_m(0, 1, 5), -1.0, 1.0e-15,
                 "horizontal drift chromatic term");
    RequireClose(jet.f3_per_m(1, 0, 0, 0), 14.0, 1.0e-14,
                 "normal octupole xxx");
    RequireClose(jet.f3_per_m(1, 0, 0, 5), 10.0 / 3.0, 1.0e-14,
                 "sextupole chromatic term");
    RequireClose(jet.f3_per_m(1, 0, 5, 5), -1.6, 1.0e-14,
                 "quadrupole second chromatic term");
    RequireClose(jet.f3_per_m(0, 1, 5, 5), 2.0, 1.0e-15,
                 "drift second chromatic term");
}

double Poisson(std::size_t row, std::size_t column) {
    if (row == 0 && column == 1) return 1.0;
    if (row == 1 && column == 0) return -1.0;
    if (row == 2 && column == 3) return 1.0;
    if (row == 3 && column == 2) return -1.0;
    if (row == 4 && column == 5) return -1.0;
    if (row == 5 && column == 4) return 1.0;
    return 0.0;
}

void TestCanonicalHamiltonianJet() {
    TransverseMagneticMultipoleExpansion expansion;
    expansion.order = 3;
    expansion.normal_t_per_m_power = {0.2, 2.4, 5.0, -7.0};
    expansion.skew_t_per_m_power = {0.0, -0.6, 1.5, 2.0};
    const double rigidity = 3.0;
    const double beta = 0.8;
    const double curvature = 0.2 / rigidity;
    const HamiltonianJet6 jet =
        radia::beam::BuildCanonicalBodyHamiltonianJet(
            expansion, rigidity, 1.0, 1.0, beta);

    RequireClose(jet.h2_per_m(1, 1), 1.0, 1.0e-15,
                 "canonical H2 px px");
    RequireClose(jet.h2_per_m(0, 0),
                 curvature * curvature + 0.8, 1.0e-15,
                 "canonical H2 x x");
    RequireClose(jet.h2_per_m(0, 5), -curvature, 1.0e-15,
                 "canonical H2 x delta");
    RequireClose(jet.h2_per_m(5, 5), 1.0 - beta * beta, 1.0e-15,
                 "canonical H2 delta delta");
    RequireClose(jet.h3_per_m(0, 0, 0), 10.0 / 3.0, 1.0e-14,
                 "canonical H3 sextupole xxx");
    RequireClose(jet.h4_per_m(0, 0, 0, 0), -14.0, 1.0e-14,
                 "canonical H4 octupole xxxx");
    RequireClose(jet.dynamics.a_per_m(1, 5), curvature, 1.0e-15,
                 "canonical A horizontal dispersion");
    RequireClose(jet.dynamics.a_per_m(4, 0), curvature, 1.0e-15,
                 "canonical A path coupling");
    RequireClose(jet.dynamics.f2_per_m(1, 0, 0), -10.0 / 3.0,
                 1.0e-14, "canonical F2 sextupole xxx");
    RequireClose(jet.dynamics.f3_per_m(1, 0, 0, 0), 14.0,
                 1.0e-14, "canonical F3 octupole xxxx");

    // A^T J + J A = 0 and -J*Fn recovers the symmetric Hamiltonian.
    double linear_defect = 0.0;
    double cubic_defect = 0.0;
    double quartic_defect = 0.0;
    for (std::size_t i = 0; i < 6; ++i)
        for (std::size_t j = 0; j < 6; ++j) {
            double residual = 0.0;
            for (std::size_t a = 0; a < 6; ++a)
                residual += jet.dynamics.a_per_m(a, i) * Poisson(a, j) +
                            Poisson(i, a) * jet.dynamics.a_per_m(a, j);
            linear_defect = std::max(linear_defect, std::abs(residual));
            for (std::size_t k = 0; k < 6; ++k) {
                double recovered_h3 = 0.0;
                for (std::size_t a = 0; a < 6; ++a)
                    recovered_h3 -=
                        Poisson(i, a) * jet.dynamics.f2_per_m(a, j, k);
                cubic_defect = std::max(
                    cubic_defect,
                    std::abs(recovered_h3 - jet.h3_per_m(i, j, k)));
                for (std::size_t l = 0; l < 6; ++l) {
                    double recovered_h4 = 0.0;
                    for (std::size_t a = 0; a < 6; ++a)
                        recovered_h4 -= Poisson(i, a) *
                            jet.dynamics.f3_per_m(a, j, k, l);
                    quartic_defect = std::max(
                        quartic_defect,
                        std::abs(recovered_h4 -
                                 jet.h4_per_m(i, j, k, l)));
                }
            }
        }
    RequireClose(linear_defect, 0.0, 2.0e-15,
                 "canonical linear Hamiltonian identity");
    RequireClose(cubic_defect, 0.0, 2.0e-15,
                 "canonical cubic generator identity");
    RequireClose(quartic_defect, 0.0, 2.0e-14,
                 "canonical quartic generator identity");
}

}  // namespace

int main() {
    try {
        TestDrift();
        TestNormalQuadrupole();
        TestNonlinearAttribution();
        TestSingleRegionCascadeAcrossSubsteps();
        TestRejectsNonsymmetricJet();
        TestMultipoleExpansionBuildsChromaticDynamicsJet();
        TestCanonicalHamiltonianJet();
        std::cout << "rad_beam_transfer: all tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "rad_beam_transfer: FAILED: " << error.what() << '\n';
        return 1;
    }
}
