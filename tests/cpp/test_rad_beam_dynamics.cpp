#include "rad_beam_dynamics.h"

#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>

namespace {

using radia::beam::Boris2;
using radia::beam::CartesianState;
using radia::beam::ClassicalRK4;
using radia::beam::IndependentVariable;
using radia::beam::LorentzEquation;
using radia::beam::ParticleSpecies;
using radia::beam::TrackPlan;
using radia::beam::Tracker;
using radia::beam::UniformField;
using radia::beam::Vec3;
using radia::beam::ZeroField;

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

void TestReferenceParticle() {
    const auto reference = radia::beam::ReferenceParticle::FromKineticEnergyEV(
        ParticleSpecies::Proton(), 220.0e6);
    RequireClose(reference.kinetic_energy_j, 220.0e6 * 1.602176634e-19,
                 1.0e-24, "reference energy");
    RequireClose(reference.momentum_kg_m_s /
                     reference.species.charge_c,
                 reference.magnetic_rigidity_t_m, 1.0e-15,
                 "magnetic rigidity");
    const auto electron = radia::beam::ReferenceParticle::FromKineticEnergyEV(
        ParticleSpecies::Electron(), 1.0e6);
    if (!(electron.magnetic_rigidity_t_m < 0.0))
        throw std::runtime_error(
            "electron magnetic rigidity must retain the charge sign");
}

void TestZeroFieldDriftRK4() {
    const ParticleSpecies species = ParticleSpecies::Proton();
    auto field = std::make_shared<ZeroField>();
    auto equation = std::make_shared<LorentzEquation>(
        species, field, IndependentVariable::time);
    auto stepper = std::make_shared<ClassicalRK4>();
    Tracker tracker(equation, stepper);
    CartesianState state;
    state.kinetic_momentum_kg_m_s = {1.0e-19, -2.0e-20, 3.0e-20};
    const auto initial = equation->Invariants(state);
    const auto trajectory = tracker.Track(
        state, TrackPlan{0.0, 2.0e-9, 1.0e-10, 100});
    RequireClose(trajectory.Samples().back().position_m.x,
                 initial.speed_m_s *
                     state.kinetic_momentum_kg_m_s.x /
                     initial.momentum_kg_m_s * 2.0e-9,
                 2.0e-16, "drift x");
    RequireClose(trajectory.Samples().back().time_s, 2.0e-9, 1.0e-24,
                 "drift time");
    RequireClose(trajectory.Summary().path_length_change_m,
                 initial.speed_m_s * 2.0e-9, 2.0e-16,
                 "drift path length");
    RequireClose(trajectory.Summary().maximum_relative_momentum_error,
                 0.0, 0.0, "drift momentum");
    if (!trajectory.Summary().momentum_conservation_applicable)
        throw std::runtime_error(
            "zero-field momentum conservation must be applicable");
}

void TestPathLengthRHS() {
    const ParticleSpecies species = ParticleSpecies::Proton();
    auto equation = std::make_shared<LorentzEquation>(
        species, std::make_shared<ZeroField>(),
        IndependentVariable::path_length);
    CartesianState state;
    state.kinetic_momentum_kg_m_s = {0.0, 3.0e-19, 4.0e-19};
    const auto rhs = equation->RHS(0.0, state);
    RequireClose(rhs.dposition_m.x, 0.0, 0.0, "path dx/ds");
    RequireClose(rhs.dposition_m.y, 0.6, 1.0e-15, "path dy/ds");
    RequireClose(rhs.dposition_m.z, 0.8, 1.0e-15, "path dz/ds");
    RequireClose(rhs.dpath_length_m, 1.0, 2.0e-15, "path ds/ds");
}

void TestAzimuthRHS() {
    const ParticleSpecies species = ParticleSpecies::Proton();
    auto equation = std::make_shared<LorentzEquation>(
        species, std::make_shared<ZeroField>(),
        IndependentVariable::azimuth);
    CartesianState state;
    state.position_m = {2.0, 0.0, 0.0};
    state.kinetic_momentum_kg_m_s = {0.0, 3.0e-19, 0.0};
    const auto rhs = equation->RHS(0.0, state);
    const auto invariants = equation->Invariants(state);
    RequireClose(rhs.dposition_m.x, 0.0, 0.0, "azimuth dx/dtheta");
    RequireClose(rhs.dposition_m.y, 2.0, 2.0e-15,
                 "azimuth dy/dtheta");
    RequireClose(rhs.dtime_s, 2.0 / invariants.speed_m_s, 2.0e-22,
                 "azimuth dt/dtheta");
    RequireClose(rhs.dpath_length_m, 2.0, 2.0e-15,
                 "azimuth ds/dtheta");
}

void TestUniformMagneticCircleBoris() {
    const ParticleSpecies species = ParticleSpecies::Proton();
    const double magnetic_t = 0.7;
    const double momentum = 3.0e-19;
    auto equation = std::make_shared<LorentzEquation>(
        species,
        std::make_shared<UniformField>(Vec3{0.0, 0.0, magnetic_t}),
        IndependentVariable::time);
    auto stepper = std::make_shared<Boris2>();
    Tracker tracker(equation, stepper);
    CartesianState state;
    state.kinetic_momentum_kg_m_s = {momentum, 0.0, 0.0};
    const auto invariants = equation->Invariants(state);
    const double omega = species.charge_c * magnetic_t /
                         (invariants.relativistic_gamma *
                          species.rest_mass_kg);
    const double period = 2.0 * 3.14159265358979323846 / omega;
    const auto trajectory = tracker.Track(
        state, TrackPlan{0.0, period, period / 8000.0, 8100});
    const auto& final = trajectory.Samples().back();
    RequireClose(radia::beam::Norm(final.kinetic_momentum_kg_m_s), momentum,
                 2.0e-32, "Boris momentum magnitude");
    RequireClose(final.position_m.x, 0.0, 3.0e-6, "Boris closed x");
    RequireClose(final.position_m.y, 0.0, 3.0e-6, "Boris closed y");
    if (trajectory.Summary().maximum_relative_momentum_error > 2.0e-13)
        throw std::runtime_error("Boris momentum conservation regressed");
}

void TestElectricFieldDisablesMomentumErrorClaim() {
    const ParticleSpecies species = ParticleSpecies::Proton();
    auto equation = std::make_shared<LorentzEquation>(
        species,
        std::make_shared<UniformField>(Vec3{}, Vec3{1.0e4, 0.0, 0.0}),
        IndependentVariable::time);
    auto stepper = std::make_shared<ClassicalRK4>();
    CartesianState state;
    state.kinetic_momentum_kg_m_s = {1.0e-19, 0.0, 0.0};
    const auto trajectory = Tracker(equation, stepper).Track(
        state, TrackPlan{0.0, 1.0e-9, 1.0e-10, 20});
    if (trajectory.Summary().momentum_conservation_applicable)
        throw std::runtime_error(
            "electric tracking must not claim momentum conservation");
    if (!std::isnan(
            trajectory.Summary().maximum_relative_momentum_error))
        throw std::runtime_error(
            "non-applicable momentum error must be NaN");
}

void TestChargeSignReversal() {
    const double magnetic_t = 0.4;
    CartesianState state;
    state.kinetic_momentum_kg_m_s = {2.0e-19, 0.0, 0.0};
    auto positive = std::make_shared<LorentzEquation>(
        ParticleSpecies::Proton(),
        std::make_shared<UniformField>(Vec3{0.0, 0.0, magnetic_t}),
        IndependentVariable::time);
    ParticleSpecies negative = ParticleSpecies::Proton();
    negative.charge_c = -negative.charge_c;
    auto reversed = std::make_shared<LorentzEquation>(
        negative,
        std::make_shared<UniformField>(Vec3{0.0, 0.0, magnetic_t}),
        IndependentVariable::time);
    const auto first = positive->RHS(0.0, state);
    const auto second = reversed->RHS(0.0, state);
    RequireClose(first.dkinetic_momentum_kg_m_s.y,
                 -second.dkinetic_momentum_kg_m_s.y, 1.0e-28,
                 "charge reversal force");
}

}  // namespace

int main() {
    try {
        TestReferenceParticle();
        TestZeroFieldDriftRK4();
        TestPathLengthRHS();
        TestAzimuthRHS();
        TestUniformMagneticCircleBoris();
        TestElectricFieldDisablesMomentumErrorClaim();
        TestChargeSignReversal();
        std::cout << "rad_beam_dynamics: all tests passed\n";
        return 0;
    } catch (const std::exception& error) {
        std::cerr << "rad_beam_dynamics: FAILED: " << error.what() << '\n';
        return 1;
    }
}
