#include "rad_beam_dynamics.h"

#include <algorithm>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <utility>

namespace radia::beam {
namespace {

constexpr double kSpeedOfLightMS = 299792458.0;
constexpr double kElementaryChargeC = 1.602176634e-19;
constexpr double kProtonMassKg = 1.67262192369e-27;
constexpr double kElectronMassKg = 9.1093837139e-31;

bool IsFinite(const Vec3& value) {
    return std::isfinite(value.x) && std::isfinite(value.y) &&
           std::isfinite(value.z);
}

void RequireFinite(const Vec3& value, const char* name) {
    if (!IsFinite(value))
        throw std::invalid_argument(std::string(name) +
                                    " must contain finite values");
}

void RequireFinite(const CartesianState& state) {
    RequireFinite(state.position_m, "state.position_m");
    RequireFinite(state.kinetic_momentum_kg_m_s,
                  "state.kinetic_momentum_kg_m_s");
    if (!std::isfinite(state.time_s) || !std::isfinite(state.path_length_m))
        throw std::invalid_argument("state time and path length must be finite");
}

void RequireSpecies(const ParticleSpecies& species) {
    if (!std::isfinite(species.charge_c))
        throw std::invalid_argument("particle charge must be finite");
    if (!std::isfinite(species.rest_mass_kg) || species.rest_mass_kg <= 0.0)
        throw std::invalid_argument(
            "particle rest mass must be finite and positive");
}

Vec3 Add(const Vec3& left, const Vec3& right) {
    return {left.x + right.x, left.y + right.y, left.z + right.z};
}

Vec3 Scale(const Vec3& value, double scale) {
    return {scale * value.x, scale * value.y, scale * value.z};
}

CartesianState AddScaled(const CartesianState& state,
                         const StateDerivative& derivative, double scale) {
    CartesianState result = state;
    result.position_m = Add(result.position_m,
                            Scale(derivative.dposition_m, scale));
    result.kinetic_momentum_kg_m_s = Add(
        result.kinetic_momentum_kg_m_s,
        Scale(derivative.dkinetic_momentum_kg_m_s, scale));
    result.time_s += scale * derivative.dtime_s;
    result.path_length_m += scale * derivative.dpath_length_m;
    return result;
}

CartesianState RK4Combine(const CartesianState& state,
                          const StateDerivative& k1,
                          const StateDerivative& k2,
                          const StateDerivative& k3,
                          const StateDerivative& k4, double step) {
    CartesianState result = state;
    const auto weighted = [&](double StateDerivative::*member) {
        return step * (k1.*member + 2.0 * k2.*member +
                       2.0 * k3.*member + k4.*member) /
               6.0;
    };
    const auto weighted_vector = [&](Vec3 StateDerivative::*member) {
        const Vec3& a = k1.*member;
        const Vec3& b = k2.*member;
        const Vec3& c = k3.*member;
        const Vec3& d = k4.*member;
        return Vec3{
            step * (a.x + 2.0 * b.x + 2.0 * c.x + d.x) / 6.0,
            step * (a.y + 2.0 * b.y + 2.0 * c.y + d.y) / 6.0,
            step * (a.z + 2.0 * b.z + 2.0 * c.z + d.z) / 6.0};
    };
    result.position_m = Add(
        result.position_m, weighted_vector(&StateDerivative::dposition_m));
    result.kinetic_momentum_kg_m_s = Add(
        result.kinetic_momentum_kg_m_s,
        weighted_vector(&StateDerivative::dkinetic_momentum_kg_m_s));
    result.time_s += weighted(&StateDerivative::dtime_s);
    result.path_length_m +=
        weighted(&StateDerivative::dpath_length_m);
    return result;
}

double RelativisticGamma(double momentum, double rest_mass) {
    const double ratio = momentum / (rest_mass * kSpeedOfLightMS);
    return std::sqrt(1.0 + ratio * ratio);
}

Vec3 Velocity(const Vec3& momentum, double rest_mass) {
    const double magnitude = Norm(momentum);
    const double gamma = RelativisticGamma(magnitude, rest_mass);
    return Scale(momentum, 1.0 / (gamma * rest_mass));
}

void RequireInside(const FieldSample& sample) {
    if (sample.domain_status != DomainStatus::inside &&
        sample.domain_status != DomainStatus::boundary)
        throw std::runtime_error("beam state is outside the field domain");
    RequireFinite(sample.electric_v_m, "electric field");
    RequireFinite(sample.magnetic_t, "magnetic field");
}

}  // namespace

double Dot(const Vec3& left, const Vec3& right) {
    return left.x * right.x + left.y * right.y + left.z * right.z;
}

Vec3 Cross(const Vec3& left, const Vec3& right) {
    return {left.y * right.z - left.z * right.y,
            left.z * right.x - left.x * right.z,
            left.x * right.y - left.y * right.x};
}

double Norm(const Vec3& value) { return std::sqrt(Dot(value, value)); }

ParticleSpecies ParticleSpecies::Proton() {
    return {kElementaryChargeC, kProtonMassKg, "proton"};
}

ParticleSpecies ParticleSpecies::Electron() {
    return {-kElementaryChargeC, kElectronMassKg, "electron"};
}

ReferenceParticle ReferenceParticle::FromKineticEnergyEV(
        ParticleSpecies species, double kinetic_energy_ev) {
    RequireSpecies(species);
    if (!std::isfinite(kinetic_energy_ev) || kinetic_energy_ev < 0.0)
        throw std::invalid_argument(
            "kinetic energy in eV must be finite and nonnegative");
    const double energy = kinetic_energy_ev * kElementaryChargeC;
    const double rest_energy = species.rest_mass_kg *
                               kSpeedOfLightMS * kSpeedOfLightMS;
    const double momentum = std::sqrt(
        energy * (energy + 2.0 * rest_energy)) / kSpeedOfLightMS;
    ReferenceParticle result;
    result.species = std::move(species);
    result.kinetic_energy_j = energy;
    result.momentum_kg_m_s = momentum;
    result.magnetic_rigidity_t_m =
        result.species.charge_c == 0.0
            ? std::numeric_limits<double>::infinity()
            : momentum / result.species.charge_c;
    return result;
}

FieldSample ZeroField::Evaluate(const Vec3& position_m, double time_s,
                                const FieldRequest&) const {
    RequireFinite(position_m, "position_m");
    if (!std::isfinite(time_s))
        throw std::invalid_argument("time_s must be finite");
    return {};
}

std::string ZeroField::TypeName() const { return "zero"; }

UniformField::UniformField(Vec3 magnetic_t, Vec3 electric_v_m)
    : magnetic_t_(magnetic_t), electric_v_m_(electric_v_m) {
    RequireFinite(magnetic_t_, "magnetic_t");
    RequireFinite(electric_v_m_, "electric_v_m");
}

FieldSample UniformField::Evaluate(const Vec3& position_m, double time_s,
                                   const FieldRequest& request) const {
    RequireFinite(position_m, "position_m");
    if (!std::isfinite(time_s))
        throw std::invalid_argument("time_s must be finite");
    FieldSample result;
    if (request.electric) result.electric_v_m = electric_v_m_;
    if (request.magnetic) result.magnetic_t = magnetic_t_;
    return result;
}

std::string UniformField::TypeName() const { return "uniform"; }
const Vec3& UniformField::MagneticFieldT() const { return magnetic_t_; }
const Vec3& UniformField::ElectricFieldVM() const { return electric_v_m_; }

LorentzEquation::LorentzEquation(
        ParticleSpecies species, std::shared_ptr<const Field> field,
        IndependentVariable independent_variable)
    : species_(std::move(species)), field_(std::move(field)),
      independent_variable_(independent_variable) {
    RequireSpecies(species_);
    if (!field_) throw std::invalid_argument("field must not be null");
}

StateDerivative LorentzEquation::RHS(
        double independent_value, const CartesianState& state) const {
    if (!std::isfinite(independent_value))
        throw std::invalid_argument("independent value must be finite");
    RequireFinite(state);
    StateDerivative result;
    result.field = field_->Evaluate(state.position_m, state.time_s);
    RequireInside(result.field);
    const Vec3 velocity = Velocity(
        state.kinetic_momentum_kg_m_s, species_.rest_mass_kg);
    const double speed = Norm(velocity);
    const Vec3 force = Scale(
        Add(result.field.electric_v_m,
            Cross(velocity, result.field.magnetic_t)),
        species_.charge_c);

    double scale = 1.0;
    if (independent_variable_ == IndependentVariable::path_length) {
        if (!(speed > 0.0))
            throw std::domain_error(
                "path-length Lorentz equation requires nonzero momentum");
        scale = 1.0 / speed;
    } else if (independent_variable_ == IndependentVariable::azimuth) {
        const double radius_squared =
            state.position_m.x * state.position_m.x +
            state.position_m.y * state.position_m.y;
        if (!(radius_squared > 0.0))
            throw std::domain_error(
                "azimuth Lorentz equation is singular on the z axis");
        const double azimuth_rate =
            (state.position_m.x * velocity.y -
             state.position_m.y * velocity.x) /
            radius_squared;
        if (std::abs(azimuth_rate) <=
            std::numeric_limits<double>::epsilon())
            throw std::domain_error(
                "azimuth Lorentz equation requires nonzero azimuth rate");
        scale = 1.0 / azimuth_rate;
    }

    result.dposition_m = Scale(velocity, scale);
    result.dkinetic_momentum_kg_m_s = Scale(force, scale);
    result.dtime_s = scale;
    result.dpath_length_m = speed * scale;
    return result;
}

InvariantReport LorentzEquation::Invariants(
        const CartesianState& state) const {
    RequireFinite(state);
    InvariantReport result;
    result.momentum_kg_m_s = Norm(state.kinetic_momentum_kg_m_s);
    result.relativistic_gamma = RelativisticGamma(
        result.momentum_kg_m_s, species_.rest_mass_kg);
    result.kinetic_energy_j =
        (result.relativistic_gamma - 1.0) * species_.rest_mass_kg *
        kSpeedOfLightMS * kSpeedOfLightMS;
    result.speed_m_s =
        result.momentum_kg_m_s /
        (result.relativistic_gamma * species_.rest_mass_kg);
    result.domain_status = field_->Evaluate(
        state.position_m, state.time_s, {}).domain_status;
    return result;
}

IndependentVariable LorentzEquation::Variable() const {
    return independent_variable_;
}
const ParticleSpecies& LorentzEquation::Species() const { return species_; }
const std::shared_ptr<const Field>& LorentzEquation::FieldObject() const {
    return field_;
}

StepResult ClassicalRK4::Step(const Equation& equation,
                              double independent_value,
                              const CartesianState& state,
                              double step) const {
    if (!std::isfinite(step) || step == 0.0)
        throw std::invalid_argument("step must be finite and nonzero");
    StepResult result;
    result.independent_before = independent_value;
    result.independent_after = independent_value + step;
    result.accepted_step = step;
    result.state_before = state;
    result.invariants_before = equation.Invariants(state);
    const StateDerivative k1 = equation.RHS(independent_value, state);
    const StateDerivative k2 = equation.RHS(
        independent_value + 0.5 * step, AddScaled(state, k1, 0.5 * step));
    const StateDerivative k3 = equation.RHS(
        independent_value + 0.5 * step, AddScaled(state, k2, 0.5 * step));
    const StateDerivative k4 = equation.RHS(
        independent_value + step, AddScaled(state, k3, step));
    result.rhs_before = k1;
    result.state_after = RK4Combine(state, k1, k2, k3, k4, step);
    RequireFinite(result.state_after);
    result.invariants_after = equation.Invariants(result.state_after);
    return result;
}

std::string ClassicalRK4::TypeName() const { return "classical-rk4"; }

StepResult Boris2::Step(const Equation& equation, double independent_value,
                        const CartesianState& state, double step) const {
    if (!std::isfinite(step) || step == 0.0)
        throw std::invalid_argument("step must be finite and nonzero");
    if (equation.Variable() != IndependentVariable::time)
        throw std::invalid_argument(
            "Boris2 currently requires time as the independent variable");
    const auto* lorentz = dynamic_cast<const LorentzEquation*>(&equation);
    if (!lorentz)
        throw std::invalid_argument("Boris2 requires a LorentzEquation");

    RequireFinite(state);
    StepResult result;
    result.independent_before = independent_value;
    result.independent_after = independent_value + step;
    result.accepted_step = step;
    result.state_before = state;
    result.invariants_before = equation.Invariants(state);
    result.rhs_before = equation.RHS(independent_value, state);

    const ParticleSpecies& species = lorentz->Species();
    const Vec3 velocity_before = Velocity(
        state.kinetic_momentum_kg_m_s, species.rest_mass_kg);
    const Vec3 midpoint = Add(state.position_m,
                              Scale(velocity_before, 0.5 * step));
    const FieldSample field = lorentz->FieldObject()->Evaluate(
        midpoint, state.time_s + 0.5 * step);
    RequireInside(field);

    const Vec3 electric_half = Scale(
        field.electric_v_m, 0.5 * species.charge_c * step);
    const Vec3 p_minus = Add(state.kinetic_momentum_kg_m_s, electric_half);
    const double gamma_minus = RelativisticGamma(
        Norm(p_minus), species.rest_mass_kg);
    const Vec3 rotation = Scale(
        field.magnetic_t,
        species.charge_c * step /
            (2.0 * gamma_minus * species.rest_mass_kg));
    const Vec3 rotation_twice = Scale(
        rotation, 2.0 / (1.0 + Dot(rotation, rotation)));
    const Vec3 p_prime = Add(p_minus, Cross(p_minus, rotation));
    const Vec3 p_plus = Add(p_minus, Cross(p_prime, rotation_twice));
    const Vec3 p_after = Add(p_plus, electric_half);
    const Vec3 velocity_after = Velocity(p_after, species.rest_mass_kg);

    result.state_after = state;
    result.state_after.position_m = Add(
        midpoint, Scale(velocity_after, 0.5 * step));
    result.state_after.kinetic_momentum_kg_m_s = p_after;
    result.state_after.time_s += step;
    result.state_after.path_length_m +=
        0.5 * step * (Norm(velocity_before) + Norm(velocity_after));
    RequireFinite(result.state_after);
    result.invariants_after = equation.Invariants(result.state_after);
    return result;
}

std::string Boris2::TypeName() const { return "boris2"; }

const std::vector<CartesianState>& Trajectory::Samples() const {
    return samples_;
}
const std::vector<StepRecord>& Trajectory::Steps() const { return steps_; }
const TrajectorySummary& Trajectory::Summary() const { return summary_; }

Tracker::Tracker(std::shared_ptr<const Equation> equation,
                 std::shared_ptr<const Stepper> stepper)
    : equation_(std::move(equation)), stepper_(std::move(stepper)) {
    if (!equation_) throw std::invalid_argument("equation must not be null");
    if (!stepper_) throw std::invalid_argument("stepper must not be null");
}

StepResult Tracker::Step(double independent_value,
                         const CartesianState& state, double step) const {
    return stepper_->Step(*equation_, independent_value, state, step);
}

Trajectory Tracker::Track(const CartesianState& initial_state,
                          const TrackPlan& plan) const {
    RequireFinite(initial_state);
    if (!std::isfinite(plan.start) || !std::isfinite(plan.stop))
        throw std::invalid_argument("track start and stop must be finite");
    if (!std::isfinite(plan.maximum_step) || plan.maximum_step <= 0.0)
        throw std::invalid_argument(
            "track maximum_step must be finite and positive");
    if (plan.maximum_steps == 0)
        throw std::invalid_argument("track maximum_steps must be positive");

    Trajectory result;
    result.samples_.push_back(initial_state);
    result.summary_.independent_start = plan.start;
    result.summary_.independent_stop = plan.start;
    if (plan.stop == plan.start) return result;

    const double direction = plan.stop > plan.start ? 1.0 : -1.0;
    const double span = std::abs(plan.stop - plan.start);
    const double terminal_tolerance =
        64.0 * std::numeric_limits<double>::epsilon() * span;
    double progress = 0.0;
    double independent = plan.start;
    CartesianState state = initial_state;
    const double initial_momentum = Norm(state.kinetic_momentum_kg_m_s);
    double maximum_error = 0.0;
    bool momentum_conservation_applicable = true;
    while (span - progress > terminal_tolerance) {
        if (result.steps_.size() >= plan.maximum_steps)
            throw std::runtime_error("beam tracker exceeded maximum_steps");
        const double remaining = span - progress;
        const bool final_step = remaining <= plan.maximum_step *
            (1.0 + 64.0 * std::numeric_limits<double>::epsilon());
        const double step = direction *
            (final_step ? remaining : plan.maximum_step);
        const StepResult accepted = Step(independent, state, step);
        StepRecord record;
        record.independent_value = independent;
        record.attempted_step = step;
        record.accepted_step = accepted.accepted_step;
        record.accepted = true;
        record.state_before = accepted.state_before;
        record.state_after = accepted.state_after;
        record.rhs_before = accepted.rhs_before;
        record.invariants_before = accepted.invariants_before;
        record.invariants_after = accepted.invariants_after;
        result.steps_.push_back(record);
        if (Norm(accepted.rhs_before.field.electric_v_m) > 0.0)
            momentum_conservation_applicable = false;
        // Reconstruct fixed-grid stations from the accepted-step count rather
        // than repeatedly adding the step. This avoids a spurious terminal
        // micro-step after long schedules while leaving the integrated state
        // and every recorded accepted step untouched.
        progress = final_step
            ? span
            : static_cast<double>(result.steps_.size()) *
                plan.maximum_step;
        independent = final_step
            ? plan.stop
            : plan.start + direction * progress;
        state = accepted.state_after;
        result.samples_.push_back(state);
        if (initial_momentum > 0.0) {
            maximum_error = std::max(
                maximum_error,
                std::abs(Norm(state.kinetic_momentum_kg_m_s) -
                         initial_momentum) /
                    initial_momentum);
        }
    }
    result.summary_.accepted_steps = result.steps_.size();
    result.summary_.independent_stop = plan.stop;
    result.summary_.path_length_change_m =
        state.path_length_m - initial_state.path_length_m;
    result.summary_.momentum_conservation_applicable =
        momentum_conservation_applicable;
    result.summary_.maximum_relative_momentum_error =
        momentum_conservation_applicable
            ? maximum_error
            : std::numeric_limits<double>::quiet_NaN();
    return result;
}

const std::shared_ptr<const Equation>& Tracker::EquationObject() const {
    return equation_;
}
const std::shared_ptr<const Stepper>& Tracker::StepperObject() const {
    return stepper_;
}

}  // namespace radia::beam
