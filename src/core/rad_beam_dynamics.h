#pragma once

#include <cstddef>
#include <memory>
#include <string>
#include <vector>

namespace radia::beam {

struct Vec3 {
    double x = 0.0;
    double y = 0.0;
    double z = 0.0;
};

double Dot(const Vec3& left, const Vec3& right);
Vec3 Cross(const Vec3& left, const Vec3& right);
double Norm(const Vec3& value);

struct ParticleSpecies {
    double charge_c = 0.0;
    double rest_mass_kg = 0.0;
    std::string name;

    static ParticleSpecies Proton();
    static ParticleSpecies Electron();
};

struct ReferenceParticle {
    ParticleSpecies species;
    double kinetic_energy_j = 0.0;
    double momentum_kg_m_s = 0.0;
    double magnetic_rigidity_t_m = 0.0;

    static ReferenceParticle FromKineticEnergyEV(
        ParticleSpecies species, double kinetic_energy_ev);
};

struct CartesianState {
    Vec3 position_m;
    Vec3 kinetic_momentum_kg_m_s;
    double time_s = 0.0;
    double path_length_m = 0.0;
};

enum class DomainStatus { inside, outside, boundary, invalid };

struct FieldRequest {
    bool electric = true;
    bool magnetic = true;
};

struct FieldSample {
    Vec3 electric_v_m;
    Vec3 magnetic_t;
    DomainStatus domain_status = DomainStatus::inside;
};

class Field {
public:
    virtual ~Field() = default;
    virtual FieldSample Evaluate(
        const Vec3& position_m, double time_s,
        const FieldRequest& request = {}) const = 0;
    virtual std::string TypeName() const = 0;
};

class ZeroField final : public Field {
public:
    FieldSample Evaluate(const Vec3& position_m, double time_s,
                         const FieldRequest& request = {}) const override;
    std::string TypeName() const override;
};

class UniformField final : public Field {
public:
    UniformField(Vec3 magnetic_t, Vec3 electric_v_m = {});

    FieldSample Evaluate(const Vec3& position_m, double time_s,
                         const FieldRequest& request = {}) const override;
    std::string TypeName() const override;
    const Vec3& MagneticFieldT() const;
    const Vec3& ElectricFieldVM() const;

private:
    Vec3 magnetic_t_;
    Vec3 electric_v_m_;
};

enum class IndependentVariable { time, path_length, azimuth };

struct StateDerivative {
    Vec3 dposition_m;
    Vec3 dkinetic_momentum_kg_m_s;
    double dtime_s = 0.0;
    double dpath_length_m = 0.0;
    FieldSample field;
};

struct InvariantReport {
    double momentum_kg_m_s = 0.0;
    double relativistic_gamma = 1.0;
    double kinetic_energy_j = 0.0;
    double speed_m_s = 0.0;
    DomainStatus domain_status = DomainStatus::inside;
};

class Equation {
public:
    virtual ~Equation() = default;
    virtual StateDerivative RHS(double independent_value,
                                const CartesianState& state) const = 0;
    virtual InvariantReport Invariants(const CartesianState& state) const = 0;
    virtual IndependentVariable Variable() const = 0;
};

class LorentzEquation final : public Equation {
public:
    LorentzEquation(ParticleSpecies species,
                    std::shared_ptr<const Field> field,
                    IndependentVariable independent_variable =
                        IndependentVariable::time);

    StateDerivative RHS(double independent_value,
                        const CartesianState& state) const override;
    InvariantReport Invariants(const CartesianState& state) const override;
    IndependentVariable Variable() const override;
    const ParticleSpecies& Species() const;
    const std::shared_ptr<const Field>& FieldObject() const;

private:
    ParticleSpecies species_;
    std::shared_ptr<const Field> field_;
    IndependentVariable independent_variable_;
};

struct StepResult {
    double independent_before = 0.0;
    double independent_after = 0.0;
    double accepted_step = 0.0;
    CartesianState state_before;
    CartesianState state_after;
    StateDerivative rhs_before;
    InvariantReport invariants_before;
    InvariantReport invariants_after;
};

class Stepper {
public:
    virtual ~Stepper() = default;
    virtual StepResult Step(const Equation& equation,
                            double independent_value,
                            const CartesianState& state,
                            double step) const = 0;
    virtual std::string TypeName() const = 0;
};

class ClassicalRK4 final : public Stepper {
public:
    StepResult Step(const Equation& equation, double independent_value,
                    const CartesianState& state,
                    double step) const override;
    std::string TypeName() const override;
};

class Boris2 final : public Stepper {
public:
    StepResult Step(const Equation& equation, double independent_value,
                    const CartesianState& state,
                    double step) const override;
    std::string TypeName() const override;
};

struct TrackPlan {
    double start = 0.0;
    double stop = 0.0;
    double maximum_step = 0.0;
    std::size_t maximum_steps = 1000000;
};

struct StepRecord {
    double independent_value = 0.0;
    double attempted_step = 0.0;
    double accepted_step = 0.0;
    bool accepted = false;
    CartesianState state_before;
    CartesianState state_after;
    StateDerivative rhs_before;
    InvariantReport invariants_before;
    InvariantReport invariants_after;
};

struct TrajectorySummary {
    std::size_t accepted_steps = 0;
    double independent_start = 0.0;
    double independent_stop = 0.0;
    double path_length_change_m = 0.0;
    bool momentum_conservation_applicable = true;
    double maximum_relative_momentum_error = 0.0;
};

class Trajectory {
public:
    const std::vector<CartesianState>& Samples() const;
    const std::vector<StepRecord>& Steps() const;
    const TrajectorySummary& Summary() const;

private:
    friend class Tracker;
    std::vector<CartesianState> samples_;
    std::vector<StepRecord> steps_;
    TrajectorySummary summary_;
};

class Tracker {
public:
    Tracker(std::shared_ptr<const Equation> equation,
            std::shared_ptr<const Stepper> stepper);

    StepResult Step(double independent_value, const CartesianState& state,
                    double step) const;
    Trajectory Track(const CartesianState& initial_state,
                     const TrackPlan& plan) const;
    const std::shared_ptr<const Equation>& EquationObject() const;
    const std::shared_ptr<const Stepper>& StepperObject() const;

private:
    std::shared_ptr<const Equation> equation_;
    std::shared_ptr<const Stepper> stepper_;
};

}  // namespace radia::beam
