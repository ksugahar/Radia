#include "rad_beam_dynamics_pybind.h"

#include "rad_beam_dynamics.h"

#include <pybind11/numpy.h>
#include <pybind11/stl.h>

#include <memory>
#include <stdexcept>
#include <string>

namespace py = pybind11;

namespace {

using radia::beam::Vec3;

Vec3 ReadVec3(const py::object& value, const char* name) {
    py::array_t<double, py::array::c_style | py::array::forcecast> array =
        py::array_t<double, py::array::c_style |
                            py::array::forcecast>::ensure(value);
    if (!array)
        throw std::invalid_argument(std::string(name) +
                                    " must be a real array");
    const auto buffer = array.request();
    if (buffer.ndim != 1 || buffer.shape[0] != 3)
        throw std::invalid_argument(std::string(name) +
                                    " must have shape (3,)");
    const double* data = static_cast<const double*>(buffer.ptr);
    return {data[0], data[1], data[2]};
}

py::array_t<double> Vec3Array(const Vec3& value) {
    py::array_t<double> result(3);
    result.mutable_data()[0] = value.x;
    result.mutable_data()[1] = value.y;
    result.mutable_data()[2] = value.z;
    return result;
}

radia::beam::IndependentVariable ReadIndependentVariable(
        const std::string& value) {
    if (value == "time") return radia::beam::IndependentVariable::time;
    if (value == "path_length")
        return radia::beam::IndependentVariable::path_length;
    if (value == "azimuth")
        return radia::beam::IndependentVariable::azimuth;
    throw std::invalid_argument(
        "independent must be 'time', 'path_length', or 'azimuth'");
}

std::string IndependentVariableName(
        radia::beam::IndependentVariable value) {
    if (value == radia::beam::IndependentVariable::time) return "time";
    if (value == radia::beam::IndependentVariable::path_length)
        return "path_length";
    return "azimuth";
}

std::string DomainStatusName(radia::beam::DomainStatus value) {
    switch (value) {
    case radia::beam::DomainStatus::inside:
        return "inside";
    case radia::beam::DomainStatus::outside:
        return "outside";
    case radia::beam::DomainStatus::boundary:
        return "boundary";
    case radia::beam::DomainStatus::invalid:
        return "invalid";
    }
    return "invalid";
}

}  // namespace

void ExportBeamDynamics(py::module_& module) {
    using namespace radia::beam;

    py::class_<ParticleSpecies>(module, "BeamParticleSpecies")
        .def(py::init([](double charge_c, double rest_mass_kg,
                         std::string name) {
                 ParticleSpecies value{charge_c, rest_mass_kg,
                                       std::move(name)};
                 // The reference factory is also the canonical species
                 // validation boundary.
                 (void)ReferenceParticle::FromKineticEnergyEV(value, 0.0);
                 return value;
             }),
             py::arg("charge_c"), py::arg("rest_mass_kg"),
             py::arg("name") = "custom")
        .def_static("proton", &ParticleSpecies::Proton)
        .def_static("electron", &ParticleSpecies::Electron)
        .def_readonly("charge_c", &ParticleSpecies::charge_c)
        .def_readonly("rest_mass_kg", &ParticleSpecies::rest_mass_kg)
        .def_readonly("name", &ParticleSpecies::name);

    py::class_<ReferenceParticle>(module, "BeamReferenceParticle")
        .def_static("from_kinetic_energy_ev",
                    &ReferenceParticle::FromKineticEnergyEV,
                    py::arg("species"), py::arg("kinetic_energy_ev"))
        .def_readonly("species", &ReferenceParticle::species)
        .def_readonly("kinetic_energy_j",
                      &ReferenceParticle::kinetic_energy_j)
        .def_readonly("momentum_kg_m_s",
                      &ReferenceParticle::momentum_kg_m_s)
        .def_readonly("magnetic_rigidity_t_m",
                      &ReferenceParticle::magnetic_rigidity_t_m);

    py::class_<CartesianState>(module, "BeamCartesianState")
        .def(py::init([](py::object position_m,
                         py::object kinetic_momentum_kg_m_s,
                         double time_s, double path_length_m) {
                 CartesianState value;
                 value.position_m = ReadVec3(position_m, "position_m");
                 value.kinetic_momentum_kg_m_s = ReadVec3(
                     kinetic_momentum_kg_m_s,
                     "kinetic_momentum_kg_m_s");
                 value.time_s = time_s;
                 value.path_length_m = path_length_m;
                 return value;
             }),
             py::arg("position_m"),
             py::arg("kinetic_momentum_kg_m_s"), py::arg("time_s") = 0.0,
             py::arg("path_length_m") = 0.0)
        .def_property(
            "position_m",
            [](const CartesianState& value) {
                return Vec3Array(value.position_m);
            },
            [](CartesianState& value, py::object input) {
                value.position_m = ReadVec3(input, "position_m");
            })
        .def_property(
            "kinetic_momentum_kg_m_s",
            [](const CartesianState& value) {
                return Vec3Array(value.kinetic_momentum_kg_m_s);
            },
            [](CartesianState& value, py::object input) {
                value.kinetic_momentum_kg_m_s = ReadVec3(
                    input, "kinetic_momentum_kg_m_s");
            })
        .def_readwrite("time_s", &CartesianState::time_s)
        .def_readwrite("path_length_m", &CartesianState::path_length_m);

    py::class_<FieldSample>(module, "BeamFieldSample")
        .def_property_readonly("electric_v_m", [](const FieldSample& value) {
            return Vec3Array(value.electric_v_m);
        })
        .def_property_readonly("magnetic_t", [](const FieldSample& value) {
            return Vec3Array(value.magnetic_t);
        })
        .def_property_readonly("domain_status", [](const FieldSample& value) {
            return DomainStatusName(value.domain_status);
        });

    py::class_<Field, std::shared_ptr<Field>>(module, "BeamField")
        .def("evaluate",
             [](const Field& field, py::object position_m, double time_s,
                bool electric, bool magnetic) {
                 return field.Evaluate(ReadVec3(position_m, "position_m"),
                                       time_s, {electric, magnetic});
             },
             py::arg("position_m"), py::arg("time_s") = 0.0,
             py::arg("electric") = true, py::arg("magnetic") = true)
        .def_property_readonly("type_name", &Field::TypeName);

    py::class_<ZeroField, Field, std::shared_ptr<ZeroField>>(
        module, "BeamZeroField")
        .def(py::init<>());

    py::class_<UniformField, Field, std::shared_ptr<UniformField>>(
        module, "BeamUniformField")
        .def(py::init([](py::object magnetic_t, py::object electric_v_m) {
                 return std::make_shared<UniformField>(
                     ReadVec3(magnetic_t, "magnetic_t"),
                     ReadVec3(electric_v_m, "electric_v_m"));
             }),
             py::arg("magnetic_t"),
             py::arg("electric_v_m") = py::make_tuple(0.0, 0.0, 0.0))
        .def_property_readonly("magnetic_t", [](const UniformField& value) {
            return Vec3Array(value.MagneticFieldT());
        })
        .def_property_readonly("electric_v_m", [](const UniformField& value) {
            return Vec3Array(value.ElectricFieldVM());
        });

    py::class_<StateDerivative>(module, "BeamStateDerivative")
        .def_property_readonly("dposition_m",
             [](const StateDerivative& value) {
                 return Vec3Array(value.dposition_m);
             })
        .def_property_readonly("dkinetic_momentum_kg_m_s",
             [](const StateDerivative& value) {
                 return Vec3Array(value.dkinetic_momentum_kg_m_s);
             })
        .def_readonly("dtime_s", &StateDerivative::dtime_s)
        .def_readonly("dpath_length_m", &StateDerivative::dpath_length_m)
        .def_readonly("field", &StateDerivative::field);

    py::class_<InvariantReport>(module, "BeamInvariantReport")
        .def_readonly("momentum_kg_m_s",
                      &InvariantReport::momentum_kg_m_s)
        .def_readonly("relativistic_gamma",
                      &InvariantReport::relativistic_gamma)
        .def_readonly("kinetic_energy_j",
                      &InvariantReport::kinetic_energy_j)
        .def_readonly("speed_m_s", &InvariantReport::speed_m_s)
        .def_property_readonly("domain_status",
             [](const InvariantReport& value) {
                 return DomainStatusName(value.domain_status);
             });

    py::class_<Equation, std::shared_ptr<Equation>>(module, "BeamEquation")
        .def("rhs", &Equation::RHS, py::arg("independent_value"),
             py::arg("state"))
        .def("invariants", &Equation::Invariants, py::arg("state"))
        .def_property_readonly("independent",
             [](const Equation& value) {
                 return IndependentVariableName(value.Variable());
             });

    py::class_<LorentzEquation, Equation,
               std::shared_ptr<LorentzEquation>>(
        module, "BeamLorentzEquation")
        .def(py::init([](ParticleSpecies species,
                         std::shared_ptr<Field> field,
                         const std::string& independent) {
                 return std::make_shared<LorentzEquation>(
                     std::move(species), std::move(field),
                     ReadIndependentVariable(independent));
             }),
             py::arg("species"), py::arg("field"),
             py::arg("independent") = "time")
        .def_property_readonly("species", &LorentzEquation::Species);

    py::class_<StepResult>(module, "BeamStepResult")
        .def_readonly("independent_before",
                      &StepResult::independent_before)
        .def_readonly("independent_after", &StepResult::independent_after)
        .def_readonly("accepted_step", &StepResult::accepted_step)
        .def_readonly("state_before", &StepResult::state_before)
        .def_readonly("state_after", &StepResult::state_after)
        .def_readonly("rhs_before", &StepResult::rhs_before)
        .def_readonly("invariants_before", &StepResult::invariants_before)
        .def_readonly("invariants_after", &StepResult::invariants_after);

    py::class_<Stepper, std::shared_ptr<Stepper>>(module, "BeamStepper")
        .def("step", &Stepper::Step, py::arg("equation"),
             py::arg("independent_value"), py::arg("state"),
             py::arg("step"))
        .def_property_readonly("type_name", &Stepper::TypeName);

    py::class_<ClassicalRK4, Stepper, std::shared_ptr<ClassicalRK4>>(
        module, "BeamClassicalRK4")
        .def(py::init<>());
    py::class_<Boris2, Stepper, std::shared_ptr<Boris2>>(
        module, "BeamBoris2")
        .def(py::init<>());

    py::class_<TrackPlan>(module, "BeamTrackPlan")
        .def(py::init<>())
        .def_readwrite("start", &TrackPlan::start)
        .def_readwrite("stop", &TrackPlan::stop)
        .def_readwrite("maximum_step", &TrackPlan::maximum_step)
        .def_readwrite("maximum_steps", &TrackPlan::maximum_steps);

    py::class_<StepRecord>(module, "BeamStepRecord")
        .def_readonly("independent_value",
                      &StepRecord::independent_value)
        .def_readonly("attempted_step", &StepRecord::attempted_step)
        .def_readonly("accepted_step", &StepRecord::accepted_step)
        .def_readonly("accepted", &StepRecord::accepted)
        .def_readonly("state_before", &StepRecord::state_before)
        .def_readonly("state_after", &StepRecord::state_after)
        .def_readonly("rhs_before", &StepRecord::rhs_before)
        .def_readonly("invariants_before", &StepRecord::invariants_before)
        .def_readonly("invariants_after", &StepRecord::invariants_after);

    py::class_<TrajectorySummary>(module, "BeamTrajectorySummary")
        .def_readonly("accepted_steps", &TrajectorySummary::accepted_steps)
        .def_readonly("independent_start",
                      &TrajectorySummary::independent_start)
        .def_readonly("independent_stop",
                      &TrajectorySummary::independent_stop)
        .def_readonly("path_length_change_m",
                      &TrajectorySummary::path_length_change_m)
        .def_readonly("momentum_conservation_applicable",
                      &TrajectorySummary::momentum_conservation_applicable)
        .def_readonly("maximum_relative_momentum_error",
                      &TrajectorySummary::maximum_relative_momentum_error);

    py::class_<Trajectory>(module, "BeamTrajectory")
        .def_property_readonly("samples", &Trajectory::Samples)
        .def_property_readonly("steps", &Trajectory::Steps)
        .def_property_readonly("summary", &Trajectory::Summary);

    py::class_<Tracker>(module, "BeamTracker")
        .def(py::init([](std::shared_ptr<Equation> equation,
                         std::shared_ptr<Stepper> stepper) {
                 return std::make_unique<Tracker>(std::move(equation),
                                                  std::move(stepper));
             }),
             py::arg("equation"), py::arg("stepper"))
        .def("step", &Tracker::Step, py::arg("independent_value"),
             py::arg("state"), py::arg("step"))
        .def("track",
             [](const Tracker& tracker, const CartesianState& initial_state,
                const TrackPlan& plan) {
                 py::gil_scoped_release release;
                 return tracker.Track(initial_state, plan);
             },
             py::arg("initial_state"), py::arg("plan"));
}
