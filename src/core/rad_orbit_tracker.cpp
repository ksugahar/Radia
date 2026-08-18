#include "rad_orbit_tracker.h"

#include "rad_hdiv_field_evaluator.h"

#include <core/taskmanager.hpp>

#include <algorithm>
#include <cmath>
#include <sstream>
#include <stdexcept>
#include <vector>

#if defined(_WIN32) && (defined(ALPHA__DLL__) || defined(MATLAB_MEX_FILE))
#define RAD_ORBIT_C_API __declspec(dllexport)
#define RAD_ORBIT_C_CALL __cdecl
#else
#define RAD_ORBIT_C_API
#define RAD_ORBIT_C_CALL
#endif
extern "C" {
RAD_ORBIT_C_API int RAD_ORBIT_C_CALL RadFldBatchSerial(
    double* flux_density, double* field_strength, int n_points,
    double* points, int object);
}
#undef RAD_ORBIT_C_API
#undef RAD_ORBIT_C_CALL

namespace rad_orbit {
namespace {

struct CompositeField {
    const rad_hdiv::HDivFieldEvaluator* iron;
    double iron_scale;
    int radia_object;
    bool mirror_z;

    // Evaluate the composite field at n points (row-major (n,3)).
    void Evaluate(const double* points, std::size_t count,
                  double* field_out) const {
        const std::size_t total = mirror_z ? 2*count : count;
        std::vector<double> stacked(3*total, 0.0);
        std::copy(points, points + 3*count, stacked.begin());
        if (mirror_z) {
            for (std::size_t index = 0; index < count; ++index) {
                stacked[3*(count+index)] = points[3*index];
                stacked[3*(count+index)+1] = points[3*index+1];
                stacked[3*(count+index)+2] = -points[3*index+2];
            }
        }
        std::vector<double> total_field(3*total, 0.0);
        if (iron) {
            std::vector<double> demag(3*total);
            iron->Evaluate(stacked.data(), total, demag.data(),
                           rad_hdiv::HDivFieldEvaluator::Algorithm::Auto);
            for (std::size_t index = 0; index < 3*total; ++index)
                total_field[index] += iron_scale * demag[index];
        }
        if (radia_object >= 0) {
            std::vector<double> flux(3*total);
            std::vector<double> strength(3*total);
            const int status = RadFldBatchSerial(
                flux.data(), strength.data(), static_cast<int>(total),
                stacked.data(), radia_object);
            if (status != 0)
                throw std::runtime_error(
                    "orbit tracker: RadFldBatchSerial failed on the Radia "
                    "object term");
            for (std::size_t index = 0; index < 3*total; ++index)
                total_field[index] += flux[index];
        }
        if (mirror_z) {
            for (std::size_t index = 0; index < count; ++index) {
                field_out[3*index] = 0.5 * (total_field[3*index]
                    - total_field[3*(count+index)]);
                field_out[3*index+1] = 0.5 * (total_field[3*index+1]
                    - total_field[3*(count+index)+1]);
                field_out[3*index+2] = 0.5 * (total_field[3*index+2]
                    + total_field[3*(count+index)+2]);
            }
        } else {
            std::copy(total_field.begin(), total_field.end(), field_out);
        }
    }
};

struct State {
    double r[3];
    double t[3];
};

void Acceleration(const CompositeField& field, const State& state,
                  double inverse_rigidity, double acceleration[3]) {
    double b[3];
    field.Evaluate(state.r, 1, b);
    acceleration[0] = (state.t[1]*b[2] - state.t[2]*b[1]) * inverse_rigidity;
    acceleration[1] = (state.t[2]*b[0] - state.t[0]*b[2]) * inverse_rigidity;
    acceleration[2] = (state.t[0]*b[1] - state.t[1]*b[0]) * inverse_rigidity;
}

// Classical RK4 step of (r' = t, t' = t x B / (B rho)); k1 is supplied so
// each stored step boundary carries its acceleration for the Hermite
// interpolation, and the k1 of the next step is evaluated fresh.
State Rk4Step(const CompositeField& field, const State& start,
              const double k1_t[3], double inverse_rigidity, double h) {
    auto shifted = [&](const State& base, const double* dr, const double* dt,
                       double scale) {
        State result;
        for (int axis = 0; axis < 3; ++axis) {
            result.r[axis] = base.r[axis] + scale * dr[axis];
            result.t[axis] = base.t[axis] + scale * dt[axis];
        }
        return result;
    };
    // k1: (t, a(start))
    const double* k1_r = start.t;
    // k2 at start + h/2 * k1
    State mid1 = shifted(start, k1_r, k1_t, 0.5*h);
    double k2_t[3];
    Acceleration(field, mid1, inverse_rigidity, k2_t);
    const double* k2_r = mid1.t;
    // k3 at start + h/2 * k2
    State mid2 = shifted(start, k2_r, k2_t, 0.5*h);
    double k3_t[3];
    Acceleration(field, mid2, inverse_rigidity, k3_t);
    const double* k3_r = mid2.t;
    // k4 at start + h * k3
    State end_stage = shifted(start, k3_r, k3_t, h);
    double k4_t[3];
    Acceleration(field, end_stage, inverse_rigidity, k4_t);
    const double* k4_r = end_stage.t;
    State result;
    for (int axis = 0; axis < 3; ++axis) {
        result.r[axis] = start.r[axis] + h/6.0 * (k1_r[axis]
            + 2.0*k2_r[axis] + 2.0*k3_r[axis] + k4_r[axis]);
        result.t[axis] = start.t[axis] + h/6.0 * (k1_t[axis]
            + 2.0*k2_t[axis] + 2.0*k3_t[axis] + k4_t[axis]);
    }
    return result;
}

// Cubic Hermite value on [0, 1] with endpoint values and h-scaled slopes.
double Hermite(double p0, double m0, double p1, double m1, double tau) {
    const double tau2 = tau*tau;
    const double tau3 = tau2*tau;
    return (2.0*tau3 - 3.0*tau2 + 1.0)*p0 + (tau3 - 2.0*tau2 + tau)*m0
        + (-2.0*tau3 + 3.0*tau2)*p1 + (tau3 - tau2)*m1;
}

struct StepRecord {
    State state;
    double acceleration[3];
};

void InterpolateState(const StepRecord& start, const StepRecord& end,
                      double h, double tau, double position[3],
                      double tangent[3]) {
    for (int axis = 0; axis < 3; ++axis) {
        position[axis] = Hermite(
            start.state.r[axis], h*start.state.t[axis],
            end.state.r[axis], h*end.state.t[axis], tau);
        tangent[axis] = Hermite(
            start.state.t[axis], h*start.acceleration[axis],
            end.state.t[axis], h*end.acceleration[axis], tau);
    }
}

}  // namespace

OrbitTrackResult TrackReferenceOrbit3D(
    const rad_hdiv::HDivFieldEvaluator* iron,
    double iron_scale,
    int radia_object,
    int mirror_z,
    double magnetic_rigidity,
    const double entrance_point[3],
    const double entrance_direction[3],
    double exit_x_m,
    double step_m,
    double maximum_path_m,
    double planarity_tolerance_m,
    std::size_t station_count,
    double* positions_out,
    double* tangents_out,
    double* stations_out,
    double* curvature_out) {
    if (!entrance_point || !entrance_direction || !positions_out
        || !tangents_out || !stations_out || !curvature_out)
        throw std::invalid_argument("orbit tracker: null array pointer");
    if (!iron && radia_object < 0)
        throw std::invalid_argument(
            "orbit tracker: at least one field source is required");
    if (!(magnetic_rigidity != 0.0) || !(step_m > 0.0)
        || !(maximum_path_m > step_m) || !(planarity_tolerance_m > 0.0)
        || station_count < 2)
        throw std::invalid_argument(
            "orbit tracker: integration controls are invalid");

    // Solve-loop self-wrap: the iron evaluator and any inner ParallelFor
    // reuse this region instead of standing up a pool per field call.
    ngcore::RegionTaskManager task_manager;

    const CompositeField field{iron, iron_scale, radia_object,
                               mirror_z != 0};
    const double inverse_rigidity = 1.0 / magnetic_rigidity;
    const long long step_cap = static_cast<long long>(
        std::ceil(maximum_path_m / step_m)) + 1;

    std::vector<StepRecord> records;
    records.reserve(static_cast<std::size_t>(step_cap) + 1);
    State current;
    for (int axis = 0; axis < 3; ++axis) {
        current.r[axis] = entrance_point[axis];
        current.t[axis] = entrance_direction[axis];
    }
    StepRecord record{current, {0.0, 0.0, 0.0}};
    Acceleration(field, current, inverse_rigidity, record.acceleration);
    records.push_back(record);

    double crossing_tau = -1.0;
    for (long long step = 0; step < step_cap; ++step) {
        const StepRecord& start = records.back();
        State next = Rk4Step(field, start.state, start.acceleration,
                             inverse_rigidity, step_m);
        StepRecord next_record{next, {0.0, 0.0, 0.0}};
        Acceleration(field, next, inverse_rigidity, next_record.acceleration);
        records.push_back(next_record);
        if (next.r[0] >= exit_x_m) {
            // Bisect the Hermite x(tau) for the exit-plane crossing.
            const StepRecord& lower = records[records.size()-2];
            const StepRecord& upper = records.back();
            double low = 0.0, high = 1.0;
            for (int iteration = 0; iteration < 80; ++iteration) {
                const double middle = 0.5*(low + high);
                const double x_value = Hermite(
                    lower.state.r[0], step_m*lower.state.t[0],
                    upper.state.r[0], step_m*upper.state.t[0], middle);
                if (x_value < exit_x_m) low = middle; else high = middle;
            }
            crossing_tau = 0.5*(low + high);
            break;
        }
    }
    if (crossing_tau < 0.0)
        throw std::runtime_error(
            "reference particle did not cross the C-magnet exit plane");

    const std::size_t full_steps = records.size() - 2;
    const double length = (static_cast<double>(full_steps)
                           + crossing_tau) * step_m;
    OrbitTrackResult result;
    result.length_m = length;

    // Measured planarity gate over every stored boundary.
    for (const StepRecord& stored : records) {
        result.out_of_plane_m = std::max(result.out_of_plane_m,
                                         std::fabs(stored.state.r[2]));
        result.out_of_plane_slope = std::max(result.out_of_plane_slope,
                                             std::fabs(stored.state.t[2]));
    }
    if (std::max(result.out_of_plane_m,
                 result.out_of_plane_slope * length)
            > planarity_tolerance_m) {
        std::ostringstream message;
        message.precision(3);
        message << std::scientific
                << "design orbit left the bend plane: max |z| = "
                << result.out_of_plane_m << " m, max |t_z| = "
                << result.out_of_plane_slope << " (gate "
                << planarity_tolerance_m << " m); the planar frame "
                << "machinery does not apply -- a 3D (Bishop-frame) orbit "
                << "chain is required for this field";
        throw std::runtime_error(message.str());
    }

    auto sample = [&](double path, double position[3], double tangent[3]) {
        double scaled = path / step_m;
        std::size_t interval = static_cast<std::size_t>(scaled);
        if (interval > records.size() - 2) interval = records.size() - 2;
        double tau = scaled - static_cast<double>(interval);
        if (tau < 0.0) tau = 0.0;
        InterpolateState(records[interval], records[interval+1], step_m,
                         tau, position, tangent);
        const double norm = std::sqrt(tangent[0]*tangent[0]
            + tangent[1]*tangent[1] + tangent[2]*tangent[2]);
        for (int axis = 0; axis < 3; ++axis) tangent[axis] /= norm;
    };
    for (std::size_t station = 0; station < station_count; ++station) {
        const double path = length * static_cast<double>(station)
            / static_cast<double>(station_count - 1);
        stations_out[station] = path;
        sample(path, positions_out + 3*station, tangents_out + 3*station);
    }

    // Collocated midpoint curvature from ONE batched field call:
    // kappa = -(B . z) / (B rho) on the measured-planar orbit.
    const std::size_t midpoint_count = station_count - 1;
    std::vector<double> midpoints(3*midpoint_count);
    double tangent_scratch[3];
    for (std::size_t midpoint = 0; midpoint < midpoint_count; ++midpoint) {
        const double path = 0.5*(stations_out[midpoint]
                                 + stations_out[midpoint+1]);
        sample(path, midpoints.data() + 3*midpoint, tangent_scratch);
    }
    std::vector<double> midpoint_field(3*midpoint_count);
    field.Evaluate(midpoints.data(), midpoint_count, midpoint_field.data());
    for (std::size_t midpoint = 0; midpoint < midpoint_count; ++midpoint)
        curvature_out[midpoint] =
            -midpoint_field[3*midpoint+2] * inverse_rigidity;
    return result;
}

}  // namespace rad_orbit
