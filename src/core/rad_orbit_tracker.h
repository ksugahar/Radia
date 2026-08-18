// Native 3D reference-orbit tracker for the beam Lie-map pipeline.
//
// Integrates dr/ds = t, dt/ds = t x B / (B rho) with fixed-step classical
// RK4 through a composite magnetic field (optional HDiv iron evaluator +
// optional Radia object, optionally z-mirror symmetrized), finds the exit
// plane crossing by cubic-Hermite interpolation of the final step, samples
// the requested stations, and evaluates the collocated midpoint curvature
// in one batched field call.  Tracking is 3D (Sugahara ruling 2026-08-18):
// all three field components drive the Lorentz force and planarity is a
// MEASURED gate, never an assumption -- a rolled/gantry-class field fails
// loudly with the message that the planar frame machinery does not apply.
//
// The ABI is plain arrays plus one optional evaluator pointer so the same
// kernel serves pybind now and a standalone MEX entry later (radia_mex.cpp
// already consumes the RadFldBatch C API this tracker uses for the Radia
// object term).

#ifndef RAD_ORBIT_TRACKER_H
#define RAD_ORBIT_TRACKER_H

#include <cstddef>

namespace rad_hdiv { class HDivFieldEvaluator; }

namespace rad_orbit {

struct OrbitTrackResult {
    double length_m = 0.0;
    double out_of_plane_m = 0.0;
    double out_of_plane_slope = 0.0;
};

// iron:          optional HDiv field evaluator (nullptr when absent); its
//                RAW output carries NO 1/(4 pi) (the evaluator contract),
//                so turning the demag solution into B takes
//                iron_scale = MU0 / (4 pi).
// radia_object:  optional Radia object key (< 0 when absent), evaluated
//                through the RadFldBatchSerial C API.
// mirror_z:      1 = z-mirror symmetrize the composite field
//                (0.5 * (B(r) + diag(-1,-1,1) B(Mr))).
// entrance/exit: start point, unit direction, and the exit plane x = exit_x.
// step_m:        fixed RK4 step; maximum_path_m caps the integration span.
// planarity_tolerance_m: measured planarity gate (max |z| and |t_z|*L).
//
// positions/tangents: (station_count, 3); stations: (station_count) path
// lengths; curvature: (station_count-1) midpoint-collocated signed
// curvature -(B . bend_axis)/(B rho) with bend_axis = +z.
//
// Throws std::invalid_argument on bad inputs, std::runtime_error when the
// exit plane is not crossed or the planarity gate trips.
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
    double* curvature_out);

}  // namespace rad_orbit

#endif  // RAD_ORBIT_TRACKER_H
