#ifndef RAD_CUBIT_PROJECTION_H
#define RAD_CUBIT_PROJECTION_H

/**
 * C++ surface projection using Cubit/ACIS geometry kernel.
 *
 * Provides fast C++ callbacks for CallbackGeometry, replacing the
 * Python callbacks that call cubit.surface(sid).u_v_from_position().
 *
 * Performance: eliminates ~1.5ms/call Python overhead.
 * For p=4 with 2790 elements: 418s (Python) -> expected <10s (C++).
 *
 * Requires: Cubit SDK (cubit_geom.lib, cubit_util.lib)
 */

#include <vector>
#include <unordered_map>
#include <tuple>

namespace radia {

/**
 * C++ surface projection cache using Cubit's RefFace API.
 *
 * Holds a mapping from Netgen surfnr (1-indexed) to Cubit surface ID,
 * and provides projection/normal functions callable from C++.
 */
class CubitProjection {
public:
    /**
     * @param surfnr_to_sid  Map from Netgen surfnr (1-indexed) to Cubit surface ID
     */
    CubitProjection(const std::unordered_map<int, int>& surfnr_to_sid);

    /**
     * Project point to surface, return (x, y, z, u, v).
     * Matches the CallbackGeometry project signature.
     */
    std::tuple<double, double, double, double, double>
    Project(int surfnr, double x, double y, double z,
            double u_hint, double v_hint) const;

    /**
     * Get surface normal at point, return (nx, ny, nz).
     * Matches the CallbackGeometry normal signature.
     */
    std::tuple<double, double, double>
    Normal(int surfnr, double x, double y, double z) const;

private:
    std::unordered_map<int, int> m_surfnr_to_sid;
};

} // namespace radia

#endif // RAD_CUBIT_PROJECTION_H
