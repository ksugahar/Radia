/**
 * C++ surface projection using Cubit/ACIS geometry kernel.
 *
 * Calls RefFace::u_v_from_position() and RefFace::normal_at() directly,
 * bypassing Python callback overhead.
 */

#include "rad_cubit_projection.h"

#include <GeometryQueryTool.hpp>
#include <RefFace.hpp>
#include <CubitVector.hpp>

namespace radia {

CubitProjection::CubitProjection(
    const std::unordered_map<int, int>& surfnr_to_sid)
    : m_surfnr_to_sid(surfnr_to_sid)
{
}

std::tuple<double, double, double, double, double>
CubitProjection::Project(int surfnr, double x, double y, double z,
                         double u_hint, double v_hint) const
{
    auto it = m_surfnr_to_sid.find(surfnr);
    if (it == m_surfnr_to_sid.end())
        return {x, y, z, u_hint, v_hint};

    RefFace* face = GeometryQueryTool::instance()->get_ref_face(it->second);
    if (!face)
        return {x, y, z, u_hint, v_hint};

    CubitVector location(x, y, z);
    double u, v;
    face->u_v_from_position(location, u, v);

    CubitVector projected = face->position_from_u_v(u, v);
    return {projected.x(), projected.y(), projected.z(), u, v};
}

std::tuple<double, double, double>
CubitProjection::Normal(int surfnr, double x, double y, double z) const
{
    auto it = m_surfnr_to_sid.find(surfnr);
    if (it == m_surfnr_to_sid.end())
        return {0.0, 0.0, 1.0};

    RefFace* face = GeometryQueryTool::instance()->get_ref_face(it->second);
    if (!face)
        return {0.0, 0.0, 1.0};

    CubitVector location(x, y, z);
    CubitVector normal = face->normal_at(location);
    return {normal.x(), normal.y(), normal.z()};
}

} // namespace radia
