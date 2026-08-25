/* *************************************************************************/
/* File:   callbackgeom.cpp                                                */
/* Author: K. Sugahara                                                     */
/* Date:   Mar. 2026                                                       */
/*                                                                         */
/* CallbackGeometry implementation — embedded in .pyd so that upstream     */
/* Netgen pip packages (without this class) can be used.                   */
/* *************************************************************************/

#include "callbackgeom.hpp"

namespace netgen
{

Vec<3> CallbackGeometry::GetNormal(int surfind, const Point<3> & p,
                                    const PointGeomInfo* gi) const
{
  if (!normal_func) return Vec<3>(0, 0, 1);
  auto [nx, ny, nz] = normal_func(surfind, p[0], p[1], p[2]);
  return Vec<3>(nx, ny, nz);
}

PointGeomInfo CallbackGeometry::ProjectPoint(int surfind, Point<3> & p) const
{
  PointGeomInfo gi;
  gi.trignum = surfind;
  if (!project_func) { gi.u = 0; gi.v = 0; return gi; }

  auto [xp, yp, zp, u, v] = project_func(surfind, p[0], p[1], p[2], 0, 0, false);
  p = Point<3>(xp, yp, zp);
  gi.u = u;
  gi.v = v;
  return gi;
}

void CallbackGeometry::ProjectPointEdge(int surfind, int surfind2, Point<3> & p,
                                         EdgePointGeomInfo* gi
#ifdef RADIA_NETGEN_EDGE_DESCRIPTOR_API
                                         , int edgenr
#endif
                                         ) const
{
#ifdef RADIA_NETGEN_EDGE_DESCRIPTOR_API
  (void)edgenr;
  const int callback_surfind = surfind;
  const int callback_surfind2 = surfind2;
#else
  const int callback_surfind = surfind + 1;
  const int callback_surfind2 = surfind2 + 1;
#endif
  if (!project_func) return;

  if (edge_project_func && surfind >= 0 && surfind2 >= 0) {
    auto [xp, yp, zp] = edge_project_func(callback_surfind, callback_surfind2,
                                           p[0], p[1], p[2]);
    p = Point<3>(xp, yp, zp);
  } else {
    // Fallback: project onto first surface only
    auto [xp, yp, zp, u, v] = project_func(callback_surfind, p[0], p[1], p[2], 0, 0, false);
    p = Point<3>(xp, yp, zp);
  }
  // Get UV on surfind (query only, do NOT move the point again)
  if (gi) {
    auto [xf, yf, zf, uf, vf] = project_func(callback_surfind, p[0], p[1], p[2], 0, 0, false);
#ifdef RADIA_NETGEN_EDGE_DESCRIPTOR_API
    gi->dist = 0;
    gi->gi.trignum = surfind;
    gi->gi.u = uf;
    gi->gi.v = vf;
#else
    gi->edgenr = 0;
    gi->body = 0;
    gi->dist = 0;
    gi->u = uf;
    gi->v = vf;
#endif
  }
}

bool CallbackGeometry::ProjectPointGI(int surfind, Point<3> & p,
                                       PointGeomInfo & gi) const
{
  if (!project_func) return false;
  auto [xp, yp, zp, u, v] = project_func(surfind, p[0], p[1], p[2],
                                           gi.u, gi.v, true);
  p = Point<3>(xp, yp, zp);
  gi.trignum = surfind;
  gi.u = u;
  gi.v = v;
  return true;
}

bool CallbackGeometry::CalcPointGeomInfo(int surfind, PointGeomInfo& gi,
                                          const Point<3> & p3) const
{
  if (!project_func) return false;
  Point<3> tmp = p3;
  auto [xp, yp, zp, u, v] = project_func(surfind, tmp[0], tmp[1], tmp[2], 0, 0, false);
  gi.trignum = surfind;
  gi.u = u;
  gi.v = v;
  return true;
}

void CallbackGeometry::PointBetweenEdge(const Point<3> & p1, const Point<3> & p2,
                                         double secpoint,
                                         int surfi1, int surfi2,
                                         const EdgePointGeomInfo & ap1,
                                         const EdgePointGeomInfo & ap2,
                                         Point<3> & newp,
                                         EdgePointGeomInfo & newgi
#ifdef RADIA_NETGEN_EDGE_DESCRIPTOR_API
                                         , int edgenr
#endif
                                         ) const
{
#ifdef RADIA_NETGEN_EDGE_DESCRIPTOR_API
  (void)edgenr;
  const int callback_surfi1 = surfi1;
  const int callback_surfi2 = surfi2;
#else
  const int callback_surfi1 = surfi1 + 1;
  const int callback_surfi2 = surfi2 + 1;
#endif
  bool placed = false;

  // Preferred path: arc-length parametric midpoint on the shared curve,
  // derived by the callback from the segment's 3D endpoints. This avoids
  // closest_point_trimmed jumping to the wrong segment on wavy curves where
  // the linear midpoint sits nearer to a different part of the curve than
  // the true arc midpoint.
  if (between_edge_func && surfi1 >= 0 && surfi2 >= 0) {
    auto [xp, yp, zp] = between_edge_func(callback_surfi1, callback_surfi2,
                                             p1[0], p1[1], p1[2],
                                             p2[0], p2[1], p2[2],
                                             secpoint);
    newp = Point<3>(xp, yp, zp);
    placed = true;
  }

  if (!placed) {
    // Linear midpoint then project (legacy path)
    newp = p1 + secpoint * (p2 - p1);
    if (edge_project_func && surfi1 >= 0 && surfi2 >= 0) {
      auto [xp, yp, zp] = edge_project_func(callback_surfi1, callback_surfi2,
                                             newp[0], newp[1], newp[2]);
      newp = Point<3>(xp, yp, zp);
    } else if (project_func && surfi1 >= 0) {
      // Fallback: project onto first surface
#ifdef RADIA_NETGEN_EDGE_DESCRIPTOR_API
      const double ap1_u = ap1.gi.u, ap1_v = ap1.gi.v;
      const double ap2_u = ap2.gi.u, ap2_v = ap2.gi.v;
#else
      const double ap1_u = ap1.u, ap1_v = ap1.v;
      const double ap2_u = ap2.u, ap2_v = ap2.v;
#endif
      auto [xp, yp, zp, u, v] = project_func(callback_surfi1, newp[0], newp[1], newp[2],
                                               0.5*(ap1_u + ap2_u),
                                               0.5*(ap1_v + ap2_v), true);
      newp = Point<3>(xp, yp, zp);
    }
  }

  // Get UV on surfi1
  if (project_func && surfi1 >= 0) {
    auto [xf, yf, zf, uf, vf] = project_func(callback_surfi1, newp[0], newp[1], newp[2], 0, 0, false);
#ifdef RADIA_NETGEN_EDGE_DESCRIPTOR_API
    newgi.gi.trignum = surfi1;
    newgi.gi.u = uf;
    newgi.gi.v = vf;
#else
    newgi.u = uf;
    newgi.v = vf;
#endif
  }
  newgi.dist = ap1.dist + secpoint * (ap2.dist - ap1.dist);
}

void CallbackGeometry::PointBetween(const Point<3> & p1, const Point<3> & p2,
                                     double secpoint,
                                     int surfi,
                                     const PointGeomInfo & gi1,
                                     const PointGeomInfo & gi2,
                                     Point<3> & newp,
                                     PointGeomInfo & newgi) const
{
  // Linear interpolation, then project onto surface
  newp = p1 + secpoint * (p2 - p1);
  if (project_func && surfi >= 0) {
    double u_hint = gi1.u + secpoint * (gi2.u - gi1.u);
    double v_hint = gi1.v + secpoint * (gi2.v - gi1.v);
    auto [xp, yp, zp, u, v] = project_func(surfi, newp[0], newp[1], newp[2],
                                             u_hint, v_hint, true);
    newp = Point<3>(xp, yp, zp);
    newgi.trignum = surfi;
    newgi.u = u;
    newgi.v = v;
  }
}

Vec<3> CallbackGeometry::GetTangent(const Point<3> & p, int surfi1,
                                     int surfi2,
                                     const EdgePointGeomInfo & egi
#ifdef RADIA_NETGEN_EDGE_DESCRIPTOR_API
                                     , int edgenr
#endif
                                     ) const
{
#ifdef RADIA_NETGEN_EDGE_DESCRIPTOR_API
  (void)edgenr;
  const int callback_surfi1 = surfi1;
  const int callback_surfi2 = surfi2;
#else
  const int callback_surfi1 = surfi1 + 1;
  const int callback_surfi2 = surfi2 + 1;
#endif
  if (tangent_func) {
    auto [tx, ty, tz] = tangent_func(callback_surfi1, callback_surfi2, p[0], p[1], p[2]);
    return Vec<3>(tx, ty, tz);
  }
  // Fallback: cross product of two surface normals
  if (normal_func && surfi1 >= 0 && surfi2 >= 0) {
    auto [n1x, n1y, n1z] = normal_func(callback_surfi1, p[0], p[1], p[2]);
    auto [n2x, n2y, n2z] = normal_func(callback_surfi2, p[0], p[1], p[2]);
    Vec<3> n1(n1x, n1y, n1z), n2(n2x, n2y, n2z);
    Vec<3> t = Cross(n1, n2);
    double len = t.Length();
    if (len > 1e-12) t /= len;
    return t;
  }
  return Vec<3>(1, 0, 0);
}

}  // namespace netgen
