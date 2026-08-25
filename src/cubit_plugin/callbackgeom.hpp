#ifndef FILE_CALLBACKGEOM
#define FILE_CALLBACKGEOM

/* *************************************************************************/
/* File:   callbackgeom.hpp                                                */
/* Author: K. Sugahara                                                     */
/* Date:   Mar. 2026                                                       */
/*                                                                         */
/* NetgenGeometry subclass that delegates surface projection to Python      */
/* callbacks. Enables mesh.Curve(order) with any CAD kernel (e.g. Cubit    */
/* ACIS) without requiring OCC or STEP files.                              */
/* *************************************************************************/

#include <meshing.hpp>
#include <functional>

namespace netgen
{

  class CallbackGeometry : public NetgenGeometry
  {
  public:
    // Callback types:
    //   project_point:     (surfnr, x, y, z, u_hint, v_hint) -> (x_proj, y_proj, z_proj, u, v)
    //   get_normal:        (surfnr, x, y, z) -> (nx, ny, nz)
    //   get_tangent:       (surfnr1, surfnr2, x, y, z) -> (tx, ty, tz)  [optional]
    //   project_edge:      (surfnr1, surfnr2, x, y, z) -> (x_proj, y_proj, z_proj) [optional]
    //   between_edge:      (surfnr1, surfnr2, dist1, dist2, secpoint) -> (x, y, z)
    //                      [optional; arc-length parametric midpoint on curve,
    //                       preferred over project_edge when distances are known]
    // NOTE: `has_hint` distinguishes "caller supplied an interpolated UV"
    // (has_hint=true, use closest_point_uv_guess) from "caller has no UV
    // and expects a global projection" (has_hint=false, use
    // closest_point_trimmed).  This is necessary because a valid UV hint
    // can legally be (0, 0) — e.g., a vertex sitting exactly at the
    // parameter-space origin — which the older (u_hint, v_hint) !=
    // (0, 0) discrimination incorrectly routed to the no-hint path.
    using ProjectFunc = std::function<std::tuple<double,double,double,double,double>
                                      (int surfnr, double x, double y, double z,
                                       double u_hint, double v_hint, bool has_hint)>;
    using NormalFunc = std::function<std::tuple<double,double,double>
                                    (int surfnr, double x, double y, double z)>;
    using TangentFunc = std::function<std::tuple<double,double,double>
                                     (int surfnr1, int surfnr2, double x, double y, double z)>;
    using EdgeProjectFunc = std::function<std::tuple<double,double,double>
                                          (int surfnr1, int surfnr2,
                                           double x, double y, double z)>;
    // Arc-length parametric midpoint on the shared curve.
    // Inputs are the 3D segment endpoints (p1 = edge start, p2 = edge end) so
    // that the callback can derive per-SEGMENT arc-length parameters via the
    // CAD kernel's u_from_position — this matches how Netgen invokes
    // PointBetweenEdge (p1/p2 are segment endpoints, and EdgePointGeomInfo.dist
    // refers to the full curve, not the segment).
    using BetweenEdgeFunc = std::function<std::tuple<double,double,double>
                                          (int surfnr1, int surfnr2,
                                           double x1, double y1, double z1,
                                           double x2, double y2, double z2,
                                           double secpoint)>;

  private:
    ProjectFunc project_func;
    NormalFunc normal_func;
    TangentFunc tangent_func;              // optional
    EdgeProjectFunc edge_project_func;     // optional
    BetweenEdgeFunc between_edge_func;     // optional (preferred for PointBetweenEdge)
    int num_surfaces;

  public:
    CallbackGeometry() : num_surfaces(0) {}

    CallbackGeometry(ProjectFunc _project, NormalFunc _normal, int _num_surfaces,
                     TangentFunc _tangent = nullptr,
                     EdgeProjectFunc _edge_project = nullptr,
                     BetweenEdgeFunc _between_edge = nullptr)
      : project_func(_project), normal_func(_normal),
        tangent_func(_tangent), edge_project_func(_edge_project),
        between_edge_func(_between_edge),
        num_surfaces(_num_surfaces) {}

    virtual Vec<3> GetNormal(int surfind, const Point<3> & p,
                             const PointGeomInfo* gi) const override;

    virtual PointGeomInfo ProjectPoint(int surfind, Point<3> & p) const override;

    virtual void ProjectPointEdge(int surfind, int surfind2, Point<3> & p,
                                  EdgePointGeomInfo* gi = nullptr
#ifdef RADIA_NETGEN_EDGE_DESCRIPTOR_API
                                  , int edgenr = -1
#endif
                                  ) const override;

    virtual bool ProjectPointGI(int surfind, Point<3> & p,
                                PointGeomInfo & gi) const override;

    virtual bool CalcPointGeomInfo(int surfind, PointGeomInfo& gi,
                                   const Point<3> & p3) const override;

    virtual void PointBetweenEdge(const Point<3> & p1, const Point<3> & p2,
                                  double secpoint,
                                  int surfi1, int surfi2,
                                  const EdgePointGeomInfo & ap1,
                                  const EdgePointGeomInfo & ap2,
                                  Point<3> & newp,
                                  EdgePointGeomInfo & newgi
#ifdef RADIA_NETGEN_EDGE_DESCRIPTOR_API
                                  , int edgenr
#endif
                                  ) const override;

    virtual void PointBetween(const Point<3> & p1, const Point<3> & p2,
                              double secpoint,
                              int surfi,
                              const PointGeomInfo & gi1,
                              const PointGeomInfo & gi2,
                              Point<3> & newp,
                              PointGeomInfo & newgi) const override;

    virtual Vec<3> GetTangent(const Point<3> & p, int surfi1,
                              int surfi2,
                              const EdgePointGeomInfo & egi
#ifdef RADIA_NETGEN_EDGE_DESCRIPTOR_API
                              , int edgenr = -1
#endif
                              ) const override;
  };

}

#endif // FILE_CALLBACKGEOM
