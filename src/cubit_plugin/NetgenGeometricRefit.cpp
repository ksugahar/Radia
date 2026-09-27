#include "NetgenGeometricRefit.hpp"

#ifdef HAVE_NETGEN

#include <meshing/meshing.hpp>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdio>
#include <cstdlib>
#include <unordered_map>
#include <vector>

namespace ng = netgen;

namespace {

// Gauss-Legendre rule on (0,1).
void gauss_rule(int n, std::vector<double> & x, std::vector<double> & w)
{
  x.resize(n);
  w.resize(n);
  for (int i = 0; i < n; i++) {
    const double pi = 3.14159265358979323846;
    double z = std::cos(pi * (i + 0.75) / (n + 0.5));
    for (int it = 0; it < 100; it++) {
      double p0 = 1, p1 = z;
      for (int k = 2; k <= n; k++) {
        double p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k;
        p0 = p1;
        p1 = p2;
      }
      double dp = n * (z * p1 - p0) / (z * z - 1);
      double dz = p1 / dp;
      z -= dz;
      if (std::fabs(dz) < 1e-16) break;
    }
    double p0 = 1, p1 = z;
    for (int k = 2; k <= n; k++) {
      double p2 = ((2 * k - 1) * z * p1 - (k - 1) * p0) / k;
      p0 = p1;
      p1 = p2;
    }
    double dp = n * (z * p1 - p0) / (z * z - 1);
    x[i] = 0.5 * (1 - z);
    w[i] = 1.0 / ((1 - z * z) * dp * dp);
  }
}

// Solve the symmetric positive semi-definite system A x = b in place
// (Gaussian elimination with partial pivoting).  Returns false if singular.
bool solve_dense(std::vector<double> & A, std::vector<double> & b, int n)
{
  for (int c = 0; c < n; c++) {
    int piv = c;
    for (int r = c + 1; r < n; r++)
      if (std::fabs(A[r * n + c]) > std::fabs(A[piv * n + c])) piv = r;
    if (std::fabs(A[piv * n + c]) < 1e-300) return false;
    if (piv != c) {
      for (int k = 0; k < n; k++) std::swap(A[c * n + k], A[piv * n + k]);
      std::swap(b[c], b[piv]);
    }
    for (int r = c + 1; r < n; r++) {
      double f = A[r * n + c] / A[c * n + c];
      if (f == 0) continue;
      for (int k = c; k < n; k++) A[r * n + k] -= f * A[c * n + k];
      b[r] -= f * b[c];
    }
  }
  for (int r = n - 1; r >= 0; r--) {
    double s = b[r];
    for (int k = r + 1; k < n; k++) s -= A[r * n + k] * b[k];
    b[r] = s / A[r * n + r];
  }
  return true;
}

// Reference coordinates of the vertices of a linear surface element, in
// Netgen's convention (TRIG: lambda = {x, y, 1-x-y}; QUAD: bilinear).
std::array<double, 2> vertex_ref(ng::ELEMENT_TYPE type, int v)
{
  static const double trig[3][2] = {{1, 0}, {0, 1}, {0, 0}};
  static const double quad[4][2] = {{0, 0}, {1, 0}, {1, 1}, {0, 1}};
  if (type == ng::TRIG) return {trig[v][0], trig[v][1]};
  return {quad[v][0], quad[v][1]};
}

// Linear shape weights of the element vertices at a reference point, used
// to interpolate the vertices' surface parameters into a projection hint.
void vertex_weights(ng::ELEMENT_TYPE type, double a, double b, double * lam)
{
  if (type == ng::TRIG) {
    lam[0] = a;
    lam[1] = b;
    lam[2] = 1 - a - b;
  } else {
    lam[0] = (1 - a) * (1 - b);
    lam[1] = a * (1 - b);
    lam[2] = a * b;
    lam[3] = (1 - a) * b;
  }
}

// One curved entity (edge or face) whose coefficients are refitted.
struct Entity
{
  ng::SurfaceElementIndex sei;   // surface element used for evaluation
  int surfnr = -1;               // surface for normal-distance entities
  bool on_curve = false;         // distance to a geometric curve instead
  int surfnr1 = -1, surfnr2 = -1;
  std::vector<std::array<double, 2>> fit_ref, check_ref;  // reference points
  std::vector<double> fit_w;
  double scale = 1;              // characteristic length (edge length)
  // Directions in which each coefficient may change.  An edge inside a
  // surface can slide laterally without changing its distance to the
  // surface; that exact null direction is excluded, leaving the chord
  // direction and the surface normal.  Empty means all three axes.
  std::vector<ng::Vec<3>> free_dirs;
  // Volume elements that share the entity; their Jacobians must stay valid.
  std::vector<ng::ElementIndex> volumes;
};

// Interior sample points of a volume element's reference domain.
const std::vector<ng::Point<3>> & volume_ref_points(ng::ELEMENT_TYPE type)
{
  static std::unordered_map<int, std::vector<ng::Point<3>>> cache;
  auto it = cache.find(int(type));
  if (it != cache.end()) return it->second;
  std::vector<ng::Point<3>> pts;
  const int m = 5;
  for (int i = 0; i < m; i++)
    for (int j = 0; j < m; j++)
      for (int k = 0; k < m; k++) {
        double a = (i + 0.5) / m, b = (j + 0.5) / m, c = (k + 0.5) / m;
        bool inside = true;
        switch (type) {
          case ng::TET: inside = a + b + c < 1; break;
          case ng::PRISM: inside = a + b < 1; break;
          case ng::PYRAMID: inside = a < 1 - c && b < 1 - c; break;
          default: break;
        }
        if (inside) pts.emplace_back(a, b, c);
      }
  return cache.emplace(int(type), std::move(pts)).first->second;
}

// Volume elements containing all the given vertices.
std::vector<ng::ElementIndex> shared_volumes(const ng::MeshTopology & top,
                                             const std::vector<ng::PointIndex> & verts)
{
  std::vector<ng::ElementIndex> out;
  if (verts.empty()) return out;
  for (auto ei : top.GetVertexElements(verts[0])) {
    bool all = true;
    for (size_t k = 1; k < verts.size() && all; k++) {
      bool found = false;
      for (auto ej : top.GetVertexElements(verts[k]))
        if (ej == ei) { found = true; break; }
      all = found;
    }
    if (all) out.push_back(ei);
  }
  return out;
}

struct Refitter
{
  ng::Mesh & mesh;
  ng::CurvedElements & curved;
  const ng::NetgenGeometry & geo;

  Refitter(ng::Mesh & m)
    : mesh(m), curved(m.GetCurvedElements()), geo(*m.GetGeometry()) {}

  ng::Point<3> eval(ng::SurfaceElementIndex sei, const std::array<double, 2> & r)
  {
    ng::Point<3> x;
    curved.CalcSurfaceTransformation(ng::Point<2>(r[0], r[1]), sei, x);
    return x;
  }

  // Closest point on the geometry and the unit directions whose distance
  // components are minimised (the surface normal, or the plane normal to
  // a curve).  Returns false when the projection is not trustworthy.
  bool project(const Entity & e, const std::array<double, 2> & r,
               const ng::Point<3> & x, ng::Point<3> & q,
               std::vector<ng::Vec<3>> & dirs)
  {
    const ng::Element2d & el = mesh[e.sei];
    q = x;
    dirs.clear();
    if (e.on_curve) {
      ng::EdgePointGeomInfo egi;
      geo.ProjectPointEdge(e.surfnr1, e.surfnr2, q, &egi);
      ng::Vec<3> t = geo.GetTangent(q, e.surfnr1, e.surfnr2, egi);
      double tl = t.Length();
      if (!(tl > 0)) return false;
      t /= tl;
      // Orthonormal basis of the plane normal to the curve tangent.
      ng::Vec<3> a = std::fabs(t(0)) < 0.9 ? ng::Vec<3>(1, 0, 0) : ng::Vec<3>(0, 1, 0);
      ng::Vec<3> n1 = a - (a * t) * t;
      n1 /= n1.Length();
      ng::Vec<3> n2 = ng::Cross(t, n1);
      dirs.push_back(n1);
      dirs.push_back(n2);
    } else {
      double lam[4];
      vertex_weights(el.GetType(), r[0], r[1], lam);
      ng::PointGeomInfo gi = el.GeomInfoPi(1);
      gi.u = gi.v = 0;
      for (int k = 0; k < el.GetNP(); k++) {
        gi.u += lam[k] * el.GeomInfoPi(k + 1).u;
        gi.v += lam[k] * el.GeomInfoPi(k + 1).v;
      }
      if (!geo.ProjectPointGI(e.surfnr, q, gi)) return false;
      ng::Vec<3> n = geo.GetNormal(e.surfnr, q, &gi);
      double nl = n.Length();
      if (!(nl > 0)) return false;
      dirs.push_back((1.0 / nl) * n);
    }
    // A projection that jumps by a sizeable fraction of the entity is a
    // wrong-branch result, not a curvature offset.
    return ng::Dist(x, q) <= 0.5 * e.scale;
  }

  // Max sampled distance to the geometry, or -1 if a projection failed.
  double max_distance(const Entity & e)
  {
    double d = 0;
    ng::Point<3> q;
    std::vector<ng::Vec<3>> dirs;
    for (const auto & r : e.check_ref) {
      ng::Point<3> x = eval(e.sei, r);
      if (!project(e, r, x, q, dirs)) return -1;
      d = std::max(d, ng::Dist(x, q));
    }
    return d;
  }

  // Surface area element |dx/dxi x dx/deta| at the entity's sample points.
  std::vector<double> area_samples(const Entity & e)
  {
    std::vector<double> a;
    for (const auto * pts : {&e.check_ref, &e.fit_ref})
      for (const auto & r : *pts) {
        ng::Point<3> x;
        ng::Mat<3, 2> jac;
        curved.CalcSurfaceTransformation(ng::Point<2>(r[0], r[1]), e.sei, x, jac);
        ng::Vec<3> c0(jac(0, 0), jac(1, 0), jac(2, 0)), c1(jac(0, 1), jac(1, 1), jac(2, 1));
        a.push_back(ng::Cross(c0, c1).Length());
      }
    return a;
  }

  // Jacobian determinants of the adjacent volume elements at interior
  // sample points (empty for a surface-only mesh).
  std::vector<double> volume_det_samples(const Entity & e)
  {
    std::vector<double> d;
    for (auto ei : e.volumes)
      for (const auto & xi : volume_ref_points(mesh[ei].GetType())) {
        ng::Point<3> x;
        ng::Mat<3, 3> jac;
        curved.CalcElementTransformation(xi, ei, x, jac);
        d.push_back(jac(0, 0) * (jac(1, 1) * jac(2, 2) - jac(1, 2) * jac(2, 1))
                  - jac(0, 1) * (jac(1, 0) * jac(2, 2) - jac(1, 2) * jac(2, 0))
                  + jac(0, 2) * (jac(1, 0) * jac(2, 1) - jac(1, 1) * jac(2, 0)));
      }
    return d;
  }

  // True when every sample keeps its sign and changes by at most a factor
  // of two.  Netgen elements may carry a negative determinant by vertex
  // order, so the comparison is pointwise against the same sample.
  static bool within_factor_two(const std::vector<double> & before,
                                const std::vector<double> & after)
  {
    for (size_t i = 0; i < before.size(); i++) {
      // A sample that was already zero or non-finite cannot be compared.
      if (before[i] == 0 || !std::isfinite(before[i])) continue;
      double r = after[i] / before[i];
      if (!(r >= 0.5 && r <= 2.0)) return false;
    }
    return true;
  }

  // Distance-minimising refit of the `n` coefficient vectors at `coef`.
  // Returns +1 accepted, 0 unchanged (no gain), -1 failed.
  int refit(const Entity & e, ng::Vec<3> * coef, int n,
            double & before, double & after)
  {
    const int np = (int)e.fit_ref.size();
    const int nu = 3 * n;
    std::vector<ng::Vec<3>> saved(coef, coef + n);

    before = max_distance(e);
    if (before < 0) return -1;
    const std::vector<double> area0 = area_samples(e);
    const std::vector<double> det0 = volume_det_samples(e);
    // Trust region: one step, and the total change, stay a fraction of the
    // entity size.  The distance cannot see a slide along the surface, so
    // without a bound an ill-conditioned step can wrap an element around it.
    const double max_step = 0.1 * e.scale, max_change = 0.25 * e.scale;

    // Basis values: the mapping is affine in the coefficients.
    std::vector<ng::Point<3>> x0(np);
    for (int j = 0; j < np; j++) x0[j] = eval(e.sei, e.fit_ref[j]);
    std::vector<double> phi(n * np);
    for (int k = 0; k < n; k++) {
      coef[k](0) += 1.0;
      for (int j = 0; j < np; j++)
        phi[k * np + j] = eval(e.sei, e.fit_ref[j])(0) - x0[j](0);
      coef[k](0) -= 1.0;
    }

    // Weighted normal-distance residual, its Gauss-Newton normal equations
    // (A, b) at the coefficient change d.  False if a projection fails.
    std::vector<ng::Vec<3>> dirs;
    ng::Point<3> q;
    auto evaluate = [&](const std::vector<double> & d, double & res,
                        std::vector<double> & A, std::vector<double> & b) {
      A.assign(nu * nu, 0.0);
      b.assign(nu, 0.0);
      res = 0;
      for (int j = 0; j < np; j++) {
        ng::Point<3> x = x0[j];
        for (int k = 0; k < n; k++)
          for (int c = 0; c < 3; c++) x(c) += phi[k * np + j] * d[3 * k + c];
        if (!project(e, e.fit_ref[j], x, q, dirs)) return false;
        for (const auto & nv : dirs) {
          double r = nv * (x - q);
          res += e.fit_w[j] * r * r;
          for (int k = 0; k < n; k++)
            for (int c = 0; c < 3; c++) {
              double jk = nv(c) * phi[k * np + j];
              b[3 * k + c] -= e.fit_w[j] * jk * r;
              for (int l = 0; l < n; l++)
                for (int cc = 0; cc < 3; cc++)
                  A[(3 * k + c) * nu + 3 * l + cc] +=
                      e.fit_w[j] * jk * nv(cc) * phi[l * np + j];
            }
        }
      }
      return true;
    };

    // Reduced unknowns: coefficient k changes by sum_m u[k*nd+m] * dir_m.
    std::vector<ng::Vec<3>> fd = e.free_dirs;
    if (fd.empty())
      fd = {ng::Vec<3>(1, 0, 0), ng::Vec<3>(0, 1, 0), ng::Vec<3>(0, 0, 1)};
    const int nd = (int)fd.size();
    const int nr = n * nd;
    auto reduce = [&](const std::vector<double> & A3, const std::vector<double> & b3,
                      std::vector<double> & Ar, std::vector<double> & br) {
      Ar.assign(nr * nr, 0.0);
      br.assign(nr, 0.0);
      for (int k = 0; k < n; k++)
        for (int m = 0; m < nd; m++) {
          int i = k * nd + m;
          for (int c = 0; c < 3; c++) br[i] += fd[m](c) * b3[3 * k + c];
          for (int l = 0; l < n; l++)
            for (int mm = 0; mm < nd; mm++) {
              double v = 0;
              for (int c = 0; c < 3; c++)
                for (int cc = 0; cc < 3; cc++)
                  v += fd[m](c) * A3[(3 * k + c) * nu + 3 * l + cc] * fd[mm](cc);
              Ar[i * nr + l * nd + mm] = v;
            }
        }
    };

    // The unknowns are the reduced changes u: coefficient k moves by
    // sum_m u[k*nd+m] * dir_m.  The objective is the weighted squared
    // distance plus a Tikhonov term lambda |u|^2.  Near-flat valleys of
    // equivalent re-parametrisations otherwise leave the minimiser
    // undetermined, and roundoff in the CAD projection then selects a
    // different point (and a different accept/reject outcome) per run.
    // lambda removes only directions whose curvature is below 1e-10 of the
    // strongest one.
    std::vector<double> d(nu, 0.0), A, b, Ar, br, M, s;
    auto to_d = [&](const std::vector<double> & uu) {
      std::vector<double> dd(nu, 0.0);
      for (int k = 0; k < n; k++)
        for (int m = 0; m < nd; m++)
          for (int c = 0; c < 3; c++) dd[3 * k + c] += uu[k * nd + m] * fd[m](c);
      return dd;
    };
    auto max_abs = [](const std::vector<double> & v) {
      double t = 0;
      for (double x : v) t = std::max(t, std::fabs(x));
      return t;
    };
    double res;
    if (!evaluate(d, res, A, b)) {
      for (int k = 0; k < n; k++) coef[k] = saved[k];
      return -1;
    }
    reduce(A, b, Ar, br);
    double diag0 = 0;
    for (int i = 0; i < nr; i++) diag0 = std::max(diag0, Ar[i * nr + i]);
    const double lambda = 1e-10 * diag0 + 1e-300;
    auto objective = [&](double r, const std::vector<double> & uu) {
      double t = r;
      for (double v : uu) t += lambda * v * v;
      return t;
    };
    // Regularised normal equations at u: (Ar + (lambda + damp) I) s = br - lambda u.
    auto step_for = [&](const std::vector<double> & uu, double damp, std::vector<double> & ss) {
      M = Ar;
      ss = br;
      for (int i = 0; i < nr; i++) {
        M[i * nr + i] += lambda + damp * std::max(Ar[i * nr + i], 1e-12 * diag0);
        ss[i] -= lambda * uu[i];
      }
      return solve_dense(M, ss, nr);
    };

    // Gauss-Newton on the full (non-monotone) path: the tangential unknowns
    // couple to the distance only through the normal's rotation, and the
    // path to the minimum typically raises the residual on its first steps,
    // where a monotone method creeps.  Steps are capped by the trust region;
    // the best iterate is kept.
    std::vector<double> u(nr, 0.0), best_u = u;
    double obj = objective(res, u), best_obj = obj;
    int lm_it = 0;
    for (; lm_it < 25; lm_it++) {
      if (!step_for(u, 0.0, s)) break;
      double step = max_abs(s);
      const double cut = step > max_step ? max_step / step : 1.0;
      for (int i = 0; i < nr; i++) u[i] += cut * s[i];
      d = to_d(u);
      if (max_abs(d) > max_change) break;    // outside the trust region
      if (!evaluate(d, res, A, b)) break;    // left the valid projection range
      reduce(A, b, Ar, br);
      obj = objective(res, u);
      if (obj < best_obj) {
        best_obj = obj;
        best_u = u;
      }
      if (step * cut < 1e-15 * e.scale) break;
    }

    // Monotone Levenberg-Marquardt polish from the best iterate.
    u = best_u;
    d = to_d(u);
    if (!evaluate(d, res, A, b)) {
      for (int k = 0; k < n; k++) coef[k] = saved[k];
      return -1;
    }
    reduce(A, b, Ar, br);
    obj = objective(res, u);
    double mu = 1e-6;
    std::vector<double> ut, dt, At, bt;
    for (int it = 0; it < 60 && mu < 1e8; it++, lm_it++) {
      if (!step_for(u, mu, s)) break;
      ut = u;
      for (int i = 0; i < nr; i++) ut[i] += s[i];
      dt = to_d(ut);
      double rt;
      if (max_abs(dt) <= max_change && evaluate(dt, rt, At, bt) &&
          objective(rt, ut) < obj) {
        u = ut;
        d = dt;
        res = rt;
        obj = objective(rt, ut);
        A = At;
        b = bt;
        reduce(A, b, Ar, br);
        mu = std::max(mu * 0.1, 1e-12);
        if (max_abs(s) < 1e-15 * e.scale) break;
      } else {
        mu *= 10;
      }
    }

    for (int k = 0; k < n; k++)
      for (int c = 0; c < 3; c++) coef[k](c) = saved[k](c) + d[3 * k + c];

    after = max_distance(e);
    // Keep the refit only if it is closer to the geometry and every surface
    // area element and adjacent volume Jacobian stays within a factor of two
    // of its previous value at the same sample: a re-parametrisation must
    // neither fold nor blow up an element.
    const std::vector<double> area1 = area_samples(e), det1 = volume_det_samples(e);
    const bool areas_ok = within_factor_two(area0, area1);
    const bool volumes_ok = within_factor_two(det0, det1);
    if (std::getenv("CUBIT_MESH_EXPORT_GEOMETRIC_REFIT_DEBUG")) {
      int nonfinite = 0;
      auto range = [&](const std::vector<double> & a, const std::vector<double> & b) {
        double lo = 1e300, hi = -1e300;
        for (size_t i = 0; i < a.size(); i++) {
          if (!std::isfinite(a[i]) || !std::isfinite(b[i])) { nonfinite++; continue; }
          if (a[i] != 0) { lo = std::min(lo, b[i] / a[i]); hi = std::max(hi, b[i] / a[i]); }
        }
        return std::pair<double, double>(lo, hi);
      };
      auto ra = range(area0, area1), rv = range(det0, det1);
      double total = 0;
      for (int i = 0; i < nu; i++) total = std::max(total, std::fabs(d[i]));
      std::fprintf(stderr, "refit n=%d on_curve=%d before=%.3e after=%.3e "
                   "fit_rms=%.3e it=%d area_ratio[%.3f,%.3f] det_ratio[%.3f,%.3f] "
                   "nonfinite=%d areas_ok=%d volumes_ok=%d change/scale=%.3e\n",
                   n, int(e.on_curve), before, after, std::sqrt(res), lm_it,
                   ra.first, ra.second, rv.first, rv.second, nonfinite,
                   int(areas_ok), int(volumes_ok), total / e.scale);
    }
    const bool failed = after < 0;
    const bool degenerate = !areas_ok || !volumes_ok;
    const bool closer = !failed && after < before;
    if (!closer || degenerate) {
      for (int k = 0; k < n; k++) coef[k] = saved[k];
      after = before;
      if (failed) return -1;
      return closer ? -2 : 0;   // -2: closer, but it would distort elements
    }
    return 1;
  }
};

}  // namespace

GeometricRefitStats geometric_refit(ng::Mesh & mesh)
{
  GeometricRefitStats st;
  if (!mesh.GetGeometry() || !mesh.GetCurvedElements().IsHighOrder()) return st;

  Refitter R(mesh);
  const ng::MeshTopology & top = mesh.GetTopology();
  const int order = R.curved.GetOrder();

  // Edges that lie on a geometric curve (mesh segments).
  std::unordered_map<int, int> curve_edge;   // edge -> segment
  for (int i = 0; i < mesh.GetNSeg(); i++)
    curve_edge[int(top.GetEdge(ng::SegmentIndex(i)))] = i;

  std::vector<double> g1, w1;
  gauss_rule(order + 6, g1, w1);
  const int ncheck = 4 * order + 9;

  // ---- Edges -------------------------------------------------------------
  std::vector<char> done(top.GetNEdges(), 0);
  std::vector<int> edge_surf(top.GetNEdges(), -2);
  for (int s = 0; s < mesh.GetNSE(); s++) {
    ng::SurfaceElementIndex sei(s);
    int sn = mesh.GetFaceDescriptor(mesh[sei].GetIndex()).SurfNr();
    for (auto e : top.GetEdges(sei)) {
      int ei = int(e);
      if (edge_surf[ei] == -2) edge_surf[ei] = sn;
      else if (edge_surf[ei] != sn) edge_surf[ei] = -1;   // two surfaces
    }
  }

  for (int s = 0; s < mesh.GetNSE(); s++) {
    ng::SurfaceElementIndex sei(s);
    const ng::Element2d & el = mesh[sei];
    auto edges = top.GetEdges(sei);
    auto local = ng::MeshTopology::GetEdges(el.GetType());
    for (int i = 0; i < edges.Size(); i++) {
      int ei = int(edges[i]);
      if (done[ei]) continue;
      done[ei] = 1;
      int n = R.curved.NumEdgeCoefficients(ei);
      if (n <= 0) continue;

      Entity E;
      E.sei = sei;
      auto it = curve_edge.find(ei);
      if (it != curve_edge.end()) {
        const ng::Segment & seg = mesh[ng::SegmentIndex(it->second)];
        E.on_curve = true;
        E.surfnr1 = seg.surfnr1;
        E.surfnr2 = seg.surfnr2;
      } else if (edge_surf[ei] >= 0) {
        E.surfnr = edge_surf[ei];
      } else {
        continue;   // an edge between surfaces without a geometric curve
      }
      auto ra = vertex_ref(el.GetType(), local[i][0]);
      auto rb = vertex_ref(el.GetType(), local[i][1]);
      E.scale = ng::Dist(mesh[el[local[i][0]]], mesh[el[local[i][1]]]);
      auto at = [&](double t) {
        return std::array<double, 2>{ra[0] + t * (rb[0] - ra[0]),
                                     ra[1] + t * (rb[1] - ra[1])};
      };
      for (size_t j = 0; j < g1.size(); j++) {
        E.fit_ref.push_back(at(g1[j]));
        E.fit_w.push_back(w1[j]);
      }
      for (int j = 1; j <= ncheck; j++) E.check_ref.push_back(at(double(j) / (ncheck + 1)));
      if (!E.on_curve && E.scale > 0) {
        ng::Vec<3> t0 = mesh[el[local[i][1]]] - mesh[el[local[i][0]]];
        t0 *= 1.0 / E.scale;
        ng::Point<3> q;
        std::vector<ng::Vec<3>> nd;
        auto rm = at(0.5);
        if (R.project(E, rm, R.eval(sei, rm), q, nd)) {
          ng::Vec<3> n0 = nd[0] - (nd[0] * t0) * t0;
          double l = n0.Length();
          if (l > 1e-8) E.free_dirs = {t0, (1.0 / l) * n0};
        }
      }

      E.volumes = shared_volumes(top, {el[local[i][0]], el[local[i][1]]});

      st.edges_tried++;
      double before = 0, after = 0;
      int r = R.refit(E, R.curved.EdgeCoefficients(ei), n, before, after);
      if (r > 0) st.edges_accepted++;
      if (r == -1) st.edges_failed++;
      if (r == -2) st.edges_rejected_distortion++;
      st.edge_dist_before = std::max(st.edge_dist_before, before);
      st.edge_dist_after = std::max(st.edge_dist_after, after);
    }
  }

  // ---- Faces (after the edges they are blended with) ---------------------
  for (int s = 0; s < mesh.GetNSE(); s++) {
    ng::SurfaceElementIndex sei(s);
    const ng::Element2d & el = mesh[sei];
    int fi = top.GetFace(sei);
    int n = R.curved.NumFaceCoefficients(fi);
    if (n <= 0) continue;

    Entity E;
    E.sei = sei;
    E.surfnr = mesh.GetFaceDescriptor(el.GetIndex()).SurfNr();
    E.scale = 0;
    for (int k = 0; k < el.GetNP(); k++)
      E.scale = std::max(E.scale, ng::Dist(mesh[el[k]], mesh[el[(k + 1) % el.GetNP()]]));
    const bool trig = el.GetType() == ng::TRIG;
    for (size_t a = 0; a < g1.size(); a++)
      for (size_t b = 0; b < g1.size(); b++) {
        if (trig) {   // collapsed (Duffy) rule on the reference triangle
          double y = g1[b], x = (1 - y) * g1[a];
          E.fit_ref.push_back({x, y});
          E.fit_w.push_back(w1[a] * w1[b] * (1 - y));
        } else {
          E.fit_ref.push_back({g1[a], g1[b]});
          E.fit_w.push_back(w1[a] * w1[b]);
        }
      }
    const int m = order + 4;
    for (int a = 1; a <= m; a++)
      for (int b = 1; b <= m; b++) {
        double x = double(a) / (m + 1), y = double(b) / (m + 1);
        if (trig && x + y >= 1) continue;
        E.check_ref.push_back({x, y});
      }

    // Faces keep all three directions: restricting them to the normal was
    // measured 2-3 orders worse, because the tangential re-parametrisation
    // is what lets the bubble degrees reduce the shape error.
    {
      std::vector<ng::PointIndex> fv;
      for (int k = 0; k < el.GetNP(); k++) fv.push_back(el[k]);
      E.volumes = shared_volumes(top, fv);
    }
    st.faces_tried++;
    double before = 0, after = 0;
    int r = R.refit(E, R.curved.FaceCoefficients(fi), n, before, after);
    if (r > 0) st.faces_accepted++;
    if (r == -1) st.faces_failed++;
    if (r == -2) st.faces_rejected_distortion++;
    st.face_dist_before = std::max(st.face_dist_before, before);
    st.face_dist_after = std::max(st.face_dist_after, after);
  }
  return st;
}

#else

GeometricRefitStats geometric_refit(netgen::Mesh &) { return {}; }

#endif
