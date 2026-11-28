// Moved from the Radia Netgen fork (CurvedElements::CopyPeriodicCoefficients);
// the logic is unchanged, the storage access goes through curved_access.hpp.
#include "periodic_copy.hpp"
#include "curved_access.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <string>

namespace radia_netgen_curving {

using netgen::Exception;

std::vector<PeriodicTransfer> copy_periodic_coefficients(
    netgen::Mesh & mesh,
    const std::vector<std::array<int, 3>> & master_tris,
    const std::vector<std::array<int, 3>> & slave_tris)
{
  CurvedAccess access = curved_access(mesh);
  if (master_tris.size() != slave_tris.size())
    throw Exception("CopyPeriodicCoefficients: master/slave triangle count mismatch");
  const int order = access.order();
  if (order < 1 || order > 3)
    throw Exception("CopyPeriodicCoefficients: supports polynomial orders 1 through 3 only");
  if (access.rational())
    throw Exception("CopyPeriodicCoefficients: rational geometry is unsupported");
  std::vector<PeriodicTransfer> plan;
  if (!access.high_order()) return plan;

  const auto & top = mesh.GetTopology();
  const size_t nedges = top.GetNEdges(), nfaces = top.GetNFaces();
  try { access.validate(nedges, nfaces); }
  catch (const std::exception & e)
    { throw Exception(std::string("CopyPeriodicCoefficients: ") + e.what()); }

  std::map<std::array<int, 2>, int> edges;
  std::map<std::array<int, 3>, int> faces;
  for (size_t e = 0; e < nedges; e++)
    {
      auto v = top.GetEdgeVertices(netgen::IndexBASE<netgen::EdgeIndex>() + e);
      std::array<int, 2> key{v[0].Nr1(), v[1].Nr1()};
      if (key[0] >= key[1] || !edges.emplace(key, int(e)).second)
        throw Exception("CopyPeriodicCoefficients: noncanonical or duplicate topology edge");
    }
  for (size_t f = 0; f < nfaces; f++)
    {
      auto v = top.GetFaceVertices(netgen::IndexBASE<netgen::FaceIndex>() + f);
      if (v.Size() != 3) continue;
      std::array<int, 3> key{v[0].Nr1(), v[1].Nr1(), v[2].Nr1()};
      std::sort(key.begin(), key.end());
      if (!faces.emplace(key, int(f)).second)
        throw Exception("CopyPeriodicCoefficients: duplicate topology face");
    }
  auto lookup = [] (auto & table, auto key, const char * kind)
  {
    auto found = table.find(key);
    if (found == table.end())
      throw Exception(std::string("CopyPeriodicCoefficients: missing ") + kind);
    return found->second;
  };

  struct Owner { int master, sign; };
  std::map<int, Owner> edge_owners, face_owners;
  auto own = [] (auto & owners, int slave, int master, int sign)
  {
    auto ins = owners.emplace(slave, Owner{master, sign});
    if (!ins.second && (ins.first->second.master != master || ins.first->second.sign != sign))
      throw Exception("CopyPeriodicCoefficients: conflicting seam coefficient ownership");
  };
  for (size_t t = 0; t < master_tris.size(); t++)
    {
      auto m = master_tris[t], sl = slave_tris[t];
      auto ms = m, ss = sl;
      std::sort(ms.begin(), ms.end()); std::sort(ss.begin(), ss.end());
      own(face_owners, lookup(faces, ss, "face"), lookup(faces, ms, "face"), 1);
      for (int a = 0; a < 3; a++)
        for (int b = a + 1; b < 3; b++)
          {
            std::array<int, 2> me{m[a], m[b]}, te{sl[a], sl[b]};
            const int sign = te[0] < te[1] ? 1 : -1;
            std::sort(me.begin(), me.end()); std::sort(te.begin(), te.end());
            // the master edge keeps its own canonical direction; the image
            // direction is compared against the master's increasing order
            const int msign = m[a] < m[b] ? 1 : -1;
            own(edge_owners, lookup(edges, te, "edge"), lookup(edges, me, "edge"), sign * msign);
          }
    }

  auto coeffs = [&] (bool face, int entity)
  { return face ? access.face_coeffs(entity) : access.edge_coeffs(entity); };
  auto stage = [&] (bool face, auto & owners)
  {
    std::set<int> sources;
    for (auto & row : owners) sources.insert(row.second.master);
    for (auto & row : owners)
      {
        const int slave = row.first, master = row.second.master;
        if (sources.count(slave))
          throw Exception("CopyPeriodicCoefficients: coefficient source and destination overlap");
        const int count = face ? access.face_count(master) : access.edge_count(master);
        const int target_count = face ? access.face_count(slave) : access.edge_count(slave);
        const int sorder = face ? access.face_order(master) : access.edge_order(master);
        const int torder = face ? access.face_order(slave) : access.edge_order(slave);
        if (count != target_count || sorder != torder)
          throw Exception("CopyPeriodicCoefficients: coefficient order or layout mismatch");
        if (!count) continue;
        const netgen::Vec<3> * src = coeffs(face, master);
        const netgen::Vec<3> * dst = coeffs(face, slave);
        for (int mode = 0; mode < count; mode++)
          {
            const int sign = (!face && mode % 2 == 1) ? row.second.sign : 1;
            plan.push_back({face, master, slave, mode, sign,
                            {dst[mode](0), dst[mode](1), dst[mode](2)},
                            {sign * src[mode](0), sign * src[mode](1), sign * src[mode](2)}});
          }
      }
  };
  stage(false, edge_owners);
  stage(true, face_owners);

  auto put = [&] (const PeriodicTransfer & tr, const std::array<double, 3> & v)
  {
    netgen::Vec<3> & c = coeffs(tr.face, tr.slave)[tr.mode];
    for (int k = 0; k < 3; k++) c(k) = v[k];
  };
  try
    {
      for (auto & tr : plan) put(tr, tr.expected);
      for (auto & tr : plan)
        for (int k = 0; k < 3; k++)
          if (coeffs(tr.face, tr.slave)[tr.mode](k) != tr.expected[k])
            throw Exception("CopyPeriodicCoefficients: paired coefficient equality failed");
    }
  catch (...)
    {
      for (auto & tr : plan) put(tr, tr.original);
      throw;
    }
  return plan;
}

} // namespace radia_netgen_curving
