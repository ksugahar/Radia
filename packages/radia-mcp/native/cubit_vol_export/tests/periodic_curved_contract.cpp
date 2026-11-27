// Guards: q2/q3 synchronization operates on real Netgen coefficient storage,
// with reversed edge orientation, shared patch edges and distinct FD/BC IDs.
// Assumes: synthetic polynomial tetrahedra, no CAD or refit operation.
// The storage-layout check inside Netgen (CurvedElements::ValidatePeriodicLayout)
// belongs to the Netgen fork; this contract drives it only through the exporter.
#include "NamedPeriodicIdentification.hpp"
#include <memory>
#include <iostream>
#include <limits>

using namespace netgen;
using namespace radia_cubit_export;

static std::unique_ptr<Mesh> fixture(int order, const std::array<int,4> &mapped = {8,7,6,5})
{
  auto mesh = std::make_unique<Mesh>();
  mesh->SetDimension(3);
  mesh->SetGeometry(std::make_shared<NetgenGeometry>());
  const std::array<std::array<double,3>,4> corners{{{0,0,0},{0,1,0},{0,0,1},{1,1,1}}};
  for (const auto &p : corners) mesh->AddPoint(Point<3>(p[0],p[1],p[2]));
  for (int point = 5; point <= 8; ++point) {
    const auto found = std::find(mapped.begin(),mapped.end(),point);
    if (found == mapped.end()) throw std::runtime_error("Fixture mapping must permute target points");
    const auto &p = corners[found-mapped.begin()]; mesh->AddPoint(Point<3>(p[0]+3,p[1],p[2]));
  }
  for (const auto &nodes : {std::array<int,4>{1,2,3,4}, mapped}) {
    Element tet(TET); tet.SetIndex(VolumeRegionIndex::FromNr1(1));
    for (int j = 0; j < 4; ++j) tet[j] = PointIndex::FromNr1(nodes[j]);
    mesh->AddVolumeElement(tet);
  }
  const std::array<std::array<int,3>,8> exterior{{{1,2,3},{2,4,3},
    {mapped[0],mapped[1],mapped[2]},{mapped[1],mapped[3],mapped[2]},
    {1,2,4},{1,3,4},{mapped[0],mapped[1],mapped[3]},
    {mapped[0],mapped[2],mapped[3]}}};
  std::set<std::array<int,3>> seeded_faces;
  for (int i = 0; i < 4; ++i) {
    auto key = exterior[i]; std::sort(key.begin(),key.end()); seeded_faces.insert(key);
  }
  int index = 0;
  for (int bc : {4,2,1,3,5,5,5,5}) {
    FaceDescriptor descriptor; descriptor.SetBCProperty(bc);
    descriptor.SetSurfNr(index);
    descriptor.SetDomainIn(1); descriptor.SetDomainOut(0);
    descriptor.SetBCName(index < 2 ? "source" : index < 4 ? "target" : "other_boundary");
    const auto fd = mesh->AddFaceDescriptor(descriptor);
    auto nodes = exterior[index];
    const std::array<int,4> tet_nodes = nodes[0] <= 4 ? std::array<int,4>{1,2,3,4} : mapped;
    int opposite = 0;
    for (int point : tet_nodes)
      if (std::find(nodes.begin(),nodes.end(),point) == nodes.end()) opposite = point;
    if (!opposite) throw std::runtime_error("Fixture boundary has no unique opposite vertex");
    const auto &a = mesh->Point(PointIndex::FromNr1(nodes[0]));
    const auto &b = mesh->Point(PointIndex::FromNr1(nodes[1]));
    const auto &c = mesh->Point(PointIndex::FromNr1(nodes[2]));
    const auto &inside = mesh->Point(PointIndex::FromNr1(opposite));
    std::array<double,3> ab, ac, inward;
    for (int axis = 0; axis < 3; ++axis) {
      ab[axis]=b(axis)-a(axis); ac[axis]=c(axis)-a(axis); inward[axis]=inside(axis)-a(axis);
    }
    const double side = (ab[1]*ac[2]-ab[2]*ac[1])*inward[0] +
      (ab[2]*ac[0]-ab[0]*ac[2])*inward[1] + (ab[0]*ac[1]-ab[1]*ac[0])*inward[2];
    if (side == 0) throw std::runtime_error("Degenerate fixture boundary face");
    if (side > 0) std::swap(nodes[1],nodes[2]); // outward from domain 1
    Element2d triangle(TRIG); triangle.SetIndex(fd);
    for (int j = 0; j < 3; ++j) triangle[j] = PointIndex::FromNr1(nodes[j]);
    mesh->AddSurfaceElement(triangle); ++index;
  }
  mesh->UpdateTopology();
  if (order > 1) {
    auto &curved = mesh->GetCurvedElements();
    curved.BuildCurvedElements(nullptr,order,false);
    curved.SetIsHighOrder(true);
    auto inactive_target_edge = std::array<int,2>{mapped[0],mapped[3]};
    std::sort(inactive_target_edge.begin(),inactive_target_edge.end());
    for (int e = 0; e < int(mesh->GetTopology().GetNEdges()); ++e) {
      const auto v = mesh->GetTopology().GetEdgeVertices(EdgeIndex::FromNr0(e));
      const std::array<int,2> key{v[0].Nr1(),v[1].Nr1()};
      const bool active = key != std::array<int,2>{1,4} && key != inactive_target_edge;
      for (int mode = 0; mode < curved.NumEdgeCoefficients(e); ++mode)
        for (int axis = 0; axis < 3; ++axis)
          curved.EdgeCoefficients(e)[mode](axis) = active ? (e+1)*(mode+1)*(axis+1)/1024.0 : 0;
    }
    for (int f = 0; f < int(mesh->GetTopology().GetNFaces()); ++f) {
      const auto v = mesh->GetTopology().GetFaceVertices(FaceIndex::FromNr0(f));
      std::array<int,3> key{v[0].Nr1(),v[1].Nr1(),v[2].Nr1()}; std::sort(key.begin(),key.end());
      for (int mode = 0; mode < curved.NumFaceCoefficients(f); ++mode)
        for (int axis = 0; axis < 3; ++axis)
          curved.FaceCoefficients(f)[mode](axis) = seeded_faces.count(key) ? (f+1)*(axis+1)/2048.0 : 0;
    }
  }
  return mesh;
}

int main()
{
  try {
    const NamedPeriodicContract contract{"kt","source","target","m",{3,0,0},0.125,{{4,1},{2,3}}};
    auto reject_unchanged = [&](auto &mesh, auto operation) {
      const auto before = periodic_geometry_serialization(*mesh);
      bool rejected = false;
      try { operation(); } catch (const std::exception &) { rejected = true; }
      if (!rejected || mesh->GetIdentifications().GetMaxNr() != 0 ||
          periodic_geometry_serialization(*mesh) != before)
        throw std::runtime_error("Rejected fixture must preserve geometry and empty IDs");
    };
    {
      // FD and BC numbers deliberately differ (FD1->BC4, FD2->BC2, FD3->BC1,
      // FD4->BC3, FD5..8->BC5): labels resolve by BC property, never by
      // treating a BC number as a descriptor position.
      auto labels = fixture(1);
      const std::map<int,std::string> expected{
        {4,"source"},{2,"source"},{1,"target"},{3,"target"},{5,"other_boundary"}};
      for (const auto &[bc, label] : expected)
        if (periodic_boundary_label(*labels,bc) != label)
          throw std::runtime_error("BC " + std::to_string(bc) + " resolves to '" +
              periodic_boundary_label(*labels,bc) + "', expected '" + label + "'");
      bool missing = false;
      try { periodic_boundary_label(*labels,6); } catch (const std::runtime_error &) { missing = true; }
      if (!missing) throw std::runtime_error("A BC with no descriptor must be rejected");
      labels->GetFaceDescriptor(periodic_face_index(8)).SetBCName("conflicting");
      bool conflicting = false;
      try { periodic_boundary_label(*labels,5); } catch (const std::runtime_error &) { conflicting = true; }
      if (!conflicting) throw std::runtime_error("Descriptors of one BC with different labels must be rejected");
      reject_unchanged(labels,[&]{
        auto relabelled = contract; relabelled.boundary_pairs = {{4,5}};
        apply_named_periodic_identification(*labels,relabelled);});
    }
    auto wrong_faces = fixture(3);
    // Same point sets, different diagonal: node coverage alone cannot detect it.
    wrong_faces->GetFaceDescriptor(periodic_face_index(2)).SetBCProperty(4);
    wrong_faces->GetFaceDescriptor(periodic_face_index(4)).SetBCProperty(1);
    auto &first = (*wrong_faces)[SurfaceElementIndex::FromNr0(2)];
    first[0] = PointIndex::FromNr1(8); first[1] = PointIndex::FromNr1(7); first[2] = PointIndex::FromNr1(5);
    auto &second = (*wrong_faces)[SurfaceElementIndex::FromNr0(3)];
    second[0] = PointIndex::FromNr1(8); second[1] = PointIndex::FromNr1(5); second[2] = PointIndex::FromNr1(6);
    auto one_patch = contract; one_patch.boundary_pairs = {{4,1}};
    reject_unchanged(wrong_faces,[&]{apply_named_periodic_identification(*wrong_faces,one_patch);});
    auto missing_face = fixture(2);
    const auto planned = periodic_triangle_plan(*missing_face,{{4,1},{2,3}},{{1,8},{2,7},{3,6},{4,5}});
    auto conflicting = planned;
    auto conflict = planned[0]; std::swap(conflict.slave[0],conflict.slave[1]);
    conflicting.push_back(conflict);
    reject_unchanged(missing_face,[&]{synchronize_periodic_curvature(*missing_face,conflicting);});
    auto unsupported = fixture(4);
    reject_unchanged(unsupported,[&]{apply_named_periodic_identification(*unsupported,contract);});
    auto rational = fixture(2);
    rational->GetCurvedElements().BuildCurvedElements(nullptr,2,true);
    rational->GetCurvedElements().SetIsHighOrder(true);
    reject_unchanged(rational,[&]{apply_named_periodic_identification(*rational,contract);});
    auto nonfinite = fixture(3);
    nonfinite->GetCurvedElements().EdgeCoefficients(0)[0](0) = std::numeric_limits<double>::quiet_NaN();
    reject_unchanged(nonfinite,[&]{apply_named_periodic_identification(*nonfinite,contract);});
    // Frozen oracle: odd bubble l_i*l_j*(l_i-l_j) at (3/4,1/4)
    // is +3/32, and exchanging its endpoints gives -3/32. The signs
    // below were tabulated for edges12,13,23,24,34, independently of
    // the producer's target-index comparison. First face covers all6 S3 cases.
    struct PermutationCase {
      std::array<int,4> mapped;
      std::array<int,5> odd_signs;
      const char *suffix;
    };
    const std::array<PermutationCase,6> cases{{
      {{8,7,6,5}, {-1,-1,-1,-1,-1}, ""},
      {{5,7,6,8}, { 1, 1,-1, 1, 1}, "-mixed"},
      {{5,6,7,8}, { 1, 1, 1, 1, 1}, "-perm123"},
      {{6,5,7,8}, {-1, 1, 1, 1, 1}, "-perm213"},
      {{6,7,5,8}, { 1,-1,-1, 1, 1}, "-perm231"},
      {{7,5,6,8}, {-1,-1, 1, 1, 1}, "-perm312"}
    }};
    for (const auto &entity : {std::string("edge [17,23]"),std::string("face [17,23,29]")}) {
      const std::map<int,int> missing;
      try {
        periodic_curved_lookup(missing,17,entity);
        throw std::runtime_error("Missing topology lookup unexpectedly succeeded");
      } catch (const std::runtime_error &error) {
        if (std::string(error.what()) != "Missing periodic " + entity) throw;
      }
    }
    for (const auto &permutation : cases) {
    const auto &mapped = permutation.mapped;
    const std::map<std::array<int,2>,int> odd_signs{
      {{1,2},permutation.odd_signs[0]}, {{1,3},permutation.odd_signs[1]},
      {{2,3},permutation.odd_signs[2]}, {{2,4},permutation.odd_signs[3]},
      {{3,4},permutation.odd_signs[4]}};
    for (int order : {1,2,3}) {
      auto mesh = fixture(order,mapped);
      const auto before = periodic_geometry_serialization(*mesh);
      // Independent fixture prediction covers every coefficient, including
      // master/unrelated entities. It does not consume the producer's plan.
      using Slot = std::array<int,3>; // face flag, native entity, mode
      std::map<Slot,std::array<double,3>> expected_coefficients;
      auto &curved = mesh->GetCurvedElements();
      const auto &top = mesh->GetTopology();
      if (order > 1) {
        std::map<std::array<int,2>,int> edge_ids;
        std::map<std::array<int,3>,int> face_ids;
        auto unpaired_edge = std::array<int,2>{mapped[0],mapped[3]};
        std::sort(unpaired_edge.begin(),unpaired_edge.end());
        auto target_a = std::array<int,3>{mapped[0],mapped[1],mapped[2]};
        auto target_b = std::array<int,3>{mapped[1],mapped[2],mapped[3]};
        std::sort(target_a.begin(),target_a.end()); std::sort(target_b.begin(),target_b.end());
        const std::set<std::array<int,3>> nonzero_faces{{1,2,3},{2,3,4},target_a,target_b};
        for (int e = 0; e < int(top.GetNEdges()); ++e) {
          auto v = top.GetEdgeVertices(EdgeIndex::FromNr0(e));
          const std::array<int,2> key{v[0].Nr1(),v[1].Nr1()};
          edge_ids[key] = e;
          for (int mode = 0; mode < curved.NumEdgeCoefficients(e); ++mode) {
            auto &expected = expected_coefficients[{0,e,mode}];
            for (int axis = 0; axis < 3; ++axis)
              expected[axis] = key == std::array<int,2>{1,4} || key == unpaired_edge ?
                0 : (e+1)*(mode+1)*(axis+1)/1024.0;
          }
        }
        for (int f = 0; f < int(top.GetNFaces()); ++f) {
          auto v = top.GetFaceVertices(FaceIndex::FromNr0(f));
          std::array<int,3> key{v[0].Nr1(),v[1].Nr1(),v[2].Nr1()};
          std::sort(key.begin(),key.end()); face_ids[key] = f;
          for (int mode = 0; mode < curved.NumFaceCoefficients(f); ++mode) {
            auto &expected = expected_coefficients[{1,f,mode}];
            for (int axis = 0; axis < 3; ++axis)
              expected[axis] = nonzero_faces.count(key) ? (f+1)*(axis+1)/2048.0 : 0;
          }
        }
        for (const auto &entry : edge_ids) {
          const auto &v = entry.first;
          if (v[1] > 4 || v == std::array<int,2>{1,4}) continue;
          std::array<int,2> image{mapped[v[0]-1],mapped[v[1]-1]};
          const int sign = odd_signs.at(v);
          std::sort(image.begin(),image.end());
          const int target = edge_ids.at(image);
          for (int mode = 0; mode < order-1; ++mode) {
            auto value = expected_coefficients.at({0,entry.second,mode});
            if (mode == 1) for (double &axis : value) axis *= sign;
            expected_coefficients[{0,target,mode}] = value;
          }
        }
        if (order == 3) for (const auto &key : {std::array<int,3>{1,2,3},std::array<int,3>{2,3,4}}) {
          const int source = face_ids.at(key);
          std::array<int,3> image{mapped[key[0]-1],mapped[key[1]-1],mapped[key[2]-1]};
          std::sort(image.begin(),image.end());
          const int target = face_ids.at(image);
          expected_coefficients[{1,target,0}] = expected_coefficients.at({1,source,0});
        }
      }
      auto result = apply_named_periodic_identification(*mesh,contract);
      for (const auto &entry : expected_coefficients) {
        const auto &key = entry.first;
        auto *actual = key[0] ? curved.FaceCoefficients(key[1]) : curved.EdgeCoefficients(key[1]);
        for (int axis = 0; axis < 3; ++axis)
          if (actual[key[2]](axis) != entry.second[axis])
            throw std::runtime_error("Coefficient differs from independently predicted complete storage");
      }
      if (result.point_pairs.size() != 4 || result.descriptor_to_boundary.at(1) != 4)
        throw std::runtime_error("Expected four mapped points and FD1->BC4");
      if (order == 1 && (!result.curved_transfers.empty() || periodic_geometry_serialization(*mesh) != before))
        throw std::runtime_error("Affine geometry must remain unchanged");
      int odd = 0, faces = 0, negative = 0;
      for (const auto &row : result.curved_transfers) {
        auto &curved = mesh->GetCurvedElements();
        auto *source = row.face ? curved.FaceCoefficients(row.master) : curved.EdgeCoefficients(row.master);
        auto *target = row.face ? curved.FaceCoefficients(row.slave) : curved.EdgeCoefficients(row.slave);
        int expected_sign = 1;
        if (!row.face && row.mode == 1) {
          const auto vertices = top.GetEdgeVertices(EdgeIndex::FromNr0(row.master));
          expected_sign = odd_signs.at({vertices[0].Nr1(),vertices[1].Nr1()});
          negative += expected_sign < 0;
        }
        if (row.sign != expected_sign) throw std::runtime_error("Wrong independently derived mode sign");
        for (int axis = 0; axis < 3; ++axis)
          if (target[row.mode](axis) != expected_sign * source[row.mode](axis))
            throw std::runtime_error("Actual slave coefficient differs from independent expected value");
        odd += !row.face && row.mode == 1; faces += row.face;
      }
      const int expected = order == 1 ? 0 : order == 2 ? 5 : 12;
      if (int(result.curved_transfers.size()) != expected ||
          (order == 3 && (odd != 5 || faces != 2 || negative != int(std::count(permutation.odd_signs.begin(),permutation.odd_signs.end(),-1)))))
        throw std::runtime_error("Expected deduplicated shared-edge transfers for declared order");
      std::ofstream saved("periodic-curved-q"+std::to_string(order)+permutation.suffix+".vol",std::ios::binary);
      mesh->Save(saved); saved.close();
      if (!saved) throw std::runtime_error("Fixture VOL write failed");
      std::cerr << "q" << order << " transfers=" << result.curved_transfers.size() << '\n';
    }
    }
    return 0;
  } catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
}
