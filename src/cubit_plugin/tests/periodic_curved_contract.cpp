// Guards: q2/q3 synchronization operates on real compact coefficient storage,
// with reversed edge orientation, shared patch edges and distinct FD/BC IDs.
// Assumes: synthetic polynomial tetrahedra, no CAD or refit operation.
#include "NamedPeriodicIdentification.hpp"
#include <memory>
#include <iostream>
#include <limits>

using namespace netgen;
using namespace cubit_mesh_export;

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
    Element tet(TET); tet.SetIndex(1);
    for (int j = 0; j < 4; ++j) tet[j] = PointIndex(nodes[j]);
    mesh->AddVolumeElement(tet);
  }
  int index = 0;
  for (int bc : {4,2,1,3}) {
    FaceDescriptor descriptor; descriptor.SetBCProperty(bc);
    descriptor.SetBCName(index < 2 ? "source" : "target");
    const auto fd = mesh->AddFaceDescriptor(descriptor);
    const std::array<std::array<int,3>,4> nodes{{{1,2,3},{2,4,3},
      {mapped[0],mapped[1],mapped[2]},{mapped[1],mapped[3],mapped[2]}}};
    Element2d triangle(TRIG); triangle.SetIndex(fd);
    for (int j = 0; j < 3; ++j) triangle[j] = PointIndex(nodes[index][j]);
    mesh->AddSurfaceElement(triangle); ++index;
  }
  mesh->SetNBCNames(4);
  for (int bc : {4,2}) mesh->SetBCName(bc-1,"source");
  for (int bc : {1,3}) mesh->SetBCName(bc-1,"target");
  mesh->GetTopology().Update();
  if (order > 1) {
    auto &curved = mesh->GetCurvedElements();
    curved.BuildCurvedElements(nullptr,order,false);
    curved.SetIsHighOrder(true);
    for (int e = 0; e < int(mesh->GetTopology().GetNEdges()); ++e)
      for (int mode = 0; mode < curved.NumEdgeCoefficients(e); ++mode)
        for (int axis = 0; axis < 3; ++axis)
          curved.EdgeCoefficients(e)[mode](axis) = (e+1)*(mode+1)*(axis+1)/1024.0;
    for (int f = 0; f < int(mesh->GetTopology().GetNFaces()); ++f)
      for (int mode = 0; mode < curved.NumFaceCoefficients(f); ++mode)
        for (int axis = 0; axis < 3; ++axis)
          curved.FaceCoefficients(f)[mode](axis) = (f+1)*(axis+1)/2048.0;
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
    auto wrong_faces = fixture(3);
    // Same point sets, different diagonal: node coverage alone cannot detect it.
    wrong_faces->GetFaceDescriptor(2).SetBCProperty(4);
    wrong_faces->GetFaceDescriptor(4).SetBCProperty(1);
    auto &first = (*wrong_faces)[SurfaceElementIndex(2)];
    first[0] = PointIndex(8); first[1] = PointIndex(7); first[2] = PointIndex(5);
    auto &second = (*wrong_faces)[SurfaceElementIndex(3)];
    second[0] = PointIndex(8); second[1] = PointIndex(5); second[2] = PointIndex(6);
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
    auto layout_guard = fixture(3);
    NgArray<int> orders, indices;
    NgArray<Vec<3>> coefficients;
    orders.SetSize(3); orders = 3;
    indices.SetSize(4); indices[0]=0; indices[1]=2; indices[2]=4; indices[3]=4;
    coefficients.SetSize(5); coefficients = Vec<3>(1,2,3);
    // The first two entities are valid; the unpaired last entity is not.
    reject_unchanged(layout_guard,[&]{CurvedElements::ValidatePeriodicLayout(orders,indices,coefficients,3,false);});
    orders[2]=1;
    CurvedElements::ValidatePeriodicLayout(orders,indices,coefficients,3,false); // finite unused tail allowed
    orders[2]=4;
    reject_unchanged(layout_guard,[&]{CurvedElements::ValidatePeriodicLayout(orders,indices,coefficients,3,false);});
    orders[2]=0;
    reject_unchanged(layout_guard,[&]{CurvedElements::ValidatePeriodicLayout(orders,indices,coefficients,3,false);});
    orders=3; indices[1]=1; indices[2]=2; indices[3]=2;
    reject_unchanged(layout_guard,[&]{CurvedElements::ValidatePeriodicLayout(orders,indices,coefficients,3,true);});
    orders[2]=1;
    CurvedElements::ValidatePeriodicLayout(orders,indices,coefficients,3,true);
    indices[3]=6;
    reject_unchanged(layout_guard,[&]{CurvedElements::ValidatePeriodicLayout(orders,indices,coefficients,3,true);});
    for (const auto &mapped : {std::array<int,4>{8,7,6,5},std::array<int,4>{5,7,6,8}}) {
    const bool mixed = mapped[0] == 5;
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
        for (int e = 0; e < int(top.GetNEdges()); ++e) {
          auto v = top.GetEdgeVertices(e);
          edge_ids[{int(v[0]),int(v[1])}] = e;
          for (int mode = 0; mode < curved.NumEdgeCoefficients(e); ++mode) {
            auto &expected = expected_coefficients[{0,e,mode}];
            for (int axis = 0; axis < 3; ++axis) expected[axis] = (e+1)*(mode+1)*(axis+1)/1024.0;
          }
        }
        for (int f = 0; f < int(top.GetNFaces()); ++f) {
          auto v = top.GetFaceVertices(f);
          std::array<int,3> key{int(v[0]),int(v[1]),int(v[2])};
          std::sort(key.begin(),key.end()); face_ids[key] = f;
          for (int mode = 0; mode < curved.NumFaceCoefficients(f); ++mode) {
            auto &expected = expected_coefficients[{1,f,mode}];
            for (int axis = 0; axis < 3; ++axis) expected[axis] = (f+1)*(axis+1)/2048.0;
          }
        }
        for (const auto &entry : edge_ids) {
          const auto &v = entry.first;
          if (v[1] > 4 || v == std::array<int,2>{1,4}) continue;
          std::array<int,2> image{mapped[v[0]-1],mapped[v[1]-1]};
          const int sign = image[0] < image[1] ? 1 : -1;
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
          const auto vertices = top.GetEdgeVertices(row.master);
          expected_sign = mapped[int(vertices[0])-1] < mapped[int(vertices[1])-1] ? 1 : -1;
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
          (order == 3 && (odd != 5 || faces != 2 || negative != (mixed ? 1 : 5))))
        throw std::runtime_error("Expected deduplicated shared-edge transfers for declared order");
      std::ofstream saved("periodic-curved-q"+std::to_string(order)+(mixed ? "-mixed" : "")+".vol",std::ios::binary);
      mesh->Save(saved); saved.close();
      if (!saved) throw std::runtime_error("Fixture VOL write failed");
      std::cerr << "q" << order << " transfers=" << result.curved_transfers.size() << '\n';
    }
    }
    return 0;
  } catch (const std::exception &error) { std::cerr << error.what() << '\n'; return 1; }
}
