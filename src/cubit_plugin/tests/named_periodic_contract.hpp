// Guards: the production writer persists a complete named periodic group;
// distinct FD/BC namespaces, seam deduplication and ambiguity fail-loud.
// Assumes: small, independently specified translational surface fixture.
#include "NamedPeriodicIdentification.hpp"
#include <memory>
#include <iostream>

inline void named_periodic_contract_test()
{
  using namespace netgen;
  using namespace cubit_mesh_export;
  auto fixture = [] {
    auto mesh = std::make_unique<Mesh>();
    mesh->SetDimension(3);
    for (int side = 0; side < 2; ++side)
      for (int y = 0; y < 3; ++y)
        for (int z = 0; z < 2; ++z)
          mesh->AddPoint(Point<3>(side == 0 ? 2.0 : 0.0, y, z));
    int index = 0;
    for (int bc : {4, 2, 1, 3}) {
      FaceDescriptor descriptor;
      descriptor.SetBCProperty(bc);
      descriptor.SetBCName(bc == 4 || bc == 2 ? "source_patch" : "target_patch");
      const auto fd = mesh->AddFaceDescriptor(descriptor);
      if (periodic_face_number(fd) != index + 1)
        throw std::runtime_error("Fixture descriptor numbering differs from explicit insertion order");
      const int base = index < 2 ? 1 + 2 * index : 7 + 2 * (index - 2);
      for (auto points : {std::array<int,3>{base, base+2, base+1},
                          std::array<int,3>{base+1, base+2, base+3}}) {
        Element2d triangle(TRIG);
        triangle.SetIndex(fd);
        for (int j = 0; j < 3; ++j) triangle[j] = periodic_point_index(points[j]);
        mesh->AddSurfaceElement(triangle);
      }
      ++index;
    }
#ifdef COMPACT_NETGEN
    mesh->SetNBCNames(4);
    for (int bc : {4, 2}) mesh->SetBCName(bc - 1, "source_patch");
    for (int bc : {1, 3}) mesh->SetBCName(bc - 1, "target_patch");
#endif
    if (mesh->GetNFD() != 4)
      throw std::runtime_error("Fixture must contain exactly four explicitly added descriptors");
    return mesh;
  };
  const NamedPeriodicContract contract{"kt", "source_patch", "target_patch", "m",
                                       {-2.0, 0.0, 0.0}, 0.125, {{4,1},{2,3}}};
  auto must_reject = [](auto action) {
    bool failed = false;
    try { action(); } catch (const std::runtime_error &) { failed = true; }
    if (!failed) throw std::runtime_error("Invalid periodic fixture was accepted");
  };
  auto missing = fixture();
  auto partial = contract; partial.boundary_pairs.pop_back();
  must_reject([&]{apply_named_periodic_identification(*missing, partial);});
  auto unmatched = fixture(); unmatched->Point(periodic_point_index(7))(0) = -1;
  must_reject([&]{apply_named_periodic_identification(*unmatched, contract);});
  auto ambiguous = fixture(); ambiguous->Point(periodic_point_index(8))(2) = 0;
  must_reject([&]{apply_named_periodic_identification(*ambiguous, contract);});
  auto duplicate = fixture(); duplicate->Point(periodic_point_index(2))(2) = 0;
  must_reject([&]{apply_named_periodic_identification(*duplicate, contract);});
  auto occupied = fixture(); occupied->GetIdentifications().Add(periodic_point_index(1),periodic_point_index(7),4);
  must_reject([&]{apply_named_periodic_identification(*occupied, contract);});
  auto extra_incidence = fixture();
  Element2d extra_triangle(TRIG);
  extra_triangle.SetIndex(periodic_face_index(4)); // FD4 is slave BC3.
  extra_triangle[0] = periodic_point_index(8);
  extra_triangle[1] = periodic_point_index(9);
  extra_triangle[2] = periodic_point_index(10);
  extra_incidence->AddSurfaceElement(extra_triangle);
  // Global point bijection and every forward master-incidence test still pass,
  // but slave BC3 has point8 without a counterpart in master BC2.
  must_reject([&]{apply_named_periodic_identification(*extra_incidence, contract);});
  if (extra_incidence->GetIdentifications().GetMaxNr() != 0)
    throw std::runtime_error("Rejected patch map mutated identification metadata");
  auto mesh = fixture();
  std::cerr << "named periodic diagnostic: GetNFD=" << mesh->GetNFD() << '\n';
  for (int fd = 1; fd <= mesh->GetNFD(); ++fd) {
    const auto &descriptor = mesh->GetFaceDescriptor(periodic_face_index(fd));
#ifdef COMPACT_NETGEN
    const auto &descriptor_label = descriptor.BCName();
#else
    const auto &descriptor_label = descriptor.GetBCName();
#endif
    std::cerr << "named periodic diagnostic: FD=" << fd
              << " BC=" << descriptor.BCProperty()
              << " descriptor_label=" << descriptor_label
              << " bc_property_label=" << periodic_boundary_label(*mesh, descriptor.BCProperty())
              << '\n';
  }
  auto result = apply_named_periodic_identification(*mesh, contract);
  std::cerr << "named periodic diagnostic: pair_count=" << result.point_pairs.size() << '\n';
  for (const auto &entry : result.descriptor_to_boundary)
    std::cerr << "named periodic diagnostic: returned FD=" << entry.first
              << " BC=" << entry.second << '\n';
  for (const auto &pair : result.point_pairs)
    std::cerr << "named periodic diagnostic: pair=" << pair.first
              << "->" << pair.second << '\n';
  if (result.point_pairs.size() != 6 || result.descriptor_to_boundary.at(1) != 4)
    throw std::runtime_error("Periodic group or FD/BC mapping differs: actual pair_count=" +
        std::to_string(result.point_pairs.size()) + ", expected=6; actual FD1 BC=" +
        std::to_string(result.descriptor_to_boundary.at(1)) + ", expected=4");
  for (int i = 0; i < 6; ++i)
    if (result.point_pairs[i] != std::make_pair(i + 1, i + 7))
      throw std::runtime_error("Periodic endpoint direction/indices differ");
  mesh->Save(std::filesystem::path("named-periodic-contract.vol"));
}
