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
    const auto &descriptor_label = descriptor.GetBCName();
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
  const std::filesystem::path output("named-periodic-contract.vol");
  const auto paths = named_periodic_output_paths(output);
  save_named_periodic_vol(*mesh, output);
  if (!std::filesystem::is_regular_file(output) || std::filesystem::exists(paths.temporary))
    throw std::runtime_error("Native periodic save did not publish the exact final path");
  auto bytes = [](const std::filesystem::path &path) {
    std::ifstream stream(path, std::ios::binary);
    if (!stream) throw std::runtime_error("Cannot read save-contract artifact");
    return std::string(std::istreambuf_iterator<char>(stream), std::istreambuf_iterator<char>());
  };
  const auto published = bytes(output);
  if (published.empty()) throw std::runtime_error("Empty published periodic VOL");
  must_reject([&]{save_named_periodic_vol(*mesh, output);});
  if (bytes(output) != published)
    throw std::runtime_error("Existing periodic output was modified");
  const std::filesystem::path stale("named-periodic-stale.vol");
  const auto stale_paths = named_periodic_output_paths(stale);
  { std::ofstream sentinel(stale_paths.temporary); sentinel << "preserve temporary"; }
  const auto sentinel = bytes(stale_paths.temporary);
  must_reject([&]{save_named_periodic_vol(*mesh, stale);});
  if (bytes(stale_paths.temporary) != sentinel || std::filesystem::exists(stale))
    throw std::runtime_error("Existing periodic temporary was changed or published");
  const std::filesystem::path blocked("named-periodic-sidecar.vol");
  const auto blocked_paths = named_periodic_output_paths(blocked);
  { std::ofstream sidecar(blocked_paths.sidecar); sidecar << "preserve sidecar"; }
  const auto existing_sidecar = bytes(blocked_paths.sidecar);
  must_reject([&]{save_named_periodic_vol(*mesh, blocked);});
  if (bytes(blocked_paths.sidecar) != existing_sidecar || std::filesystem::exists(blocked))
    throw std::runtime_error("Existing periodic sidecar was changed or published");
  must_reject([&]{save_named_periodic_vol(*mesh, "unsupported-periodic-output.partial");});
  const std::filesystem::path missing_parent("named-periodic-no-parent/output.vol");
  if (std::filesystem::exists(missing_parent.parent_path()))
    throw std::runtime_error("Save-error fixture parent must be absent");
  must_reject([&]{save_named_periodic_vol(*mesh, missing_parent);});
  if (std::filesystem::exists(missing_parent))
    throw std::runtime_error("Failed periodic writer published an output");
  const auto rename_error = named_periodic_output_paths("named-periodic-rename-error.vol");
  if (std::filesystem::exists(rename_error.temporary) || std::filesystem::exists(rename_error.final))
    throw std::runtime_error("Rename-error fixture paths must be absent");
  must_reject([&]{publish_named_periodic_output(rename_error);});
  if (std::filesystem::exists(rename_error.final))
    throw std::runtime_error("Failed native rename published an output");
  // Exercise the OS no-replace primitive directly, including the state that
  // can arise between a caller's freshness check and publication.
  const std::filesystem::path race_source("named-periodic-race-source.vol");
  const std::filesystem::path race_target("named-periodic-race-target.vol");
  { std::ofstream file(race_source); file << "owned source"; }
  { std::ofstream file(race_target); file << "other target"; }
  must_reject([&]{move_named_periodic_no_replace(race_source, race_target);});
  if (bytes(race_source) != "owned source" || bytes(race_target) != "other target")
    throw std::runtime_error("No-replace publication modified an existing target");

  const std::filesystem::path record_path("named-periodic-mapping.json");
  write_named_periodic_sidecar(record_path, [&](std::ostream &record) {
    record << "{\"type\":" << static_cast<int>(Identifications::PERIODIC) << "}";
  });
  const auto record_bytes = bytes(record_path);
  if (record_bytes != "{\"type\":" + std::to_string(static_cast<int>(Identifications::PERIODIC)) + "}")
    throw std::runtime_error("Periodic mapping record incomplete");
  must_reject([&]{write_named_periodic_sidecar(record_path, [](std::ostream &record){record << "overwrite";});});
  if (bytes(record_path) != record_bytes)
    throw std::runtime_error("Exclusive sidecar writer overwrote existing data");
  const std::filesystem::path failed_record("named-periodic-failed-record.json");
  must_reject([&]{write_named_periodic_sidecar(failed_record, [](std::ostream &record){record.setstate(std::ios::badbit);});});
  if (std::filesystem::exists(failed_record))
    throw std::runtime_error("Failed sidecar serialization published a mapping record");
  must_reject([&]{write_named_periodic_sidecar("named-periodic-no-parent/mapping.json", [](std::ostream &record){record << "mapping";});});
  const auto bundle = named_periodic_output_paths("named-periodic-bundle.vol");
  save_named_periodic_bundle(*mesh, bundle.final, [&](std::ostream &record) {
    if (std::filesystem::exists(bundle.final) || !std::filesystem::is_regular_file(bundle.temporary))
      throw std::runtime_error("Final VOL must not publish before its mandatory record");
    record << "{\"export_passed\":true}";
  });
  if (!std::filesystem::is_regular_file(bundle.final) ||
      bytes(bundle.sidecar) != "{\"export_passed\":true}" || std::filesystem::exists(bundle.temporary))
    throw std::runtime_error("Complete named bundle did not publish after its record");
  const auto failed_bundle = named_periodic_output_paths("named-periodic-bundle-failed.vol");
  must_reject([&]{save_named_periodic_bundle(*mesh, failed_bundle.final,
      [](std::ostream &record){record.setstate(std::ios::badbit);});});
  if (std::filesystem::exists(failed_bundle.final) || std::filesystem::exists(failed_bundle.sidecar))
    throw std::runtime_error("Failed mandatory record left a published VOL or success record");
  const auto blocked_bundle = named_periodic_output_paths("named-periodic-bundle-blocked.vol");
  must_reject([&]{save_named_periodic_bundle(*mesh, blocked_bundle.final,
      [&](std::ostream &record) {
        std::ofstream concurrent(blocked_bundle.sidecar); concurrent << "preserve other record";
        record << "mapping";
      });});
  if (std::filesystem::exists(blocked_bundle.final) || bytes(blocked_bundle.sidecar) != "preserve other record")
    throw std::runtime_error("Exclusive sidecar failure published VOL or overwrote another record");
  const auto guidance = named_periodic_retry_guidance(paths);
  for (const auto &path : {paths.final, paths.temporary, paths.sidecar})
    if (guidance.find(path.u8string()) == std::string::npos)
      throw std::runtime_error("Failed export guidance omitted an exact leftover path");
  std::cerr << "named periodic save contract: native Save/rename passed; existing final/temp/sidecar preserved\n";
}
