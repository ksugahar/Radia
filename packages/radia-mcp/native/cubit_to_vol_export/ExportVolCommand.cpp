#include "ExportVolCommand.hpp"
#include "MeshData.hpp"
#include "LearnEditionMessageFilter.hpp"
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <filesystem>
#include <sstream>

#include "NetgenCurver.hpp"
#include "NamedPeriodicIdentification.hpp"
#include "LiveVolMesh.hpp"
#include "netgen_meshing.hpp"

// Cubit geometry for CAD reference values (companion JSON)
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "RefVolume.hpp"
#include "GeometryQueryTool.hpp"

#ifdef _WIN32
#include <windows.h>
#endif

#include <fstream>
#include <map>
#include <set>
#include "utf8_path.hpp"

static std::string explicit_surface_label(int surface_id)
{
  std::string label = CubitInterface::get_entity_name("surface", surface_id);
  // Cubit reports its generated display name even when the user did not name
  // the entity.  That is not an authored boundary contract.
  if (label == "Surface " + std::to_string(surface_id) ||
      label == "Surface_" + std::to_string(surface_id) ||
      label == "surface_" + std::to_string(surface_id))
    return {};
  return label;
}

static double triangle_area(const netgen::Point<3> &p0,
                            const netgen::Point<3> &p1,
                            const netgen::Point<3> &p2)
{
  double ux = p1(0) - p0(0), uy = p1(1) - p0(1), uz = p1(2) - p0(2);
  double vx = p2(0) - p0(0), vy = p2(1) - p0(1), vz = p2(2) - p0(2);
  double cx = uy * vz - uz * vy;
  double cy = uz * vx - ux * vz;
  double cz = ux * vy - uy * vx;
  return 0.5 * std::sqrt(cx * cx + cy * cy + cz * cz);
}

static double bilinear_quad_area(const netgen::Point<3> &p0,
                                 const netgen::Point<3> &p1,
                                 const netgen::Point<3> &p2,
                                 const netgen::Point<3> &p3)
{
  // Integrate the bilinear surface Jacobian with an 8x8 Gauss rule.  The
  // Jacobian norm is not polynomial for a warped quad, so the usual 2x2 rule
  // leaves a small but visible discrepancy against NGSolve's area integral.
  const double q[8] = {
      0.0198550717512319,
      0.1016667612931866,
      0.2372337950418355,
      0.4082826787521751,
      0.5917173212478249,
      0.7627662049581645,
      0.8983332387068134,
      0.9801449282487681,
  };
  const double w[8] = {
      0.0506142681451881,
      0.1111905172266872,
      0.1568533229389437,
      0.1813418916891810,
      0.1813418916891810,
      0.1568533229389437,
      0.1111905172266872,
      0.0506142681451881,
  };
  double area = 0.0;
  for (int iu = 0; iu < 8; iu++) {
    for (int iv = 0; iv < 8; iv++) {
      double u = q[iu];
      double v = q[iv];
      double du[3], dv[3];
      for (int k = 0; k < 3; k++) {
        du[k] = -(1.0 - v) * p0(k) + (1.0 - v) * p1(k)
              + v * p2(k) - v * p3(k);
        dv[k] = -(1.0 - u) * p0(k) - u * p1(k)
              + u * p2(k) + (1.0 - u) * p3(k);
      }
      double cx = du[1] * dv[2] - du[2] * dv[1];
      double cy = du[2] * dv[0] - du[0] * dv[2];
      double cz = du[0] * dv[1] - du[1] * dv[0];
      area += w[iu] * w[iv] * std::sqrt(cx * cx + cy * cy + cz * cz);
    }
  }
  return area;
}

static double surface_descriptor_mesh_area(netgen::Mesh &mesh, int descriptor)
{
  double area = 0.0;
  for (int sei = 1; sei <= mesh.GetNSE(); sei++) {
    const auto element = mesh.SurfaceElement(netgen::SurfaceElementIndex::FromNr1(sei));
    if (element.GetIndex().Nr1() != descriptor)
      continue;
    int np = element.GetNP();
    if (np != 3 && np != 4)
      continue;
    const auto p0 = mesh.Point(netgen::PointIndex(element[0]));
    const auto p1 = mesh.Point(netgen::PointIndex(element[1]));
    const auto p2 = mesh.Point(netgen::PointIndex(element[2]));
    if (np == 3) {
      area += triangle_area(p0, p1, p2);
    } else {
      const auto p3 = mesh.Point(netgen::PointIndex(element[3]));
      area += bilinear_quad_area(p0, p1, p2, p3);
    }
  }
  return area;
}

ExportVolCommand::ExportVolCommand() {}
ExportVolCommand::~ExportVolCommand() {}

std::vector<std::string> ExportVolCommand::get_syntax()
{
  std::vector<std::string> syntax_list;
  syntax_list.push_back(
    "export vol <string:label='filename',help='<filename>'> "
    "[order <value:label='order',help='<1-5>'>] "
    "[overwrite] "
    "[periodic_name <string:label='periodic_name'>] "
    "[periodic_master <string:label='periodic_master'>] "
    "[periodic_slave <string:label='periodic_slave'>] "
    "[periodic_unit <string:label='periodic_unit'>] "
    "[periodic_patches <string:label='periodic_patches'>] "
    "[periodic_tx <value:label='periodic_tx'>] "
    "[periodic_ty <value:label='periodic_ty'>] "
    "[periodic_tz <value:label='periodic_tz'>] "
    "[periodic_match <value:label='periodic_match'>] "
    "[kelvin_tx <value:label='kelvin_tx'>] "
    "[kelvin_ty <value:label='kelvin_ty'>] "
    "[kelvin_tz <value:label='kelvin_tz'>] "
    "[kelvin_match <value:label='kelvin_match'>]"
  );
  return syntax_list;
}

std::vector<std::string> ExportVolCommand::get_syntax_help()
{
  std::vector<std::string> help;
  help.push_back(
    "export vol \"filename.vol\" [order {1|2|3|4|5}] [overwrite]\n"
  );
  return help;
}

std::vector<std::string> ExportVolCommand::get_help()
{
  std::vector<std::string> help;
  help.push_back(
    "Export mesh as Netgen .vol with high-order curving and labels.\n"
    "Uses C++ NetgenCurver (no Python subprocess, no .cub5 needed).\n"
    "Produces .vol file + companion .json with CAD reference values.\n"
    "Order 1-5 supported (default: 2). Requires Netgen.\n"
    "Named periodic metadata requires all fields explicitly: periodic_name kt,\n"
    "  periodic_master/slave <BC label>, periodic_unit m,\n"
    "  periodic_patches <masterBC:slaveBC,...>, periodic_tx/ty/tz <translation>,\n"
    "  periodic_match <positive geometric distance>. Fresh output required.\n"
    "Kelvin meshes (kelvin_int/kelvin_ext labels) require kelvin_tx/ty/tz\n"
    "  <Kelvin sphere offset> and kelvin_match <positive distance>; the export\n"
    "  writes the 'kelvin' periodic group. Fresh output required.\n"
    "One named periodic group per export (kt or kelvin).\n"
    "No inferred pairs or consumer repair; q3 geometry is preserved.\n"
    "\n"
  );
  return help;
}

// What the callbacks of an in-memory mesh query: re-curving refuses a changed model.
static std::string cubit_model_fingerprint()
{
  std::ostringstream text;
  for (const char *kind : {"volume", "surface", "curve"}) {
    const auto ids = CubitInterface::get_entities(kind);
    long long sum = 0;
    for (int id : ids) sum += id;
    text << kind << ' ' << ids.size() << ':' << sum << ' ';
  }
  text << "nodes " << CubitInterface::get_node_count()
       << " tets " << CubitInterface::get_tet_count()
       << " hexes " << CubitInterface::get_hex_count()
       << " wedges " << CubitInterface::get_wedge_count()
       << " pyramids " << CubitInterface::get_pyramid_count()
       << " tris " << CubitInterface::get_tri_count()
       << " quads " << CubitInterface::get_quad_count();
  return text.str();
}

bool ExportVolCommand::execute(CubitCommandData &data)
{
  VolExportRequest request;
  data.get_string("filename", request.filename);
  if (request.filename.empty()) {
    PRINT_ERROR("Filename required.\n");
    return false;
  }
  data.get_value("order", request.order);

  const bool periodic_fields[] = {
      data.get_string("periodic_name", request.periodic.name),
      data.get_string("periodic_master", request.periodic.master_label),
      data.get_string("periodic_slave", request.periodic.slave_label),
      data.get_string("periodic_unit", request.periodic.coordinate_unit),
      data.get_string("periodic_patches", request.periodic_patch_text),
      data.get_value("periodic_tx", request.periodic.translation[0]),
      data.get_value("periodic_ty", request.periodic.translation[1]),
      data.get_value("periodic_tz", request.periodic.translation[2]),
      data.get_value("periodic_match", request.periodic.matching_distance)};
  int periodic_present = 0;
  for (bool present : periodic_fields) periodic_present += present ? 1 : 0;
  request.has_periodic = periodic_present != 0;
  if (request.has_periodic && periodic_present != 9) {
    PRINT_ERROR("Named periodic export requires all nine explicit periodic fields.\n");
    return false;
  }
  const bool kelvin_fields[] = {
      data.get_value("kelvin_tx", request.kelvin_translation[0]),
      data.get_value("kelvin_ty", request.kelvin_translation[1]),
      data.get_value("kelvin_tz", request.kelvin_translation[2]),
      data.get_value("kelvin_match", request.kelvin_match)};
  int kelvin_present = 0;
  for (bool present : kelvin_fields) kelvin_present += present ? 1 : 0;
  request.has_kelvin = kelvin_present != 0;
  if (request.has_kelvin && kelvin_present != 4) {
    PRINT_ERROR("Kelvin export requires all four fields: kelvin_tx, kelvin_ty, kelvin_tz, kelvin_match.\n");
    return false;
  }
  return run(request, nullptr);
}

bool ExportVolCommand::run(const VolExportRequest &request, std::shared_ptr<netgen::Mesh> *live)
{
  // Suppress Cubit Learn Edition's 50k-cap ERROR during the export.
  // The exporter bypasses the cap and completes successfully, so
  // the ERROR line is misleading noise confusing users in logs.
  cubit_to_vol_export::ScopedLearnEditionFilter _lef_guard;

  auto t_start = std::chrono::high_resolution_clock::now();

  const std::string &filename = request.filename;
  const int order = request.order;
  if (order < 1 || order > 5) {
    PRINT_ERROR("Order must be 1-5.\n");
    return false;
  }

  cubit_to_vol_export::NamedPeriodicContract periodic = request.periodic;
  const std::string &periodic_patch_text = request.periodic_patch_text;
  const bool has_periodic = request.has_periodic;
  const std::array<double, 3> kelvin_translation = request.kelvin_translation;
  const double kelvin_match = request.kelvin_match;
  const bool has_kelvin = request.has_kelvin;
  if (has_kelvin && has_periodic) {
    PRINT_ERROR("One named periodic group per export: kt and kelvin pairs cannot be combined.\n");
    return false;
  }
  // Both named groups publish through the same fresh-output bundle.
  const bool named_output = has_periodic || has_kelvin;
  if (has_kelvin && !live) {
    try {
      cubit_to_vol_export::require_fresh_named_periodic_output(
          cubit_to_vol_export::named_periodic_output_paths(u8_string_to_path(filename)));
    } catch (const std::exception &error) {
      PRINT_ERROR("Kelvin export rejected: %s\n", error.what());
      return false;
    }
  }
  if (has_periodic) {
    try {
      if (!live)
        cubit_to_vol_export::require_fresh_named_periodic_output(
            cubit_to_vol_export::named_periodic_output_paths(u8_string_to_path(filename)));
      std::istringstream patches(periodic_patch_text);
      std::string patch;
      while (std::getline(patches, patch, ',')) {
        const auto colon = patch.find(':');
        if (colon == std::string::npos || colon == 0 || colon + 1 == patch.size() ||
            patch.find(':', colon + 1) != std::string::npos)
          throw std::runtime_error("Expected BC patch list master:slave,...");
        for (size_t j = 0; j < patch.size(); ++j)
          if (j != colon && (patch[j] < '0' || patch[j] > '9'))
            throw std::runtime_error("BC patch IDs must contain digits only");
        periodic.boundary_pairs.emplace_back(std::stoi(patch.substr(0, colon)),
                                             std::stoi(patch.substr(colon + 1)));
      }
      if (periodic_patch_text.empty() || periodic_patch_text.back() == ',')
        throw std::runtime_error("Empty BC patch entry");
    } catch (const std::exception &error) {
      PRINT_ERROR("Named periodic contract rejected: %s\n", error.what());
      return false;
    }
  }

  cubit_to_vol_export::NamedPeriodicResult periodic_result;

  // Extract linear mesh from Cubit (MeshData)
  MeshData md;
  md.require_volume_elements = true;  // .vol has no surface-only mode
  int build_order = (order >= 2) ? order : 2;
  if (!md.extract(build_order)) {
    PRINT_ERROR("Mesh extraction failed.\n");
    auto failed = md.get_netgen_curver();
    if (failed && !named_output && !live) {
      std::ofstream failure(u8_string_to_path(filename + ".json"));
      if (failure) failure << "{\"export_passed\":false,\"curving\":"
                           << failed->edge_diagnostics_json() << "}\n";
    }
    if (named_output && failed)
      PRINT_ERROR("Named periodic curving failed: %s\n", failed->edge_diagnostics_json().c_str());
    return false;
  }

  // Get the NetgenCurver that was used during extraction
  auto nc = md.get_netgen_curver();
  if (!nc) {
    PRINT_ERROR("NetgenCurver not available. Is Netgen installed?\n");
    return false;
  }

  auto ng_mesh = nc->get_ng_mesh();
  if (!ng_mesh) {
    PRINT_ERROR("Netgen mesh is null.\n");
    return false;
  }

  // ---- Set material labels (block name -> domain index) ----
  // Build domain_index map: block_id -> 1-based domain index
  std::map<int, int> block_to_domain;
  int dom_idx = 1;
  for (int bid : md.block_ids) {
    block_to_domain[bid] = dom_idx++;
  }
  int ndomains = (int)md.block_ids.size();

  // Set material names from Cubit block names
  for (int bid : md.block_ids) {
    std::string bname = CubitInterface::get_block_name(bid);
    if (bname.empty())
      bname = "volume_" + std::to_string(bid);
    ng_mesh->SetMaterial(block_to_domain[bid], bname);
  }

  // Update volume element domain indices (was hardcoded to 1)
  for (int ve_idx = 1; ve_idx <= ng_mesh->GetNE(); ve_idx++) {
    int elem_i = ve_idx - 1;  // 0-based into md.elements
    if (elem_i < (int)md.elements.size()) {
      auto it = block_to_domain.find(md.elements[elem_i].group_id);
      int di = (it != block_to_domain.end()) ? it->second : 1;
      ng_mesh->VolumeElement(netgen::ElementIndex::FromNr1(ve_idx))
          .SetIndex(netgen::VolumeRegionIndex::FromNr1(di));
    }
  }

  // Record user-authored sideset labels before same-material seam cleanup.
  // A same-domain internal surface is not a valid Netgen boundary, but an
  // explicitly labelled one must not disappear silently: it may represent a
  // user-intended interface or cut that needs distinct domain topology.
  std::map<int, std::string> surf_to_ssname;
  for (auto &sg : md.sidesets) {
    std::string ssname = sg.name.empty()
        ? "sideset_" + std::to_string(sg.id) : sg.name;
    surf_to_ssname[-sg.id] = ssname;
    for (int sid : CubitInterface::get_sideset_surfaces(sg.id))
      surf_to_ssname[sid] = ssname;
  }

  // ---- Fix FaceDescriptor DomainIn/DomainOut ----
  // NetgenCurver hardcodes DomainIn=1, DomainOut=0 for all faces.
  // Fix: use Cubit surface-volume topology to set correct domain indices.
  // This is required for periodic identification to work correctly.
  std::set<int> same_material_face_descriptors;
  std::map<int, std::string> labelled_same_material_surfaces;
  std::set<int> mesh_support_surface_ids;
  {
    // Build volume_id -> domain_index map
    std::map<int, int> vol_to_domain;
    for (int bid : md.block_ids) {
      std::vector<int> bvols = CubitInterface::get_block_volumes(bid);
      int di = block_to_domain[bid];
      for (int vid : bvols)
        vol_to_domain[vid] = di;
    }

    int nfd_fix = ng_mesh->GetNFD();
    for (int fi = 1; fi <= nfd_fix; fi++) {
      // BCProperty still holds original Cubit surface ID at this point
      int cubit_sid = ng_mesh->GetFaceDescriptor(netgen::FaceRegionIndex::FromNr1(fi)).BCProperty();
      if (cubit_sid > 0)
        mesh_support_surface_ids.insert(cubit_sid);
      if (cubit_sid < 0) {
        // NetgenCurver recovered this direct/free face's owner from adjacent
        // volume elements.  Validate that handoff instead of overwriting every
        // synthetic descriptor with domain 1.
        int domin = ng_mesh->GetFaceDescriptor(netgen::FaceRegionIndex::FromNr1(fi)).DomainIn();
        int domout = ng_mesh->GetFaceDescriptor(netgen::FaceRegionIndex::FromNr1(fi)).DomainOut();
        bool valid_domin = domin >= 1 && domin <= ndomains;
        bool valid_domout = domout == 0 || (domout >= 1 && domout <= ndomains);
        if (!valid_domin || !valid_domout || domin == domout) {
          PRINT_ERROR("Free sideset descriptor %d has invalid domain ownership "
                      "(%d -> %d) for a %d-domain mesh.\n",
                      -cubit_sid, domin, domout, ndomains);
          return false;
        }
        continue;
      }
      if (cubit_sid == 0) continue;

      // Get parent volumes of this surface
      std::vector<int> parent_vols =
          CubitInterface::get_relatives("surface", cubit_sid, "volume");

      int domin = 0, domout = 0;
      if (parent_vols.size() >= 1) {
        auto it = vol_to_domain.find(parent_vols[0]);
        domin = (it != vol_to_domain.end()) ? it->second : 0;
      }
      if (parent_vols.size() >= 2) {
        auto it = vol_to_domain.find(parent_vols[1]);
        domout = (it != vol_to_domain.end()) ? it->second : 0;
      }
      ng_mesh->GetFaceDescriptor(netgen::FaceRegionIndex::FromNr1(fi)).SetDomainIn(domin);
      ng_mesh->GetFaceDescriptor(netgen::FaceRegionIndex::FromNr1(fi)).SetDomainOut(domout);
      if (domin > 0 && domin == domout) {
        std::string explicit_label;
        auto sideset = surf_to_ssname.find(cubit_sid);
        if (sideset != surf_to_ssname.end())
          explicit_label = sideset->second;
        else
          explicit_label = explicit_surface_label(cubit_sid);
        if (explicit_label.empty())
          same_material_face_descriptors.insert(fi);
        else
          labelled_same_material_surfaces[cubit_sid] = explicit_label;
      }
    }
  }

  if (!labelled_same_material_surfaces.empty()) {
    for (const auto &[surface_id, label] : labelled_same_material_surfaces) {
      PRINT_ERROR("Same-material internal surface %d is explicitly labelled "
                  "'%s'. Netgen cannot encode it as a boundary with identical "
                  "DomainIn/DomainOut. Remove that surface from the sideset, "
                  "or model the intended interface/cut with distinct domains.\n",
                  surface_id, label.c_str());
    }
    return false;
  }

  // A Cubit webcut may split one material volume only to provide CAD curves
  // for a conforming mesh copy (the full-sphere Kelvin route does this at
  // z=0). Those conforming faces are not FE boundaries: both adjacent volume
  // elements carry the same material/domain index. Keeping their triangles in
  // the .vol creates invalid DomainIn == DomainOut surface rows and makes
  // NGSolve treat a meshing seam as physics. Remove only the surface elements;
  // shared volume nodes/edges remain untouched and conforming.
  if (!same_material_face_descriptors.empty()) {
    int removed_same_material_surfaces = 0;
    for (int sei = ng_mesh->GetNSE() - 1; sei >= 0; --sei) {
      const auto sei_index = netgen::SurfaceElementIndex::FromNr0(sei);
      auto surface_element = ng_mesh->SurfaceElements()[sei_index];
      if (same_material_face_descriptors.count(surface_element.GetIndex().Nr1())) {
        ng_mesh->SurfaceElements().DeleteElement(sei_index);
        removed_same_material_surfaces++;
      }
    }

    const int old_nfd = ng_mesh->GetNFD();
    std::vector<int> old_to_new(old_nfd + 1, 0);
    std::vector<netgen::FaceDescriptor> retained_descriptors;
    retained_descriptors.reserve(
        old_nfd - static_cast<int>(same_material_face_descriptors.size()));
    for (int fi = 1; fi <= old_nfd; ++fi) {
      if (same_material_face_descriptors.count(fi))
        continue;
      const auto &old_fd = ng_mesh->GetFaceDescriptor(netgen::FaceRegionIndex::FromNr1(fi));
      netgen::FaceDescriptor new_fd(
          old_fd.SurfNr(), old_fd.DomainIn(), old_fd.DomainOut(),
          old_fd.TLOSurface());
      new_fd.SetBCProperty(old_fd.BCProperty());
      new_fd.SetSurfColour(old_fd.SurfColour());
      new_fd.SetDomainInSingular(old_fd.DomainInSingular());
      new_fd.SetDomainOutSingular(old_fd.DomainOutSingular());
      retained_descriptors.push_back(new_fd);
      old_to_new[fi] = static_cast<int>(retained_descriptors.size());
    }
    for (auto &&surface_element : ng_mesh->SurfaceElements()) {
      const int old_index = surface_element.GetIndex().Nr1();
      if (old_index <= 0 || old_index > old_nfd || old_to_new[old_index] == 0) {
        PRINT_ERROR("Cannot remap retained surface element descriptor %d.\n",
                    old_index);
        return false;
      }
      surface_element.SetIndex(
          netgen::FaceRegionIndex::FromNr1(old_to_new[old_index]));
    }
    ng_mesh->ClearFaceDescriptors();
    for (const auto &descriptor : retained_descriptors)
      ng_mesh->AddFaceDescriptor(descriptor);

    PRINT_INFO("Removed %d same-material internal seam surface elements and "
               "%d unused face descriptors.\n",
               removed_same_material_surfaces,
               static_cast<int>(same_material_face_descriptors.size()));
    same_material_face_descriptors.clear();
  }

  // ---- Set boundary labels (sideset name -> bc number) ----
  // FaceDescriptors have BCProperty set to Cubit surface ID, which can be
  // larger than the number of FaceDescriptors.  NGSolve expects BCProperty
  // to be a contiguous 1-based index into the bcnames array.  Remap here.
  int nfd = ng_mesh->GetNFD();

  // Save original Cubit surface IDs before remapping (for companion JSON)
  std::vector<int> orig_surf_ids(nfd);
  for (int fi = 1; fi <= nfd; fi++)
    orig_surf_ids[fi - 1] = ng_mesh->GetFaceDescriptor(netgen::FaceRegionIndex::FromNr1(fi)).BCProperty();

  // Set BCName and remap BCProperty to contiguous 1-based index
  ng_mesh->SetNBCNames(nfd);
  for (int fi = 1; fi <= nfd; fi++) {
    int bc_prop = orig_surf_ids[fi - 1];

    // Determine label
    std::string label;
    auto it = surf_to_ssname.find(bc_prop);
    if (it != surf_to_ssname.end()) {
      label = it->second;
    } else {
      label = CubitInterface::get_entity_name("surface", bc_prop);
      if (label.empty())
        label = "surface_" + std::to_string(bc_prop);
    }
    for (auto &ch : label) if (ch == ' ') ch = '_';

    // Remap BCProperty: Cubit surface ID -> contiguous index fi
    ng_mesh->GetFaceDescriptor(netgen::FaceRegionIndex::FromNr1(fi)).SetBCProperty(fi);
    ng_mesh->SetBCName(fi - 1, label);
  }

  // For order=1: reset curving
  if (order == 1) {
    ng_mesh->BuildCurvedElements(1);
  }

  // ---- Kelvin auto-detection ----
  // If blocks named "air" and "kelvin" exist, automatically label the
  // shared air|kelvin interface as "kelvin_int" and the outer kelvin
  // boundary (DomainOut=0) as "kelvin_ext".  This saves the user from
  // having to add manual sidesets in the .jou file.
  //
  // IMPORTANT (2026-04-25): a face descriptor that ALREADY has a
  // sideset-derived name MUST NOT be overridden here.  In reduction
  // mode, `add_kelvin_cubit(reduction=...)` creates sym_bn=0_<axis>
  // / sym_ht=0_<axis> sidesets on the Kelvin's flat cut faces, which
  // also border domain 0 -- matching the "kelvin_ext" condition.
  // Overwriting them corrupts the label map.  `surf_to_ssname` (set
  // up earlier from md.sidesets) is the authoritative source for
  // user-chosen labels; the auto-detect only fills gaps.
  {
    int air_dom = 0, kelvin_dom = 0;
    for (int bid : md.block_ids) {
      std::string bname = CubitInterface::get_block_name(bid);
      for (auto &ch : bname) ch = (char)tolower((unsigned char)ch);
      if (bname == "air")    air_dom = block_to_domain[bid];
      if (bname == "kelvin") kelvin_dom = block_to_domain[bid];
    }

    if (air_dom > 0 && kelvin_dom > 0) {
      int n_int = 0, n_ext = 0, n_skipped = 0;
      for (int fi = 1; fi <= ng_mesh->GetNFD(); fi++) {
        int di = ng_mesh->GetFaceDescriptor(netgen::FaceRegionIndex::FromNr1(fi)).DomainIn();
        int dout = ng_mesh->GetFaceDescriptor(netgen::FaceRegionIndex::FromNr1(fi)).DomainOut();
        // Skip faces the user has already named via a sideset.  Lookup
        // by original Cubit surface ID (orig_surf_ids[fi-1]) matches
        // the sideset's surface-list exactly; the remapped BCProperty
        // (= fi) is the contiguous index, not the Cubit surface ID.
        if (surf_to_ssname.count(orig_surf_ids[fi - 1])) {
          n_skipped++;
          continue;
        }
        // Shared face between air and kelvin = kelvin_int
        if ((di == air_dom && dout == kelvin_dom) ||
            (di == kelvin_dom && dout == air_dom)) {
          ng_mesh->SetBCName(fi - 1, "kelvin_int");
          n_int++;
        }
        // Outer boundary of kelvin shell (one side is kelvin, other is 0)
        else if ((di == kelvin_dom && dout == 0) ||
                 (di == 0 && dout == kelvin_dom)) {
          ng_mesh->SetBCName(fi - 1, "kelvin_ext");
          n_ext++;
        }
      }
      if (n_int > 0 && n_ext > 0) {
        PRINT_INFO("Kelvin auto-detect: %d inner + %d outer face "
                   "descriptors (%d skipped due to user sideset)\n",
                   n_int, n_ext, n_skipped);
      } else if (n_int > 0 || n_ext > 0) {
        PRINT_WARNING("Kelvin auto-detect: found %d inner, %d outer "
                      "(need both > 0; %d skipped due to user sideset)\n",
                      n_int, n_ext, n_skipped);
      }
    }
  }

  // ---- Cubit curve sidesets -> Netgen CD2 (BBND edge) names ----
  // CD2 names use the 1-based EdgeDescriptor index, not the Cubit curve ID
  // or the adjacent face descriptor. Each represented curve owns one label.
  {
    const auto &curve_edges = nc->get_cubit_curve_to_ng_edge();
    std::map<int, std::string> edge_labels;
    for (const auto &sg : md.sidesets) {
      auto curves = CubitInterface::get_sideset_curves(sg.id);
      if (curves.empty()) continue;
      std::string label = sg.name.empty() ? "sideset_" + std::to_string(sg.id) : sg.name;
      for (auto &ch : label) if (ch == ' ') ch = '_';
      if (label == "default") {
        PRINT_ERROR("Curve sideset %d uses reserved BBND name 'default'.\n", sg.id);
        return false;
      }
      for (int cid : curves) {
        auto edge = curve_edges.find(cid);
        if (edge == curve_edges.end()) {
          PRINT_ERROR("Curve sideset %d '%s': curve %d cannot be represented "
                      "as a complete BBND edge in the exported volume mesh.\n",
                      sg.id, label.c_str(), cid);
          return false;
        }
        auto previous = edge_labels.find(edge->second);
        if (previous != edge_labels.end() && previous->second != label) {
          PRINT_ERROR("Curve %d has conflicting BBND sideset names '%s' and '%s'.\n",
                      cid, previous->second.c_str(), label.c_str());
          return false;
        }
        edge_labels[edge->second] = label;
      }
    }
    for (const auto &entry : edge_labels)
      ng_mesh->SetCD2Name(entry.first, entry.second);
  }

  // ---- Cubit nodesets -> Netgen CD3 (BBBND vertex) names ----
  // Each Cubit nodeset (e.g. "GND" at the Kelvin sphere centre, used as
  // the Dirichlet anchor for Omega-reduced FEM) becomes a named Element0d
  // in the Netgen mesh's `pointelements` array.  NGSolve exposes these
  // via `mesh.GetBBBoundaries()`, and calc_accel_magnet.py looks up the
  // vertex by that name to place the Dirichlet BC.
  {
    const auto &nid_to_pi = nc->get_cubit_nid_to_ng_pi();
    int cd3_names_written = 0;
    int point_elements_written = 0;
    for (const auto &ns : md.nodesets) {
      if (ns.name.empty()) continue;

      // Collect the mesh-node ids actually in this nodeset.  The
      // MeshExportInterface path (ns.node_ids) works for node-based
      // nodesets but returns empty for nodesets defined by geometric
      // entities (vertex/curve/surface/volume).  Fall back to
      // CubitInterface::get_nodeset_children and expand any vertex
      // entries via parse_cubit_list("node", "in vertex N").
      std::vector<int> mesh_nodes(ns.node_ids.begin(), ns.node_ids.end());
      if (mesh_nodes.empty()) {
        std::vector<int> nl, vl, sl, cl, xl;
        CubitInterface::get_nodeset_children(ns.id, nl, vl, sl, cl, xl);
        mesh_nodes.insert(mesh_nodes.end(), nl.begin(), nl.end());
        // Expand vertices to mesh nodes
        for (int vid : xl) {
          std::string spec = "in vertex " + std::to_string(vid);
          auto nids = CubitInterface::parse_cubit_list("node", spec);
          mesh_nodes.insert(mesh_nodes.end(), nids.begin(), nids.end());
        }
        // Expand curves/surfaces/volumes similarly (unusual for GND-style
        // anchors but cheap).
        for (int eid : cl) {
          auto nids = CubitInterface::parse_cubit_list(
              "node", "in curve " + std::to_string(eid));
          mesh_nodes.insert(mesh_nodes.end(), nids.begin(), nids.end());
        }
        for (int eid : sl) {
          auto nids = CubitInterface::parse_cubit_list(
              "node", "in surface " + std::to_string(eid));
          mesh_nodes.insert(mesh_nodes.end(), nids.begin(), nids.end());
        }
        for (int eid : vl) {
          auto nids = CubitInterface::parse_cubit_list(
              "node", "in volume " + std::to_string(eid));
          mesh_nodes.insert(mesh_nodes.end(), nids.begin(), nids.end());
        }

        // Free-floating vertex: not merged into any meshed volume, so it
        // has no mesh node of its own.  Fall back to nearest-neighbor in
        // md.nodes.  This matches the use case where add_kelvin_cubit
        // creates a bare `create vertex X Y Z` and puts it in the "GND"
        // nodeset -- the nearest mesh node to the Kelvin centre is the
        // right anchor for the Omega-reduced Dirichlet.
        if (mesh_nodes.empty() && !xl.empty()) {
          // Restrict the search to ORIGINAL Cubit nodes (indices
          // 0 .. num_original_nodes).  High-order (curved) nodes
          // generated by NetgenCurver after AddPoint() are NOT in
          // nid_to_pi and would silently be dropped by the Element0d
          // append below.
          int n_orig = md.num_original_nodes > 0
                       ? md.num_original_nodes
                       : (int)md.nodes.size();
          for (int vid : xl) {
            auto xyz = CubitInterface::get_center_point("vertex", vid);
            int best_nid = -1;
            double best_d2 = 1e300;
            for (int i = 0; i < n_orig; i++) {
              const auto &nd = md.nodes[i];
              double dx = nd.x - xyz[0];
              double dy = nd.y - xyz[1];
              double dz = nd.z - xyz[2];
              double d2 = dx*dx + dy*dy + dz*dz;
              if (d2 < best_d2) { best_d2 = d2; best_nid = nd.id; }
            }
            if (best_nid > 0) {
              mesh_nodes.push_back(best_nid);
              PRINT_INFO("Nodeset export: vertex %d (%.4g,%.4g,%.4g) -> "
                         "nearest original mesh node %d "
                         "(dist=%.3e m, searched %d nodes)\n",
                         vid, xyz[0], xyz[1], xyz[2],
                         best_nid, std::sqrt(best_d2), n_orig);
            }
          }
        }
      }

      if (mesh_nodes.empty()) {
        PRINT_WARNING("Nodeset export: nodeset '%s' (id=%d) has no mesh nodes; "
                      "CD3 name will be unused.\n", ns.name.c_str(), ns.id);
        continue;
      }

      // AddCD3Name returns the 0-based slot it assigned.  Element0d's
      // `index` field is what Netgen's .vol writer emits as the "1 GND"
      // line in the cd3names section (1-based), so bump to match.
      int cd3_slot = ng_mesh->AddCD3Name(ns.name);
      int cd3_idx = cd3_slot + 1;  // 1-based index for Element0d
      cd3_names_written++;
      PRINT_INFO("Nodeset export: AddCD3Name('%s') -> slot %d, element index %d\n",
                 ns.name.c_str(), cd3_slot, cd3_idx);
      for (int nid : mesh_nodes) {
        auto it = nid_to_pi.find(nid);
        if (it == nid_to_pi.end()) continue;
        const auto pi = netgen::PointIndex::FromNr1(it->second);
        ng_mesh->pointelements.Append(netgen::Element0d(
            pi, netgen::VertexRegionIndex::FromNr1(cd3_idx)));
        point_elements_written++;
      }
    }
    if (cd3_names_written > 0) {
      PRINT_INFO("Nodeset export: %d CD3 name(s), %d Element0d point(s)\n",
                 cd3_names_written, point_elements_written);
    }
  }

  // Explicit metadata only, after curving and final descriptor/BC labeling.
  // No centroid/tolerance inference and no solver-side repair fallback.
  // A Kelvin mesh is exported only with its pairs, so the .vol is self-contained.
  bool kelvin_labelled = false;
  for (int fi = 1; fi <= ng_mesh->GetNFD(); fi++) {
    const auto &name = ng_mesh->GetBCName(fi - 1);
    if (name == "kelvin_int" || name == "kelvin_ext") kelvin_labelled = true;
  }
  if (kelvin_labelled && !has_kelvin) {
    PRINT_ERROR("Kelvin boundaries (kelvin_int/kelvin_ext) need their pairs: pass kelvin_tx, "
                "kelvin_ty, kelvin_tz (the Kelvin sphere offset; add_kelvin_cubit returns it "
                "as 'center') and kelvin_match.\n");
    return false;
  }
  if (has_kelvin && !kelvin_labelled) {
    PRINT_ERROR("kelvin_* fields given, but the mesh has no kelvin_int/kelvin_ext boundaries.\n");
    return false;
  }
  if (has_kelvin) {
    try {
      periodic = cubit_to_vol_export::kelvin_periodic_contract(
          *ng_mesh, kelvin_translation, kelvin_match);
    } catch (const std::exception &error) {
      PRINT_ERROR("Kelvin pairing rejected: %s\n", error.what());
      return false;
    }
  }
  if (named_output) {
    try {
      periodic_result = cubit_to_vol_export::apply_named_periodic_identification(*ng_mesh, periodic);
    } catch (const std::exception &error) {
      PRINT_ERROR("Named periodic export rejected: %s\n", error.what());
      return false;
    }
    PRINT_INFO("Periodic group '%s': %d point pairs, max distance %.3e m\n",
               periodic.name.c_str(), (int)periodic_result.point_pairs.size(),
               periodic_result.maximum_distance);
  }

  // In-memory mesh: keep the live geometry, owning the curver its callbacks
  // read, and write nothing.
  if (live) {
    auto geometry = std::dynamic_pointer_cast<cubit_to_vol_export::LiveCubitGeometry>(
        ng_mesh->GetGeometry());
    if (!geometry) {
      PRINT_ERROR("The curved mesh carries no live Cubit geometry.\n");
      return false;
    }
    geometry->owner = nc;
    geometry->periodic_triangles = periodic_result.triangles;
    geometry->model_fingerprint = cubit_model_fingerprint;
    geometry->built_fingerprint = cubit_model_fingerprint();
    nc->release_mesh();
    *live = ng_mesh;
    PRINT_INFO("Built in-memory Netgen mesh (order %d): %d nodes, %d elements\n",
               order, (int)ng_mesh->GetNP(), (int)ng_mesh->GetNE());
    return true;
  }

  // Detach CallbackGeometry — it has function pointers, not serializable.
  // Curving data is preserved in CurvedElements (written as "curvedelements"
  // section in .vol text format, upstream Netgen master feature).
  ng_mesh->SetGeometry(nullptr);

  // Save .vol text format (includes curvedelements section for order >= 2).
  // NGSolve reads this with Mesh("file.vol") — no STEP/Cubit needed.
  // Pass as std::filesystem::path built from UTF-8 so Unicode / Japanese
  // directories work on Windows (cp932 narrow API path would otherwise
  // throw "No mapping for the Unicode character").
  if (!named_output) {
    cubit_to_vol_export::save_native_mesh(*ng_mesh, u8_string_to_path(filename));
  }

  int ne = ng_mesh->GetNE();
  int np = ng_mesh->GetNP();


  // The named export is incomplete without its companion mapping record.
  auto write_companion = [&](std::ostream &jf) {
      jf << std::setprecision(16);
      jf << "{\n";
      jf << "  \"export_passed\": true,\n";
      if (named_output) {
        jf << "  \"periodic_identification\": {\"name\":\"" << periodic.name
           << "\",\"nr\":1,\"type\":"
           << static_cast<int>(netgen::Identifications::PERIODIC) << ","
           << "\"pair_count\":" << periodic_result.point_pairs.size()
           << ",\"max_distance\":" << std::setprecision(17) << periodic_result.maximum_distance
           << ",\"descriptor_mapping\":[";
        bool first_mapping = true;
        for (const auto &mapping : periodic_result.descriptor_to_boundary) {
          if (!first_mapping) jf << ",";
          first_mapping = false;
          jf << "{\"fd\":" << mapping.first << ",\"bc\":" << mapping.second
             << ",\"cubit_surface\":" << orig_surf_ids.at(mapping.first - 1) << "}";
        }
        jf << "],\"curved_transfers\":[";
        bool first_transfer = true;
        for (const auto &transfer : periodic_result.curved_transfers) {
          if (!first_transfer) jf << ",";
          first_transfer = false;
          jf << "{\"entity\":\"" << (transfer.face ? "face" : "edge")
             << "\",\"master_topology_index\":" << transfer.master
             << ",\"slave_topology_index\":" << transfer.slave
             << ",\"mode\":" << transfer.mode << ",\"sign\":" << transfer.sign
             << ",\"original\":[";
          for (int axis = 0; axis < 3; ++axis) {
            if (axis) jf << ",";
            jf << transfer.original[axis];
          }
          jf << "],\"expected\":[";
          for (int axis = 0; axis < 3; ++axis) {
            if (axis) jf << ",";
            jf << transfer.expected[axis];
          }
          jf << "]}";
        }
        jf << "]},\n";
      }
      jf << "  \"curving\": " << nc->edge_diagnostics_json() << ",\n";

      // Materials (per-block volume)
      // MeshExportInterface may return a default block (all elements) whose ID
      // does not exist as a user-defined Cubit block. Skip such blocks to avoid
      // "No block with ID N was found" errors from parse_cubit_list.
      jf << "  \"materials\": {";
      bool first = true;
      std::vector<std::string> mesh_only_materials;
      std::vector<int> cubit_block_ids = CubitInterface::parse_cubit_list("block", "all");
      std::set<int> cubit_block_set(cubit_block_ids.begin(), cubit_block_ids.end());
      for (int bid : md.block_ids) {
        if (cubit_block_set.find(bid) == cubit_block_set.end())
          continue;  // skip MeshExportInterface default block
        std::string bname = CubitInterface::get_block_name(bid);
        if (bname.empty()) bname = "volume_" + std::to_string(bid);
        // Get volumes in this block
        double total_vol = 0.0;
        std::vector<int> vols_in_block = CubitInterface::parse_cubit_list(
            "volume", "in block " + std::to_string(bid));
        if (vols_in_block.empty()) {
          mesh_only_materials.push_back(bname);
          continue;
        }
        for (int vid : vols_in_block) {
          // 2026-05-25: Cubit 2025.12 removed GeometryQueryTool::instance().
          // Use CubitInterface::get_volume_volume() instead -- returns the
          // same CAD volume measure RefVolume::measure() does, and works
          // identically on a missing vid (returns 0).
          total_vol += CubitInterface::get_volume_volume(vid);
        }
        if (!first) jf << ",";
        jf << "\n    \"" << bname << "\": " << std::scientific << total_vol;
        first = false;
      }
      jf << "\n  },\n";

      // Sculpt and imported Exodus meshes can have material blocks without
      // owning CAD volumes.  Distinguish that state from a real zero-volume
      // CAD body so downstream gates do not silently accept a false reference.
      jf << "  \"mesh_only_materials\": [";
      for (size_t i = 0; i < mesh_only_materials.size(); i++) {
        if (i > 0) jf << ", ";
        jf << "\"" << mesh_only_materials[i] << "\"";
      }
      jf << "],\n";

      // Boundaries: sum CAD area per unique bcname.  Multiple face
      // descriptors can share a bcname (e.g. sym_ht=0_y for both the
      // air's and the Kelvin's y=0 cut faces in 1/8 reduction); the
      // JSON should report the TOTAL area, not just the last entry's,
      // otherwise key collision in JSON drops the others.
      jf << "  \"boundaries\": {";
      int nfd_json = ng_mesh->GetNFD();
      std::map<std::string, double> bname_to_area;
      for (int fi = 1; fi <= nfd_json; fi++) {
        int cubit_surf_id = orig_surf_ids[fi - 1];
        std::string bname = ng_mesh->GetBCName(fi - 1);
        // 2026-05-25 Cubit 2025.12: GeometryQueryTool::instance() removed.
        double cad_area = 0.0;
        if (cubit_surf_id < 0) {
          // A Skin-generated sideset can overlap an existing exterior
          // sideset.  Surface connectivity is deduplicated during export, so
          // the raw Cubit sideset area would count faces that are not present
          // under this descriptor.  Measure the actual exported linear map.
          cad_area = surface_descriptor_mesh_area(*ng_mesh, fi);
        } else {
          cad_area = CubitInterface::get_surface_area(cubit_surf_id);
        }
        bname_to_area[bname] += cad_area;
      }
      first = true;
      for (auto& kv : bname_to_area) {
        if (!first) jf << ",";
        jf << "\n    \"" << kv.first << "\": " << std::scientific << kv.second;
        first = false;
      }
      jf << "\n  },\n";

      // Edges (per-curve length, curves represented by Netgen segments).
      // Imported STL/Sculpt geometry can retain CAD curves even when no mesh
      // edge belongs to them.  Do not publish those curves as BBND reference
      // data: the checker would otherwise compare a non-existent Netgen edge
      // set against stale CAD topology.
      jf << "  \"edges\": {";
      first = true;
      // Include curves on same-material webcut surfaces too. Those surfaces
      // are removed as FE boundaries, but their curves still partition and
      // support the conforming volume mesh and therefore remain as BBND
      // segments in the .vol file.
      const std::set<int> &exported_surface_ids = mesh_support_surface_ids;
      std::vector<int> curve_ids = CubitInterface::parse_cubit_list("curve", "all");
      for (int cid : curve_ids) {
        std::vector<int> parent_surfs = CubitInterface::parse_cubit_list(
            "surface", "in curve " + std::to_string(cid));
        if (parent_surfs.size() < 2) continue;
        if (exported_surface_ids.count(parent_surfs[0]) == 0 ||
            exported_surface_ids.count(parent_surfs[1]) == 0)
          continue;
        std::vector<int> edge_ids = CubitInterface::parse_cubit_list(
            "edge", "in curve " + std::to_string(cid));
        if (edge_ids.empty()) continue;
        // 2026-05-25 Cubit 2025.12: GeometryQueryTool::instance() removed.
        double cad_len = CubitInterface::get_curve_length(cid);
        std::string ename = "curve_" + std::to_string(cid);
        if (!first) jf << ",";
        jf << "\n    \"" << ename << "\": " << std::scientific << cad_len;
        first = false;
      }
      jf << "\n  },\n";

      auto t_end = std::chrono::high_resolution_clock::now();
      double t_sec = std::chrono::duration<double>(t_end - t_start).count();

      jf << "  \"n_elements\": " << ne << ",\n";
      jf << "  \"n_points\": " << np << ",\n";
      jf << "  \"order\": " << order << ",\n";
      jf << "  \"export_time_s\": " << std::fixed << std::setprecision(3) << t_sec << "\n";
      jf << "}\n";
  };
  if (named_output) {
    try {
      cubit_to_vol_export::save_named_periodic_bundle(
          *ng_mesh, u8_string_to_path(filename), write_companion);
    } catch (const std::exception &error) {
      PRINT_ERROR("Named periodic bundle publication failed: %s\n", error.what());
      PRINT_ERROR("%s\n", cubit_to_vol_export::named_periodic_retry_guidance(
          cubit_to_vol_export::named_periodic_output_paths(u8_string_to_path(filename))).c_str());
      return false;
    }
  } else {
    std::ofstream jf(u8_string_to_path(filename + ".json"));
    if (jf.is_open()) { write_companion(jf); jf.close(); }
  }

  PRINT_INFO("Exported Netgen Vol (order %d): %s (%d nodes, %d elements)\n",
             order, filename.c_str(), np, ne);

  return true;
}
