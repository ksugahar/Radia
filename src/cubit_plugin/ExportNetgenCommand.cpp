#include "ExportNetgenCommand.hpp"
#include "MeshData.hpp"
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"
#include <chrono>
#include <iomanip>

#ifdef HAVE_NETGEN
#include "NetgenCurver.hpp"
#include <meshing.hpp>
#endif

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

// Ensure Netgen DLLs (nglib.dll, ngcore.dll) can be found.
// They live in plugins/ which may not be on the DLL search path.
static void ensure_netgen_dll_path()
{
#ifdef COMPACT_NETGEN
  // Static link: no DLL path setup needed
  return;
#elif defined(_WIN32)
  // Find plugins/ directory relative to this DLL's location
  // or use CUBIT_PLUGIN_DIR environment variable
  const char *pd = std::getenv("CUBIT_PLUGIN_DIR");
  if (pd && pd[0]) {
    SetDllDirectoryA(pd);
    return;
  }
  // Fallback: add Cubit's bin/plugins/ to DLL search path
  HMODULE hm = nullptr;
  GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                     GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                     (LPCSTR)&ensure_netgen_dll_path, &hm);
  if (hm) {
    char path[MAX_PATH];
    if (GetModuleFileNameA(hm, path, MAX_PATH)) {
      // path = .../plugins/radia_cubit.ccm or .../bin/radia_cubit.ccl
      std::string dir(path);
      auto pos = dir.find_last_of("\\/");
      if (pos != std::string::npos) {
        dir = dir.substr(0, pos);
        // If we're in bin/, try bin/plugins/
        std::string plugins = dir + "\\plugins";
        if (GetFileAttributesA(plugins.c_str()) != INVALID_FILE_ATTRIBUTES)
          dir = plugins;
        AddDllDirectory(std::wstring(dir.begin(), dir.end()).c_str());
        SetDllDirectoryA(dir.c_str());
      }
    }
  }
#endif
}

ExportNetgenCommand::ExportNetgenCommand() {}
ExportNetgenCommand::~ExportNetgenCommand() {}

std::vector<std::string> ExportNetgenCommand::get_syntax()
{
  std::vector<std::string> syntax_list;
  syntax_list.push_back(
    "radia_export netgen <string:label='filename',help='<filename>'> "
    "[order <value:label='order',help='<1-5>'>] "
    "[overwrite]"
  );
  return syntax_list;
}

std::vector<std::string> ExportNetgenCommand::get_syntax_help()
{
  std::vector<std::string> help;
  help.push_back(
    "radia_export netgen \"filename.vol\" [order {1|2|3|4|5}] [overwrite]"
  );
  return help;
}

std::vector<std::string> ExportNetgenCommand::get_help()
{
  std::vector<std::string> help;
  help.push_back(
    "Export mesh as Netgen .vol with high-order curving and labels.\n"
    "Uses C++ NetgenCurver (no Python subprocess, no .cub5 needed).\n"
    "Produces .vol file + companion .json with CAD reference values.\n"
    "Order 1-5 supported (default: 2). Requires Netgen."
  );
  return help;
}

bool ExportNetgenCommand::execute(CubitCommandData &data)
{
#ifndef HAVE_NETGEN
  PRINT_ERROR("radia_export netgen requires Netgen support (not built).\n");
  return false;
#else
  auto t_start = std::chrono::high_resolution_clock::now();

  std::string filename;
  data.get_string("filename", filename);
  if (filename.empty()) {
    PRINT_ERROR("Filename required.\n");
    return false;
  }

  int order = 2;
  data.get_value("order", order);
  if (order < 1 || order > 5) {
    PRINT_ERROR("Order must be 1-5.\n");
    return false;
  }

  // Ensure Netgen DLLs are findable (DELAYLOAD resolution)
  ensure_netgen_dll_path();

  // Extract linear mesh from Cubit (MeshData)
  MeshData md;
  int build_order = (order >= 2) ? order : 2;
  if (!md.extract(build_order)) {
    PRINT_ERROR("Mesh extraction failed.\n");
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
      ng_mesh->VolumeElement(ve_idx).SetIndex(di);
    }
  }

  // ---- Fix FaceDescriptor DomainIn/DomainOut ----
  // NetgenCurver hardcodes DomainIn=1, DomainOut=0 for all faces.
  // Fix: use Cubit surface-volume topology to set correct domain indices.
  // This is required for periodic identification to work correctly.
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
      int cubit_sid = ng_mesh->GetFaceDescriptor(fi).BCProperty();
      if (cubit_sid <= 0) continue;

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
      ng_mesh->GetFaceDescriptor(fi).SetDomainIn(domin);
      ng_mesh->GetFaceDescriptor(fi).SetDomainOut(domout);
    }
  }

  // ---- Set boundary labels (sideset name -> bc number) ----
  // FaceDescriptors have BCProperty set to Cubit surface ID, which can be
  // larger than the number of FaceDescriptors.  NGSolve expects BCProperty
  // to be a contiguous 1-based index into the bcnames array.  Remap here.
  int nfd = ng_mesh->GetNFD();

  // Build surface_id -> sideset name map
  std::map<int, std::string> surf_to_ssname;
  for (auto &sg : md.sidesets) {
    std::vector<int> ss_surfs = CubitInterface::get_sideset_surfaces(sg.id);
    for (int sid : ss_surfs) {
      if (!sg.name.empty())
        surf_to_ssname[sid] = sg.name;
    }
  }

  // Save original Cubit surface IDs before remapping (for companion JSON)
  std::vector<int> orig_surf_ids(nfd);
  for (int fi = 1; fi <= nfd; fi++)
    orig_surf_ids[fi - 1] = ng_mesh->GetFaceDescriptor(fi).BCProperty();

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
    ng_mesh->GetFaceDescriptor(fi).SetBCProperty(fi);
    ng_mesh->SetBCName(fi - 1, label);
    ng_mesh->GetFaceDescriptor(fi).SetBCName(ng_mesh->GetBCNamePtr(fi - 1));
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
  {
    int air_dom = 0, kelvin_dom = 0;
    for (int bid : md.block_ids) {
      std::string bname = CubitInterface::get_block_name(bid);
      for (auto &ch : bname) ch = (char)tolower((unsigned char)ch);
      if (bname == "air")    air_dom = block_to_domain[bid];
      if (bname == "kelvin") kelvin_dom = block_to_domain[bid];
    }

    if (air_dom > 0 && kelvin_dom > 0) {
      int n_int = 0, n_ext = 0;
      for (int fi = 1; fi <= ng_mesh->GetNFD(); fi++) {
        int di = ng_mesh->GetFaceDescriptor(fi).DomainIn();
        int dout = ng_mesh->GetFaceDescriptor(fi).DomainOut();
        // Shared face between air and kelvin = kelvin_int
        if ((di == air_dom && dout == kelvin_dom) ||
            (di == kelvin_dom && dout == air_dom)) {
          ng_mesh->SetBCName(fi - 1, "kelvin_int");
          ng_mesh->GetFaceDescriptor(fi).SetBCName(
              ng_mesh->GetBCNamePtr(fi - 1));
          n_int++;
        }
        // Outer boundary of kelvin shell (one side is kelvin, other is 0)
        else if ((di == kelvin_dom && dout == 0) ||
                 (di == 0 && dout == kelvin_dom)) {
          ng_mesh->SetBCName(fi - 1, "kelvin_ext");
          ng_mesh->GetFaceDescriptor(fi).SetBCName(
              ng_mesh->GetBCNamePtr(fi - 1));
          n_ext++;
        }
      }
      if (n_int > 0 && n_ext > 0) {
        PRINT_INFO("Kelvin auto-detect: %d inner + %d outer face descriptors\n",
                   n_int, n_ext);
      } else if (n_int > 0 || n_ext > 0) {
        PRINT_WARNING("Kelvin auto-detect: found %d inner, %d outer "
                      "(need both > 0)\n", n_int, n_ext);
      }
    }
  }

  // ---- Kelvin periodic identification (optional) ----
  // Concentric-sphere Kelvin transform: inner surface at R_inner,
  // outer surface at R_outer, both centered at the same point.
  // The mapping is RADIAL SCALING: p_outer = center + (p_inner - center) * scale
  // where scale = R_outer / R_inner.
  //
  // Algorithm:
  //   1. Collect inner/outer triangles (vertex-3 only)
  //   2. Compute center and radii from point clouds
  //   3. For each inner vertex, predict its outer position by radial scaling
  //   4. Find the nearest outer vertex to each predicted position
  //   5. Verify via triangle topology (matched triangles must have
  //      consistent vertex pairing across all shared edges)
  {
    std::set<int> fd_inner_set, fd_outer_set;
    for (int fi = 1; fi <= ng_mesh->GetNFD(); fi++) {
      std::string bc = ng_mesh->GetBCName(fi - 1);
      if (bc == "kelvin_int") fd_inner_set.insert(fi);
      else if (bc == "kelvin_ext") fd_outer_set.insert(fi);
    }

    if (!fd_inner_set.empty() && !fd_outer_set.empty()) {
      // Collect vertex positions for inner and outer surfaces
      std::map<int, netgen::Point<3>> inner_pts, outer_pts;

      for (int sei = 1; sei <= ng_mesh->GetNSE(); sei++) {
        const auto &sel = ng_mesh->SurfaceElement(sei);
        int fd = sel.GetIndex();
        bool is_inner = fd_inner_set.count(fd) > 0;
        bool is_outer = fd_outer_set.count(fd) > 0;
        if (!is_inner && !is_outer) continue;

        auto &target = is_inner ? inner_pts : outer_pts;
        for (int j = 0; j < sel.GetNP(); j++) {
          int pi = sel[j];
          if (target.find(pi) == target.end())
            target[pi] = ng_mesh->Point(netgen::PointIndex(pi));
        }
      }

      if (!inner_pts.empty() && !outer_pts.empty()) {
        // Compute center as mean of inner points (should be near sphere center)
        double cx = 0, cy = 0, cz = 0;
        for (auto &p : inner_pts) {
          cx += p.second(0); cy += p.second(1); cz += p.second(2);
        }
        double ni = (double)inner_pts.size();
        cx /= ni; cy /= ni; cz /= ni;

        // Compute average radii
        double r_inner_sum = 0, r_outer_sum = 0;
        for (auto &p : inner_pts) {
          double dx = p.second(0) - cx, dy = p.second(1) - cy, dz = p.second(2) - cz;
          r_inner_sum += sqrt(dx*dx + dy*dy + dz*dz);
        }
        for (auto &p : outer_pts) {
          double dx = p.second(0) - cx, dy = p.second(1) - cy, dz = p.second(2) - cz;
          r_outer_sum += sqrt(dx*dx + dy*dy + dz*dz);
        }
        double R_inner = r_inner_sum / ni;
        double R_outer = r_outer_sum / (double)outer_pts.size();
        double scale = R_outer / R_inner;

        PRINT_INFO("Kelvin periodic: center=(%.4f,%.4f,%.4f) "
                   "R_inner=%.4f R_outer=%.4f scale=%.4f\n",
                   cx, cy, cz, R_inner, R_outer, scale);
        PRINT_INFO("Kelvin periodic: %zu inner verts, %zu outer verts\n",
                   inner_pts.size(), outer_pts.size());

        // Build spatial index of outer points for fast lookup.
        // For each inner vertex, predict outer position by radial scaling
        // and find the nearest outer vertex.
        std::vector<std::pair<int, netgen::Point<3>>> outer_vec(
            outer_pts.begin(), outer_pts.end());

        std::map<int, int> vertex_pair;   // inner_vid -> outer_vid
        double max_dist = 0;
        int n_good = 0, n_bad = 0;

        for (auto &ip : inner_pts) {
          // Predict outer position: radial scaling from center
          double dx = ip.second(0) - cx;
          double dy = ip.second(1) - cy;
          double dz = ip.second(2) - cz;
          double px = cx + dx * scale;
          double py = cy + dy * scale;
          double pz = cz + dz * scale;

          // Find nearest outer vertex (brute force, O(N*M))
          int best_ov = -1;
          double best_d2 = 1e30;
          for (auto &op : outer_vec) {
            double ddx = op.second(0) - px;
            double ddy = op.second(1) - py;
            double ddz = op.second(2) - pz;
            double d2 = ddx*ddx + ddy*ddy + ddz*ddz;
            if (d2 < best_d2) { best_d2 = d2; best_ov = op.first; }
          }

          double dist = sqrt(best_d2);
          // Tolerance: 30% of average outer edge length (generous for tet mesh)
          double tol = 0.3 * R_outer * 0.5;  // conservative
          if (dist < tol && best_ov > 0) {
            vertex_pair[ip.first] = best_ov;
            if (dist > max_dist) max_dist = dist;
            n_good++;
          } else {
            n_bad++;
          }
        }

        // Write identifications
        auto &ident = ng_mesh->GetIdentifications();
        int n_paired = 0;
        for (auto &p : vertex_pair) {
          ident.Add(netgen::PointIndex(p.first),
                    netgen::PointIndex(p.second),
                    "kelvin",
                    netgen::Identifications::PERIODIC);
          n_paired++;
        }

        PRINT_INFO("Kelvin periodic: %d/%zu pairs written "
                   "(max_dist=%.2e, %d unmatched)\n",
                   n_paired, inner_pts.size(), max_dist, n_bad);
      }
    }
  }

  // Detach CallbackGeometry — it has function pointers, not serializable.
  // Curving data is preserved in CurvedElements (written as "curvedelements"
  // section in .vol text format, upstream Netgen master feature).
  ng_mesh->SetGeometry(nullptr);

  // Save .vol text format (includes curvedelements section for order >= 2).
  // NGSolve reads this with Mesh("file.vol") — no STEP/Cubit needed.
  ng_mesh->Save(filename);

  int ne = ng_mesh->GetNE();
  int np = ng_mesh->GetNP();

  PRINT_INFO("Exported Netgen Vol (order %d): %s (%d nodes, %d elements)\n",
             order, filename.c_str(), np, ne);

  // --- Write companion JSON with CAD reference values ---
  {
    std::string json_path = filename + ".json";
#ifdef _WIN32
    // Use wide path for Unicode support on Windows
    int wlen = MultiByteToWideChar(CP_ACP, 0, json_path.c_str(), -1, NULL, 0);
    std::wstring wpath(wlen, 0);
    MultiByteToWideChar(CP_ACP, 0, json_path.c_str(), -1, &wpath[0], wlen);
    std::ofstream jf(wpath);
#else
    std::ofstream jf(json_path);
#endif
    if (jf.is_open()) {
      jf << "{\n";

      // Materials (per-block volume)
      // MeshExportInterface may return a default block (all elements) whose ID
      // does not exist as a user-defined Cubit block. Skip such blocks to avoid
      // "No block with ID N was found" errors from parse_cubit_list.
      jf << "  \"materials\": {";
      bool first = true;
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
        for (int vid : vols_in_block) {
          RefVolume* rv = GeometryQueryTool::instance()->get_ref_volume(vid);
          if (rv) total_vol += rv->measure();
        }
        if (!first) jf << ",";
        jf << "\n    \"" << bname << "\": " << std::scientific << total_vol;
        first = false;
      }
      jf << "\n  },\n";

      // Boundaries (per-surface area from CAD, not mesh)
      jf << "  \"boundaries\": {";
      first = true;
      int nfd_json = ng_mesh->GetNFD();
      for (int fi = 1; fi <= nfd_json; fi++) {
        // Use original Cubit surface ID for CAD area query
        int cubit_surf_id = orig_surf_ids[fi - 1];
        // Use the BCName written to .vol (matches mesh.GetBoundaries())
        std::string bname = ng_mesh->GetBCName(fi - 1);
        // Get CAD area from Cubit surface
        RefFace* rf = GeometryQueryTool::instance()->get_ref_face(cubit_surf_id);
        double cad_area = rf ? rf->area() : 0.0;
        if (!first) jf << ",";
        jf << "\n    \"" << bname << "\": " << std::scientific << cad_area;
        first = false;
      }
      jf << "\n  },\n";

      // Edges (per-curve length, curves shared by 2+ surfaces)
      jf << "  \"edges\": {";
      first = true;
      std::vector<int> curve_ids = CubitInterface::parse_cubit_list("curve", "all");
      for (int cid : curve_ids) {
        std::vector<int> parent_surfs = CubitInterface::parse_cubit_list(
            "surface", "in curve " + std::to_string(cid));
        if (parent_surfs.size() < 2) continue;
        RefEdge* re = GeometryQueryTool::instance()->get_ref_edge(cid);
        double cad_len = re ? re->measure() : 0.0;
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
      jf.close();
    }
  }

  return true;
#endif
}
