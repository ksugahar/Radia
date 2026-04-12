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

  // ---- Kelvin periodic identification (optional) ----
  // Two-sphere approach: interior sphere (physical) + exterior sphere (Kelvin image).
  // Both have the same radius R, offset in space.
  //
  // Pairing must preserve TRIANGLE TOPOLOGY: every inner triangle must have a
  // matching outer triangle when its vertices are remapped via the pair table.
  // A naive nearest-neighbor vertex match does not guarantee this — small
  // numerical noise can pair vertices in a way that breaks triangle adjacency.
  //
  // Algorithm:
  //   1. Collect inner/outer triangles (3 vertex IDs each)
  //   2. Estimate translation T = mean(outer) - mean(inner)
  //   3. For each inner triangle, find the outer triangle with closest centroid
  //      after translation
  //   4. For each matched pair, try the 3 cyclic rotations of the inner vertices
  //      against the outer vertices (orientation may flip, so also try reverse)
  //      and pick the rotation with minimum total vertex distance
  //   5. Record (inner_vid -> outer_vid) from the chosen rotation
  //   6. Verify consistency: every inner vertex must always map to the same
  //      outer vertex across all triangles it belongs to
  {
    std::set<int> fd_inner_set, fd_outer_set;
    for (int fi = 1; fi <= ng_mesh->GetNFD(); fi++) {
      std::string bc = ng_mesh->GetBCName(fi - 1);
      if (bc == "kelvin_int") fd_inner_set.insert(fi);
      else if (bc == "kelvin_ext") fd_outer_set.insert(fi);
    }

    if (!fd_inner_set.empty() && !fd_outer_set.empty()) {
      struct Tri { int v[3]; double cx, cy, cz; };
      std::vector<Tri> inner_tris, outer_tris;
      std::map<int, netgen::Point<3>> inner_pts, outer_pts;

      for (int sei = 1; sei <= ng_mesh->GetNSE(); sei++) {
        const auto &sel = ng_mesh->SurfaceElement(sei);
        int fd = sel.GetIndex();
        bool is_inner = fd_inner_set.count(fd) > 0;
        bool is_outer = fd_outer_set.count(fd) > 0;
        if (!is_inner && !is_outer) continue;
        if (sel.GetNP() < 3) continue;  // Skip non-triangle elements

        Tri t;
        t.cx = t.cy = t.cz = 0;
        for (int j = 0; j < 3; j++) {
          int pi = sel[j];
          t.v[j] = pi;
          auto &pt = ng_mesh->Point(netgen::PointIndex(pi));
          t.cx += pt(0); t.cy += pt(1); t.cz += pt(2);
          auto &target = is_inner ? inner_pts : outer_pts;
          if (target.find(pi) == target.end()) target[pi] = pt;
        }
        t.cx /= 3; t.cy /= 3; t.cz /= 3;
        (is_inner ? inner_tris : outer_tris).push_back(t);
      }

      if (!inner_tris.empty() && !outer_tris.empty()) {
        // Translation: centroid(outer) - centroid(inner)
        double mean_ix = 0, mean_iy = 0, mean_iz = 0;
        double mean_ox = 0, mean_oy = 0, mean_oz = 0;
        for (auto &p : inner_pts) {
          mean_ix += p.second(0); mean_iy += p.second(1); mean_iz += p.second(2);
        }
        for (auto &p : outer_pts) {
          mean_ox += p.second(0); mean_oy += p.second(1); mean_oz += p.second(2);
        }
        double ni = inner_pts.size(), no = outer_pts.size();
        double tx = mean_ox/no - mean_ix/ni;
        double ty = mean_oy/no - mean_iy/ni;
        double tz = mean_oz/no - mean_iz/ni;

        // Step 1: For each inner triangle, find matching outer triangle by centroid
        // (linear search is O(N^2); fine for typical Kelvin sphere sizes)
        std::map<int, int> vertex_pair;          // inner_vid -> outer_vid
        std::map<int, int> vertex_pair_count;    // for consistency check
        std::map<int, std::map<int,int>> vertex_pair_votes;  // inner -> {outer -> count}
        double max_vertex_dist = 0;
        int n_tri_matched = 0;
        int n_tri_unmatched = 0;

        for (auto &t_in : inner_tris) {
          double ex = t_in.cx + tx, ey = t_in.cy + ty, ez = t_in.cz + tz;

          // Find nearest outer triangle by centroid
          int best_oi = -1;
          double best_d2 = 1e30;
          for (size_t k = 0; k < outer_tris.size(); k++) {
            auto &t_out = outer_tris[k];
            double dx = t_out.cx - ex, dy = t_out.cy - ey, dz = t_out.cz - ez;
            double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 < best_d2) { best_d2 = d2; best_oi = (int)k; }
          }
          if (best_oi < 0) { n_tri_unmatched++; continue; }
          n_tri_matched++;
          auto &t_out = outer_tris[best_oi];

          // Try 6 permutations: 3 cyclic rotations × 2 mirror (orientation flip)
          // Cycle: (0,1,2), (1,2,0), (2,0,1)
          // Reverse: (0,2,1), (2,1,0), (1,0,2)
          int perms[6][3] = {
            {0,1,2}, {1,2,0}, {2,0,1},
            {0,2,1}, {2,1,0}, {1,0,2}
          };
          double best_perm_d2 = 1e30;
          int best_perm = -1;
          for (int p = 0; p < 6; p++) {
            double sum = 0;
            for (int j = 0; j < 3; j++) {
              auto &pi = inner_pts[t_in.v[j]];
              auto &po = outer_pts[t_out.v[perms[p][j]]];
              double dx = po(0) - (pi(0)+tx);
              double dy = po(1) - (pi(1)+ty);
              double dz = po(2) - (pi(2)+tz);
              sum += dx*dx + dy*dy + dz*dz;
            }
            if (sum < best_perm_d2) { best_perm_d2 = sum; best_perm = p; }
          }

          // Vote: each triangle proposes a vertex pairing
          for (int j = 0; j < 3; j++) {
            int vi = t_in.v[j];
            int vo = t_out.v[perms[best_perm][j]];
            vertex_pair_votes[vi][vo]++;
          }
          double d = sqrt(best_perm_d2 / 3.0);
          if (d > max_vertex_dist) max_vertex_dist = d;
        }

        // Step 2: Resolve final pairing by majority vote per inner vertex.
        // Each vertex appears in multiple triangles; the correct outer pair
        // is the one most-voted (should be unanimous if mesh is conformal).
        int n_consistent = 0;
        int n_conflicts = 0;
        for (auto &kv : vertex_pair_votes) {
          int vi = kv.first;
          int best_vo = -1, best_votes = 0, total = 0;
          for (auto &vo_count : kv.second) {
            total += vo_count.second;
            if (vo_count.second > best_votes) {
              best_votes = vo_count.second;
              best_vo = vo_count.first;
            }
          }
          if (best_vo > 0) {
            vertex_pair[vi] = best_vo;
            if (best_votes < total) n_conflicts++;
            else n_consistent++;
          }
        }

        // Step 3: Write identifications
        auto &ident = ng_mesh->GetIdentifications();
        int n_paired = 0;
        for (auto &p : vertex_pair) {
          ident.Add(netgen::PointIndex(p.first),
                    netgen::PointIndex(p.second),
                    "kelvin",
                    netgen::Identifications::PERIODIC);
          n_paired++;
        }

        PRINT_INFO("Kelvin periodic: %d vertex pairs from %d/%d matched triangles "
                   "(%zu inner pts, %zu outer pts, offset=(%.4f,%.4f,%.4f), "
                   "max_vert_dist=%.2e, conflicts=%d/%d)\n",
                   n_paired, n_tri_matched, (int)inner_tris.size(),
                   inner_pts.size(), outer_pts.size(), tx, ty, tz,
                   max_vertex_dist, n_conflicts, n_consistent + n_conflicts);
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
