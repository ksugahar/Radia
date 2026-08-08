#include "ExportNetgenCommand.hpp"
#include "MeshData.hpp"
#include "RadiaMessageFilter.hpp"
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"
#include <chrono>
#include <cmath>
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
#include "utf8_path.hpp"

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
      // path = .../plugins/cubit_mesh_export.ccm or .../bin/cubit_mesh_export.ccl
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

#ifdef HAVE_NETGEN
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
    const auto &element = mesh.SurfaceElement(sei);
    if (element.GetIndex() != descriptor)
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
#endif

ExportNetgenCommand::ExportNetgenCommand() {}
ExportNetgenCommand::~ExportNetgenCommand() {}

std::vector<std::string> ExportNetgenCommand::get_syntax()
{
  std::vector<std::string> syntax_list;
  // Kelvin / symmetry options were APREPRO-promoted on 2026-05-05.  The
  // .vol is the only mesh format that consumes Kelvin (open-boundary
  // FEM), so these knobs live on the netgen subcommand instead of as a
  // global GUI launcher.  Each kelvin_sym_<axis> takes a string in
  // {off, bn, ht}; "off" is the default and means "no reduction".
  syntax_list.push_back(
    "export netgen <string:label='filename',help='<filename>'> "
    "[order <value:label='order',help='<1-5>'>] "
    "[overwrite] "
    "[add_kelvin] "
    "[kelvin_air <string:label='kelvin_air',help='<block name, default \"air\">'>] "
    "[kelvin_block <string:label='kelvin_block',help='<block name, default \"kelvin\">'>] "
    "[kelvin_mesh <value:label='kelvin_mesh',help='<size in m, blank=auto>'>] "
    "[kelvin_sym_x <string:label='kelvin_sym_x',help='<off|bn|ht>'>] "
    "[kelvin_sym_y <string:label='kelvin_sym_y',help='<off|bn|ht>'>] "
    "[kelvin_sym_z <string:label='kelvin_sym_z',help='<off|bn|ht>'>]"
  );
  return syntax_list;
}

std::vector<std::string> ExportNetgenCommand::get_syntax_help()
{
  std::vector<std::string> help;
  help.push_back(
    "export netgen \"filename.vol\" [order {1|2|3|4|5}] [overwrite]\n"
    "                  [add_kelvin]\n"
    "                  [kelvin_air \"air\"] [kelvin_block \"kelvin\"]\n"
    "                  [kelvin_mesh <size_m>]\n"
    "                  [kelvin_sym_x {off|bn|ht}]\n"
    "                  [kelvin_sym_y {off|bn|ht}]\n"
    "                  [kelvin_sym_z {off|bn|ht}]"
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
    "Order 1-5 supported (default: 2). Requires Netgen.\n"
    "\n"
    "Kelvin open-boundary (optional, .vol only):\n"
    "  add_kelvin                   Auto-create the exterior Kelvin sphere\n"
    "                               (idempotent: skipped if a 'kelvin' block\n"
    "                               already exists; needs an 'air' block).\n"
    "  kelvin_mesh <size_m>         Tet edge length on the Kelvin shell;\n"
    "                               omit to inherit from the air surface.\n"
    "  kelvin_sym_<axis> bn|ht      Per-axis symmetry-plane BC label on a\n"
    "                               domain-reduced (1/2 or 1/4) model.\n"
    "                               'bn' = B.n=0 (flux parallel, Radia '+').\n"
    "                               'ht' = HxN=0 (flux perp,  Radia '-').\n"
    "                               'off' = no reduction (default)."
  );
  return help;
}

// ----------------------------------------------------------------
// find_plugin_dir() -- return the directory containing this .ccm
// (and, by deployment convention, the cubit_helpers/ subdirectory
// holding add_kelvin.py and auto_kelvin_entry.py).  Mirrors the
// search used by ensure_netgen_dll_path.
// ----------------------------------------------------------------
static std::string find_plugin_dir()
{
#ifdef _WIN32
  const char *pd = std::getenv("CUBIT_PLUGIN_DIR");
  if (pd && pd[0])
    return std::string(pd);
  HMODULE hm = nullptr;
  GetModuleHandleExA(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS |
                     GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
                     (LPCSTR)&find_plugin_dir, &hm);
  if (hm) {
    char path[MAX_PATH];
    if (GetModuleFileNameA(hm, path, MAX_PATH)) {
      std::string dir(path);
      auto pos = dir.find_last_of("\\/");
      if (pos != std::string::npos)
        return dir.substr(0, pos);
    }
  }
#endif
  return std::string();
}

// ----------------------------------------------------------------
// run_auto_kelvin -- write a JSON config and `play` the Cubit-side
// auto_kelvin_entry.py, which calls add_kelvin_cubit() in Cubit's
// embedded Python.  Runs BEFORE mesh extract so the new Kelvin
// volume / sidesets are present when extract() runs.
//
// add_kelvin == false short-circuits without writing the JSON.
// Failures are logged but do NOT abort the export -- user gets
// Dirichlet truncation on the air outer surface as a fallback.
// ----------------------------------------------------------------
static void run_auto_kelvin(bool add_kelvin,
                            const std::string &air_block,
                            const std::string &kelvin_block,
                            bool has_mesh_size, double mesh_size,
                            const std::string &sym_x,
                            const std::string &sym_y,
                            const std::string &sym_z,
                            const std::string &filename)
{
  if (!add_kelvin) return;

  std::string plugin_dir = find_plugin_dir();
  if (plugin_dir.empty()) {
    PRINT_WARNING("Auto-Kelvin skipped: cannot locate plugin directory.\n");
    return;
  }
  std::string helpers_dir = plugin_dir + "/cubit_helpers";
  std::string entry = helpers_dir + "/auto_kelvin_entry.py";
  // Normalize separators
  for (auto &c : entry)       if (c == '\\') c = '/';
  for (auto &c : helpers_dir) if (c == '\\') c = '/';

  std::ifstream probe(entry);
  if (!probe.good()) {
    PRINT_WARNING("Auto-Kelvin skipped: %s not found.  Re-run "
                  "cubit-plugin-install to deploy cubit_helpers/.\n",
                  entry.c_str());
    return;
  }
  probe.close();

  // Place the config JSON next to the output .vol so the user can see
  // exactly what was passed.
  std::string out_dir = filename;
  auto pos = out_dir.find_last_of("\\/");
  if (pos != std::string::npos) out_dir = out_dir.substr(0, pos);
  else                          out_dir = ".";
  std::string cfg_path = out_dir + "/radia_kelvin_config.json";
  for (auto &c : cfg_path) if (c == '\\') c = '/';

  // Build the JSON config.  Schema matches auto_kelvin_entry.py.
  auto bc_token = [](const std::string &s) -> std::string {
    if (s == "bn" || s == "bn=0") return "bn=0";
    if (s == "ht" || s == "ht=0") return "ht=0";
    return std::string();   // "off" / "" -> no reduction on this axis
  };
  std::string rx = bc_token(sym_x);
  std::string ry = bc_token(sym_y);
  std::string rz = bc_token(sym_z);
  bool has_reduction = !rx.empty() || !ry.empty() || !rz.empty();

  {
    std::ofstream jf(cfg_path);
    if (!jf.good()) {
      PRINT_WARNING("Auto-Kelvin: cannot write %s; aborting Kelvin step.\n",
                    cfg_path.c_str());
      return;
    }
    jf << "{";
    jf << "\"add_kelvin\": true";
    jf << ", \"kelvin_air_block\": \""  << air_block    << "\"";
    jf << ", \"kelvin_block_name\": \"" << kelvin_block << "\"";
    if (has_mesh_size && mesh_size > 0.0)
      jf << ", \"kelvin_mesh_size\": " << std::scientific << mesh_size;
    else
      jf << ", \"kelvin_mesh_size\": null";
    if (has_reduction) {
      jf << ", \"kelvin_reduction\": {";
      bool first = true;
      // Cannot name this lambda `emit` -- the .ccl build path links
      // Qt5 and Qt's `emit` keyword is a preprocessor macro defined
      // to nothing, which would expand `auto emit = ...` to
      // `auto = ...` and fail with C2513.
      auto write_axis = [&](const char *axis, const std::string &bc) {
        if (bc.empty()) return;
        if (!first) jf << ", ";
        jf << "\"" << axis << "\": \"" << bc << "\"";
        first = false;
      };
      write_axis("x", rx);
      write_axis("y", ry);
      write_axis("z", rz);
      jf << "}";
    } else {
      jf << ", \"kelvin_reduction\": null";
    }
    jf << "}\n";
  }

  PRINT_INFO("Auto-Kelvin: config -> %s\n", cfg_path.c_str());

#ifdef _WIN32
  // Wide-char env vars survive non-ASCII paths cleanly.
  std::wstring wcfg(cfg_path.begin(), cfg_path.end());
  std::wstring whlp(helpers_dir.begin(), helpers_dir.end());
  SetEnvironmentVariableW(L"RADIA_LAUNCHER_CONFIG", wcfg.c_str());
  SetEnvironmentVariableW(L"CUBIT_HELPERS_DIR",     whlp.c_str());
#else
  setenv("RADIA_LAUNCHER_CONFIG", cfg_path.c_str(),    1);
  setenv("CUBIT_HELPERS_DIR",     helpers_dir.c_str(), 1);
#endif

  std::string play_cmd = "play \"" + entry + "\"";
  PRINT_INFO("Auto-Kelvin: %s\n", play_cmd.c_str());
  CubitInterface::cmd(play_cmd.c_str());

#ifdef _WIN32
  SetEnvironmentVariableW(L"RADIA_LAUNCHER_CONFIG", nullptr);
  SetEnvironmentVariableW(L"CUBIT_HELPERS_DIR",     nullptr);
#else
  unsetenv("RADIA_LAUNCHER_CONFIG");
  unsetenv("CUBIT_HELPERS_DIR");
#endif
}

bool ExportNetgenCommand::execute(CubitCommandData &data)
{
#ifndef HAVE_NETGEN
  PRINT_ERROR("export netgen requires Netgen support (not built).\n");
  return false;
#else
  // Suppress Cubit Learn Edition's 50k-cap ERROR during the export.
  // Radia bypasses the cap and exports successfully regardless, so
  // the ERROR line is misleading noise confusing users in logs.
  radia::ScopedLearnEditionFilter _lef_guard;

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

  // ---- Auto-Kelvin (optional) ----------------------------------
  // `add_kelvin` is a bare keyword flag (matches `overwrite` /
  // `nopyramid` in the other export commands).  find_keyword is the
  // Cubit SDK probe for "did the user supply this token".
  bool add_kelvin = data.find_keyword("add_kelvin");
  std::string kelvin_air = "air";
  std::string kelvin_block = "kelvin";
  data.get_string("kelvin_air", kelvin_air);
  data.get_string("kelvin_block", kelvin_block);
  if (kelvin_air.empty())   kelvin_air   = "air";
  if (kelvin_block.empty()) kelvin_block = "kelvin";

  // get_value leaves the variable untouched when the keyword is
  // absent.  Detect "user supplied a positive size" by initialising
  // to 0 and re-checking.
  double kelvin_mesh = 0.0;
  data.get_value("kelvin_mesh", kelvin_mesh);
  bool has_mesh_size = (kelvin_mesh > 0.0);

  std::string sym_x = "off", sym_y = "off", sym_z = "off";
  data.get_string("kelvin_sym_x", sym_x);
  data.get_string("kelvin_sym_y", sym_y);
  data.get_string("kelvin_sym_z", sym_z);

  if (add_kelvin) {
    run_auto_kelvin(true, kelvin_air, kelvin_block,
                    has_mesh_size, kelvin_mesh,
                    sym_x, sym_y, sym_z, filename);
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
      if (cubit_sid < 0) {
        // NetgenCurver recovered this direct/free face's owner from adjacent
        // volume elements.  Validate that handoff instead of overwriting every
        // synthetic descriptor with domain 1.
        int domin = ng_mesh->GetFaceDescriptor(fi).DomainIn();
        int domout = ng_mesh->GetFaceDescriptor(fi).DomainOut();
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
    std::string ssname = sg.name.empty()
        ? "sideset_" + std::to_string(sg.id) : sg.name;
    surf_to_ssname[-sg.id] = ssname;
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
        int di = ng_mesh->GetFaceDescriptor(fi).DomainIn();
        int dout = ng_mesh->GetFaceDescriptor(fi).DomainOut();
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
        netgen::PointIndex pi(it->second);
        ng_mesh->pointelements.Append(netgen::Element0d(pi, cd3_idx));
        point_elements_written++;
      }
    }
    if (cd3_names_written > 0) {
      PRINT_INFO("Nodeset export: %d CD3 name(s), %d Element0d point(s)\n",
                 cd3_names_written, point_elements_written);
    }
  }

  // ---- Kelvin periodic identification ----
  // Kelvin = two identical spheres, offset in space.
  // kelvin_int / kelvin_ext labels come from sidesets (set in .py) or
  // auto-detect (concentric shell only).
  // Identification = TRANSLATION: offset = mean(outer) - mean(inner).
  // With copy mesh, inner and outer have 1:1 vertex correspondence.
  //
  // Algorithm: translation-based nearest-neighbor vertex matching.
  //   1. Collect inner/outer vertex positions from bc names
  //   2. Offset = mean(outer) - mean(inner)
  //   3. For each inner vertex, find nearest unused outer vertex at (inner + offset)
  //   4. Write identification pairs via ident.Add()
  //
  // IMPORTANT (2026-04-25): match ONLY face descriptors with the exact
  // bcnames kelvin_int / kelvin_ext -- not a substring test.  In
  // reduction mode the sym_*_* cut faces border the Kelvin domain
  // too; matching them as "kelvin" candidates would corrupt the
  // pair set.
  {
    std::set<int> fd_inner_set, fd_outer_set;
    for (int fi = 1; fi <= ng_mesh->GetNFD(); fi++) {
      std::string bc = ng_mesh->GetBCName(fi - 1);
      if (bc == "kelvin_int") fd_inner_set.insert(fi);
      else if (bc == "kelvin_ext") fd_outer_set.insert(fi);
    }

    if (!fd_inner_set.empty() && !fd_outer_set.empty()) {
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
        // Translation offset = mean(outer) - mean(inner)
        double mix = 0, miy = 0, miz = 0;
        for (auto &p : inner_pts) {
          mix += p.second(0); miy += p.second(1); miz += p.second(2);
        }
        double ni = (double)inner_pts.size();
        mix /= ni; miy /= ni; miz /= ni;

        double mox = 0, moy = 0, moz = 0;
        for (auto &p : outer_pts) {
          mox += p.second(0); moy += p.second(1); moz += p.second(2);
        }
        double no = (double)outer_pts.size();
        mox /= no; moy /= no; moz /= no;

        double tx = mox - mix, ty = moy - miy, tz = moz - miz;

        // Pairing tolerance: was hard-coded at 1e-2 m which is too
        // tight for accelerator-scale models (R~0.4 m, mesh edge
        // ~20 mm).  Scale to 5% of the largest mesh edge length so
        // copy-mesh imperfections at corner vertices still match.
        // Cap at 5% of |offset| to avoid pairing wrong vertices for
        // models with very fine meshes near the curved cap boundary.
        double offset_len = sqrt(tx*tx + ty*ty + tz*tz);
        // Estimate inner-side characteristic mesh edge length from
        // the bounding-box span of the inner vertex set.
        double ix_lo = 1e30, ix_hi = -1e30, iy_lo = 1e30, iy_hi = -1e30,
               iz_lo = 1e30, iz_hi = -1e30;
        for (auto &ip : inner_pts) {
          double x = ip.second(0), y = ip.second(1), z = ip.second(2);
          if (x < ix_lo) ix_lo = x;  if (x > ix_hi) ix_hi = x;
          if (y < iy_lo) iy_lo = y;  if (y > iy_hi) iy_hi = y;
          if (z < iz_lo) iz_lo = z;  if (z > iz_hi) iz_hi = z;
        }
        double inner_diag = sqrt((ix_hi-ix_lo)*(ix_hi-ix_lo) +
                                 (iy_hi-iy_lo)*(iy_hi-iy_lo) +
                                 (iz_hi-iz_lo)*(iz_hi-iz_lo));
        // Heuristic edge length: diagonal / sqrt(N_pts) for a 2D
        // surface mesh of N points.  (std::min/max wrapped in parens
        // to dodge the Windows.h min/max macros.)
        size_t n_pts = inner_pts.size();
        if (n_pts < 4) n_pts = 4;
        double mesh_edge = inner_diag / sqrt((double)n_pts);
        double cap_offset = 0.05 * offset_len;
        double cap_edge   = 0.5  * mesh_edge;
        double tol_upper  = (cap_offset < cap_edge) ? cap_offset : cap_edge;
        double match_tol  = (1e-3 > tol_upper) ? 1e-3 : tol_upper;

        PRINT_INFO("Kelvin periodic: %zu inner, %zu outer verts, "
                   "offset=(%.4f,%.4f,%.4f), match_tol=%.3e\n",
                   inner_pts.size(), outer_pts.size(), tx, ty, tz,
                   match_tol);

        // Match: for each inner vertex, find nearest unused outer at (inner + offset)
        std::vector<std::pair<int, netgen::Point<3>>> outer_vec(
            outer_pts.begin(), outer_pts.end());
        std::map<int, int> vertex_pair;
        std::set<int> used_outer;
        double max_dist = 0;
        int n_bad = 0;

        for (auto &ip : inner_pts) {
          double px = ip.second(0) + tx;
          double py = ip.second(1) + ty;
          double pz = ip.second(2) + tz;

          int best_ov = -1;
          double best_d2 = 1e30;
          for (auto &op : outer_vec) {
            if (used_outer.count(op.first)) continue;
            double dx = op.second(0) - px;
            double dy = op.second(1) - py;
            double dz = op.second(2) - pz;
            double d2 = dx*dx + dy*dy + dz*dz;
            if (d2 < best_d2) { best_d2 = d2; best_ov = op.first; }
          }

          double dist = sqrt(best_d2);
          if (dist < match_tol && best_ov > 0) {
            vertex_pair[ip.first] = best_ov;
            used_outer.insert(best_ov);
            if (dist > max_dist) max_dist = dist;
          } else {
            n_bad++;
          }
        }

        // All-or-nothing policy: only write identification pairs if
        // every inner vertex matched.  A partial identification leaves
        // 4-out-of-N vertices "almost paired" -- NGSolve's Mesh()
        // wrapper rejects such .vol files with NgException("Ask for
        // unused hash-value") because it expects either a complete
        // 1:1 pairing or none at all.  When there are unmatched
        // vertices (typically corner / edge vertices where copy-mesh
        // can leave sub-pixel mismatches at large geometry scale),
        // we skip the C++ side entirely and rely on NGSolve's
        // `add_periodic_kelvin` (in calc_common.py) to set up the
        // identification at solve time using bcname-driven matching.
        if (n_bad == 0) {
          auto &ident = ng_mesh->GetIdentifications();
          for (auto &p : vertex_pair) {
            ident.Add(netgen::PointIndex(p.first),
                      netgen::PointIndex(p.second),
                      "kelvin",
                      netgen::Identifications::PERIODIC);
          }
        } else {
          PRINT_WARNING("Kelvin periodic: %d vertex pair(s) unmatched "
                        "(max_dist=%.2e, tol=%.2e); skipping C++ "
                        "identification entirely.  NGSolve's "
                        "add_periodic_kelvin will build identification "
                        "at solve time.\n",
                        n_bad, max_dist, match_tol);
        }

        PRINT_INFO("Kelvin periodic: %d/%zu pairs (max_dist=%.2e, %d unmatched)\n",
                   (int)vertex_pair.size(), inner_pts.size(), max_dist, n_bad);
      }
    }
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
  ng_mesh->Save(u8_string_to_path(filename));

  int ne = ng_mesh->GetNE();
  int np = ng_mesh->GetNP();

  PRINT_INFO("Exported Netgen Vol (order %d): %s (%d nodes, %d elements)\n",
             order, filename.c_str(), np, ne);

  // --- Write companion JSON with CAD reference values ---
  {
    std::string json_path = filename + ".json";
    // ofstream on MSVC accepts const wchar_t* / std::filesystem::path for
    // Unicode filenames.  Build via u8_string_to_path so UTF-8 path survives.
    std::ofstream jf(u8_string_to_path(json_path));
    if (jf.is_open()) {
      jf << std::setprecision(16);
      jf << "{\n";

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
      std::set<int> exported_surface_ids;
      for (int surface_id : orig_surf_ids) {
        if (surface_id > 0)
          exported_surface_ids.insert(surface_id);
      }
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
      jf.close();
    }
  }

  return true;
#endif
}
