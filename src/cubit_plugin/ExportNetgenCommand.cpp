#include "ExportNetgenCommand.hpp"
#include "MeshData.hpp"
#include "LearnEditionMessageFilter.hpp"
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"
#include <chrono>
#include <cmath>
#include <iomanip>
#include <filesystem>

#ifdef HAVE_NETGEN
#include "NetgenCurver.hpp"
#include "NamedPeriodicIdentification.hpp"
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
    "[periodic_name <string:label='periodic_name'>] "
    "[periodic_master <string:label='periodic_master'>] "
    "[periodic_slave <string:label='periodic_slave'>] "
    "[periodic_unit <string:label='periodic_unit'>] "
    "[periodic_patches <string:label='periodic_patches'>] "
    "[periodic_tx <value:label='periodic_tx'>] "
    "[periodic_ty <value:label='periodic_ty'>] "
    "[periodic_tz <value:label='periodic_tz'>] "
    "[periodic_match <value:label='periodic_match'>] "
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
    "Named periodic metadata requires all fields explicitly: periodic_name kt,\n"
    "  periodic_master/slave <BC label>, periodic_unit m,\n"
    "  periodic_patches <masterBC:slaveBC,...>, periodic_tx/ty/tz <translation>,\n"
    "  periodic_match <positive geometric distance>. Fresh output required.\n"
    "No inferred pairs or consumer repair; q3 geometry is preserved.\n"
    "\n"
    "Kelvin open-boundary (optional, .vol only):\n"
    "  add_kelvin                   Auto-create the exterior Kelvin sphere\n"
    "                               (idempotent: skipped if a 'kelvin' block\n"
    "                               already exists; needs an 'air' block).\n"
    "  kelvin_mesh <size_m>         Tet edge length on the Kelvin shell;\n"
    "                               omit to inherit from the air surface.\n"
    "  kelvin_sym_<axis> bn|ht      Per-axis symmetry-plane BC label on a\n"
    "                               domain-reduced (1/2 or 1/4) model.\n"
    "                               'bn' = B.n=0 (flux parallel).\n"
    "                               'ht' = HxN=0 (flux perpendicular).\n"
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
  std::string cfg_path = out_dir + "/cubit_mesh_export_kelvin_config.json";
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
  SetEnvironmentVariableW(L"CUBIT_MESH_EXPORT_KELVIN_CONFIG", wcfg.c_str());
  SetEnvironmentVariableW(L"CUBIT_HELPERS_DIR",     whlp.c_str());
#else
  setenv("CUBIT_MESH_EXPORT_KELVIN_CONFIG", cfg_path.c_str(), 1);
  setenv("CUBIT_HELPERS_DIR",     helpers_dir.c_str(), 1);
#endif

  std::string play_cmd = "play \"" + entry + "\"";
  PRINT_INFO("Auto-Kelvin: %s\n", play_cmd.c_str());
  CubitInterface::cmd(play_cmd.c_str());

#ifdef _WIN32
  SetEnvironmentVariableW(L"CUBIT_MESH_EXPORT_KELVIN_CONFIG", nullptr);
  SetEnvironmentVariableW(L"CUBIT_HELPERS_DIR",     nullptr);
#else
  unsetenv("CUBIT_MESH_EXPORT_KELVIN_CONFIG");
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
  // The exporter bypasses the cap and completes successfully, so
  // the ERROR line is misleading noise confusing users in logs.
  cubit_mesh_export::ScopedLearnEditionFilter _lef_guard;

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

  cubit_mesh_export::NamedPeriodicContract periodic;
  std::string periodic_patch_text;
  const bool periodic_fields[] = {
      data.get_string("periodic_name", periodic.name),
      data.get_string("periodic_master", periodic.master_label),
      data.get_string("periodic_slave", periodic.slave_label),
      data.get_string("periodic_unit", periodic.coordinate_unit),
      data.get_string("periodic_patches", periodic_patch_text),
      data.get_value("periodic_tx", periodic.translation[0]),
      data.get_value("periodic_ty", periodic.translation[1]),
      data.get_value("periodic_tz", periodic.translation[2]),
      data.get_value("periodic_match", periodic.matching_distance)};
  int periodic_present = 0;
  for (bool present : periodic_fields) periodic_present += present ? 1 : 0;
  const bool has_periodic = periodic_present != 0;
  if (has_periodic && periodic_present != 9) {
    PRINT_ERROR("Named periodic export requires all nine explicit periodic fields.\n");
    return false;
  }
  if (has_periodic) {
    try {
      cubit_mesh_export::require_fresh_named_periodic_output(
          cubit_mesh_export::named_periodic_output_paths(u8_string_to_path(filename)));
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
  cubit_mesh_export::NamedPeriodicResult periodic_result;

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
  md.require_volume_elements = true;  // .vol has no surface-only mode
  int build_order = (order >= 2) ? order : 2;
  if (!md.extract(build_order)) {
    PRINT_ERROR("Mesh extraction failed.\n");
    auto failed = md.get_netgen_curver();
    if (failed && !has_periodic) {
      std::ofstream failure(u8_string_to_path(filename + ".json"));
      if (failure) failure << "{\"export_passed\":false,\"curving\":"
                           << failed->edge_diagnostics_json() << "}\n";
    }
    if (has_periodic && failed)
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
      ng_mesh->VolumeElement(ve_idx).SetIndex(di);
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
      int cubit_sid = ng_mesh->GetFaceDescriptor(fi).BCProperty();
      if (cubit_sid > 0)
        mesh_support_surface_ids.insert(cubit_sid);
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
      auto &surface_element = ng_mesh->SurfaceElements()[sei];
      if (same_material_face_descriptors.count(surface_element.GetIndex())) {
        ng_mesh->SurfaceElements().DeleteElement(sei);
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
      const auto &old_fd = ng_mesh->GetFaceDescriptor(fi);
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
    for (auto &surface_element : ng_mesh->SurfaceElements()) {
      const int old_index = surface_element.GetIndex();
      if (old_index <= 0 || old_index > old_nfd || old_to_new[old_index] == 0) {
        PRINT_ERROR("Cannot remap retained surface element descriptor %d.\n",
                    old_index);
        return false;
      }
      surface_element.SetIndex(old_to_new[old_index]);
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

  // Explicit metadata only, after curving and final descriptor/BC labeling.
  // No centroid/tolerance inference and no solver-side repair fallback.
  if (has_periodic) {
    try {
      periodic_result = cubit_mesh_export::apply_named_periodic_identification(*ng_mesh, periodic);
    } catch (const std::exception &error) {
      PRINT_ERROR("Named periodic export rejected: %s\n", error.what());
      return false;
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
  if (has_periodic) {
    try {
      cubit_mesh_export::save_named_periodic_vol(*ng_mesh, u8_string_to_path(filename));
    } catch (const std::exception &error) {
      PRINT_ERROR("Named periodic save failed: %s\n", error.what());
      PRINT_ERROR("%s\n", cubit_mesh_export::named_periodic_retry_guidance(
          cubit_mesh_export::named_periodic_output_paths(u8_string_to_path(filename))).c_str());
      return false;
    }
  } else {
    ng_mesh->Save(u8_string_to_path(filename));
  }

  int ne = ng_mesh->GetNE();
  int np = ng_mesh->GetNP();


  // The named export is incomplete without its companion mapping record.
  auto write_companion = [&](std::ostream &jf) {
      jf << std::setprecision(16);
      jf << "{\n";
      jf << "  \"export_passed\": true,\n";
      if (has_periodic) {
        jf << "  \"periodic_identification\": {\"name\":\"kt\",\"nr\":1,\"type\":"
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
  if (has_periodic) {
    try {
      cubit_mesh_export::write_named_periodic_sidecar(
          u8_string_to_path(filename + ".json"), write_companion);
    } catch (const std::exception &error) {
      PRINT_ERROR("Named periodic sidecar failed: %s\n", error.what());
      PRINT_ERROR("%s\n", cubit_mesh_export::named_periodic_retry_guidance(
          cubit_mesh_export::named_periodic_output_paths(u8_string_to_path(filename))).c_str());
      return false;
    }
  } else {
    std::ofstream jf(u8_string_to_path(filename + ".json"));
    if (jf.is_open()) { write_companion(jf); jf.close(); }
  }

  PRINT_INFO("Exported Netgen Vol (order %d): %s (%d nodes, %d elements)\n",
             order, filename.c_str(), np, ne);

  return true;
#endif
}
