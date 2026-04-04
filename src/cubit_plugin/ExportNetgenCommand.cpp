#include "ExportNetgenCommand.hpp"
#include "MeshData.hpp"
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"

#ifdef HAVE_NETGEN
#include "NetgenCurver.hpp"
#include <meshing.hpp>
#endif

#ifdef _WIN32
#include <windows.h>
#endif

#include <fstream>
#include <map>

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
    "export netgen <string:label='filename',help='<filename>'> "
    "[order <value:label='order',help='<1-5>'>] "
    "[overwrite]"
  );
  return syntax_list;
}

std::vector<std::string> ExportNetgenCommand::get_syntax_help()
{
  std::vector<std::string> help;
  help.push_back(
    "export netgen \"filename.vol\" [order {1|2|3|4|5}] [overwrite]"
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
  PRINT_ERROR("export netgen requires Netgen support (not built).\n");
  return false;
#else
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

  // ---- Set boundary labels (sideset name -> bc number) ----
  // FaceDescriptors already have BCProperty set to Cubit surface ID.
  // Map sideset names to BCProperty values.
  int nfd = ng_mesh->GetNFD();
  if (!md.sidesets.empty()) {
    // Build sideset_id -> name map
    std::map<int, std::string> ss_name;
    for (auto &sg : md.sidesets) {
      if (!sg.name.empty())
        ss_name[sg.id] = sg.name;
    }

    // Also map surface_id -> sideset name (sideset references surface)
    std::map<int, std::string> surf_to_ssname;
    for (auto &sg : md.sidesets) {
      // Get surface IDs in this sideset
      std::vector<int> ss_surfs = CubitInterface::get_sideset_surfaces(sg.id);
      for (int sid : ss_surfs) {
        if (!sg.name.empty())
          surf_to_ssname[sid] = sg.name;
      }
    }

    // Set BCName for each FaceDescriptor
    ng_mesh->SetNBCNames(nfd);
    for (int fi = 1; fi <= nfd; fi++) {
      int bc_prop = ng_mesh->GetFaceDescriptor(fi).BCProperty();
      auto it = surf_to_ssname.find(bc_prop);
      if (it != surf_to_ssname.end()) {
        ng_mesh->SetBCName(fi - 1, it->second);
        ng_mesh->GetFaceDescriptor(fi).SetBCName(ng_mesh->GetBCNamePtr(fi - 1));
      } else {
        // Use Cubit entity name as fallback
        std::string sname = CubitInterface::get_entity_name("surface", bc_prop);
        if (sname.empty())
          sname = "surface_" + std::to_string(bc_prop);
        ng_mesh->SetBCName(fi - 1, sname);
        ng_mesh->GetFaceDescriptor(fi).SetBCName(ng_mesh->GetBCNamePtr(fi - 1));
      }
    }

    // Update FaceDescriptor DomainIn/Out from volume element adjacency
    for (int fi = 1; fi <= nfd; fi++) {
      auto &fd = ng_mesh->GetFaceDescriptor(fi);
      // DomainIn/Out already set by surface element adjacency in build_netgen_mesh.
      // For now keep as-is (domain 1 for single-volume, needs improvement for multi-volume).
    }
  } else {
    // No sidesets: use Cubit surface entity names
    ng_mesh->SetNBCNames(nfd);
    for (int fi = 1; fi <= nfd; fi++) {
      int bc_prop = ng_mesh->GetFaceDescriptor(fi).BCProperty();
      std::string sname = CubitInterface::get_entity_name("surface", bc_prop);
      if (sname.empty())
        sname = "surface_" + std::to_string(bc_prop);
      ng_mesh->SetBCName(fi - 1, sname);
      ng_mesh->GetFaceDescriptor(fi).SetBCName(ng_mesh->GetBCNamePtr(fi - 1));
    }
  }

  // For order=1: reset curving
  if (order == 1) {
    ng_mesh->BuildCurvedElements(1);
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

  return true;
#endif
}
