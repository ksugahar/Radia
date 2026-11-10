// Neutral, linear volume-mesh acquisition. No Python or downstream runtime.
#include "CubitCommandInterface.hpp"
#include "CubitPluginExport.hpp"
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"
#include "RefEntityFactory.hpp"
#include "CGMApp.hpp"
#include "RefFace.hpp"
#include "RefVolume.hpp"
#include <algorithm>
#include <array>
#include <cmath>
#include <filesystem>
#include <fstream>
#include <iomanip>
#include <map>
#include <set>
#include <sstream>
#include <stdexcept>
#include <string>
#include <tuple>
#include <vector>
#ifdef _WIN32
#include <windows.h>
#endif

namespace {
namespace fs = std::filesystem;
std::string quoted(const std::string& text) {
  std::ostringstream out;
  out << '"';
  for (unsigned char c : text) {
    if (c == '"' || c == '\\') out << '\\' << c;
    else if (c < 32) out << "\\u" << std::hex << std::setw(4) << std::setfill('0') << int(c) << std::dec;
    else out << c;
  }
  out << '"';
  return out.str();
}
template<class T> void numbers(std::ostream& out, const T& values) {
  out << '['; bool first = true;
  for (auto v : values) { if (!first) out << ','; first = false; out << v; }
  out << ']';
}
void require(bool value, const std::string& message) {
  if (!value) throw std::runtime_error(message);
}
struct Cell { std::string type; int id, region; std::vector<int> nodes; };
struct Face { int surface, id, volume; std::string sense; std::vector<int> nodes; };
struct Region { int id; std::string name; };
struct SideSet { int id; std::string name; std::vector<Face> faces; std::set<int> nodes; };
struct NodeSet { int id; std::string name; std::vector<int> nodes; };
using Point = std::array<double,3>;
struct Capture {
  std::map<int,Point> nodes;
  std::vector<Cell> cells;
  std::vector<Region> regions;
  std::vector<Face> faces;
  std::vector<SideSet> sidesets;
  std::vector<NodeSet> nodesets;
};
std::string module_path() {
#ifdef _WIN32
  HMODULE module = nullptr;
  require(GetModuleHandleExW(GET_MODULE_HANDLE_EX_FLAG_FROM_ADDRESS | GET_MODULE_HANDLE_EX_FLAG_UNCHANGED_REFCOUNT,
    reinterpret_cast<LPCWSTR>(&module_path), &module) != 0, "Cannot identify exporter module");
  std::vector<wchar_t> buffer(32768);
  DWORD length = GetModuleFileNameW(module, buffer.data(), DWORD(buffer.size()));
  require(length > 0 && length < buffer.size(), "Exporter module path unavailable/truncated");
  return fs::path(std::wstring(buffer.data(),length)).u8string();
#else
  throw std::runtime_error("This acquisition module requires qualified Windows module identity");
#endif
}
Capture acquire(bool labels_only) {
  Capture result;
  require(labels_only || (CubitInterface::get_sideset_id_list().empty() && CubitInterface::get_nodeset_id_list().empty()),
          "Boundary sets require a separately qualified marker profile; none may be silently dropped");
  std::set<std::pair<std::string,int>> assigned;
  for (int bid : CubitInterface::get_block_id_list()) {
    require(bid > 0, "Nonpositive block ID");
    auto type = CubitInterface::get_block_element_type(bid);
    require(type == "HEX8" || type == "TETRA4" || type == "TET4", "Only explicit HEX8/TETRA4 blocks supported: " + type);
    std::string kind = type == "HEX8" ? "hex" : "tet";
    // Expand geometric block membership through Cubit's native list parser.
    // Direct get_block_hexes/tets queries may omit a block's volume members.
    auto ids = CubitInterface::parse_cubit_list(kind, "in block " + std::to_string(bid));
    require(!ids.empty(), "Empty block " + std::to_string(bid));
    result.regions.push_back({bid, CubitInterface::get_exodus_entity_name("block",bid)});
    for (int id : ids) {
      require(id > 0 && assigned.insert({kind,id}).second, "Element belongs to multiple blocks");
      auto nodes = CubitInterface::get_expanded_connectivity(kind,id);
      require(nodes.size() == (kind == "hex" ? 8u : 4u), "Higher-order or unsupported cell; no order reduction");
      require(std::set<int>(nodes.begin(),nodes.end()).size() == nodes.size(), "Repeated cell node");
      result.cells.push_back({kind == "hex" ? "HEX8" : "TET4",id,bid,nodes});
      for (int nid : nodes) {
        require(nid > 0, "Nonpositive node ID");
        auto p = CubitInterface::get_nodal_coordinates(nid);
        require(std::all_of(p.begin(),p.end(),[](double x){return std::isfinite(x);}), "Nonfinite coordinate");
        result.nodes[nid] = p;
      }
    }
  }
  std::set<std::pair<std::string,int>> all;
  for (const std::string kind : {"hex","tet","wedge","pyramid"})
    for (int id : CubitInterface::get_entities(kind)) all.insert({kind,id});
  require(!all.empty() && assigned == all, "Every volume element must belong to exactly one supported block");
  for (int sid : CubitInterface::get_entities("surface")) {
    auto* app = CubitInterface::cgm();
    require(app != nullptr && app->factory() != nullptr, "Cubit geometry factory unavailable");
    auto* face = app->factory()->get_ref_face(sid);
    require(face != nullptr, "Surface lookup failed");
    DLIList<RefVolume*> volumes;
    face->ref_volumes(volumes);
    if (volumes.size() != 1) continue;
    auto* volume = volumes.get();
    auto sense = face->sense(volume);
    require(sense == CUBIT_FORWARD || sense == CUBIT_REVERSED, "Unknown exterior surface/volume sense");
    for (const std::string kind : {"tri","quad"}) {
      auto ids = kind == "tri" ? CubitInterface::get_surface_tris(sid) : CubitInterface::get_surface_quads(sid);
      for (int id : ids) {
        auto nodes = CubitInterface::get_expanded_connectivity(kind,id);
        require(nodes.size() == (kind == "tri" ? 3u : 4u), "High-order skin unsupported");
        for (int n : nodes) require(result.nodes.count(n) != 0, "Skin node outside captured volume mesh");
        Point center{0,0,0};
        for (int n : nodes) for (int k=0;k<3;++k) center[k] += result.nodes.at(n)[k]/double(nodes.size());
        CubitVector normal = face->normal_at(CubitVector(center[0],center[1],center[2]),volume);
        require(std::isfinite(normal.length_squared()) && normal.length_squared() > 0, "Outward CAD normal unavailable");
        auto a = result.nodes.at(nodes[0]), b = result.nodes.at(nodes[1]), c = result.nodes.at(nodes[2]);
        CubitVector ab(b[0]-a[0],b[1]-a[1],b[2]-a[2]), ac(c[0]-a[0],c[1]-a[1],c[2]-a[2]);
        double dot = (ab * ac) % normal;
        require(std::isfinite(dot) && dot != 0, "Degenerate/ambiguous skin orientation");
        if (dot < 0) std::reverse(nodes.begin(),nodes.end());
        result.faces.push_back({sid,id,volume->id(),sense == CUBIT_FORWARD ? "Forward" : "Reversed",nodes});
      }
    }
  }
  require(!result.faces.empty(), "Geometry-backed independent outward skin required");
  for (int id : CubitInterface::get_sideset_id_list()) {
    std::vector<int> direct_faces, surfaces, curves;
    CubitInterface::get_sideset_children(id,direct_faces,surfaces,curves);
    require(curves.empty(), "3D label capture does not support curve sidesets");
    SideSet set{id,CubitInterface::get_exodus_entity_name("sideset",id),{}, {}};
    require(!set.name.empty(),"An explicit sideset name is required");
    std::set<int> requested_faces(direct_faces.begin(),direct_faces.end()), requested_surfaces(surfaces.begin(),surfaces.end());
    std::set<int> matched_faces, matched_surfaces;
    for (const auto& face : result.faces) {
      if (!requested_faces.count(face.id) && !requested_surfaces.count(face.surface)) continue;
      set.faces.push_back(face);
      set.nodes.insert(face.nodes.begin(),face.nodes.end());
      if (requested_faces.count(face.id)) matched_faces.insert(face.id);
      if (requested_surfaces.count(face.surface)) matched_surfaces.insert(face.surface);
    }
    require(!set.faces.empty() && matched_faces==requested_faces && matched_surfaces==requested_surfaces,
            "Sideset contains unsupported/internal/unmeshed entities: " + std::to_string(id));
    result.sidesets.push_back(std::move(set));
  }
  for (int id : CubitInterface::get_nodeset_id_list()) {
    NodeSet set{id,CubitInterface::get_exodus_entity_name("nodeset",id),
                CubitInterface::parse_cubit_list("node","in nodeset " + std::to_string(id))};
    require(!set.name.empty() && !set.nodes.empty(),"Named nonempty nodeset required");
    std::sort(set.nodes.begin(),set.nodes.end());
    require(std::adjacent_find(set.nodes.begin(),set.nodes.end())==set.nodes.end(),"Repeated nodeset node");
    for(int node:set.nodes) require(result.nodes.count(node)!=0,"Nodeset node outside volume mesh");
    result.nodesets.push_back(std::move(set));
  }
  std::sort(result.cells.begin(),result.cells.end(),[](const Cell& a,const Cell& b){return std::tie(a.type,a.id)<std::tie(b.type,b.id);});
  return result;
}
void write_capture(const Capture& mesh, const fs::path& directory, const std::string& unit) {
  std::ofstream out(directory/"native-capture.json",std::ios::binary);
  require(bool(out),"Cannot create capture");
  out << std::setprecision(17) << "{\"schema\":\"cubit.native-mesh-capture.v1\",\"units\":" << quoted(unit)
      << ",\"cubit_version\":" << quoted(CubitInterface::get_version()) << ",\"exporter_module_path\":" << quoted(module_path()) << ",\"nodes\":[";
  bool first = true;
  for (const auto& node : mesh.nodes) { if (!first) out << ','; first=false; out << "{\"id\":" << node.first << ",\"xyz\":";numbers(out,node.second);out << '}'; }
  out << "],\"elements\":["; first=true;
  for (const auto& cell : mesh.cells) { if(!first)out << ',';first=false;out << "{\"source_id\":" << cell.id << ",\"type\":" << quoted(cell.type) << ",\"region_id\":" << cell.region << ",\"node_ids\":";numbers(out,cell.nodes);out << '}'; }
  out << "],\"regions\":[";first=true;
  for(const auto& region:mesh.regions){if(!first)out << ',';first=false;out << "{\"id\":" << region.id << ",\"name\":" << quoted(region.name) << '}';}
  out << "],\"faces\":[";first=true;
  for(const auto& face:mesh.faces){if(!first)out << ',';first=false;out << "{\"source_surface_id\":" << face.surface << ",\"source_face_id\":" << face.id << ",\"source_volume_id\":" << face.volume << ",\"source_volume_sense\":" << quoted(face.sense) << ",\"node_ids\":";numbers(out,face.nodes);out << '}';}
  out << "],\"sets\":{\"sidesets\":[";first=true;
  for(const auto& set:mesh.sidesets){
    if(!first)out << ',';first=false;
    out << "{\"id\":" << set.id << ",\"name\":" << quoted(set.name) << ",\"kind\":\"face\",\"edges\":[],\"node_ids\":";
    numbers(out,set.nodes);out << ",\"faces\":[";bool first_face=true;
    for(const auto& f:set.faces){if(!first_face)out << ',';first_face=false;out << "{\"source_surface_id\":" << f.surface << ",\"source_face_id\":" << f.id << ",\"source_volume_id\":" << f.volume << ",\"source_volume_sense\":" << quoted(f.sense) << ",\"node_ids\":";numbers(out,f.nodes);out << '}';}
    out << "]}";
  }
  out << "],\"nodesets\":[";first=true;
  for(const auto& set:mesh.nodesets){if(!first)out << ',';first=false;out << "{\"id\":" << set.id << ",\"name\":" << quoted(set.name) << ",\"node_ids\":";numbers(out,set.nodes);out << '}';}
  out << "]}}\n";out.close();require(bool(out),"Capture write failed");
}
class ExportNeutralMesh : public CubitCommand {
public:
  std::vector<std::string> get_syntax() override { return {"export neutral_mesh <string:label='directory'> units <string:label='units'> [labels_only]"}; }
  std::vector<std::string> get_syntax_help() override { return {"export neutral_mesh \"new-directory\" units \"m|mm\" [labels_only]"}; }
  std::vector<std::string> get_help() override { return {"Acquire linear volume mesh, independent oriented CAD skin and native alternate export; external bundle validation required. No physical conditions or downstream operations."}; }
  bool execute(CubitCommandData& data) override {
    try {
      std::string directory,unit;
      require(data.get_string("directory",directory) && data.get_string("units",unit), "Directory and length unit required");
      require(unit=="m" || unit=="mm", "Explicit length unit must be m or mm");
      require(directory.find_first_of("\"\r\n") == std::string::npos, "Unsafe output path");
      auto mesh=acquire(data.find_keyword("labels_only")); auto path=fs::u8path(directory);
      require(!fs::exists(path),"Output directory exists; no overwrite");
      require(fs::create_directory(path),"Cannot exclusively create output directory");
      write_capture(mesh,path,unit);
      std::string command="export mesh \""+(path/"alternate-mesh.e").generic_u8string()+"\" dimension 3";
      const int before=CubitInterface::get_error_count();
      require(CubitInterface::cmd(command.c_str()) && CubitInterface::get_error_count()==before,"Native alternate export failed");
      require(fs::is_regular_file(path/"alternate-mesh.e") && fs::file_size(path/"alternate-mesh.e")>0,"Missing native alternate export");
      PRINT_INFO("Neutral capture written; independent validation required before use.\n");
      return true;
    } catch (const std::exception& e) { PRINT_ERROR("Neutral mesh export failed: %s\n",e.what());return false; }
      catch (...) { PRINT_ERROR("Neutral mesh export failed: unknown exception\n");return false; }
  }
};
}
class NeutralMeshPlugin : public CubitCommandInterface {
public:
  std::vector<std::string> get_keys() override { return {"ExportNeutralMesh"}; }
  CubitCommand* create_command(const std::string& key) override { return key=="ExportNeutralMesh" ? new ExportNeutralMesh() : nullptr; }
};
CUBIT_PLUGIN(NeutralMeshPlugin)
