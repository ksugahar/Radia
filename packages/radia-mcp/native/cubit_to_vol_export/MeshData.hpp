#ifndef MESH_DATA_HPP
#define MESH_DATA_HPP

// ============================================================
// MeshData
//
// Shared mesh data extraction from Cubit's MeshExportInterface.
// Extracts blocks, sidesets, nodesets in one pass.  High-order
// node generation (NetgenCurver) is applied uniformly to both
// block elements and sideset faces.
//
// Each exporter reads from MeshData and writes format-specific
// output; no exporter touches MeshExportInterface directly.
// ============================================================

#include <vector>
#include <array>
#include <string>
#include <unordered_map>

class MeshExportInterface;
class NetgenCurver;

#include <memory>  // for shared_ptr

// Pull in Cubit's ElementType enum
#include "ElementType.h"

// --- Edge tables for high-order connectivity ---
// Internal order matches VTK convention (used by build_ho_conn_nc).
//
// Internal (= VTK) edge order:
//   HEX:   bottom(4) + top(4) + vertical(4)
//   WEDGE: bottom_tri(3) + top_tri(3) + vertical(3)
//
namespace EdgeTables {
  // --- Internal (= VTK) edge tables ---
  inline constexpr int tet[][2] = {
    {0,1},{1,2},{0,2},{0,3},{1,3},{2,3}
  };
  inline constexpr int hex[][2] = {
    {0,1},{1,2},{2,3},{3,0},     // bottom face
    {4,5},{5,6},{6,7},{7,4},     // top face
    {0,4},{1,5},{2,6},{3,7}      // vertical
  };
  inline constexpr int wedge[][2] = {
    {0,1},{1,2},{2,0},           // bottom tri
    {3,4},{4,5},{5,3},           // top tri
    {0,3},{1,4},{2,5}            // vertical
  };
  inline constexpr int pyramid[][2] = {
    {0,1},{1,2},{2,3},{3,0},
    {0,4},{1,4},{2,4},{3,4}
  };
  inline constexpr int tri[][2] = {
    {0,1},{1,2},{2,0}
  };
  inline constexpr int quad[][2] = {
    {0,1},{1,2},{2,3},{3,0}
  };

}

// --- Data types ---

struct NodeCoord {
  int    id;
  double x, y, z;
};

struct MeshElement {
  int         group_id;     // block ID or sideset ID
  ElementType type;         // Cubit element type (TETRA4, HEX8, TRI3, etc.)
  int         nv;           // vertex count (linear)
  std::vector<int> conn;    // linear vertex node IDs
};

struct SidesetGroup {
  int         id;
  std::string name;
  std::vector<MeshElement> faces;  // extracted face elements (TRI3/QUAD4 etc.)
};

struct NodesetGroup {
  int         id;
  std::string name;
  std::vector<int> node_ids;
};

// --- Main class ---

class MeshData {
public:
  // Configuration
  int  order = 1;
  // Fail extract() when no block holds a volume element (volume-only formats)
  bool require_volume_elements = false;

  // Access the NetgenCurver used during extraction (null if order < 2 or failed)
  std::shared_ptr<NetgenCurver> get_netgen_curver() const { return netgen_curver_; }

  // Nodes (linear Cubit mesh)
  std::vector<NodeCoord> nodes;
  int num_original_nodes  = 0;
  int max_original_node_id = 0;

  // Block elements
  std::vector<MeshElement> elements;
  std::vector<int> block_ids;  // unique block IDs in encounter order

  // Sidesets
  std::vector<SidesetGroup> sidesets;

  // Nodesets
  std::vector<NodesetGroup> nodesets;

  // Lookup: Cubit node ID -> 0-based index into nodes[]
  std::unordered_map<int, int> node_id_to_index;

  // Get node coordinates by Cubit node ID
  std::array<double,3> get_node_coords(int node_id) const {
    auto it = node_id_to_index.find(node_id);
    if (it != node_id_to_index.end()) {
      auto &n = nodes[it->second];
      return {n.x, n.y, n.z};
    }
    return {0, 0, 0};
  }

  // ---- Main entry point ----
  // Acquires MeshExportInterface, extracts everything, releases it.
  // Returns false on failure.
  bool extract(int order = 1);

  // ---- Convenience ----
  int total_node_count() const { return (int)nodes.size(); }
  int total_element_count() const;  // block elements + all sideset faces

private:
  void extract_nodes(MeshExportInterface *iface);
  void extract_elements(MeshExportInterface *iface);
  bool check_block_coverage() const;
  void extract_sidesets(MeshExportInterface *iface);
  void extract_nodesets(MeshExportInterface *iface);

  // Preserved NetgenCurver (keeps ng_mesh_ alive for the .vol writer)
  std::shared_ptr<NetgenCurver> netgen_curver_;
};

#endif // MESH_DATA_HPP
