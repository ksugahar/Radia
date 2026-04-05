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
#ifdef HAVE_NETGEN
class NetgenCurver;
#endif

#include <memory>  // for shared_ptr

// Pull in Cubit's ElementType enum
#include "ElementType.h"

// --- Edge tables for high-order connectivity ---
// Vertex-index pairs for each element type, in Nastran/Gmsh/ExodusII order.
namespace EdgeTables {
  // TET4: 6 edges
  inline constexpr int tet[][2] = {
    {0,1},{1,2},{0,2},{0,3},{1,3},{2,3}
  };
  // HEX8: 12 edges
  inline constexpr int hex[][2] = {
    {0,1},{1,2},{2,3},{3,0},
    {4,5},{5,6},{6,7},{7,4},
    {0,4},{1,5},{2,6},{3,7}
  };
  // WEDGE6: 9 edges
  inline constexpr int wedge[][2] = {
    {0,1},{1,2},{2,0},
    {3,4},{4,5},{5,3},
    {0,3},{1,4},{2,5}
  };
  // PYRAMID5: 8 edges
  inline constexpr int pyramid[][2] = {
    {0,1},{1,2},{2,3},{3,0},
    {0,4},{1,4},{2,4},{3,4}
  };
  // TRI3: 3 edges
  inline constexpr int tri[][2] = {
    {0,1},{1,2},{2,0}
  };
  // QUAD4: 4 edges
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
  std::vector<int> ho_conn; // high-order full connectivity (empty if order 1)
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
  bool has_netgen = false;  // true if NetgenCurver was used

#ifdef HAVE_NETGEN
  // Access the NetgenCurver used during extraction (null if order < 2 or failed)
  std::shared_ptr<NetgenCurver> get_netgen_curver() const { return netgen_curver_; }
#endif

  // Nodes (original + HO generated)
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
  void extract_sidesets(MeshExportInterface *iface);
  void extract_nodesets(MeshExportInterface *iface);

  // Apply HO connectivity to elements and sideset faces
#ifdef HAVE_NETGEN
  void apply_ho(const NetgenCurver &nc);

  // Build ho_conn for a single element using NetgenCurver edge nodes
  static std::vector<int> build_ho_conn_nc(
      const std::vector<int> &conn, int nv, ElementType type,
      const NetgenCurver &nc);

  // Preserved NetgenCurver (keeps ng_mesh_ alive for pybind11 export)
  std::shared_ptr<NetgenCurver> netgen_curver_;
#endif
};

#endif // MESH_DATA_HPP
