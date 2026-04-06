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
// Internal order matches VTK convention (used by build_ho_conn_nc).
// GMSH and Nastran use DIFFERENT orders — reorder tables below.
//
// Internal (= VTK) edge order:
//   HEX:   bottom(4) + top(4) + vertical(4)
//   WEDGE: bottom_tri(3) + top_tri(3) + vertical(3)
// Nastran edge order:
//   CHEXA: bottom(4) + vertical(4) + top(4)
//   CPENTA: bottom_tri(3) + vertical(3) + top_tri(3)
// GMSH edge order: vertex-grouped (see gmsh.info Section 10.2)
//
// Reorder semantics: out[nv + reorder[i]] = conn[nv + i]
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

  // --- GMSH reorder: internal pos i -> GMSH file pos reorder[i] ---
  // Verified: gmsh.info Section 10.2 + gmsh.model.mesh.getElementProperties()
  inline constexpr int tet_gmsh_reorder[6] = {0,1,2,3,5,4};
  inline constexpr int hex_gmsh_reorder[12] = {0,3,5,1,8,10,11,9,2,4,6,7};
  inline constexpr int wedge_gmsh_reorder[9] = {0,3,1,6,8,7,2,4,5};
  inline constexpr int pyramid_gmsh_reorder[8] = {0,3,5,1,2,4,6,7};
  inline constexpr int tri_gmsh_reorder[3] = {0,1,2};
  inline constexpr int quad_gmsh_reorder[4] = {0,1,2,3};

  // --- Nastran BDF reorder: internal pos i -> Nastran file pos reorder[i] ---
  // Derived from GMSH source (MHexahedron.h getVertexBDF, MPrism.h getVertexBDF)
  // Nastran CHEXA: bottom(4) + vertical(4) + top(4) (swaps top/vertical vs internal)
  // Nastran CPENTA: bottom(3) + vertical(3) + top(3) (swaps top/vertical vs internal)
  inline constexpr int tet_bdf_reorder[6] = {0,1,2,3,4,5};  // identity
  inline constexpr int hex_bdf_reorder[12] = {0,1,2,3,8,9,10,11,4,5,6,7};
  inline constexpr int wedge_bdf_reorder[9] = {0,1,2,6,7,8,3,4,5};
  inline constexpr int pyramid_bdf_reorder[8] = {0,1,2,3,4,5,6,7};  // identity
  inline constexpr int tri_bdf_reorder[3] = {0,1,2};  // identity
  inline constexpr int quad_bdf_reorder[4] = {0,1,2,3};  // identity

  // --- VTK reorder: internal pos i -> VTK file pos reorder[i] ---
  // Derived from GMSH source (GModelIO_VTK.cpp reader + MHexahedron.h getVertexVTK)
  inline constexpr int tet_vtk_reorder[6] = {0,1,2,3,4,5};  // identity
  inline constexpr int hex_vtk_reorder[12] = {0,1,10,3,2,6,7,4,5,8,11,9};
  inline constexpr int wedge_vtk_reorder[9] = {0,6,3,2,5,4,1,8,7};
  inline constexpr int pyramid_vtk_reorder[8] = {0,1,2,3,4,5,6,7};  // identity (GMSH VTK reader doesn't handle type 27)
  inline constexpr int tri_vtk_reorder[3] = {0,1,2};  // identity
  inline constexpr int quad_vtk_reorder[4] = {0,1,2,3};  // identity
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
