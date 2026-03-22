#ifndef RAD_HIGHORDER_NODES_H
#define RAD_HIGHORDER_NODES_H

/**
 * High-order node computation for GMSH export (C++ accelerated).
 *
 * Replaces the Python loop in gmsh_post_export.py that calls
 * mesh.GetTrafo(el) for each BND element. For p=4 with 2790 elements,
 * Python takes ~288 seconds; C++ takes <1 second.
 *
 * GMSH Lagrange triangle node ordering:
 *   - 3 corners: (0,0), (1,0), (0,1)
 *   - 3*(p-1) edge nodes: equidistant on edges 0-1, 1-2, 2-0
 *   - (p-1)*(p-2)/2 interior nodes: row-by-row equidistant
 *   Total: (p+1)*(p+2)/2 nodes per element
 */

#include <vector>
#include <array>
#include <string>
#include <memory>

namespace ngcomp { class MeshAccess; }

namespace radia {

struct HighOrderBndResult {
    // All node positions: [0..n_vertices-1] = mesh vertices,
    // [n_vertices..] = high-order nodes (edge + interior)
    std::vector<std::array<double, 3>> nodes;

    // Per-element: GMSH connectivity (all node indices in GMSH order)
    // elem_conn[i] has (p+1)*(p+2)/2 entries for TRIG
    std::vector<std::vector<int>> elem_conn;

    // Per-element: material name (from mesh boundary label)
    std::vector<std::string> elem_materials;

    // Per-element: GMSH element type code
    std::vector<int> elem_gmsh_types;

    // Per-element: original NGSolve element index
    std::vector<int> elem_orig_idx;

    // Number of vertex nodes
    size_t n_vertices;
};

/**
 * Compute high-order BND nodes for GMSH export.
 *
 * @param ma  NGSolve MeshAccess (must have Curve(p) applied)
 * @param curve_order  Polynomial order (1..5)
 * @return  HighOrderBndResult with nodes, connectivity, materials
 */
HighOrderBndResult ComputeHighOrderBndNodes(
    std::shared_ptr<ngcomp::MeshAccess> ma,
    int curve_order);

} // namespace radia

#endif // RAD_HIGHORDER_NODES_H
