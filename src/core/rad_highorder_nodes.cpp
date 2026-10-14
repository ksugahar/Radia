/**
 * High-order node computation for GMSH export.
 * C++ accelerated replacement for Python _compute_highorder_nodes().
 */

#include "rad_highorder_nodes.h"

#include <comp.hpp>    // MeshAccess, ElementId
#include <fem.hpp>     // IntegrationPoint, ElementTransformation

#include <unordered_map>
#include <algorithm>
#include <cmath>

namespace radia {

// GMSH element type codes for Lagrange triangles
static const int GMSH_TRIG_TYPE[] = {
    0,  // p=0 (invalid)
    2,  // p=1: Tri3
    9,  // p=2: Tri6
    21, // p=3: Tri10
    23, // p=4: Tri15
    25, // p=5: Tri21
};

/**
 * Compute equidistant GMSH reference points for Lagrange triangle of order p.
 * Node ordering follows GMSH convention:
 *   corners -> edge 0-1 -> edge 1-2 -> edge 2-0 -> interior (row-by-row)
 */
static std::vector<std::pair<double, double>> GetGmshTrigRefPoints(int p)
{
    std::vector<std::pair<double, double>> pts;
    int n_total = (p + 1) * (p + 2) / 2;
    pts.reserve(n_total);

    // Corners
    pts.push_back({0.0, 0.0});
    pts.push_back({1.0, 0.0});
    pts.push_back({0.0, 1.0});

    // Edge 0->1: (k/p, 0)
    for (int k = 1; k < p; k++)
        pts.push_back({(double)k / p, 0.0});

    // Edge 1->2: (1-k/p, k/p)
    for (int k = 1; k < p; k++)
        pts.push_back({1.0 - (double)k / p, (double)k / p});

    // Edge 2->0: (0, 1-k/p)
    for (int k = 1; k < p; k++)
        pts.push_back({0.0, 1.0 - (double)k / p});

    // Interior: row-by-row equidistant
    for (int j = 1; j < p; j++)
        for (int i = 1; i < p - j; i++)
            pts.push_back({(double)i / p, (double)j / p});

    return pts;
}

/**
 * Evaluate ElementTransformation at a reference point (u, v).
 * Returns the 3D mapped point.
 */
static inline std::array<double, 3> EvalTrafo(
    const ngfem::ElementTransformation& trafo,
    double u, double v)
{
    ngfem::IntegrationPoint ip(u, v, 0.0, 1.0);
    ngfem::MappedIntegrationPoint<2, 3> mip(ip, trafo);
    auto& pt = mip.GetPoint();
    return {pt(0), pt(1), pt(2)};
}

/**
 * Squared distance between two 3D points.
 */
static inline double Dist2(const std::array<double, 3>& a,
                           const std::array<double, 3>& b)
{
    double dx = a[0] - b[0], dy = a[1] - b[1], dz = a[2] - b[2];
    return dx * dx + dy * dy + dz * dz;
}

/**
 * Hash key for edge cache: combine two vertex indices into uint64_t.
 */
static inline uint64_t EdgeKey(int va, int vb)
{
    int lo = std::min(va, vb);
    int hi = std::max(va, vb);
    return ((uint64_t)lo << 32) | (uint64_t)hi;
}

HighOrderBndResult ComputeHighOrderBndNodes(
    std::shared_ptr<ngcomp::MeshAccess> ma,
    int curve_order)
{
    using namespace ngcomp;
    using namespace ngfem;

    HighOrderBndResult result;
    int p = curve_order;

    // Vertex nodes
    size_t nv = ma->GetNV();
    result.n_vertices = nv;
    result.nodes.resize(nv);
    for (size_t i = 0; i < nv; i++) {
        auto pt = ma->GetPoint<3>(i);
        result.nodes[i] = {pt(0), pt(1), pt(2)};
    }

    if (p < 2) {
        // p=1: no high-order nodes, just build connectivity from corners
        auto gmsh_type = GMSH_TRIG_TYPE[1];
        auto bnd_labels = ma->GetMaterials(BND);
        size_t nse = ma->GetNSE();

        for (size_t i = 0; i < nse; i++) {
            ElementId eid(BND, i);
            auto el = ma->GetElement(eid);
            auto et = el.GetType();
            if (et != ET_TRIG) continue;

            auto verts = el.Vertices();
            std::string mat = std::string(bnd_labels[el.GetIndex()]);

            std::vector<int> conn = {(int)verts[0], (int)verts[1], (int)verts[2]};
            result.elem_conn.push_back(std::move(conn));
            result.elem_materials.push_back(mat);
            result.elem_gmsh_types.push_back(gmsh_type);
            result.elem_orig_idx.push_back((int)i);
        }
        return result;
    }

    // High-order (p >= 2)
    auto gmsh_ref = GetGmshTrigRefPoints(p);
    int n_total = (int)gmsh_ref.size();  // (p+1)(p+2)/2
    int n_edge_per = p - 1;
    int gmsh_type = GMSH_TRIG_TYPE[p];

    // Edge cache: EdgeKey -> (start_vertex, vector of node indices)
    std::unordered_map<uint64_t, std::pair<int, std::vector<int>>> edge_cache;

    auto bnd_labels = ma->GetMaterials(BND);
    size_t nse = ma->GetNSE();

    ngstd::LocalHeap lh(1000000, "ho_nodes");

    for (size_t i = 0; i < nse; i++) {
        ngstd::HeapReset hr(lh);
        ElementId eid(BND, i);
        auto el = ma->GetElement(eid);
        auto et = el.GetType();
        if (et != ET_TRIG) continue;

        auto verts = el.Vertices();
        int v0 = verts[0], v1 = verts[1], v2 = verts[2];

        auto& trafo = ma->GetTrafo(eid, lh);

        // Match ref corners to physical vertices
        std::array<double, 3> corner_pts[3];
        for (int ci = 0; ci < 3; ci++)
            corner_pts[ci] = EvalTrafo(trafo, gmsh_ref[ci].first, gmsh_ref[ci].second);

        // Find permutation: ref_corner[ci] -> vertex index
        int vert_arr[3] = {v0, v1, v2};
        std::array<double, 3> vert_pts[3] = {
            result.nodes[v0], result.nodes[v1], result.nodes[v2]
        };

        int ref_to_vert[3];
        for (int ci = 0; ci < 3; ci++) {
            double best_d2 = 1e30;
            int best_vi = 0;
            for (int vi = 0; vi < 3; vi++) {
                double d2 = Dist2(corner_pts[ci], vert_pts[vi]);
                if (d2 < best_d2) { best_d2 = d2; best_vi = vi; }
            }
            ref_to_vert[ci] = vert_arr[best_vi];
        }

        // Build full connectivity: corners + edge nodes + interior nodes
        std::vector<int> conn;
        conn.reserve(n_total);

        // Corners (in GMSH ref order -> mapped to vertex indices)
        conn.push_back(ref_to_vert[0]);
        conn.push_back(ref_to_vert[1]);
        conn.push_back(ref_to_vert[2]);

        // 3 GMSH edges: (ref0,ref1), (ref1,ref2), (ref2,ref0)
        int gmsh_edges[3][2] = {{0, 1}, {1, 2}, {2, 0}};

        for (int ei = 0; ei < 3; ei++) {
            int va = ref_to_vert[gmsh_edges[ei][0]];
            int vb = ref_to_vert[gmsh_edges[ei][1]];
            uint64_t ekey = EdgeKey(va, vb);
            int start_idx = 3 + ei * n_edge_per;

            auto it = edge_cache.find(ekey);
            if (it != edge_cache.end()) {
                // Cached: check direction
                auto& [cached_start, cached_nodes] = it->second;
                if (cached_start == va) {
                    for (int k = 0; k < n_edge_per; k++)
                        conn.push_back(cached_nodes[k]);
                } else {
                    for (int k = n_edge_per - 1; k >= 0; k--)
                        conn.push_back(cached_nodes[k]);
                }
            } else {
                // Compute and cache
                std::vector<int> edge_nodes;
                edge_nodes.reserve(n_edge_per);
                for (int k = 0; k < n_edge_per; k++) {
                    int ri = start_idx + k;
                    auto pt = EvalTrafo(trafo, gmsh_ref[ri].first, gmsh_ref[ri].second);
                    int nidx = (int)result.nodes.size();
                    result.nodes.push_back(pt);
                    edge_nodes.push_back(nidx);
                }
                edge_cache[ekey] = {va, edge_nodes};
                for (int k = 0; k < n_edge_per; k++)
                    conn.push_back(edge_nodes[k]);
            }
        }

        // Interior nodes
        int int_start = 3 + 3 * n_edge_per;
        for (int ri = int_start; ri < n_total; ri++) {
            auto pt = EvalTrafo(trafo, gmsh_ref[ri].first, gmsh_ref[ri].second);
            int nidx = (int)result.nodes.size();
            result.nodes.push_back(pt);
            conn.push_back(nidx);
        }

        std::string mat = std::string(bnd_labels[el.GetIndex()]);
        result.elem_conn.push_back(std::move(conn));
        result.elem_materials.push_back(mat);
        result.elem_gmsh_types.push_back(gmsh_type);
        result.elem_orig_idx.push_back((int)i);
    }

    return result;
}

} // namespace radia
