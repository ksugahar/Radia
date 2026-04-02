/* ===================================================================
 * NetgenCurverPure — Cubit-free implementation
 *
 * Phase 1: Build netgen::Mesh from passed-in data (no Cubit API)
 * Phase 2: Attach CallbackGeometry with Python callbacks
 * Phase 3: BuildCurvedElements (Netgen C++)
 * =================================================================== */

#include "NetgenCurverPure.hpp"

#include <meshing.hpp>
#include "callbackgeom.hpp"

#include <cstdio>
#include <stdexcept>

namespace ng = netgen;

// ============================================================
// build() — main entry point
// ============================================================
std::shared_ptr<ng::Mesh> NetgenCurverPure::build(
    const std::vector<double>& node_coords,
    const std::vector<int>&    node_ids,
    const std::vector<VolumeElement>& vol_elems,
    const std::vector<SurfaceInfo>& surfaces,
    ProjectFunc project_func,
    NormalFunc  normal_func,
    int order,
    EdgeProjectFunc edge_project_func,
    const std::vector<EdgeInfo>& edge_elements)
{
  if (order < 2)
    throw std::runtime_error("NetgenCurverPure: order must be >= 2");

  if (!build_netgen_mesh(node_coords, node_ids, vol_elems, surfaces, edge_elements))
    throw std::runtime_error("NetgenCurverPure: failed to build Netgen mesh");

  if (!attach_callback_geometry(project_func, normal_func, (int)surfaces.size(),
                                 edge_project_func))
    throw std::runtime_error("NetgenCurverPure: failed to attach geometry");

  if (!curve(order))
    throw std::runtime_error("NetgenCurverPure: BuildCurvedElements failed");

  return ng_mesh_;
}

// ============================================================
// Phase 1: Build netgen::Mesh from external data
// ============================================================
bool NetgenCurverPure::build_netgen_mesh(
    const std::vector<double>& node_coords,
    const std::vector<int>&    node_ids,
    const std::vector<VolumeElement>& vol_elems,
    const std::vector<SurfaceInfo>& surfaces,
    const std::vector<EdgeInfo>& edge_elements)
{
  ng_mesh_ = std::make_shared<ng::Mesh>();
  ng_mesh_->SetDimension(3);

  int num_nodes = (int)node_ids.size();
  if ((int)node_coords.size() != num_nodes * 3)
    return false;

  // Add nodes
  for (int i = 0; i < num_nodes; i++) {
    ng::PointIndex pi = ng_mesh_->AddPoint(
        ng::MeshPoint(ng::Point<3>(node_coords[3*i],
                                    node_coords[3*i+1],
                                    node_coords[3*i+2])));
    cubit_nid_to_ng_pi_[node_ids[i]] = pi;
  }

  // Build UV lookup: (surface_fd_index, cubit_node_id) -> (u, v)
  std::unordered_map<int, std::unordered_map<int, std::array<double,2>>> uv_map;

  // Build FaceDescriptors and add surface elements
  int fd_idx = 1;
  for (auto& surf : surfaces) {
    ng::FaceDescriptor fd(fd_idx, 1, 0, fd_idx);
    fd.SetBCProperty(surf.id);
    ng_mesh_->AddFaceDescriptor(fd);

    // UV map for this surface
    for (auto& [nid, uv] : surf.uvs) {
      uv_map[fd_idx][nid] = uv;
    }

    // Add triangles
    for (auto& tri : surf.tris) {
      ng::Element2d el(ng::TRIG);
      el.SetIndex(fd_idx);
      for (int k = 0; k < 3; k++) {
        auto it = cubit_nid_to_ng_pi_.find(tri[k]);
        if (it == cubit_nid_to_ng_pi_.end()) return false;
        el[k] = it->second;

        // Set GeomInfo (UV)
        ng::PointGeomInfo gi;
        gi.trignum = fd_idx;
        auto uv_it = uv_map[fd_idx].find(tri[k]);
        if (uv_it != uv_map[fd_idx].end()) {
          gi.u = uv_it->second[0];
          gi.v = uv_it->second[1];
        } else {
          gi.u = 0; gi.v = 0;
        }
        el.GeomInfo()[k] = gi;
      }
      ng_mesh_->AddSurfaceElement(el);
    }

    // Add quads
    for (auto& quad : surf.quads) {
      ng::Element2d el(ng::QUAD);
      el.SetIndex(fd_idx);
      for (int k = 0; k < 4; k++) {
        auto it = cubit_nid_to_ng_pi_.find(quad[k]);
        if (it == cubit_nid_to_ng_pi_.end()) return false;
        el[k] = it->second;

        ng::PointGeomInfo gi;
        gi.trignum = fd_idx;
        auto uv_it = uv_map[fd_idx].find(quad[k]);
        if (uv_it != uv_map[fd_idx].end()) {
          gi.u = uv_it->second[0];
          gi.v = uv_it->second[1];
        } else {
          gi.u = 0; gi.v = 0;
        }
        el.GeomInfo()[k] = gi;
      }
      ng_mesh_->AddSurfaceElement(el);
    }

    fd_idx++;
  }

  // Add volume elements
  for (auto& ve : vol_elems) {
    switch (ve.type) {
      case PURE_TET: {
        if (ve.conn.size() < 4) continue;
        ng::Element el(ng::TET);
        for (int k = 0; k < 4; k++)
          el[k] = cubit_nid_to_ng_pi_[ve.conn[k]];
        el.SetIndex(1);
        ng_mesh_->AddVolumeElement(el);
        break;
      }
      case PURE_HEX: {
        if (ve.conn.size() < 8) continue;
        ng::Element el(ng::HEX);
        for (int k = 0; k < 8; k++)
          el[k] = cubit_nid_to_ng_pi_[ve.conn[k]];
        el.SetIndex(1);
        ng_mesh_->AddVolumeElement(el);
        break;
      }
      case PURE_PRISM: {
        if (ve.conn.size() < 6) continue;
        ng::Element el(ng::PRISM);
        for (int k = 0; k < 6; k++)
          el[k] = cubit_nid_to_ng_pi_[ve.conn[k]];
        el.SetIndex(1);
        ng_mesh_->AddVolumeElement(el);
        break;
      }
      case PURE_PYRAMID: {
        if (ve.conn.size() < 5) continue;
        ng::Element el(ng::PYRAMID);
        for (int k = 0; k < 5; k++)
          el[k] = cubit_nid_to_ng_pi_[ve.conn[k]];
        el.SetIndex(1);
        ng_mesh_->AddVolumeElement(el);
        break;
      }
    }
  }

  // Add 1D segment elements on geometry edges (curves).
  // Required for BuildCurvedElements to call PointBetweenEdge for edge snapping.
  //
  // Key fields in ng::Segment:
  //   [0], [1]     : PointIndex of endpoints
  //   si           : FaceDescriptor index (1-based) of one adjacent surface
  //   edgenr       : unique edge identifier (1-based, same for all segs on one curve)
  //   surfnr1/2    : 0-based indices of the two adjacent surfaces
  //   epgeominfo[] : EdgePointGeomInfo with edgenr and u/v (set to 0 for now)
  //
  int edgenr = 1;
  for (auto& ei : edge_elements) {
    for (auto& seg : ei.segments) {
      auto it0 = cubit_nid_to_ng_pi_.find(seg[0]);
      auto it1 = cubit_nid_to_ng_pi_.find(seg[1]);
      if (it0 == cubit_nid_to_ng_pi_.end() || it1 == cubit_nid_to_ng_pi_.end())
        continue;

      ng::Segment nseg;
      nseg[0] = it0->second;
      nseg[1] = it1->second;
      nseg.si = ei.surfnr1;              // 1-based FD index
      nseg.edgenr = edgenr;
      nseg.surfnr1 = ei.surfnr1 - 1;    // 0-based
      nseg.surfnr2 = ei.surfnr2 - 1;    // 0-based

      // EdgePointGeomInfo: edgenr must match, u/v = 0 (not used for projection)
      nseg.epgeominfo[0].edgenr = edgenr;
      nseg.epgeominfo[0].dist = 0;
      nseg.epgeominfo[0].u = 0;
      nseg.epgeominfo[0].v = 0;
      nseg.epgeominfo[1].edgenr = edgenr;
      nseg.epgeominfo[1].dist = 0;
      nseg.epgeominfo[1].u = 0;
      nseg.epgeominfo[1].v = 0;

      ng_mesh_->AddSegment(nseg);
    }
    edgenr++;
  }

  ng_mesh_->UpdateTopology();
  return true;
}

// ============================================================
// Phase 2: Attach CallbackGeometry
// ============================================================
bool NetgenCurverPure::attach_callback_geometry(
    ProjectFunc proj, NormalFunc norm, int num_surfaces,
    EdgeProjectFunc edge_proj)
{
  auto geom = std::make_shared<ng::CallbackGeometry>(
      proj, norm, num_surfaces, nullptr, edge_proj);
  ng_mesh_->SetGeometry(geom);
  return true;
}

// ============================================================
// Phase 3: BuildCurvedElements
// ============================================================
bool NetgenCurverPure::curve(int order)
{
  ng_mesh_->BuildCurvedElements(order);

  auto& curved = ng_mesh_->GetCurvedElements();
  if (!curved.IsHighOrder()) {
    return false;
  }
  return true;
}
