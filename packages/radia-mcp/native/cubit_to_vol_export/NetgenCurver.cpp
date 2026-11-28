#include "NetgenCurver.hpp"


// Include Netgen BEFORE Cubit to manage namespace conflicts
// Netgen defines: netgen::QUAD, netgen::HEX, netgen::PYRAMID, netgen::Surface
#include "netgen_meshing.hpp"
#include "CurveIdentity.hpp"

#include "MeshData.hpp"
#include "callbackgeom.hpp"
#include "netgen_curving/geometric_refit.hpp"
#include "netgen_curving/curved_access.hpp"

// Avoid Cubit/Netgen name clashes by aliasing Netgen element types
namespace ng = netgen;

// Evaluate a volume element through Netgen's multi-point transformation:
// NGSolve assembles with that path, and the scalar CalcElementTransformation
// ignores the quad-face coefficients of hexes (its derivative rows for them
// are even left uninitialised).
static void map_volume(ng::CurvedElements & ce, ng::ElementIndex ei,
                       const ng::Point<3> & ref, ng::Point<3> & phys)
{
  double xi[3] = {ref(0), ref(1), ref(2)}, x[3];
  ce.CalcMultiPointElementTransformation(ei, 1, xi, 3, x, 3, (double *)nullptr, 0);
  phys = ng::Point<3>(x[0], x[1], x[2]);
}

// Now include Cubit headers
// Cubit defines: ::QUAD, ::HEX, ::PYRAMID, ::Surface (class)
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"

// cubit_geom headers for ACIS projection
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "Body.hpp"
#include "DLIList.hpp"
#include "Surface.hpp"
#include "Curve.hpp"
#include "CubitVector.hpp"
// 2026-05-25: Cubit 2025.12 removed GeometryQueryTool::instance() static
// accessor.  CubitRefAccessor wraps the new path
//   CubitCoreModel::instance()->cgmapp()->factory()->get_ref_{face,edge,volume}()
// so the 6 lookup sites below stay one-liners.
#include "CubitRefAccessor.hpp"
#include "CurveIdentity.hpp"

// Use explicit namespace prefixes to disambiguate
// Cubit ElementType enums: ::TETRA4, ::TETRA, ::HEX8, ::HEX, etc.
// Netgen ELEMENT_TYPE enums: netgen::TET, netgen::HEX, netgen::PRISM, etc.

#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <fstream>
#include <limits>
#include <map>
#include <set>
#include <stdexcept>
#include <sstream>
#include <cstdio>

NetgenCurver::NetgenCurver() {}
NetgenCurver::~NetgenCurver() {}

bool NetgenCurver::build(const MeshData &md, int order) try
{
  last_error_.clear();
  if (order < 2) return false;
  order_ = order;

  PRINT_INFO("NetgenCurver: Phase 1 - build_netgen_mesh (MeshData)...\n");
  if (!build_netgen_mesh(md)) return false;
  PRINT_INFO("NetgenCurver: Phase 2 - attach_callback_geometry...\n");
  if (!attach_callback_geometry()) return false;
  PRINT_INFO("NetgenCurver: Phase 3 - curve_and_extract(order=%d)...\n", order);
  if (!curve_and_extract(order)) return false;

  if (!diag_edge_worst_.empty()) PRINT_INFO("%s", diag_edge_worst_.c_str());
  // Diagnostic report: projection reject stats
  PRINT_INFO("NetgenCurver: projection stats — "
             "project=%ld calls, %ld rejects (>%g), max_disp=%g (surf %d) [model length units]\n",
             diag_project_calls_, diag_project_rejects_,
             project_reject_threshold_,
             diag_project_max_disp_, diag_project_max_disp_surf_);
  PRINT_INFO("NetgenCurver: edge_project stats — %ld calls, %ld rejects, max_disp=%g [model length units]\n",
             diag_edge_calls_, diag_edge_rejects_, diag_edge_max_disp_);
  return true;
}
catch (const std::exception & e) {
  last_error_ = e.what();
  PRINT_ERROR("NetgenCurver: %s\n", e.what());
  if (!diag_edge_worst_.empty()) PRINT_ERROR("%s", diag_edge_worst_.c_str());
  PRINT_ERROR("NetgenCurver export diagnostics: %s\n", edge_diagnostics_json().c_str());
  return false;
}
catch (...) {
  last_error_ = "unknown native curving exception";
  PRINT_ERROR("NetgenCurver: %s\n", last_error_.c_str());
  return false;
}

std::string NetgenCurver::edge_diagnostics_json() const
{
  std::string escaped;
  for (unsigned char c : last_error_) {
    if (c == '"' || c == '\\') escaped += '\\';
    if (c < 32) { escaped += ' '; continue; }
    escaped += static_cast<char>(c);
  }
  std::ostringstream out;
  out << "{\"edge_calls\":" << diag_edge_calls_
      << ",\"edge_rejects\":" << diag_edge_rejects_
      << ",\"edge_failures\":" << diag_edge_failures_
      << ",\"error\":\"" << escaped << "\"}";
  return out.str();
}


// ============================================================
// Phase 1: Convert linear mesh (from MeshData) to netgen::Mesh
// ============================================================
bool NetgenCurver::build_netgen_mesh(const MeshData &md)
{
  PRINT_INFO("  Phase1a: creating ng::Mesh...\n");
  try {
    ng_mesh_ = std::make_shared<ng::Mesh>();
  } catch (const std::exception &e) {
    PRINT_ERROR("  Phase1a FAILED (exception): %s\n", e.what());
    return false;
  } catch (...) {
    PRINT_ERROR("  Phase1a FAILED (unknown exception)\n");
    return false;
  }
  PRINT_INFO("  Phase1a: ng::Mesh created OK\n");
  PRINT_INFO("  Phase1b: SetDimension(3)\n");
  ng_mesh_->SetDimension(3);
  // Mark the mesh as ACIS-derived so Netgen's Save() writes
  // `surfaceelementsuv` (not plain `surfaceelements`) -- the UV
  // parameters are already populated per-vertex via
  // `el.GeomInfo()[k].u/.v` in Phase1e.  NGSolve uses these UVs to
  // resolve high-order DOF coupling at periodic boundaries; without
  // the UVs (geomtype == NO_GEOM, the constructor default), order=2
  // FEM Omega + Periodic Kelvin produces wrong field magnitudes
  // (verified 2026-04-25 on em_elf_quarter_xz: order=2 gave -138 mT
  // instead of expected -228 mT).
  ng_mesh_->geomtype = ng::Mesh::GEOM_ACIS;

  int num_nodes = (int)md.nodes.size();
  PRINT_INFO("  Phase1c: AddPoint for %d nodes\n", num_nodes);

  for (int i = 0; i < num_nodes; i++) {
    auto &nd = md.nodes[i];
    if (i == 0) { PRINT_INFO("  Phase1c: first AddPoint (id=%d, %.4f, %.4f, %.4f)\n", nd.id, nd.x, nd.y, nd.z); }
    ng::PointIndex pi;
    try {
      pi = ng_mesh_->AddPoint(ng::MeshPoint(ng::Point<3>(nd.x, nd.y, nd.z)));
    } catch (const std::exception &e) {
      PRINT_ERROR("  Phase1c: AddPoint FAILED at node %d: %s\n", i, e.what());
      return false;
    } catch (...) {
      PRINT_ERROR("  Phase1c: AddPoint FAILED at node %d (unknown)\n", i);
      return false;
    }
    cubit_nid_to_ng_pi_[nd.id] = pi.Nr1();
    ng_pi_to_cubit_nid_[pi.Nr1()] = nd.id;
    if (i == 0 || i == num_nodes - 1)
      PRINT_INFO("  Phase1c: node %d/%d OK\n", i + 1, num_nodes);
  }


  PRINT_INFO("  Phase1d: FaceDescriptors\n");
  // --- Build FaceDescriptors for each Cubit surface ---
  PRINT_INFO("  Phase1d: calling get_entities('surface')...\n");
  std::vector<int> surface_ids = CubitInterface::get_entities("surface");
  PRINT_INFO("  Phase1d: got %d surfaces\n", (int)surface_ids.size());
  int fd_idx = 1;  // Netgen FD indices are 1-based
  for (int sid : surface_ids) {
    std::string surface_query = "in surface " + std::to_string(sid);
    std::vector<int> tri_ids =
        CubitInterface::parse_cubit_list("tri", surface_query);
    std::vector<int> quad_ids =
        CubitInterface::parse_cubit_list("face", surface_query);
    if (tri_ids.empty() && quad_ids.empty())
      continue;
    ng::FaceDescriptor fd(fd_idx, 1, 0, fd_idx);
    fd.SetBCProperty(sid);
    ng_mesh_->AddFaceDescriptor(fd);
    cubit_sid_to_ng_fd_[sid] = fd_idx;
    ng_fd_to_cubit_sid_[fd_idx] = sid;
    fd_idx++;
  }

  // --- Phase1e: Add surface elements via parse_cubit_list + get_connectivity ---
  // get_surface_tris() returns UNSHARED node IDs (rendering mesh).
  // Instead, use parse_cubit_list("tri"/"face", "in surface N") which returns
  // mesh entity IDs, then get_connectivity() for shared node IDs.
  PRINT_INFO("  Phase1e: Add surface elements (parse_cubit_list)\n");
  {
    auto canonical_face = [](const std::vector<int> &conn, int nv) {
      std::vector<int> key(conn.begin(), conn.begin() + nv);
      std::sort(key.begin(), key.end());
      return key;
    };
    std::set<std::vector<int>> added_surface_faces;

    // Recover domain ownership from volume-element adjacency.  Free Sculpt
    // sidesets are not attached to CAD surfaces, so Cubit's surface-volume
    // topology cannot tell the exporter which material owns a face.  One
    // sideset can also span several blocks; grouping by the adjacent-domain
    // pair prevents every free face from silently inheriting domain 1.
    std::unordered_map<int, int> block_to_domain;
    for (size_t i = 0; i < md.block_ids.size(); i++)
      block_to_domain[md.block_ids[i]] = (int)i + 1;

    std::map<std::vector<int>, std::set<int>> face_adjacent_domains;
    auto add_volume_face = [&](const MeshElement &elem,
                               std::initializer_list<int> local_nodes,
                               int domain) {
      std::vector<int> key;
      key.reserve(local_nodes.size());
      for (int local : local_nodes) {
        if (local < 0 || local >= elem.nv || local >= (int)elem.conn.size())
          return;
        key.push_back(elem.conn[local]);
      }
      std::sort(key.begin(), key.end());
      face_adjacent_domains[key].insert(domain);
    };
    for (const auto &elem : md.elements) {
      auto domain_it = block_to_domain.find(elem.group_id);
      if (domain_it == block_to_domain.end()) {
        PRINT_ERROR("Volume element belongs to unknown block %d while recovering free sideset ownership.\n",
                    elem.group_id);
        return false;
      }
      int domain = domain_it->second;
      if (elem.nv == 4 && (elem.type == ::TETRA4 || elem.type == ::TETRA)) {
        add_volume_face(elem, {0, 1, 2}, domain);
        add_volume_face(elem, {0, 1, 3}, domain);
        add_volume_face(elem, {1, 2, 3}, domain);
        add_volume_face(elem, {2, 0, 3}, domain);
      } else if (elem.nv == 8 && (elem.type == ::HEX8 || elem.type == ::HEX)) {
        add_volume_face(elem, {0, 1, 2, 3}, domain);
        add_volume_face(elem, {4, 5, 6, 7}, domain);
        add_volume_face(elem, {0, 1, 5, 4}, domain);
        add_volume_face(elem, {1, 2, 6, 5}, domain);
        add_volume_face(elem, {2, 3, 7, 6}, domain);
        add_volume_face(elem, {3, 0, 4, 7}, domain);
      } else if (elem.nv == 6 &&
                 (elem.type == ::WEDGE6 || elem.type == ::WEDGE)) {
        add_volume_face(elem, {0, 1, 2}, domain);
        add_volume_face(elem, {3, 4, 5}, domain);
        add_volume_face(elem, {0, 1, 4, 3}, domain);
        add_volume_face(elem, {1, 2, 5, 4}, domain);
        add_volume_face(elem, {2, 0, 3, 5}, domain);
      } else if (elem.nv == 5 &&
                 (elem.type == ::PYRAMID5 || elem.type == ::PYRAMID)) {
        add_volume_face(elem, {0, 1, 2, 3}, domain);
        add_volume_face(elem, {0, 1, 4}, domain);
        add_volume_face(elem, {1, 2, 4}, domain);
        add_volume_face(elem, {2, 3, 4}, domain);
        add_volume_face(elem, {3, 0, 4}, domain);
      }
    }

    // Cache RefFace pointers per surface
    std::unordered_map<int, RefFace*> rf_cache;
    for (int sid : surface_ids)
      rf_cache[sid] = cubit_to_vol_export::cubit_get_ref_face(sid);

    // Cache per-surface periodicity. Seam shift (below) is applied only to
    // surfaces that ACIS flags as actually periodic in that parameter
    // direction — for non-periodic (open) NURBS, `u_v_from_position` returns
    // a unique UV in the parameter range, so there is no "seam" and shifting
    // by the parameter range would push the UV outside the valid interval.
    // period == 0 means "don't apply seam shift in this direction".
    std::unordered_map<int, std::array<double,2>> period_cache;
    for (int sid : surface_ids) {
      RefFace* rf = rf_cache[sid];
      double up = 0, vp = 0;
      if (rf) {
        double tmp = 0;
        if (rf->is_periodic_in_U(tmp)) up = tmp;
        if (rf->is_periodic_in_V(tmp)) vp = tmp;
      }
      period_cache[sid] = {up, vp};
    }

    // Detect polar-singular planar disks (wire end-cap faces produced by
    // `unite volume all` on lofted coil assemblies). These have at least one
    // singular pole on the parametric surface AND trim to a flat plane, so
    // the radial center maps to a single 3D point with arbitrary UV angle.
    // `closest_point_uv_guess` on such surfaces produces spikes on the order
    // of the face radius (~wire diameter). Since a planar face requires no
    // high-order curving at all, we short-circuit projection to identity in
    // the callback.  Non-planar polar surfaces (sphere, cone apex) are NOT
    // flagged — their curving paths work fine as long as meshes avoid the
    // exact pole, which is the common case.
    polar_surfaces_.clear();
    int diag_n_poled = 0, diag_n_planar = 0, diag_n_both = 0;
    for (int sid : surface_ids) {
      RefFace* rf = rf_cache[sid];
      if (!rf) continue;
      int n_poles = rf->num_poles();
      bool planar = (rf->is_planar() == CUBIT_TRUE);
      if (n_poles > 0) diag_n_poled++;
      if (planar) diag_n_planar++;
      if (n_poles > 0 && planar) {
        diag_n_both++;
        polar_surfaces_.insert(sid);
      }
    }
    PRINT_INFO("  Phase1: geometry inspection - %d surfaces, "
               "%d with poles (num_poles>0), %d planar, %d both (=polar disk)\n",
               (int)surface_ids.size(), diag_n_poled, diag_n_planar, diag_n_both);
    if (!polar_surfaces_.empty()) {
      PRINT_INFO("  Phase1: polar planar-disk surfaces detected = %zu "
                 "(HO projection will be skipped)\n",
                 polar_surfaces_.size());
      // Dump their IDs for debugging
      std::string ids_str;
      for (int sid : polar_surfaces_) {
        ids_str += std::to_string(sid) + " ";
      }
      PRINT_INFO("    polar surface ids: %s\n", ids_str.c_str());
    }

    // Helper: shift u1 by ±period to minimize |u1 - u0|. Returns adjusted u1.
    auto seam_consistent = [](double u_ref, double u, double period) -> double {
      if (period <= 0.0) return u;
      double du = u - u_ref;
      while (du > 0.5 * period) { u -= period; du -= period; }
      while (du < -0.5 * period) { u += period; du += period; }
      return u;
    };

    // Pre-pass: build per-node canonical UV per surface.
    // Key = (int64)sid << 32 | cubit_nid. Value = (u, v) from first
    // u_v_from_position call. For non-periodic surfaces this is the
    // unique UV; for periodic surfaces it resolves the seam ambiguity
    // deterministically (first-seen wins) so all triangles sharing that
    // node get the same UV before per-triangle seam shift.
    // In parallel we also populate surf_vertex_uvs_[sid] — a per-surface
    // table of (x, y, z, u, v) rows used by the project() callback to
    // derive a UV hint from the nearest stored vertex on its no-hint
    // first call.
    std::unordered_map<int64_t, std::pair<double,double>> node_canon_uv;
    surf_vertex_uvs_.clear();
    auto collect_uvs = [&](int sid, const std::vector<int> &conn) {
      RefFace* rf = rf_cache[sid];
      if (!rf) return;
      for (int nid : conn) {
        int64_t key = ((int64_t)sid << 32) | (uint32_t)nid;
        if (node_canon_uv.count(key)) continue;
        auto coords = CubitInterface::get_nodal_coordinates(nid);
        CubitVector loc(coords[0], coords[1], coords[2]);
        double u = 0, v = 0;
        rf->u_v_from_position(loc, u, v);
        node_canon_uv[key] = {u, v};
        surf_vertex_uvs_[sid].push_back(
            {coords[0], coords[1], coords[2], u, v});
      }
    };
    for (int sid : surface_ids) {
      RefFace* rf = rf_cache[sid];
      if (!rf) continue;
      std::vector<int> tri_ids = CubitInterface::parse_cubit_list("tri",
          "in surface " + std::to_string(sid));
      for (int tid : tri_ids) {
        std::vector<int> conn = CubitInterface::get_connectivity("tri", tid);
        collect_uvs(sid, conn);
      }
      std::vector<int> face_ids = CubitInterface::parse_cubit_list("face",
          "in surface " + std::to_string(sid));
      for (int fid : face_ids) {
        std::vector<int> conn = CubitInterface::get_connectivity("face", fid);
        collect_uvs(sid, conn);
      }
    }

    int nse_added = 0;
    long diag_shift_count = 0;
    for (int sid : surface_ids) {
      auto ng_fd_it = cubit_sid_to_ng_fd_.find(sid);
      if (ng_fd_it == cubit_sid_to_ng_fd_.end())
        continue;
      int ng_fd = ng_fd_it->second;
      RefFace* rf = rf_cache[sid];
      double u_period = period_cache[sid][0];
      double v_period = period_cache[sid][1];

      // Triangles on this surface
      std::vector<int> tri_ids = CubitInterface::parse_cubit_list("tri",
          "in surface " + std::to_string(sid));
      for (int tid : tri_ids) {
        std::vector<int> conn = CubitInterface::get_connectivity("tri", tid);
        if ((int)conn.size() < 3) continue;

        ng::Element2d el(ng::TRIG);
        el.SetIndex(ng::FaceRegionIndex::FromNr1(ng_fd));
        // Pass 1: collect UVs from per-node canonical UV map (computed above)
        double uu[3], vv[3];
        int pis[3];
        bool uv_ok = true;
        for (int k = 0; k < 3; k++) {
          auto it = cubit_nid_to_ng_pi_.find(conn[k]);
          if (it == cubit_nid_to_ng_pi_.end()) { uv_ok = false; goto skip_tri; }
          pis[k] = it->second;
          el[k] = ng::PointIndex::FromNr1(it->second);
          int64_t ckey = ((int64_t)sid << 32) | (uint32_t)conn[k];
          auto cit = node_canon_uv.find(ckey);
          if (cit != node_canon_uv.end()) {
            uu[k] = cit->second.first;
            vv[k] = cit->second.second;
          } else {
            // Fallback: recompute (should not happen since pre-pass covered all tris)
            auto coords = CubitInterface::get_nodal_coordinates(conn[k]);
            CubitVector loc(coords[0], coords[1], coords[2]);
            double u = 0, v = 0;
            if (rf) rf->u_v_from_position(loc, u, v);
            uu[k] = u; vv[k] = v;
          }
        }
        // Seam-consistent UV pass: use vertex 0 as reference, shift v1, v2 by
        // ±period so their UV is closest to vertex 0 in the periodic wraparound
        // sense. Prevents Netgen from averaging (u=0.9, u=0.1) into hint=0.5
        // which would land on the opposite side of a periodic NURBS seam.
        if (u_period > 0.0) {
          uu[1] = seam_consistent(uu[0], uu[1], u_period);
          uu[2] = seam_consistent(uu[0], uu[2], u_period);
        }
        if (v_period > 0.0) {
          vv[1] = seam_consistent(vv[0], vv[1], v_period);
          vv[2] = seam_consistent(vv[0], vv[2], v_period);
        }
        // Record UV spread AFTER seam shift (should be small)
        {
          double du01 = std::fabs(uu[1]-uu[0]);
          double du12 = std::fabs(uu[2]-uu[1]);
          double du02 = std::fabs(uu[2]-uu[0]);
          double dv01 = std::fabs(vv[1]-vv[0]);
          double dv12 = std::fabs(vv[2]-vv[1]);
          double dv02 = std::fabs(vv[2]-vv[0]);
          double max_du = std::max({du01, du12, du02});
          double max_dv = std::max({dv01, dv12, dv02});
          if (u_period > 0 && max_du > 0.5 * u_period) diag_shift_count++;
          if (v_period > 0 && max_dv > 0.5 * v_period) diag_shift_count++;
        }
        // Pass 2: emit GeomInfo with seam-consistent UVs
        for (int k = 0; k < 3; k++) {
          ng::PointGeomInfo gi;
          gi.trignum = ng_fd;
          gi.u = uu[k];
          gi.v = vv[k];
          el.GeomInfo()[k] = gi;
        }
        ng_mesh_->AddSurfaceElement(el);
        added_surface_faces.insert(canonical_face(conn, 3));
        nse_added++;
        skip_tri:;
      }

      // Quads on this surface (hex meshes)
      std::vector<int> face_ids = CubitInterface::parse_cubit_list("face",
          "in surface " + std::to_string(sid));
      for (int fid : face_ids) {
        std::vector<int> conn = CubitInterface::get_connectivity("face", fid);
        if ((int)conn.size() < 4) continue;

        ng::Element2d el(ng::QUAD);
        el.SetIndex(ng::FaceRegionIndex::FromNr1(ng_fd));
        double qu[4], qv[4];
        for (int k = 0; k < 4; k++) {
          auto it = cubit_nid_to_ng_pi_.find(conn[k]);
          if (it == cubit_nid_to_ng_pi_.end()) goto skip_quad;
          el[k] = ng::PointIndex::FromNr1(it->second);

          int64_t ckey = ((int64_t)sid << 32) | (uint32_t)conn[k];
          auto cit = node_canon_uv.find(ckey);
          if (cit != node_canon_uv.end()) {
            qu[k] = cit->second.first;
            qv[k] = cit->second.second;
          } else {
            auto coords = CubitInterface::get_nodal_coordinates(conn[k]);
            CubitVector loc(coords[0], coords[1], coords[2]);
            double u = 0, v = 0;
            if (rf) rf->u_v_from_position(loc, u, v);
            qu[k] = u; qv[k] = v;
          }
        }
        // Seam-consistent UV (anchor at vertex 0) for periodic surfaces
        if (u_period > 0.0) {
          qu[1] = seam_consistent(qu[0], qu[1], u_period);
          qu[2] = seam_consistent(qu[0], qu[2], u_period);
          qu[3] = seam_consistent(qu[0], qu[3], u_period);
        }
        if (v_period > 0.0) {
          qv[1] = seam_consistent(qv[0], qv[1], v_period);
          qv[2] = seam_consistent(qv[0], qv[2], v_period);
          qv[3] = seam_consistent(qv[0], qv[3], v_period);
        }
        for (int k = 0; k < 4; k++) {
          ng::PointGeomInfo gi;
          gi.trignum = ng_fd;
          gi.u = qu[k];
          gi.v = qv[k];
          el.GeomInfo()[k] = gi;
        }
        ng_mesh_->AddSurfaceElement(el);
        added_surface_faces.insert(canonical_face(conn, 4));
        nse_added++;
        skip_quad:;
      }
    }

    // Sculpt and imported Exodus meshes can carry boundary faces directly in
    // a sideset without associating them with a geometric surface.  Add those
    // free faces under one synthetic FaceDescriptor per sideset.  Negative
    // BCProperty values are an internal handoff to ExportVolCommand, which
    // restores the sideset name after all descriptors are made contiguous.
    int free_face_count = 0;
    int free_descriptor_count = 0;
    for (const auto &sg : md.sidesets) {
      std::map<std::pair<int, int>, std::vector<const MeshElement*>> pending;
      for (const auto &face : sg.faces) {
        if (face.nv != 3 && face.nv != 4) continue;
        auto key = canonical_face(face.conn, face.nv);
        if (added_surface_faces.count(key)) continue;

        auto adjacent_it = face_adjacent_domains.find(key);
        if (adjacent_it == face_adjacent_domains.end() ||
            adjacent_it->second.empty()) {
          PRINT_ERROR("Free sideset %d contains a face with no adjacent volume element.\n",
                      sg.id);
          return false;
        }
        if (adjacent_it->second.size() > 2) {
          PRINT_ERROR("Free sideset %d contains a non-manifold face adjacent to %d domains.\n",
                      sg.id, (int)adjacent_it->second.size());
          return false;
        }
        auto domain_it = adjacent_it->second.begin();
        int domain_in = *domain_it++;
        int domain_out = domain_it == adjacent_it->second.end() ? 0 : *domain_it;
        pending[{domain_in, domain_out}].push_back(&face);
      }
      if (pending.empty()) continue;

      for (const auto &domain_faces : pending) {
        int domain_in = domain_faces.first.first;
        int domain_out = domain_faces.first.second;
        int ng_fd = fd_idx++;
        ng::FaceDescriptor fd(ng_fd, domain_in, domain_out, ng_fd);
        fd.SetBCProperty(-sg.id);
        ng_mesh_->AddFaceDescriptor(fd);
        free_descriptor_count++;

        for (const MeshElement *face : domain_faces.second) {
          ng::ELEMENT_TYPE etype = face->nv == 3 ? ng::TRIG : ng::QUAD;
          ng::Element2d el(etype);
          el.SetIndex(ng::FaceRegionIndex::FromNr1(ng_fd));
          bool complete = true;
          for (int k = 0; k < face->nv; k++) {
            auto it = cubit_nid_to_ng_pi_.find(face->conn[k]);
            if (it == cubit_nid_to_ng_pi_.end()) {
              complete = false;
              break;
            }
            el[k] = ng::PointIndex::FromNr1(it->second);
            ng::PointGeomInfo gi;
            gi.trignum = ng_fd;
            gi.u = 0.0;
            gi.v = 0.0;
            el.GeomInfo()[k] = gi;
          }
          if (!complete) continue;
          ng_mesh_->AddSurfaceElement(el);
          added_surface_faces.insert(canonical_face(face->conn, face->nv));
          nse_added++;
          free_face_count++;
        }
      }
    }
    if (free_face_count > 0) {
      PRINT_INFO("  Phase1e: added %d free sideset faces under %d domain-aware synthetic descriptors\n",
                 free_face_count, free_descriptor_count);
    }
    PRINT_INFO("  Phase1e: %d surface elements added (%ld seam-consistent UV shifts applied)\n",
               nse_added, diag_shift_count);
  }

  // --- Phase1e2: Add segment elements on curves (inter-surface edges) ---
  PRINT_INFO("  Phase1e2: Add segments on curves\n");
  {
    std::vector<int> curve_ids = CubitInterface::parse_cubit_list("curve", "all");
    cubit_curve_to_ng_edge_.clear();
    int edgenr = 1;
    for (int cid : curve_ids) {
      // Get parent surfaces of this curve
      std::vector<int> parent_surfs = CubitInterface::parse_cubit_list(
          "surface", "in curve " + std::to_string(cid));
      if (parent_surfs.size() < 2) continue;

      int sid1 = parent_surfs[0], sid2 = parent_surfs[1];
      auto it1 = cubit_sid_to_ng_fd_.find(sid1);
      auto it2 = cubit_sid_to_ng_fd_.find(sid2);
      if (it1 == cubit_sid_to_ng_fd_.end() || it2 == cubit_sid_to_ng_fd_.end())
        continue;
      int fd1 = it1->second;  // 1-based
      int fd2 = it2->second;

      // Get edge (1D) elements on this curve
      std::vector<int> edge_ids = CubitInterface::parse_cubit_list(
          "edge", "in curve " + std::to_string(cid));

      // Netgen >= 6.2.2605 keeps geometric edge metadata in a 1-based
      // EdgeDescriptor table.  Its surface numbers are 1-based face
      // descriptor indices, passed to CallbackGeometry as-is.
      ng::EdgeDescriptor ed(edgenr, fd1, fd2);
      ed.SetIndex(ng::FaceRegionIndex::FromNr1(fd2));
      const auto edge_descriptor_index = ng_mesh_->AddEdgeDescriptor(ed);

      // Get curve length for dist normalization
      RefEdge* re = cubit_to_vol_export::cubit_get_ref_edge(cid);
      double crv_length = re ? re->measure() : 1.0;
      double u_start = re ? re->start_param() : 0.0;

      int exported_segments = 0;
      for (int eid : edge_ids) {
        std::vector<int> conn = CubitInterface::get_connectivity("edge", eid);
        if ((int)conn.size() < 2) continue;
        auto itn0 = cubit_nid_to_ng_pi_.find(conn[0]);
        auto itn1 = cubit_nid_to_ng_pi_.find(conn[1]);
        if (itn0 == cubit_nid_to_ng_pi_.end() || itn1 == cubit_nid_to_ng_pi_.end())
          continue;

        // Compute normalized arc-length parameter for each endpoint
        double dist0 = 0, dist1 = 0;
        if (re && crv_length > 0) {
          auto c0 = CubitInterface::get_nodal_coordinates(conn[0]);
          auto c1 = CubitInterface::get_nodal_coordinates(conn[1]);
          CubitVector cv0(c0[0], c0[1], c0[2]);
          CubitVector cv1(c1[0], c1[1], c1[2]);
          double u0 = re->u_from_position(cv0);
          double u1 = re->u_from_position(cv1);
          dist0 = re->length_from_u(u_start, u0) / crv_length;
          dist1 = re->length_from_u(u_start, u1) / crv_length;
        }

        ng::Segment nseg;
        nseg[0] = ng::PointIndex::FromNr1(itn0->second);
        nseg[1] = ng::PointIndex::FromNr1(itn1->second);
        nseg.SetIndex(edge_descriptor_index);
        nseg.EPGeomInfo(0).dist = dist0;
        nseg.EPGeomInfo(1).dist = dist1;
        ng_mesh_->AddSegment(nseg);
        exported_segments++;
      }
      if (exported_segments != static_cast<int>(edge_ids.size()))
        throw std::invalid_argument("Cubit curve " + std::to_string(cid) +
            " exported only " + std::to_string(exported_segments) + "/" +
            std::to_string(edge_ids.size()) + " segments");
      if (exported_segments > 0)
        cubit_curve_to_ng_edge_[cid] = edge_descriptor_index.Nr1();
      edgenr++;
    }
    PRINT_INFO("  Phase1e2: %d edge groups added\n", edgenr - 1);
  }

  PRINT_INFO("  Phase1f: Add volume elements (%d)\n", (int)md.elements.size());
  int vol_count = 0;

  // --- Vertex reordering: GMSH/Cubit -> Netgen convention ---
  // From netgen/read_gmsh.py: ordering[ng_idx] = gmsh_idx
  // TET: identity (same vertex numbering)
  // HEX: Netgen swaps vertices 2<->5, 3<->4 relative to GMSH
  //   GMSH:   0=(0,0,0) 1=(1,0,0) 2=(1,1,0) 3=(0,1,0) 4=(0,0,1) 5=(1,0,1) 6=(1,1,1) 7=(0,1,1)
  //   Netgen: 0=(0,0,0) 1=(1,0,0) 2=(1,0,1) 3=(0,0,1) 4=(0,1,0) 5=(1,1,0) 6=(1,1,1) 7=(0,1,1)
  static const int hex_gmsh_to_ng[8]  = {0, 1, 5, 4, 3, 2, 6, 7};
  // PRISM: Netgen reverses triangle winding (0,2,1 instead of 0,1,2)
  static const int prism_gmsh_to_ng[6] = {0, 2, 1, 3, 5, 4};
  // PYRAMID: Netgen reverses base quad winding
  static const int pyr_gmsh_to_ng[5]  = {3, 2, 1, 0, 4};

  // --- Add volume elements from MeshData ---
  for (auto &elem : md.elements) {
    int nn = elem.nv;
    ElementType et = elem.type;

    if (nn == 4 && (et == ::TETRA4 || et == ::TETRA)) {
      ng::Element el(ng::TET);
      for (int k = 0; k < 4; k++)
        el[k] = ng::PointIndex::FromNr1(cubit_nid_to_ng_pi_[elem.conn[k]]);
      el.SetIndex(ng::VolumeRegionIndex::FromNr1(1));
      ng_mesh_->AddVolumeElement(el);
      vol_count++;
    }
    else if (nn == 8 && (et == ::HEX8 || et == ::HEX)) {
      ng::Element el(ng::HEX);
      // Keep GMSH/Cubit vertex ordering — NOT reordered to Netgen convention.
      // Vertex reordering would fix .vol export but breaks CalcElementTransformation
      // for HO node computation (face-surface mapping inconsistency in BuildCurvedElements).
      // Interior volume edges use Cubit coordinate linear interpolation instead.
      for (int k = 0; k < 8; k++)
        el[k] = ng::PointIndex::FromNr1(cubit_nid_to_ng_pi_[elem.conn[k]]);
      el.SetIndex(ng::VolumeRegionIndex::FromNr1(1));
      ng_mesh_->AddVolumeElement(el);
      vol_count++;
    }
    else if (nn == 6 && (et == ::WEDGE6 || et == ::WEDGE)) {
      ng::Element el(ng::PRISM);
      for (int k = 0; k < 6; k++)
        el[k] = ng::PointIndex::FromNr1(cubit_nid_to_ng_pi_[elem.conn[k]]);
      el.SetIndex(ng::VolumeRegionIndex::FromNr1(1));
      ng_mesh_->AddVolumeElement(el);
      vol_count++;
    }
    else if (nn == 5 && (et == ::PYRAMID5 || et == ::PYRAMID)) {
      ng::Element el(ng::PYRAMID);
      for (int k = 0; k < 5; k++)
        el[k] = ng::PointIndex::FromNr1(cubit_nid_to_ng_pi_[elem.conn[k]]);
      el.SetIndex(ng::VolumeRegionIndex::FromNr1(1));
      ng_mesh_->AddVolumeElement(el);
      vol_count++;
    }
  }

  PRINT_INFO("  Phase1f: %d volume elements added\n", vol_count);

  PRINT_INFO("  Phase1g: UpdateTopology...\n");
  ng_mesh_->UpdateTopology();
  PRINT_INFO("  Phase1g: UpdateTopology done\n");

  return true;
}

// ============================================================
// Phase 2: Attach CallbackGeometry with ACIS projection
// ============================================================
bool NetgenCurver::attach_callback_geometry()
{
  // Includes synthetic free-mesh descriptors. Their callbacks intentionally
  // fall back to identity because they have no CAD RefFace.
  int num_surfaces = ng_mesh_ ? ng_mesh_->GetNFD() : 0;

  // Reset diagnostic counters for this curving run
  diag_project_calls_ = 0;
  diag_project_rejects_ = 0;
  diag_project_max_disp_ = 0.0;
  diag_project_max_disp_surf_ = -1;
  diag_edge_calls_ = 0;
  diag_edge_rejects_ = 0;
  diag_edge_max_disp_ = 0.0;
  diag_edge_failures_ = 0;
  diag_edge_worst_.clear();
  diag_polar_skips_ = 0;

  // Build per-surface bounding boxes from the Phase1e vertex-UV table.
  // `surf_bbox_pad_` is set to 1.5 x the largest mesh-edge length we
  // observe — large enough to accommodate a legitimate curvature offset
  // on a curved surface while still catching "interior midpoint assigned
  // to a distant terminal face" spikes (tens of mm).
  surf_bbox_.clear();
  double global_diam = 0.0;
  for (const auto &kv : surf_vertex_uvs_) {
    if (kv.second.empty()) continue;
    BBox bb{kv.second[0][0], kv.second[0][1], kv.second[0][2],
            kv.second[0][0], kv.second[0][1], kv.second[0][2]};
    for (const auto &row : kv.second) {
      if (row[0] < bb.xmin) bb.xmin = row[0];
      if (row[0] > bb.xmax) bb.xmax = row[0];
      if (row[1] < bb.ymin) bb.ymin = row[1];
      if (row[1] > bb.ymax) bb.ymax = row[1];
      if (row[2] < bb.zmin) bb.zmin = row[2];
      if (row[2] > bb.zmax) bb.zmax = row[2];
    }
    double dx = bb.xmax - bb.xmin;
    double dy = bb.ymax - bb.ymin;
    double dz = bb.zmax - bb.zmin;
    double diam = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (diam > global_diam) global_diam = diam;
    surf_bbox_[kv.first] = bb;
  }
  // Pad = 5% of the global diameter (a fraction of the whole body's size).
  // This keeps legitimate HO curving offsets (normally small relative to
  // the body) admissible while catching projections of points that sit a
  // substantial fraction of the body away from a candidate surface.
  surf_bbox_pad_ = 0.05 * global_diam;
  PRINT_INFO("NetgenCurver: built bbox for %d surfaces (pad=%g, global_diam=%g [model length units])\n",
             (int)surf_bbox_.size(), surf_bbox_pad_, global_diam);

  // Project callback: (surfnr, x,y,z, u_hint, v_hint, has_hint) -> (xp,yp,zp, u,v)
  auto project = [this](int surfnr, double x, double y, double z,
                         double u_hint, double v_hint, bool has_hint)
    -> std::tuple<double,double,double,double,double>
  {
    diag_project_calls_++;

    auto it = ng_fd_to_cubit_sid_.find(surfnr);
    if (it == ng_fd_to_cubit_sid_.end())
      return {x, y, z, u_hint, v_hint};

    int cubit_sid = it->second;

    // Polar planar-disk short-circuit: for faces flagged in Phase1 as
    // polar-singular planar surfaces (e.g. wire terminal end-caps produced
    // by `unite volume all` on lofted coils), there is no well-defined
    // projection via UV — the radial-center pole maps all U values to one
    // 3D point, and `closest_point_uv_guess` spikes on the order of the
    // face radius (61 mm observed on 3turnCoil ø6.3 mm wire after unite).
    // Planar faces require no high-order curving, so return identity.
    if (polar_surfaces_.count(cubit_sid)) {
      diag_polar_skips_++;
      return {x, y, z, u_hint, v_hint};
    }

    RefFace* rf = cubit_to_vol_export::cubit_get_ref_face(cubit_sid);
    if (!rf) return {x, y, z, u_hint, v_hint};

    // Early reject: if the query point lies noticeably outside this
    // surface's own bounding box, return identity. This catches the
    // dominant failure mode on `unite volume all` coil bodies, where the
    // curving pipeline asks us to project a coil-interior midpoint onto
    // a distant terminal face and the ACIS projection clamps to that
    // face's boundary at tens of mm displacement.
    auto bb_it = surf_bbox_.find(cubit_sid);
    if (bb_it != surf_bbox_.end()) {
      const BBox &bb = bb_it->second;
      double pad = surf_bbox_pad_;
      if (x < bb.xmin - pad || x > bb.xmax + pad ||
          y < bb.ymin - pad || y > bb.ymax + pad ||
          z < bb.zmin - pad || z > bb.zmax + pad) {
        return {x, y, z, u_hint, v_hint};
      }
    }

    ::Surface* surf = rf->get_surface_ptr();
    CubitVector loc(x, y, z);
    CubitVector closest;

    // Path: 0=uv_direct_eval, 1=uv_guess_fallback_trimmed,
    //       2=nearest-vertex-hint_then_uv_guess, 3=no_hint_table_trimmed
    int path = -1;

    // Netgen's curving loop calls this with u_hint, v_hint set to the
    // UV-interpolated position between two surface vertices.  Use them to
    // seed a UV-local Newton iteration (`closest_point_uv_guess`) — this
    // is the OCC equivalent of "project with hint": the hint fixes the
    // branch of the surface, and the Newton step recovers the chord-
    // perpendicular midpoint that Netgen's polynomial fit expects.
    //
    // Experiments showed plain `position_from_u_v(hint)` (strictly UV-
    // midpoint without projection) is slightly WORSE on 3turnCoil because
    // Netgen's LSQ curving assumes chord-perpendicular offsets, not UV-
    // midpoint offsets.  Keep the closest_point search but trust the
    // UV-hinted result unconditionally.
    double u, v;
    if (has_hint) {
      u = u_hint;
      v = v_hint;
      CubitStatus stat = surf->closest_point_uv_guess(loc, u, v, &closest, nullptr);
      if (stat != CUBIT_SUCCESS) {
        surf->closest_point_trimmed(loc, closest);
        rf->u_v_from_position(closest, u, v);
        path = 1;
      } else {
        path = 0;
      }
    } else {
      // First call: no UV hint from Netgen, but we have a per-surface UV
      // table from Phase1e. Find the nearest stored vertex on this surface
      // and reuse its canonical UV as a hint for closest_point_uv_guess —
      // this converts the global search into a local Newton iteration,
      // keeping the projection on the correct branch of thin-cross-section
      // lofts where closest_point_trimmed tends to jump to the opposite
      // wire face.
      // Netgen called ProjectPoint without a UV hint — this happens for
      // first-time projections (no GeomInfo yet).  Find the nearest stored
      // surface vertex and use its UV as hint for closest_point_uv_guess;
      // this is the OCC-equivalent of "project with hint" and keeps the
      // Newton iteration on the correct branch of thin-cross-section lofts.
      auto sit = surf_vertex_uvs_.find(cubit_sid);
      bool guess_ok = false;
      if (sit != surf_vertex_uvs_.end() && !sit->second.empty()) {
        const auto &verts = sit->second;
        double best_d2 = std::numeric_limits<double>::infinity();
        double uh = 0, vh = 0;
        for (const auto &row : verts) {
          double dx = row[0] - x;
          double dy = row[1] - y;
          double dz = row[2] - z;
          double d2 = dx*dx + dy*dy + dz*dz;
          if (d2 < best_d2) { best_d2 = d2; uh = row[3]; vh = row[4]; }
        }
        u = uh; v = vh;
        CubitStatus stat = surf->closest_point_uv_guess(loc, u, v, &closest, nullptr);
        if (stat == CUBIT_SUCCESS) {
          path = 2;
          guess_ok = true;
        }
      }
      if (!guess_ok) {
        // No hint table — fall back to global search.
        surf->closest_point_trimmed(loc, closest);
        rf->u_v_from_position(closest, u, v);
        path = 3;
      }
    }

    // Diagnostic + defensive reject: if ACIS projection displaces the point
    // beyond project_reject_threshold_, the projection is almost certainly
    // a spike (seam jump or wrong local minimum on loft NURBS). Fall back
    // to identity (returned point = original), which Phase 2 then handles
    // via its linear-interpolation fill pass.
    //
    // path=0 is UV-direct eval (position_from_u_v), which is mathematically
    // exact for the supplied (u_hint, v_hint) — it involves no closest-point
    // search, so the displacement has no "wrong-branch" meaning.  Trust it
    // unconditionally.
    double dx = closest.x() - x;
    double dy = closest.y() - y;
    double dz = closest.z() - z;
    double disp = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (disp > diag_project_max_disp_) {
      diag_project_max_disp_ = disp;
      diag_project_max_disp_surf_ = cubit_sid;
    }
    // path 0 (UV hint supplied by Netgen, eval via position_from_u_v) is
    // trusted unconditionally — the hint came from interpolating two
    // surface vertices' UVs, which is mathematically well-defined on a
    // regular NURBS.  For paths 1, 2, 3 we synthesized or recovered the
    // UV ourselves, so a displacement well beyond the expected curvature
    // offset (~ surface thickness / mesh size) signals a degenerate
    // parametrization — typically a polar/pole singularity on a disk
    // end-cap where u_v_from_position returns arbitrary angle UVs.
    // Reject those with a generous threshold so normal curvature offsets
    // on curved lofts pass (~1-3 mm on 3turnCoil wire) but polar spikes
    // (~wire diameter or more) don't.
    // path 0 (UV-hint from Netgen, direct surface eval): trust unconditionally.
    // paths 1, 2, 3 (closest-point searches): reject jumps beyond
    // project_reject_threshold_.
    double path_threshold = (path == 0) ? 1e30 : project_reject_threshold_;
    if (disp > path_threshold) {
      diag_project_rejects_++;
      // Return identity for the 3D point (HO node stays on the linear
      // midpoint) but return a UV derived from the *input* point rather
      // than the boundary-clamped projection result.  This stops bogus
      // boundary UVs from propagating to downstream PointBetween calls
      // as their UV hint, which would otherwise compound the error.
      double u_clean = 0, v_clean = 0;
      rf->u_v_from_position(loc, u_clean, v_clean);
      return {x, y, z, u_clean, v_clean};
    }

    return {closest.x(), closest.y(), closest.z(), u, v};
  };

  // Normal callback: (surfnr, x,y,z) -> (nx,ny,nz)
  auto normal = [this](int surfnr, double x, double y, double z)
    -> std::tuple<double,double,double>
  {
    auto it = ng_fd_to_cubit_sid_.find(surfnr);
    if (it == ng_fd_to_cubit_sid_.end())
      return {0.0, 0.0, 1.0};

    int cubit_sid = it->second;
    RefFace* rf = cubit_to_vol_export::cubit_get_ref_face(cubit_sid);
    if (!rf) return {0.0, 0.0, 1.0};

    CubitVector loc(x, y, z);
    CubitVector n = rf->normal_at(loc);
    return {n.x(), n.y(), n.z()};
  };

  // Netgen passes EdgeDescriptor::EdgeNr(), NOT the descriptor table index.
  // Two distinct curves may bound the same pair of split smooth surfaces.
  auto edge_bindings = curve_identities(*ng_mesh_, cubit_curve_to_ng_edge_);
  auto curve_for_edge = [edge_bindings](int edge, int face1, int face2) -> RefEdge* {
    const auto & b = edge_bindings.lookup(edge, face1, face2);
    RefEdge* re = cubit_to_vol_export::cubit_get_ref_edge(b.curve);
    if (!re || !re->get_curve_ptr())
      throw std::invalid_argument("Missing Cubit curve " + std::to_string(b.curve) +
                                  " for Netgen geometric edge " + std::to_string(edge));
    return re;
  };

  auto edge_project = [this, curve_for_edge, edge_bindings](int edgenr, int surfnr1, int surfnr2,
                                                   double x, double y, double z)
    -> std::tuple<double,double,double>
  {
    diag_edge_calls_++;
    RefEdge* re = curve_for_edge(edgenr, surfnr1, surfnr2);
    CubitVector loc(x, y, z);
    CubitVector closest;
    if (re->get_curve_ptr()->closest_point_trimmed(loc, closest) != CUBIT_SUCCESS) {
      ++diag_edge_failures_;
      throw std::runtime_error("closest_point_trimmed failed for Cubit curve " +
          std::to_string(re->id()) + " edge " + std::to_string(edgenr));
    }
    double dx = closest.x() - x;
    double dy = closest.y() - y;
    double dz = closest.z() - z;
    double disp = std::sqrt(dx*dx + dy*dy + dz*dz);
    if (!std::isfinite(disp)) {
      ++diag_edge_failures_;
      throw std::invalid_argument("Non-finite projection for Netgen geometric edge " + std::to_string(edgenr));
    }
    if (disp > diag_edge_max_disp_) {
      diag_edge_max_disp_ = disp;
      const auto & b = edge_bindings.lookup(edgenr, surfnr1, surfnr2);
      std::string bodies;
      DLIList<Body*> parents;
      re->bodies(parents);
      for (int i = 0; i < parents.size(); ++i) {
        Body* body = parents.get_and_step();
        if (body) bodies += (bodies.empty() ? "" : ",") + std::to_string(body->id());
      }
      char worst[2048];
      std::snprintf(worst, sizeof(worst), "NetgenCurver: worst edge projection edgenr=%d descriptor=%d curve=%d "
                 "surfaces=%d,%d bodies=[%s] in=(%.17g,%.17g,%.17g) "
                 "out=(%.17g,%.17g,%.17g) displacement=%.17g [model length units]\n",
                 edgenr, b.descriptor, b.curve,
                 ng_fd_to_cubit_sid_.at(surfnr1), ng_fd_to_cubit_sid_.at(surfnr2), bodies.c_str(),
                 x, y, z, closest.x(), closest.y(), closest.z(), disp);
      diag_edge_worst_ = worst;
    }
    if (disp > edge_reject_threshold_) {
      diag_edge_rejects_++;
      ++diag_edge_failures_;
      throw std::runtime_error("Curve projection displacement exceeds unchanged limit for Cubit curve " +
          std::to_string(re->id()) + " edge " + std::to_string(edgenr));
    }
    return {closest.x(), closest.y(), closest.z()};
  };

  // Segment-bounded interpolation on the shared curve.
  // Periodic curves require special handling: Cubit normalizes u values to a
  // fundamental period, so an edge crossing the seam can otherwise look like
  // the almost-complete curve. Work on the nearest unwrapped branch and use a
  // signed, endpoint-relative arc length for the fallback.
  auto between_edge = [this, curve_for_edge, edge_bindings](int edgenr, int surfnr1, int surfnr2,
                                                   double x1, double y1, double z1,
                                                   double x2, double y2, double z2,
                                                   double secpoint)
    -> std::tuple<double,double,double>
  {
    auto fail_curve = [&](const char* reason) -> std::tuple<double,double,double> {
      ++diag_edge_failures_;
      const auto & b = edge_bindings.lookup(edgenr, surfnr1, surfnr2);
      throw std::runtime_error(std::string(reason) + " for Cubit curve " +
          std::to_string(b.curve) + " edge " + std::to_string(edgenr));
    };

    RefEdge* re = curve_for_edge(edgenr, surfnr1, surfnr2);

    double crv_length = re->measure();
    if (!std::isfinite(crv_length) || crv_length <= 0.0) return fail_curve("Invalid curve length");

    CubitVector cv1(x1, y1, z1), cv2(x2, y2, z2);
    double u1 = re->u_from_position(cv1);
    double u2 = re->u_from_position(cv2);

    double period = 0.0;
    const bool periodic =
        re->is_periodic(period) == CUBIT_TRUE && period > 0.0;
    auto unwrap_near = [period](double reference, double u) {
      return u + period * std::round((reference - u) / period);
    };
    double u2_local = periodic ? unwrap_near(u1, u2) : u2;

    // Primary path: chord-perpendicular projection on the curve. This
    // matches Netgen's uniform-t LSQ curving convention. closest_point_trimmed
    // can jump to the wrong segment of a wavy curve, so we validate the
    // result by checking that the returned point's curve parameter falls
    // inside this segment on the same local parameter branch. If yes, use it.
    // Otherwise fall back to endpoint-relative arc-length interpolation.
    double lmx = x1 + secpoint * (x2 - x1);
    double lmy = y1 + secpoint * (y2 - y1);
    double lmz = z1 + secpoint * (z2 - z1);
    CubitVector loc(lmx, lmy, lmz);
    CubitVector cp;
    if (re->get_curve_ptr()->closest_point_trimmed(loc, cp) != CUBIT_SUCCESS)
      return fail_curve("closest_point_trimmed interpolation failed");
    double u_cp = re->u_from_position(cp);
    const double u_expected = u1 + secpoint * (u2_local - u1);
    double u_cp_local = periodic ? unwrap_near(u_expected, u_cp) : u_cp;
    double umin = std::min(u1, u2_local);
    double umax = std::max(u1, u2_local);
    // Allow a tiny tolerance (1% of range) to cover floating-point noise
    // at the exact endpoint.
    double tol = 0.01 * std::fabs(umax - umin);
    double chord = std::sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) +
                             (z2-z1)*(z2-z1));
    double cp1 = std::sqrt((cp.x()-x1)*(cp.x()-x1) +
                           (cp.y()-y1)*(cp.y()-y1) +
                           (cp.z()-z1)*(cp.z()-z1));
    double cp2 = std::sqrt((cp.x()-x2)*(cp.x()-x2) +
                           (cp.y()-y2)*(cp.y()-y2) +
                           (cp.z()-z2)*(cp.z()-z2));
    const double endpoint_tol = 1.0e-8 * chord + 1.0e-15;
    const bool stuck_at_endpoint =
        secpoint > 1.0e-8 && secpoint < 1.0 - 1.0e-8 &&
        (cp1 <= endpoint_tol || cp2 <= endpoint_tol);
    if (!stuck_at_endpoint &&
        u_cp_local >= umin - tol && u_cp_local <= umax + tol) {
      // chord-perpendicular projection lies inside the segment — use it
      return {cp.x(), cp.y(), cp.z()};
    }

    CubitVector out;
    CubitStatus st;
    if (periodic) {
      // RefEdge::point_from_arc_length normalizes the periodic root before
      // applying the distance and can select the long trimmed branch at the
      // seam. position_from_u explicitly accepts an unwrapped parameter and
      // normalizes only after the local branch has been selected.
      st = re->position_from_u(u_expected, out);
    } else {
      double segment_length = re->length_from_u(u1, u2_local);
      st = re->point_from_arc_length(
          u1, secpoint * segment_length, out);
    }
    if (st != CUBIT_SUCCESS) return fail_curve("Arc interpolation failed");

    // Sanity: the arc-midpoint should not lie farther from the linear
    // midpoint than the segment's chord length.
    double dd = std::sqrt((out.x()-lmx)*(out.x()-lmx) +
                          (out.y()-lmy)*(out.y()-lmy) +
                          (out.z()-lmz)*(out.z()-lmz));
    if (!std::isfinite(dd) || dd > chord) return fail_curve("Arc interpolation exceeds unchanged chord limit");

    return {out.x(), out.y(), out.z()};
  };

  auto tangent = [this, curve_for_edge](int edgenr, int face1, int face2,
                                  double x, double y, double z)
      -> std::tuple<double,double,double> {
    RefEdge* re = curve_for_edge(edgenr, face1, face2);
    CubitVector t;
    re->tangent(CubitVector(x, y, z), t);
    const double length = std::sqrt(t.x()*t.x() + t.y()*t.y() + t.z()*t.z());
    if (!std::isfinite(length) || !(length > 0.0)) {
      ++diag_edge_failures_;
      throw std::invalid_argument("Invalid tangent for Netgen geometric edge " + std::to_string(edgenr));
    }
    return {t.x()/length, t.y()/length, t.z()/length};
  };

  // between_edge: hybrid chord-perpendicular (primary) + arc-length
  // midpoint (fallback when closest-point jumps to wrong curve segment).
  // Using between_edge consistently gives Netgen segment-bounded curve
  // projections, which is what its uniform-t LSQ curving expects on
  // per-segment curving.
  auto geom = std::make_shared<ng::CallbackGeometry>(
      project, normal, num_surfaces, tangent, edge_project, between_edge);
  ng_mesh_->SetGeometry(geom);

  return true;
}

// ============================================================
// Phase 3: Curve and extract high-order node positions
// ============================================================
bool NetgenCurver::curve_and_extract(int order)
{
  // Curve the mesh
  ng_mesh_->BuildCurvedElements(order);

  auto & curved = ng_mesh_->GetCurvedElements();
  if (!curved.IsHighOrder()) {
    PRINT_WARNING("NetgenCurver: BuildCurvedElements did not produce high-order data.\n");
    return false;
  }

  radia_netgen_curving::RefitStats rs = radia_netgen_curving::geometric_refit(*ng_mesh_);
  PRINT_INFO("NetgenCurver: geometric refit - edges %d/%d accepted "
             "(%d failed, %d kept to avoid distortion), max distance %.3e -> %.3e; "
             "faces %d/%d accepted (%d failed, %d kept to avoid distortion), "
             "max distance %.3e -> %.3e\n",
             rs.edges_accepted, rs.edges_tried, rs.edges_failed,
             rs.edges_rejected_distortion,
             rs.edge_dist_before, rs.edge_dist_after,
             rs.faces_accepted, rs.faces_tried, rs.faces_failed,
             rs.faces_rejected_distortion,
             rs.face_dist_before, rs.face_dist_after);

  // Pyramids are evaluated correctly only by a Netgen 6.2.2608 runtime with
  // the apex clamp and the base-edge blend in the vectorised map; this
  // snapshot carries both, the stock wheel neither.  check-vol rejects
  // pyramid meshes on an unpatched runtime, so say so at export time.
  {
    const auto & top = ng_mesh_->GetTopology();
    int curved_pyramids = 0;
    for (int e = 0; e < ng_mesh_->GetNE(); e++) {
      const auto ei = ng::ElementIndex::FromNr0(e);
      if ((*ng_mesh_)[ei].GetType() != ng::PYRAMID) continue;
      for (auto edge : top.GetEdges(ei))
        if (radia_netgen_curving::CurvedAccess{curved}.edge_count(edge.Nr0()) > 0) { curved_pyramids++; break; }
    }
    if (curved_pyramids > 0)
      PRINT_WARNING("NetgenCurver: %d pyramids have curved edges. Read this "
                    ".vol only with the Netgen 6.2.2608 runtime that carries "
                    "the pyramid apex clamp and base-edge blend; check-vol "
                    "rejects it on any other runtime.\n",
                    curved_pyramids);
  }

  return true;
}


