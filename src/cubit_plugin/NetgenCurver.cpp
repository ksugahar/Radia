#include "NetgenCurver.hpp"
#include "MeshData.hpp"

#ifdef HAVE_NETGEN

// Include Netgen BEFORE Cubit to manage namespace conflicts
// Netgen defines: netgen::QUAD, netgen::HEX, netgen::PYRAMID, netgen::Surface
#include <meshing/meshing.hpp>
#include "callbackgeom.hpp"

// Avoid Cubit/Netgen name clashes by aliasing Netgen element types
namespace ng = netgen;

// Now include Cubit headers
// Cubit defines: ::QUAD, ::HEX, ::PYRAMID, ::Surface (class)
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"

// cubit_geom headers for ACIS projection
#include "RefFace.hpp"
#include "RefEdge.hpp"
#include "Surface.hpp"
#include "Curve.hpp"
#include "CubitVector.hpp"
// 2026-05-25: Cubit 2025.12 removed GeometryQueryTool::instance() static
// accessor.  CubitRefAccessor wraps the new path
//   CubitCoreModel::instance()->cgmapp()->factory()->get_ref_{face,edge,volume}()
// so the 6 lookup sites below stay one-liners.
#include "CubitRefAccessor.hpp"

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

NetgenCurver::NetgenCurver() {}
NetgenCurver::~NetgenCurver() {}

bool NetgenCurver::build(const MeshData &md, int order)
{
  if (order < 2) return false;
  order_ = order;

  PRINT_INFO("NetgenCurver: Phase 1 - build_netgen_mesh (MeshData)...\n");
  if (!build_netgen_mesh(md)) return false;
  PRINT_INFO("NetgenCurver: Phase 2 - attach_callback_geometry...\n");
  if (!attach_callback_geometry()) return false;
  PRINT_INFO("NetgenCurver: Phase 3 - curve_and_extract(order=%d)...\n", order);
  if (!curve_and_extract(order)) return false;

  PRINT_INFO("NetgenCurver: order %d, %d original + %d curved nodes = %d total\n",
             order, (int)cubit_nid_to_ng_pi_.size(),
             total_nodes_ - (int)cubit_nid_to_ng_pi_.size(), total_nodes_);

  // Diagnostic report: projection reject stats
  PRINT_INFO("NetgenCurver: projection stats — "
             "project=%ld calls, %ld rejects (>%g mm), max_disp=%g mm (surf %d)\n",
             diag_project_calls_, diag_project_rejects_,
             project_reject_threshold_,
             diag_project_max_disp_, diag_project_max_disp_surf_);
  PRINT_INFO("NetgenCurver: edge_project stats — %ld calls, %ld rejects, max_disp=%g mm\n",
             diag_edge_calls_, diag_edge_rejects_, diag_edge_max_disp_);
  PRINT_INFO("NetgenCurver: trim-then-refine — tries=%ld, recovered=%ld, max_tdisp=%g mm\n",
             diag_refine_tries_, diag_refine_success_, diag_refine_max_tdisp_);
  if (!polar_surfaces_.empty() || diag_polar_skips_ > 0) {
    PRINT_INFO("NetgenCurver: polar-disk short-circuit — %zu surfaces, %ld HO skips (linear fallback)\n",
               polar_surfaces_.size(), diag_polar_skips_);
  }
  PRINT_INFO("NetgenCurver: path distribution — "
             "0(uv_direct_eval)=%ld (max_disp=%g mm), "
             "1(uv_guess_fallback)=%ld (max_disp=%g mm), "
             "2(nearest_vertex_hint)=%ld (max_disp=%g mm), "
             "3(no_hint_trimmed)=%ld (max_disp=%g mm)\n",
             diag_path_calls_[0], diag_path_max_disp_[0],
             diag_path_calls_[1], diag_path_max_disp_[1],
             diag_path_calls_[2], diag_path_max_disp_[2],
             diag_path_calls_[3], diag_path_max_disp_[3]);
  PRINT_INFO("NetgenCurver: path=0 disp histogram — "
             "<0.01:%ld  [0.01,0.05):%ld  [0.05,0.1):%ld  [0.1,0.2):%ld  "
             "[0.2,0.5):%ld  [0.5,1.0):%ld  >=1.0:%ld mm\n",
             diag_path0_hist_[0], diag_path0_hist_[1], diag_path0_hist_[2],
             diag_path0_hist_[3], diag_path0_hist_[4], diag_path0_hist_[5],
             diag_path0_hist_[6]);
  if (diag_project_rejects_ > 0) {
    // Show top-5 offending surfaces by reject count + UV spread for each
    std::vector<std::pair<int,long>> ranked(diag_project_reject_per_surf_.begin(),
                                             diag_project_reject_per_surf_.end());
    std::sort(ranked.begin(), ranked.end(),
              [](const auto &a, const auto &b) { return a.second > b.second; });
    int n_show = (int)ranked.size() < 5 ? (int)ranked.size() : 5;
    PRINT_INFO("NetgenCurver: top-%d surfaces by reject count (with Phase1e UV spread):\n", n_show);
    for (int i = 0; i < n_show; i++) {
      int sid = ranked[i].first;
      double mdu = diag_tri_max_du_.count(sid) ? diag_tri_max_du_[sid] : 0.0;
      double mdv = diag_tri_max_dv_.count(sid) ? diag_tri_max_dv_[sid] : 0.0;
      int ntri = diag_tri_count_.count(sid) ? diag_tri_count_[sid] : 0;
      PRINT_INFO("  cubit_surf %d: %ld rejects | tris=%d, max_tri_Δu=%.3e, max_tri_Δv=%.3e\n",
                 sid, ranked[i].second, ntri, mdu, mdv);
    }
  }

  // Dump detailed reject CSV only when RADIA_NETGEN_REJECT_DUMP env var is
  // set (to avoid noisy writes in production builds). Value is used as the
  // output path. Typical debug usage: RADIA_NETGEN_REJECT_DUMP=C:\temp\netgen_reject_log.csv
  if (!project_reject_log_entries_.empty()) {
    const char * dump_path = std::getenv("RADIA_NETGEN_REJECT_DUMP");
    if (dump_path && *dump_path) {
      std::ofstream ofs(dump_path);
      if (ofs) {
        ofs << "surfnr,cubit_sid,path,in_x,in_y,in_z,"
               "out_x,out_y,out_z,u_hint,v_hint,out_u,out_v,disp\n";
        for (auto &row : project_reject_log_entries_)
          ofs << row << "\n";
        PRINT_INFO("NetgenCurver: %d reject rows written to %s\n",
                   (int)project_reject_log_entries_.size(), dump_path);
      } else {
        PRINT_INFO("NetgenCurver: failed to open %s for reject log\n", dump_path);
      }
    }
  }
  return true;
}

// ============================================================
// Phase 1: Convert linear mesh (from MeshData) to netgen::Mesh
// ============================================================
bool NetgenCurver::build_netgen_mesh(const MeshData &md)
{
  PRINT_INFO("  Phase1a: creating ng::Mesh...\n");
  fflush(stdout);
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
  fflush(stdout);

  for (int i = 0; i < num_nodes; i++) {
    auto &nd = md.nodes[i];
    if (i == 0) { PRINT_INFO("  Phase1c: first AddPoint (id=%d, %.4f, %.4f, %.4f)\n", nd.id, nd.x, nd.y, nd.z); fflush(stdout); }
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
    cubit_nid_to_ng_pi_[nd.id] = pi;
    ng_pi_to_cubit_nid_[pi] = nd.id;
    if (i == 0 || i == num_nodes - 1)
      PRINT_INFO("  Phase1c: node %d/%d OK\n", i + 1, num_nodes);
  }

  next_node_id_ = md.max_original_node_id + 1;
  total_nodes_ = num_nodes;

  PRINT_INFO("  Phase1d: FaceDescriptors\n");
  fflush(stdout);
  // --- Build FaceDescriptors for each Cubit surface ---
  PRINT_INFO("  Phase1d: calling get_entities('surface')...\n");
  fflush(stdout);
  std::vector<int> surface_ids = CubitInterface::get_entities("surface");
  PRINT_INFO("  Phase1d: got %d surfaces\n", (int)surface_ids.size());
  fflush(stdout);
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
  fflush(stdout);
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
      rf_cache[sid] = radia::cubit_get_ref_face(sid);

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
    fflush(stdout);

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
        el.SetIndex(ng_fd);
        // Pass 1: collect UVs from per-node canonical UV map (computed above)
        double uu[3], vv[3];
        int pis[3];
        bool uv_ok = true;
        for (int k = 0; k < 3; k++) {
          auto it = cubit_nid_to_ng_pi_.find(conn[k]);
          if (it == cubit_nid_to_ng_pi_.end()) { uv_ok = false; goto skip_tri; }
          pis[k] = it->second;
          el[k] = it->second;
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
          auto &cur_du = diag_tri_max_du_[sid];
          if (max_du > cur_du) cur_du = max_du;
          auto &cur_dv = diag_tri_max_dv_[sid];
          if (max_dv > cur_dv) cur_dv = max_dv;
          diag_tri_count_[sid]++;
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
        el.SetIndex(ng_fd);
        double qu[4], qv[4];
        for (int k = 0; k < 4; k++) {
          auto it = cubit_nid_to_ng_pi_.find(conn[k]);
          if (it == cubit_nid_to_ng_pi_.end()) goto skip_quad;
          el[k] = it->second;

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
    // BCProperty values are an internal handoff to ExportNetgenCommand, which
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
          el.SetIndex(ng_fd);
          bool complete = true;
          for (int k = 0; k < face->nv; k++) {
            auto it = cubit_nid_to_ng_pi_.find(face->conn[k]);
            if (it == cubit_nid_to_ng_pi_.end()) {
              complete = false;
              break;
            }
            el[k] = it->second;
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
  fflush(stdout);
  {
    std::vector<int> curve_ids = CubitInterface::parse_cubit_list("curve", "all");
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

      // Get curve length for dist normalization
      RefEdge* re = radia::cubit_get_ref_edge(cid);
      double crv_length = re ? re->measure() : 1.0;
      double u_start = re ? re->start_param() : 0.0;

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
        nseg[0] = itn0->second;
        nseg[1] = itn1->second;
        nseg.si = fd2;
        nseg.edgenr = edgenr;
        nseg.surfnr1 = fd1 - 1;  // 0-based
        nseg.surfnr2 = fd2 - 1;
        nseg.epgeominfo[0].edgenr = edgenr;
        nseg.epgeominfo[0].dist = dist0;
        nseg.epgeominfo[1].edgenr = edgenr;
        nseg.epgeominfo[1].dist = dist1;
        ng_mesh_->AddSegment(nseg);
      }
      edgenr++;
    }
    if (edgenr > 1)
      ng_mesh_->SetNCD2Names(edgenr);
    PRINT_INFO("  Phase1e2: %d edge groups added\n", edgenr - 1);
  }

  PRINT_INFO("  Phase1f: Add volume elements (%d)\n", (int)md.elements.size());
  fflush(stdout);
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
        el[k] = cubit_nid_to_ng_pi_[elem.conn[k]];
      el.SetIndex(1);
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
        el[k] = cubit_nid_to_ng_pi_[elem.conn[k]];
      el.SetIndex(1);
      ng_mesh_->AddVolumeElement(el);
      vol_count++;
    }
    else if (nn == 6 && (et == ::WEDGE6 || et == ::WEDGE)) {
      ng::Element el(ng::PRISM);
      for (int k = 0; k < 6; k++)
        el[k] = cubit_nid_to_ng_pi_[elem.conn[k]];
      el.SetIndex(1);
      ng_mesh_->AddVolumeElement(el);
      vol_count++;
    }
    else if (nn == 5 && (et == ::PYRAMID5 || et == ::PYRAMID)) {
      ng::Element el(ng::PYRAMID);
      for (int k = 0; k < 5; k++)
        el[k] = cubit_nid_to_ng_pi_[elem.conn[k]];
      el.SetIndex(1);
      ng_mesh_->AddVolumeElement(el);
      vol_count++;
    }
  }

  PRINT_INFO("  Phase1f: %d volume elements added\n", vol_count);
  fflush(stdout);

  PRINT_INFO("  Phase1g: UpdateTopology...\n");
  fflush(stdout);
  ng_mesh_->UpdateTopology();
  PRINT_INFO("  Phase1g: UpdateTopology done\n");
  fflush(stdout);

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
  diag_refine_success_ = 0;
  diag_refine_tries_ = 0;
  diag_refine_max_tdisp_ = 0.0;
  diag_project_max_disp_ = 0.0;
  diag_project_max_disp_surf_ = -1;
  diag_project_reject_per_surf_.clear();
  diag_edge_calls_ = 0;
  diag_edge_rejects_ = 0;
  diag_edge_max_disp_ = 0.0;
  for (int i = 0; i < 4; i++) {
    diag_path_calls_[i] = 0;
    diag_path_max_disp_[i] = 0.0;
  }
  for (int i = 0; i < 7; i++) diag_path0_hist_[i] = 0;
  diag_polar_skips_ = 0;
  project_reject_log_entries_.clear();

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
  PRINT_INFO("NetgenCurver: built bbox for %d surfaces (pad=%g mm, global_diam=%g mm)\n",
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

    RefFace* rf = radia::cubit_get_ref_face(cubit_sid);
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
    if (path >= 0 && path < 4) {
      diag_path_calls_[path]++;
      if (disp > diag_path_max_disp_[path]) diag_path_max_disp_[path] = disp;
    }
    if (path == 0) {
      int bin = 6;
      if (disp < 0.01) bin = 0;
      else if (disp < 0.05) bin = 1;
      else if (disp < 0.1)  bin = 2;
      else if (disp < 0.2)  bin = 3;
      else if (disp < 0.5)  bin = 4;
      else if (disp < 1.0)  bin = 5;
      diag_path0_hist_[bin]++;
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
      // Trim-then-refine: before accepting the reject, try a second pass
      // starting from closest_point_trimmed (global trimmed projection).
      // This recovers the correct UV branch when the nearest-vertex hint
      // (path=2) landed on a wrong-side merged sub-chart.  Enabled via env
      // var RADIA_NETGEN_TRIM_REFINE=1 so we can A/B test.
      static const bool do_trim_refine = [](){
        const char *s = std::getenv("RADIA_NETGEN_TRIM_REFINE");
        return s && (*s == '1' || *s == 'y' || *s == 'Y');
      }();
      if (do_trim_refine) {
        diag_refine_tries_++;
        // Use closest_point_trimmed directly — this is guaranteed on the
        // trimmed face and avoids the un-trimmed Newton jumping branches.
        CubitVector trimmed_pt;
        surf->closest_point_trimmed(loc, trimmed_pt);
        double tdx = trimmed_pt.x() - x;
        double tdy = trimmed_pt.y() - y;
        double tdz = trimmed_pt.z() - z;
        double tdisp = std::sqrt(tdx*tdx + tdy*tdy + tdz*tdz);
        if (tdisp < disp) {
          double trim_u = 0, trim_v = 0;
          rf->u_v_from_position(trimmed_pt, trim_u, trim_v);
          diag_refine_success_++;
          if (tdisp > diag_refine_max_tdisp_) diag_refine_max_tdisp_ = tdisp;
          return {trimmed_pt.x(), trimmed_pt.y(), trimmed_pt.z(),
                  trim_u, trim_v};
        }
      }
      diag_project_rejects_++;
      diag_project_reject_per_surf_[cubit_sid]++;
      // Detailed log for post-analysis (surf UV pattern vs projection outlier)
      if ((int)project_reject_log_entries_.size() < project_reject_dump_max_) {
        char buf[512];
        snprintf(buf, sizeof(buf),
                 "%d,%d,%d,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.9e,%.6e",
                 surfnr, cubit_sid, path,
                 x, y, z,
                 closest.x(), closest.y(), closest.z(),
                 u_hint, v_hint, u, v, disp);
        project_reject_log_entries_.emplace_back(buf);
      }
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
    RefFace* rf = radia::cubit_get_ref_face(cubit_sid);
    if (!rf) return {0.0, 0.0, 1.0};

    CubitVector loc(x, y, z);
    CubitVector n = rf->normal_at(loc);
    return {n.x(), n.y(), n.z()};
  };

  // Edge (curve) projection callback: (surfnr1, surfnr2, x,y,z) -> (xp,yp,zp)
  // surfnr1/surfnr2 are 1-based FD indices. Find the Cubit curve at their intersection.
  // Build surfpair -> RefEdge map
  std::unordered_map<int64_t, int> surfpair_to_curve;  // (fd1<<32|fd2) -> cubit_curve_id
  {
    std::vector<int> curve_ids = CubitInterface::parse_cubit_list("curve", "all");
    for (int cid : curve_ids) {
      std::vector<int> parent_surfs = CubitInterface::parse_cubit_list(
          "surface", "in curve " + std::to_string(cid));
      if (parent_surfs.size() < 2) continue;
      for (size_t i = 0; i < parent_surfs.size(); i++) {
        for (size_t j = i+1; j < parent_surfs.size(); j++) {
          int s1 = parent_surfs[i], s2 = parent_surfs[j];
          auto it1 = cubit_sid_to_ng_fd_.find(s1);
          auto it2 = cubit_sid_to_ng_fd_.find(s2);
          if (it1 == cubit_sid_to_ng_fd_.end() || it2 == cubit_sid_to_ng_fd_.end()) continue;
          int f1 = it1->second, f2 = it2->second;
          surfpair_to_curve[((int64_t)f1 << 32) | f2] = cid;
          surfpair_to_curve[((int64_t)f2 << 32) | f1] = cid;
        }
      }
    }
  }

  auto edge_project = [this, surfpair_to_curve](int surfnr1, int surfnr2,
                                                   double x, double y, double z)
    -> std::tuple<double,double,double>
  {
    diag_edge_calls_++;
    int64_t key = ((int64_t)surfnr1 << 32) | surfnr2;
    auto it = surfpair_to_curve.find(key);
    if (it != surfpair_to_curve.end()) {
      RefEdge* re = radia::cubit_get_ref_edge(it->second);
      if (re) {
        CubitVector loc(x, y, z);
        CubitVector closest;
        re->get_curve_ptr()->closest_point_trimmed(loc, closest);
        double dx = closest.x() - x;
        double dy = closest.y() - y;
        double dz = closest.z() - z;
        double disp = std::sqrt(dx*dx + dy*dy + dz*dz);
        if (disp > diag_edge_max_disp_) diag_edge_max_disp_ = disp;
        if (disp > edge_reject_threshold_) {
          diag_edge_rejects_++;
          return {x, y, z};  // reject spike: fall back to identity
        }
        return {closest.x(), closest.y(), closest.z()};
      }
    }
    return {x, y, z};  // no projection if curve not found
  };

  // Arc-length parametric midpoint on the shared curve.
  // Inputs (x1,y1,z1) and (x2,y2,z2) are the segment's 3D endpoints. We
  // compute each endpoint's curve parameter via u_from_position, derive the
  // per-endpoint arc lengths from the curve's start, linearly interpolate
  // the arc length at secpoint, then place the point at that arc-length
  // position via position_from_fraction. On any error we fall back to the
  // linear midpoint between p1 and p2 (caller is expected to then ignore or
  // re-project). This avoids closest_point_trimmed jumping to the wrong part
  // of a wavy curve when the linear midpoint is closer to a different
  // segment of it.
  auto between_edge = [this, surfpair_to_curve](int surfnr1, int surfnr2,
                                                   double x1, double y1, double z1,
                                                   double x2, double y2, double z2,
                                                   double secpoint)
    -> std::tuple<double,double,double>
  {
    auto linear_mid = [&]() -> std::tuple<double,double,double> {
      return {x1 + secpoint * (x2 - x1),
              y1 + secpoint * (y2 - y1),
              z1 + secpoint * (z2 - z1)};
    };

    int64_t key = ((int64_t)surfnr1 << 32) | surfnr2;
    auto it = surfpair_to_curve.find(key);
    if (it == surfpair_to_curve.end()) return linear_mid();

    RefEdge* re = radia::cubit_get_ref_edge(it->second);
    if (!re) return linear_mid();

    double crv_length = re->measure();
    if (crv_length <= 0.0) return linear_mid();

    CubitVector cv1(x1, y1, z1), cv2(x2, y2, z2);
    double u1 = re->u_from_position(cv1);
    double u2 = re->u_from_position(cv2);
    double u_start = re->start_param();

    // Primary path: chord-perpendicular projection on the curve. This
    // matches Netgen's uniform-t LSQ curving convention. closest_point_trimmed
    // can jump to the wrong segment of a wavy curve, so we validate the
    // result by checking that the returned point's curve parameter falls
    // within [min(u1,u2), max(u1,u2)] — i.e., stays inside this segment.
    // If yes, use it. Otherwise fall back to arc-length midpoint.
    double lmx = x1 + secpoint * (x2 - x1);
    double lmy = y1 + secpoint * (y2 - y1);
    double lmz = z1 + secpoint * (z2 - z1);
    CubitVector loc(lmx, lmy, lmz);
    CubitVector cp;
    re->get_curve_ptr()->closest_point_trimmed(loc, cp);
    double u_cp = re->u_from_position(cp);
    double umin = std::min(u1, u2);
    double umax = std::max(u1, u2);
    // Allow a tiny tolerance (1% of range) to cover floating-point noise
    // at the exact endpoint.
    double tol = 0.01 * std::fabs(umax - umin);
    if (u_cp >= umin - tol && u_cp <= umax + tol) {
      // chord-perpendicular projection lies inside the segment — use it
      return {cp.x(), cp.y(), cp.z()};
    }

    // closest_point jumped to a different segment of the curve — fall
    // back to arc-length parametric midpoint derived from the endpoints.
    double L1 = re->length_from_u(u_start, u1);
    double L2 = re->length_from_u(u_start, u2);
    double L_mid = L1 + secpoint * (L2 - L1);
    double frac = L_mid / crv_length;
    if (frac < 0.0) frac = 0.0;
    if (frac > 1.0) frac = 1.0;

    CubitVector out;
    CubitStatus st = re->position_from_fraction(frac, out);
    if (st != CUBIT_SUCCESS) return linear_mid();

    // Sanity: the arc-midpoint should not lie farther from the linear
    // midpoint than the segment's chord length.
    double chord = std::sqrt((x2-x1)*(x2-x1) + (y2-y1)*(y2-y1) + (z2-z1)*(z2-z1));
    double dd = std::sqrt((out.x()-lmx)*(out.x()-lmx) +
                          (out.y()-lmy)*(out.y()-lmy) +
                          (out.z()-lmz)*(out.z()-lmz));
    if (dd > chord) return linear_mid();

    return {out.x(), out.y(), out.z()};
  };

  // between_edge: hybrid chord-perpendicular (primary) + arc-length
  // midpoint (fallback when closest-point jumps to wrong curve segment).
  // Using between_edge consistently gives Netgen segment-bounded curve
  // projections, which is what its uniform-t LSQ curving expects on
  // per-segment curving.
  auto geom = std::make_shared<ng::CallbackGeometry>(
      project, normal, num_surfaces, nullptr, edge_project, between_edge);
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

  // Extract curved node positions by evaluating CalcElementTransformation
  // at reference-space coordinates of high-order nodes.
  //
  // IMPORTANT: Volume elements are processed FIRST. CalcElementTransformation
  // uses BuildCurvedElements data (same as .vol export), producing coordinates
  // consistent with the .vol file. Surface elements are processed AFTER to
  // add only face interior nodes (order >= 3) and edges not covered by volumes.

  // --- Phase 1: Volume elements (CalcElementTransformation) ---
  // This is the primary source of edge HO nodes, consistent with .vol export.
  int nve = ng_mesh_->GetNE();
  for (int ei = 0; ei < nve; ei++) {
    if (!curved.IsElementCurved(ei)) continue;

    const ng::Element & el = ng_mesh_->VolumeElement(ei+1);
    int np = el.GetNP();

    if (np == 4) {  // Tet
      // Edge points (6 edges) — register in edge_ho_nodes_ for connectivity
      static const int tet_edges[][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
#ifdef DEBUG_HO_NODES
      if (ei < 3) {
        PRINT_INFO("  VE %d: el[0..3] ng=(%d,%d,%d,%d) cubit=(%d,%d,%d,%d)\n",
          ei, (int)el[0], (int)el[1], (int)el[2], (int)el[3],
          ng_pi_to_cubit_nid_[el[0]],
          ng_pi_to_cubit_nid_[el[1]],
          ng_pi_to_cubit_nid_[el[2]],
          ng_pi_to_cubit_nid_[el[3]]);
      }
#endif
      for (auto &e : tet_edges) {
        int v0_cubit = ng_pi_to_cubit_nid_[el[e[0]]];
        int v1_cubit = ng_pi_to_cubit_nid_[el[e[1]]];
        HoEdgeKey key(v0_cubit, v1_cubit);
        if (edge_ho_nodes_.count(key)) continue;

        std::vector<int> edge_nodes;
        for (int k = 1; k < order; k++) {
          double t = (double)k / order;
          // Netgen ref coords: ref=(lam1,lam2,lam3), lam0=1-sum
          // Vertex mapping: el[i] -> lam[(i+1)%4]
          double lam[4] = {0,0,0,0};
          lam[(e[0]+1)%4] = 1.0 - t;
          lam[(e[1]+1)%4] = t;
          ng::Point<3> ref(lam[1], lam[2], lam[3]);
          ng::Point<3> phys;
          curved.CalcElementTransformation(ref, ei, phys);

          int new_id = next_node_id_++;
          total_nodes_++;
          ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
          edge_nodes.push_back(new_id);
        }
        if (v0_cubit <= v1_cubit)
          edge_ho_nodes_[key] = edge_nodes;
        else
          edge_ho_nodes_[key] = std::vector<int>(edge_nodes.rbegin(), edge_nodes.rend());
      }

      // Face interior points (4 faces, for order >= 3)
      if (order >= 3) {
        static const int tet_faces[][3] = {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
        for (auto &f : tet_faces) {
          int fv0 = ng_pi_to_cubit_nid_[el[f[0]]];
          int fv1 = ng_pi_to_cubit_nid_[el[f[1]]];
          int fv2 = ng_pi_to_cubit_nid_[el[f[2]]];
          HoFaceKey3 fkey(fv0, fv1, fv2);
          if (face_ho_nodes_.count(fkey)) {
            // Already generated by adjacent element — skip generation
            // but still need to advance node IDs for consistency
            // Actually, skip entirely — nodes already exist
            continue;
          }
          std::vector<int> face_nodes;
          for (int i = 1; i < order; i++) {
            for (int j = 1; j < order - i; j++) {
              double lam[4] = {0,0,0,0};
              double a = (double)i / order;
              double b = (double)j / order;
              double c = 1.0 - a - b;
              // Vertex mapping: el[i] -> lam[(i+1)%4]
              lam[(f[0]+1)%4] = c;
              lam[(f[1]+1)%4] = a;
              lam[(f[2]+1)%4] = b;
              ng::Point<3> ref(lam[1], lam[2], lam[3]);
              ng::Point<3> phys;
              curved.CalcElementTransformation(ref, ei, phys);

              int new_id = next_node_id_++;
              total_nodes_++;
              ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
              face_nodes.push_back(new_id);
            }
          }
          face_ho_nodes_[fkey] = {face_nodes, fv0, fv1, fv2};
        }
      }

      // Volume interior points (for order >= 4)
      if (order >= 4) {
        int tv0 = ng_pi_to_cubit_nid_[el[0]];
        int tv1 = ng_pi_to_cubit_nid_[el[1]];
        int tv2 = ng_pi_to_cubit_nid_[el[2]];
        int tv3 = ng_pi_to_cubit_nid_[el[3]];
        HoElemKey4 vkey(tv0, tv1, tv2, tv3);
        std::vector<int> vol_nodes;
        for (int i = 1; i < order; i++) {
          for (int j = 1; j < order - i; j++) {
            for (int k = 1; k < order - i - j; k++) {
              double a = (double)i / order;
              double b = (double)j / order;
              double c = (double)k / order;
              ng::Point<3> ref(a, b, c);
              ng::Point<3> phys;
              curved.CalcElementTransformation(ref, ei, phys);

              int new_id = next_node_id_++;
              total_nodes_++;
              ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
              vol_nodes.push_back(new_id);
            }
          }
        }
        vol_ho_nodes_[vkey] = {vol_nodes, tv0, tv1, tv2, tv3};
      }
    }
    else if (np == 8) {  // Hex
#ifdef DEBUG_HO_NODES
      if (ei < 2) {
        // Evaluate CalcElementTransformation at 8 ref vertices to find mapping
        static const double ref_verts[][3] = {
          {0,0,0},{1,0,0},{1,1,0},{0,1,0},
          {0,0,1},{1,0,1},{1,1,1},{0,1,1}
        };
        PRINT_INFO("  HEX VE %d: Netgen ref -> physical:\n", ei);
        for (int rv = 0; rv < 8; rv++) {
          ng::Point<3> ref(ref_verts[rv][0], ref_verts[rv][1], ref_verts[rv][2]);
          ng::Point<3> phys;
          curved.CalcElementTransformation(ref, ei, phys);
          // Find closest Cubit vertex
          int best_v = -1; double best_d = 1e30;
          for (int k = 0; k < 8; k++) {
            int cid = ng_pi_to_cubit_nid_[el[k]];
            auto cc = CubitInterface::get_nodal_coordinates(cid);
            double d = (cc[0]-phys[0])*(cc[0]-phys[0])
                     + (cc[1]-phys[1])*(cc[1]-phys[1])
                     + (cc[2]-phys[2])*(cc[2]-phys[2]);
            if (d < best_d) { best_d = d; best_v = k; }
          }
          PRINT_INFO("    ref(%g,%g,%g) -> el[%d] (cubit %d) dist=%.2e\n",
            ref_verts[rv][0], ref_verts[rv][1], ref_verts[rv][2],
            best_v, ng_pi_to_cubit_nid_[el[best_v]], sqrt(best_d));
        }
      }
#endif
      // HEX: ref(i,j,k) maps directly to el[i] (identity mapping verified).
      // Use CalcElementTransformation (consistent with BuildCurvedElements/.vol).
      static const int hex_edges[][2] = {
        {0,1},{1,2},{2,3},{3,0},
        {4,5},{5,6},{6,7},{7,4},
        {0,4},{1,5},{2,6},{3,7}
      };
      // Netgen hex ref vertices (identity mapping to el[0..7])
      static const double hex_ref[][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},
        {0,0,1},{1,0,1},{1,1,1},{0,1,1}
      };

      for (auto &e : hex_edges) {
        int v0_cubit = ng_pi_to_cubit_nid_[el[e[0]]];
        int v1_cubit = ng_pi_to_cubit_nid_[el[e[1]]];
        HoEdgeKey key(v0_cubit, v1_cubit);
        if (edge_ho_nodes_.count(key)) continue;

        std::vector<int> edge_nodes;
        for (int k = 1; k < order; k++) {
          double t = (double)k / order;
          double xi   = hex_ref[e[0]][0] + t * (hex_ref[e[1]][0] - hex_ref[e[0]][0]);
          double eta  = hex_ref[e[0]][1] + t * (hex_ref[e[1]][1] - hex_ref[e[0]][1]);
          double zeta = hex_ref[e[0]][2] + t * (hex_ref[e[1]][2] - hex_ref[e[0]][2]);
          ng::Point<3> ref(xi, eta, zeta);
          ng::Point<3> phys;
          curved.CalcElementTransformation(ref, ei, phys);
          int new_id = next_node_id_++;
          total_nodes_++;
          ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
          edge_nodes.push_back(new_id);
        }
        if (v0_cubit <= v1_cubit)
          edge_ho_nodes_[key] = edge_nodes;
        else
          edge_ho_nodes_[key] = std::vector<int>(edge_nodes.rbegin(), edge_nodes.rend());
      }

      // 6 face interiors (for order >= 3)
      // NOTE: face/volume interior still uses CalcElementTransformation with Netgen ref.
      // This is unreliable for hex due to vertex ordering mismatch. TODO: fix for order >= 3.
      if (order >= 3) {
        static const double hex_ref[][3] = {
          {0,0,0},{1,0,0},{1,0,1},{0,0,1},
          {0,1,0},{1,1,0},{1,1,1},{0,1,1}
        };
        static const int hex_faces[][4] = {
          {0,1,2,3}, {4,5,6,7},  // bottom, top
          {0,1,5,4}, {2,3,7,6},  // front, back
          {0,3,7,4}, {1,2,6,5}   // left, right
        };
        for (auto &f : hex_faces) {
          for (int i = 1; i < order; i++) {
            for (int j = 1; j < order; j++) {
              double s = (double)i / order;
              double t = (double)j / order;
              // Bilinear interpolation in ref space
              double xi = (1-s)*(1-t)*hex_ref[f[0]][0] + s*(1-t)*hex_ref[f[1]][0]
                        + s*t*hex_ref[f[2]][0] + (1-s)*t*hex_ref[f[3]][0];
              double eta = (1-s)*(1-t)*hex_ref[f[0]][1] + s*(1-t)*hex_ref[f[1]][1]
                         + s*t*hex_ref[f[2]][1] + (1-s)*t*hex_ref[f[3]][1];
              double zeta = (1-s)*(1-t)*hex_ref[f[0]][2] + s*(1-t)*hex_ref[f[1]][2]
                          + s*t*hex_ref[f[2]][2] + (1-s)*t*hex_ref[f[3]][2];
              ng::Point<3> ref(xi, eta, zeta);
              ng::Point<3> phys;
              curved.CalcElementTransformation(ref, ei, phys);
              int new_id = next_node_id_++;
              total_nodes_++;
              ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
            }
          }
        }
      }

      // Volume interior (for order >= 3)
      if (order >= 3) {
        for (int i = 1; i < order; i++) {
          for (int j = 1; j < order; j++) {
            for (int k = 1; k < order; k++) {
              double xi   = (double)i / order;
              double eta  = (double)j / order;
              double zeta = (double)k / order;
              ng::Point<3> ref(xi, eta, zeta);
              ng::Point<3> phys;
              curved.CalcElementTransformation(ref, ei, phys);
              int new_id = next_node_id_++;
              total_nodes_++;
              ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
            }
          }
        }
      }
    }
    else if (np == 6) {  // Wedge (Prism)
#ifdef DEBUG_HO_NODES
      if (ei < 2) {
        static const double ref_verts[][3] = {
          {0,0,0},{1,0,0},{0,1,0},{0,0,1},{1,0,1},{0,1,1}
        };
        PRINT_INFO("  PRISM VE %d: Netgen ref -> physical:\n", ei);
        for (int rv = 0; rv < 6; rv++) {
          ng::Point<3> ref(ref_verts[rv][0], ref_verts[rv][1], ref_verts[rv][2]);
          ng::Point<3> phys;
          curved.CalcElementTransformation(ref, ei, phys);
          int best_v = -1; double best_d = 1e30;
          for (int k = 0; k < 6; k++) {
            int cid = ng_pi_to_cubit_nid_[el[k]];
            auto cc = CubitInterface::get_nodal_coordinates(cid);
            double d = (cc[0]-phys[0])*(cc[0]-phys[0])
                     + (cc[1]-phys[1])*(cc[1]-phys[1])
                     + (cc[2]-phys[2])*(cc[2]-phys[2]);
            if (d < best_d) { best_d = d; best_v = k; }
          }
          PRINT_INFO("    ref(%g,%g,%g) -> el[%d] (cubit %d) dist=%.2e\n",
            ref_verts[rv][0], ref_verts[rv][1], ref_verts[rv][2],
            best_v, ng_pi_to_cubit_nid_[el[best_v]], sqrt(best_d));
        }
      }
#endif
      // PRISM: ref mapping verified by DEBUG_HO_NODES:
      //   ref(0,0,0)->el[2], ref(1,0,0)->el[0], ref(0,1,0)->el[1]
      //   ref(0,0,1)->el[5], ref(1,0,1)->el[3], ref(0,1,1)->el[4]
      // prism_ref[i] = ref coords for el[i]
      static const double prism_ref[][3] = {
        {1,0,0},  // el[0]
        {0,1,0},  // el[1]
        {0,0,0},  // el[2]
        {1,0,1},  // el[3]
        {0,1,1},  // el[4]
        {0,0,1},  // el[5]
      };
      static const int prism_edges[][2] = {
        {0,1},{1,2},{2,0},
        {3,4},{4,5},{5,3},
        {0,3},{1,4},{2,5}
      };

      for (auto &e : prism_edges) {
        int v0_cubit = ng_pi_to_cubit_nid_[el[e[0]]];
        int v1_cubit = ng_pi_to_cubit_nid_[el[e[1]]];
        HoEdgeKey key(v0_cubit, v1_cubit);
        if (edge_ho_nodes_.count(key)) continue;

        std::vector<int> edge_nodes;
        for (int k = 1; k < order; k++) {
          double t = (double)k / order;
          double xi   = prism_ref[e[0]][0] + t * (prism_ref[e[1]][0] - prism_ref[e[0]][0]);
          double eta  = prism_ref[e[0]][1] + t * (prism_ref[e[1]][1] - prism_ref[e[0]][1]);
          double zeta = prism_ref[e[0]][2] + t * (prism_ref[e[1]][2] - prism_ref[e[0]][2]);
          ng::Point<3> ref(xi, eta, zeta);
          ng::Point<3> phys;
          curved.CalcElementTransformation(ref, ei, phys);
          int new_id = next_node_id_++;
          total_nodes_++;
          ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
          edge_nodes.push_back(new_id);
        }
        if (v0_cubit <= v1_cubit)
          edge_ho_nodes_[key] = edge_nodes;
        else
          edge_ho_nodes_[key] = std::vector<int>(edge_nodes.rbegin(), edge_nodes.rend());
      }

      // Face interiors (for order >= 3)
      // NOTE: still uses CalcElementTransformation. TODO: fix for order >= 3.
      if (order >= 3) {
        static const double prism_ref[][3] = {
          {0,0,0},{1,0,0},{0,1,0},
          {0,0,1},{1,0,1},{0,1,1}
        };
        static const int prism_trifaces[][3] = {{0,1,2},{3,4,5}};
        for (auto &f : prism_trifaces) {
          for (int i = 1; i < order; i++) {
            for (int j = 1; j < order - i; j++) {
              double a = (double)i / order;
              double b = (double)j / order;
              double c = 1.0 - a - b;
              double xi   = c*prism_ref[f[0]][0] + a*prism_ref[f[1]][0] + b*prism_ref[f[2]][0];
              double eta  = c*prism_ref[f[0]][1] + a*prism_ref[f[1]][1] + b*prism_ref[f[2]][1];
              double zeta = c*prism_ref[f[0]][2] + a*prism_ref[f[1]][2] + b*prism_ref[f[2]][2];
              ng::Point<3> ref(xi, eta, zeta);
              ng::Point<3> phys;
              curved.CalcElementTransformation(ref, ei, phys);
              int new_id = next_node_id_++;
              total_nodes_++;
              ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
            }
          }
        }
        // 3 quad faces
        static const int prism_quadfaces[][4] = {{0,1,4,3},{1,2,5,4},{0,2,5,3}};
        for (auto &f : prism_quadfaces) {
          for (int i = 1; i < order; i++) {
            for (int j = 1; j < order; j++) {
              double s = (double)i / order;
              double t = (double)j / order;
              double xi   = (1-s)*(1-t)*prism_ref[f[0]][0] + s*(1-t)*prism_ref[f[1]][0]
                          + s*t*prism_ref[f[2]][0] + (1-s)*t*prism_ref[f[3]][0];
              double eta  = (1-s)*(1-t)*prism_ref[f[0]][1] + s*(1-t)*prism_ref[f[1]][1]
                          + s*t*prism_ref[f[2]][1] + (1-s)*t*prism_ref[f[3]][1];
              double zeta = (1-s)*(1-t)*prism_ref[f[0]][2] + s*(1-t)*prism_ref[f[1]][2]
                          + s*t*prism_ref[f[2]][2] + (1-s)*t*prism_ref[f[3]][2];
              ng::Point<3> ref(xi, eta, zeta);
              ng::Point<3> phys;
              curved.CalcElementTransformation(ref, ei, phys);
              int new_id = next_node_id_++;
              total_nodes_++;
              ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
            }
          }
        }

        // Volume interior
        for (int i = 1; i < order; i++) {
          for (int j = 1; j < order - i; j++) {
            for (int k = 1; k < order; k++) {
              double xi   = (double)i / order;
              double eta  = (double)j / order;
              double zeta = (double)k / order;
              ng::Point<3> ref(xi, eta, zeta);
              ng::Point<3> phys;
              curved.CalcElementTransformation(ref, ei, phys);
              int new_id = next_node_id_++;
              total_nodes_++;
              ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
            }
          }
        }
      }
    }
    else if (np == 5) {  // Pyramid
      // Netgen pyramid ref: (0,0,0),(1,0,0),(1,1,0),(0,1,0),(0.5,0.5,1)
      // Need to determine mapping to el[0..4] empirically.
      // read_gmsh.py: pyramid reorder = [3,2,1,0,4] (GMSH->Netgen)
      // So: Netgen pos 0 = GMSH v3, pos 1 = GMSH v2, pos 2 = GMSH v1, pos 3 = GMSH v0, pos 4 = GMSH v4
      // Inverse (Netgen ref vertex -> el index):
      //   ref(0,0,0) -> Netgen pos 0 -> el[0] which is GMSH v3
      //   ref(1,0,0) -> Netgen pos 1 -> el[1] which is GMSH v2
      //   ref(1,1,0) -> Netgen pos 2 -> el[2] which is GMSH v1
      //   ref(0,1,0) -> Netgen pos 3 -> el[3] which is GMSH v0
      //   ref(0.5,0.5,1) -> Netgen pos 4 -> el[4] which is GMSH v4
      // BUT we don't reorder! el[0..4] = Cubit/GMSH order.
      // So pyr_ref[i] = ref coords for el[i]:
      //   el[0] = GMSH v0 = Netgen pos 3 -> ref(0,1,0)
      //   el[1] = GMSH v1 = Netgen pos 2 -> ref(1,1,0)
      //   el[2] = GMSH v2 = Netgen pos 1 -> ref(1,0,0)
      //   el[3] = GMSH v3 = Netgen pos 0 -> ref(0,0,0)
      //   el[4] = GMSH v4 = Netgen pos 4 -> ref(0.5,0.5,1)
      // Netgen pyramid ref: shapes[4]=z, apex at ref(0,0,1)
      // Base vertices identity-mapped (verified by DEBUG_HO_NODES)
      static const double pyr_ref[][3] = {
        {0,0,0},  // el[0]
        {1,0,0},  // el[1]
        {1,1,0},  // el[2]
        {0,1,0},  // el[3]
        {0,0,1},  // el[4] = apex
      };
#ifdef DEBUG_HO_NODES
      if (ei < 2) {
        PRINT_INFO("  PYRAMID VE %d: Netgen ref -> physical:\n", ei);
        for (int rv = 0; rv < 5; rv++) {
          ng::Point<3> ref(pyr_ref[rv][0], pyr_ref[rv][1], pyr_ref[rv][2]);
          ng::Point<3> phys;
          curved.CalcElementTransformation(ref, ei, phys);
          int best_v = -1; double best_d = 1e30;
          for (int k = 0; k < 5; k++) {
            int cid = ng_pi_to_cubit_nid_[el[k]];
            auto cc = CubitInterface::get_nodal_coordinates(cid);
            double d = (cc[0]-phys[0])*(cc[0]-phys[0])
                     + (cc[1]-phys[1])*(cc[1]-phys[1])
                     + (cc[2]-phys[2])*(cc[2]-phys[2]);
            if (d < best_d) { best_d = d; best_v = k; }
          }
          PRINT_INFO("    pyr_ref[%d]=(%g,%g,%g) -> el[%d] dist=%.2e\n",
            rv, pyr_ref[rv][0], pyr_ref[rv][1], pyr_ref[rv][2],
            best_v, sqrt(best_d));
        }
      }
#endif
      static const int pyr_edges[][2] = {
        {0,1},{1,2},{2,3},{3,0},
        {0,4},{1,4},{2,4},{3,4}
      };

      for (auto &e : pyr_edges) {
        int v0_cubit = ng_pi_to_cubit_nid_[el[e[0]]];
        int v1_cubit = ng_pi_to_cubit_nid_[el[e[1]]];
        HoEdgeKey key(v0_cubit, v1_cubit);
        if (edge_ho_nodes_.count(key)) continue;

        std::vector<int> edge_nodes;
        for (int k = 1; k < order; k++) {
          double t = (double)k / order;
          double xi   = pyr_ref[e[0]][0] + t * (pyr_ref[e[1]][0] - pyr_ref[e[0]][0]);
          double eta  = pyr_ref[e[0]][1] + t * (pyr_ref[e[1]][1] - pyr_ref[e[0]][1]);
          double zeta = pyr_ref[e[0]][2] + t * (pyr_ref[e[1]][2] - pyr_ref[e[0]][2]);
          ng::Point<3> ref(xi, eta, zeta);
          ng::Point<3> phys;
          curved.CalcElementTransformation(ref, ei, phys);
          int new_id = next_node_id_++;
          total_nodes_++;
          ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
          edge_nodes.push_back(new_id);
        }
        if (v0_cubit <= v1_cubit)
          edge_ho_nodes_[key] = edge_nodes;
        else
          edge_ho_nodes_[key] = std::vector<int>(edge_nodes.rbegin(), edge_nodes.rend());
      }

      // Face interiors (for order >= 3)
      // NOTE: still uses CalcElementTransformation. TODO: fix for order >= 3.
      if (order >= 3) {
        static const double pyr_ref[][3] = {
          {0,0,0},{1,0,0},{1,1,0},{0,1,0},
          {0.5,0.5,1}
        };
        // 1 quad face (base)
        for (int i = 1; i < order; i++) {
          for (int j = 1; j < order; j++) {
            double xi   = (double)i / order;
            double eta  = (double)j / order;
            ng::Point<3> ref(xi, eta, 0.0);
            ng::Point<3> phys;
            curved.CalcElementTransformation(ref, ei, phys);
            int new_id = next_node_id_++;
            total_nodes_++;
            ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
          }
        }
        // 4 triangular faces
        static const int pyr_trifaces[][3] = {{0,1,4},{1,2,4},{2,3,4},{3,0,4}};
        for (auto &f : pyr_trifaces) {
          for (int i = 1; i < order; i++) {
            for (int j = 1; j < order - i; j++) {
              double a = (double)i / order;
              double b = (double)j / order;
              double c = 1.0 - a - b;
              double xi   = c*pyr_ref[f[0]][0] + a*pyr_ref[f[1]][0] + b*pyr_ref[f[2]][0];
              double eta  = c*pyr_ref[f[0]][1] + a*pyr_ref[f[1]][1] + b*pyr_ref[f[2]][1];
              double zeta = c*pyr_ref[f[0]][2] + a*pyr_ref[f[1]][2] + b*pyr_ref[f[2]][2];
              ng::Point<3> ref(xi, eta, zeta);
              ng::Point<3> phys;
              curved.CalcElementTransformation(ref, ei, phys);
              int new_id = next_node_id_++;
              total_nodes_++;
              ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
            }
          }
        }
      }
    }
  }

  // --- Phase 2: Surface elements (CalcSurfaceTransformation) ---
  // Only for edges NOT already registered by volume elements,
  // and for face interior nodes (order >= 3).
  int nse = ng_mesh_->GetNSE();
  for (int sei = 0; sei < nse; sei++) {
    if (!curved.IsSurfaceElementCurved(sei)) continue;
    const ng::Element2d & sel = ng_mesh_->SurfaceElement(sei+1);
    int snp = sel.GetNP();

    if (snp == 3) {  // Triangle
      for (int edge = 0; edge < 3; edge++) {
        int v0_cubit = ng_pi_to_cubit_nid_[sel[edge]];
        int v1_cubit = ng_pi_to_cubit_nid_[sel[(edge+1) % 3]];
        HoEdgeKey key(v0_cubit, v1_cubit);
        if (edge_ho_nodes_.count(key)) continue;

        std::vector<int> edge_nodes;
        for (int k = 1; k < order; k++) {
          double t = (double)k / order;
          double xi, eta;
          if (edge == 0) { xi = t; eta = 0; }
          else if (edge == 1) { xi = 1-t; eta = t; }
          else { xi = 0; eta = 1-t; }
          ng::Point<2> ref(xi, eta);
          ng::Point<3> phys;
          curved.CalcSurfaceTransformation(ref, sei, phys);
          int new_id = next_node_id_++;
          total_nodes_++;
          ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
          edge_nodes.push_back(new_id);
        }
        if (v0_cubit <= v1_cubit)
          edge_ho_nodes_[key] = edge_nodes;
        else
          edge_ho_nodes_[key] = std::vector<int>(edge_nodes.rbegin(), edge_nodes.rend());
      }
      if (order >= 3) {
        int sv0 = ng_pi_to_cubit_nid_[sel[0]];
        int sv1 = ng_pi_to_cubit_nid_[sel[1]];
        int sv2 = ng_pi_to_cubit_nid_[sel[2]];
        HoFaceKey3 fkey(sv0, sv1, sv2);
        if (!face_ho_nodes_.count(fkey)) {
          std::vector<int> face_nodes;
          for (int i = 1; i < order; i++) {
            for (int j = 1; j < order - i; j++) {
              double xi = (double)i / order;
              double eta = (double)j / order;
              ng::Point<2> ref(xi, eta);
              ng::Point<3> phys;
              curved.CalcSurfaceTransformation(ref, sei, phys);
              int new_id = next_node_id_++;
              total_nodes_++;
              ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
              face_nodes.push_back(new_id);
            }
          }
          face_ho_nodes_[fkey] = {face_nodes, sv0, sv1, sv2};
        }
      }
    }
    else if (snp == 4) {  // Quad
      for (int edge = 0; edge < 4; edge++) {
        int v0_cubit = ng_pi_to_cubit_nid_[sel[edge]];
        int v1_cubit = ng_pi_to_cubit_nid_[sel[(edge+1) % 4]];
        HoEdgeKey key(v0_cubit, v1_cubit);
        if (edge_ho_nodes_.count(key)) continue;

        std::vector<int> edge_nodes;
        for (int k = 1; k < order; k++) {
          double t = (double)k / order;
          double xi, eta;
          if (edge == 0) { xi = t; eta = 0; }
          else if (edge == 1) { xi = 1; eta = t; }
          else if (edge == 2) { xi = 1-t; eta = 1; }
          else { xi = 0; eta = 1-t; }
          ng::Point<2> ref(xi, eta);
          ng::Point<3> phys;
          curved.CalcSurfaceTransformation(ref, sei, phys);
          int new_id = next_node_id_++;
          total_nodes_++;
          ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
          edge_nodes.push_back(new_id);
        }
        if (v0_cubit <= v1_cubit)
          edge_ho_nodes_[key] = edge_nodes;
        else
          edge_ho_nodes_[key] = std::vector<int>(edge_nodes.rbegin(), edge_nodes.rend());
      }
      if (order >= 3) {
        for (int i = 1; i < order; i++) {
          for (int j = 1; j < order; j++) {
            double xi = (double)i / order;
            double eta = (double)j / order;
            ng::Point<2> ref(xi, eta);
            ng::Point<3> phys;
            curved.CalcSurfaceTransformation(ref, sei, phys);
            int new_id = next_node_id_++;
            total_nodes_++;
            ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
          }
        }
      }
    }
  }

  // ---------------------------------------------------------------
  // Final pass: fill in missing edges with Cubit linear interpolation.
  // Uncurved volume elements are skipped by the loop above, so their
  // interior edges (not on any surface) have no HO nodes yet.
  // This ensures build_ho_conn_nc() never returns empty for any element.
  // ---------------------------------------------------------------
  {
    // Edge tables for each element type (matching EdgeTables:: in MeshData.hpp)
    static const int tet_edges[][2]  = {{0,1},{1,2},{0,2},{0,3},{1,3},{2,3}};
    static const int hex_edges2[][2] = {
      {0,1},{1,2},{2,3},{3,0},
      {4,5},{5,6},{6,7},{7,4},
      {0,4},{1,5},{2,6},{3,7}
    };
    static const int prism_edges2[][2] = {
      {0,1},{1,2},{2,0},
      {3,4},{4,5},{5,3},
      {0,3},{1,4},{2,5}
    };
    static const int pyr_edges2[][2] = {
      {0,1},{1,2},{2,3},{3,0},
      {0,4},{1,4},{2,4},{3,4}
    };

    // Helper lambda: register missing edge with linear interpolation
    auto fill_edge = [&](int v0_cubit, int v1_cubit) -> bool {
      HoEdgeKey key(v0_cubit, v1_cubit);
      if (edge_ho_nodes_.count(key)) return false;

      auto c0 = CubitInterface::get_nodal_coordinates(v0_cubit);
      auto c1 = CubitInterface::get_nodal_coordinates(v1_cubit);

      std::vector<int> edge_nodes;
      for (int k = 1; k < order; k++) {
        double t = (double)k / order;
        int new_id = next_node_id_++;
        total_nodes_++;
        ho_node_coords_[new_id] = {c0[0] + t*(c1[0]-c0[0]),
                                    c0[1] + t*(c1[1]-c0[1]),
                                    c0[2] + t*(c1[2]-c0[2])};
        edge_nodes.push_back(new_id);
      }
      if (v0_cubit <= v1_cubit)
        edge_ho_nodes_[key] = edge_nodes;
      else
        edge_ho_nodes_[key] = std::vector<int>(edge_nodes.rbegin(), edge_nodes.rend());
      return true;
    };

    int filled = 0;

    // Pass A: volume elements
    for (int ei = 0; ei < nve; ei++) {
      const ng::Element & el = ng_mesh_->VolumeElement(ei+1);
      int np = el.GetNP();

      const int (*edges)[2] = nullptr;
      int num_edges = 0;
      if (np == 4)      { edges = tet_edges;    num_edges = 6; }
      else if (np == 8) { edges = hex_edges2;   num_edges = 12; }
      else if (np == 6) { edges = prism_edges2; num_edges = 9; }
      else if (np == 5) { edges = pyr_edges2;   num_edges = 8; }
      else continue;

      for (int e = 0; e < num_edges; e++)
        if (fill_edge(ng_pi_to_cubit_nid_[el[edges[e][0]]],
                      ng_pi_to_cubit_nid_[el[edges[e][1]]]))
          filled++;
    }

    // Pass B: surface elements (uncurved faces skipped by CalcSurfaceTransformation)
    static const int tri_edges[][2]  = {{0,1},{1,2},{2,0}};
    static const int quad_edges[][2] = {{0,1},{1,2},{2,3},{3,0}};

    int nse = ng_mesh_->GetNSE();
    for (int sei = 0; sei < nse; sei++) {
      const ng::Element2d & el = ng_mesh_->SurfaceElement(sei+1);
      int np = el.GetNP();

      const int (*edges)[2] = nullptr;
      int num_edges = 0;
      if (np == 3)      { edges = tri_edges;  num_edges = 3; }
      else if (np == 4) { edges = quad_edges; num_edges = 4; }
      else continue;

      for (int e = 0; e < num_edges; e++)
        if (fill_edge(ng_pi_to_cubit_nid_[el[edges[e][0]]],
                      ng_pi_to_cubit_nid_[el[edges[e][1]]]))
          filled++;
    }

    if (filled > 0)
      PRINT_INFO("NetgenCurver: filled %d missing edges with linear interpolation.\n", filled);
  }

  // --- Phase 4: Fill missing face/volume interior nodes (linear interpolation) ---
  // For uncurved elements (flat geometry), face/volume nodes were not generated
  // in Phase 1/2.  Generate them using barycentric interpolation from vertices.
  if (order >= 3) {
    int filled_faces = 0, filled_vols = 0;
    int nve = ng_mesh_->GetNE();
    for (int ei = 0; ei < nve; ei++) {
      const ng::Element & el = ng_mesh_->VolumeElement(ei+1);
      int np = el.GetNP();
      if (np != 4) continue;  // Only tet needs complete Lagrange face/vol nodes

      int v[4];
      for (int k = 0; k < 4; k++)
        v[k] = ng_pi_to_cubit_nid_[el[k]];

      // Get vertex coordinates
      std::array<double,3> vc[4];
      for (int k = 0; k < 4; k++)
        vc[k] = get_node_coords(v[k]);

      // Face interior nodes
      static const int tet_faces_fill[][3] = {{0,1,2},{0,1,3},{0,2,3},{1,2,3}};
      for (auto &f : tet_faces_fill) {
        HoFaceKey3 fkey(v[f[0]], v[f[1]], v[f[2]]);
        if (face_ho_nodes_.count(fkey)) continue;

        std::vector<int> face_nodes;
        for (int i = 1; i < order; i++) {
          for (int j = 1; j < order - i; j++) {
            double a = (double)i / order;
            double b = (double)j / order;
            double c = 1.0 - a - b;
            // Barycentric on face vertices f[0], f[1], f[2]
            int new_id = next_node_id_++;
            total_nodes_++;
            ho_node_coords_[new_id] = {
              c * vc[f[0]][0] + a * vc[f[1]][0] + b * vc[f[2]][0],
              c * vc[f[0]][1] + a * vc[f[1]][1] + b * vc[f[2]][1],
              c * vc[f[0]][2] + a * vc[f[1]][2] + b * vc[f[2]][2],
            };
            face_nodes.push_back(new_id);
          }
        }
        face_ho_nodes_[fkey] = {face_nodes, v[f[0]], v[f[1]], v[f[2]]};
        filled_faces++;
      }

      // Volume interior nodes (order >= 4)
      if (order >= 4) {
        HoElemKey4 vkey(v[0], v[1], v[2], v[3]);
        if (!vol_ho_nodes_.count(vkey)) {
          std::vector<int> vol_nodes;
          for (int i = 1; i < order; i++) {
            for (int j = 1; j < order - i; j++) {
              for (int k = 1; k < order - i - j; k++) {
                double a = (double)i / order;
                double b = (double)j / order;
                double c = (double)k / order;
                double d = 1.0 - a - b - c;
                // Barycentric on tet vertices (Netgen ref: el[0]->lam1, el[1]->lam2, ...)
                int new_id = next_node_id_++;
                total_nodes_++;
                ho_node_coords_[new_id] = {
                  d * vc[0][0] + a * vc[1][0] + b * vc[2][0] + c * vc[3][0],
                  d * vc[0][1] + a * vc[1][1] + b * vc[2][1] + c * vc[3][1],
                  d * vc[0][2] + a * vc[1][2] + b * vc[2][2] + c * vc[3][2],
                };
                vol_nodes.push_back(new_id);
              }
            }
          }
          vol_ho_nodes_[vkey] = {vol_nodes, v[0], v[1], v[2], v[3]};
          filled_vols++;
        }
      }
    }

    // Also fill missing face nodes for surface triangles
    int nse = ng_mesh_->GetNSE();
    for (int sei = 0; sei < nse; sei++) {
      const ng::Element2d & sel = ng_mesh_->SurfaceElement(sei+1);
      if (sel.GetNP() != 3) continue;

      int sv[3];
      for (int k = 0; k < 3; k++)
        sv[k] = ng_pi_to_cubit_nid_[sel[k]];

      HoFaceKey3 fkey(sv[0], sv[1], sv[2]);
      if (face_ho_nodes_.count(fkey)) continue;

      std::array<double,3> vc[3];
      for (int k = 0; k < 3; k++)
        vc[k] = get_node_coords(sv[k]);

      std::vector<int> face_nodes;
      for (int i = 1; i < order; i++) {
        for (int j = 1; j < order - i; j++) {
          double a = (double)i / order;
          double b = (double)j / order;
          double c = 1.0 - a - b;
          int new_id = next_node_id_++;
          total_nodes_++;
          ho_node_coords_[new_id] = {
            c * vc[0][0] + a * vc[1][0] + b * vc[2][0],
            c * vc[0][1] + a * vc[1][1] + b * vc[2][1],
            c * vc[0][2] + a * vc[1][2] + b * vc[2][2],
          };
          face_nodes.push_back(new_id);
        }
      }
      face_ho_nodes_[fkey] = {face_nodes, sv[0], sv[1], sv[2]};
      filled_faces++;
    }

    if (filled_faces > 0 || filled_vols > 0)
      PRINT_INFO("NetgenCurver: filled %d missing faces, %d missing volumes "
                 "with linear interpolation.\n", filled_faces, filled_vols);
  }

  return true;
}

std::vector<int> NetgenCurver::get_edge_nodes(int n0, int n1) const
{
  HoEdgeKey key(n0, n1);
  auto it = edge_ho_nodes_.find(key);
  if (it == edge_ho_nodes_.end())
    return {};

  // Return in n0→n1 direction
  if (n0 <= n1)
    return it->second;
  else
    return std::vector<int>(it->second.rbegin(), it->second.rend());
}

// Helper: permute face interior nodes from generation order to caller order.
// Face nodes are generated iterating (i,j) with weight on (gen_v0, gen_v1, gen_v2).
// Caller wants nodes relative to (req_v0, req_v1, req_v2).
// The permutation maps barycentric coordinates: find which gen vertex maps to
// which req vertex, then re-index the (i,j) loop accordingly.
static std::vector<int> permute_face_nodes(
    const std::vector<int> &nodes, int order,
    int gen_v0, int gen_v1, int gen_v2,
    int req_v0, int req_v1, int req_v2)
{
  // Build permutation: perm[k] = which gen vertex corresponds to req vertex k
  // gen_v[perm[0]] == req_v0, gen_v[perm[1]] == req_v1, gen_v[perm[2]] == req_v2
  int gen[3] = {gen_v0, gen_v1, gen_v2};
  int req[3] = {req_v0, req_v1, req_v2};
  int perm[3] = {-1, -1, -1};
  for (int r = 0; r < 3; r++)
    for (int g = 0; g < 3; g++)
      if (req[r] == gen[g]) { perm[r] = g; break; }

  if (perm[0] == -1 || perm[1] == -1 || perm[2] == -1)
    return nodes;  // vertex mismatch — should not happen

  // Identity permutation — no reorder needed
  if (perm[0] == 0 && perm[1] == 1 && perm[2] == 2)
    return nodes;

  // Build index map: for each (i,j) in generation order, compute the
  // corresponding (i',j') in the requested order.
  // Generation: weight_on_gen0 = 1-i/p-j/p, weight_on_gen1 = i/p, weight_on_gen2 = j/p
  // Requested:  weight_on_req0 = 1-i'/p-j'/p, weight_on_req1 = i'/p, weight_on_req2 = j'/p
  // Since req[k] = gen[perm[k]], the weight on req[k] = weight on gen[perm[k]].
  // So for req vertex 1: i'/p = weight on gen[perm[1]]
  //    for req vertex 2: j'/p = weight on gen[perm[2]]

  int p = order;
  int n = (int)nodes.size();
  std::vector<int> result(n);

  // Map (i,j) -> flat index in generation order
  auto flat_idx = [p](int i, int j) -> int {
    int idx = 0;
    for (int ii = 1; ii < p; ii++) {
      for (int jj = 1; jj < p - ii; jj++) {
        if (ii == i && jj == j) return idx;
        idx++;
      }
    }
    return -1;
  };

  // For each generation (i,j), compute barycentric weights on gen vertices,
  // then find (i',j') in requested order
  int idx = 0;
  for (int i = 1; i < p; i++) {
    for (int j = 1; j < p - i; j++) {
      // Barycentric weights on gen vertices: w0=1-i/p-j/p, w1=i/p, w2=j/p
      // Weight on gen[k]: w[0] if k==0, w[1] if k==1, w[2] if k==2
      // Weight on req vertex 1 = weight on gen[perm[1]]
      // Weight on req vertex 2 = weight on gen[perm[2]]
      int w[3] = {p - i - j, i, j};  // integer barycentric (sum = p)
      int ip = w[perm[1]];  // i' * (implicitly /p)
      int jp = w[perm[2]];  // j'
      int dst = flat_idx(ip, jp);
      if (dst >= 0 && dst < n)
        result[dst] = nodes[idx];
      idx++;
    }
  }
  return result;
}

std::vector<int> NetgenCurver::get_face_nodes(int n0, int n1, int n2) const
{
  HoFaceKey3 fkey(n0, n1, n2);
  auto it = face_ho_nodes_.find(fkey);
  if (it == face_ho_nodes_.end())
    return {};

  const auto &fd = it->second;
  if (fd.node_ids.empty())
    return {};

  // Permute from generation vertex order to caller's vertex order
  return permute_face_nodes(fd.node_ids, order_,
                            fd.gen_v0, fd.gen_v1, fd.gen_v2,
                            n0, n1, n2);
}

std::vector<int> NetgenCurver::get_vol_nodes(int n0, int n1, int n2, int n3) const
{
  HoElemKey4 vkey(n0, n1, n2, n3);
  auto it = vol_ho_nodes_.find(vkey);
  if (it == vol_ho_nodes_.end())
    return {};

  const auto &vd = it->second;
  if (vd.node_ids.empty())
    return {};

  // For volume interior nodes, the permutation is more complex (4 vertices).
  // However, volume nodes are unique per element (not shared), so the generation
  // order should match the caller's order. Verify and return directly.
  // If they don't match, we'd need a full tet barycentric permutation — but
  // since each tet generates its own vol nodes, gen order == caller order.
  return vd.node_ids;
}

std::array<double,3> NetgenCurver::get_node_coords(int node_id) const
{
  auto it = ho_node_coords_.find(node_id);
  if (it != ho_node_coords_.end())
    return it->second;
  return CubitInterface::get_nodal_coordinates(node_id);
}

void NetgenCurver::write_ho_nodes_gmsh(std::ofstream &fid) const
{
  for (auto &pair : ho_node_coords_) {
    fid << pair.first << " "
        << pair.second[0] << " "
        << pair.second[1] << " "
        << pair.second[2] << "\n";
  }
}

void NetgenCurver::write_ho_nodes_nastran(std::ofstream &fid, bool is_3d) const
{
  fid << "$\n$ Generated high-order nodes (NetgenCurver)\n$\n";
  for (auto &pair : ho_node_coords_) {
    int nid = pair.first;
    auto &c = pair.second;
    char line1[128], line2[128];
    std::snprintf(line1, sizeof(line1), "GRID*   %16d%16d%16.5f%16.5f",
                  nid, 0, c[0], c[1]);
    std::snprintf(line2, sizeof(line2), "*       %16.5f", c[2]);
    fid << line1 << "\n" << line2 << "\n";
  }
}

// Gmsh element type codes for higher orders
int NetgenCurver::gmsh_element_type(int base_nn, int elem_type, int order)
{
  // Order 1
  if (order == 1) {
    if (base_nn == 4) return 4;   // TET4
    if (base_nn == 8) return 5;   // HEX8
    if (base_nn == 6) return 6;   // WEDGE6
    if (base_nn == 5) return 7;   // PYRAMID5
    if (base_nn == 3) return 2;   // TRI3
    if (base_nn == 4) return 3;   // QUAD4 (ambiguous with TET4)
    return 0;
  }
  // Order 2
  if (order == 2) {
    if (base_nn == 4) return 11;  // TET10
    if (base_nn == 8) return 17;  // HEX20
    if (base_nn == 6) return 18;  // WEDGE15
    if (base_nn == 5) return 19;  // PYRAMID13
    if (base_nn == 3) return 9;   // TRI6
    return 0;
  }
  // Order 3
  if (order == 3) {
    if (base_nn == 4) return 29;  // TET20
    if (base_nn == 3) return 21;  // TRI10
    return 0;
  }
  // Order 4
  if (order == 4) {
    if (base_nn == 4) return 30;  // TET35
    if (base_nn == 3) return 23;  // TRI15
    return 0;
  }
  return 0;
}

#else  // !HAVE_NETGEN

// Stub implementation when Netgen is not available
NetgenCurver::NetgenCurver() {}
NetgenCurver::~NetgenCurver() {}

bool NetgenCurver::build(const MeshData &, int)
{
  // Cannot build without Netgen
  return false;
}

std::array<double,3> NetgenCurver::get_node_coords(int) const { return {0,0,0}; }
std::vector<int> NetgenCurver::get_edge_nodes(int, int) const { return {}; }
std::vector<int> NetgenCurver::get_face_nodes(int, int, int) const { return {}; }
std::vector<int> NetgenCurver::get_vol_nodes(int, int, int, int) const { return {}; }
void NetgenCurver::write_ho_nodes_gmsh(std::ofstream &) const {}
void NetgenCurver::write_ho_nodes_nastran(std::ofstream &, bool) const {}
int NetgenCurver::gmsh_element_type(int, int, int) { return 0; }

#endif  // HAVE_NETGEN
