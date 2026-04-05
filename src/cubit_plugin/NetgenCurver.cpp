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
#include "GeometryQueryTool.hpp"
#include "CubitVector.hpp"

// Use explicit namespace prefixes to disambiguate
// Cubit ElementType enums: ::TETRA4, ::TETRA, ::HEX8, ::HEX, etc.
// Netgen ELEMENT_TYPE enums: netgen::TET, netgen::HEX, netgen::PRISM, etc.

#include <algorithm>
#include <cmath>
#include <fstream>
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
    // Cache RefFace pointers per surface
    std::unordered_map<int, RefFace*> rf_cache;
    for (int sid : surface_ids)
      rf_cache[sid] = GeometryQueryTool::instance()->get_ref_face(sid);

    int nse_added = 0;
    for (int sid : surface_ids) {
      int ng_fd = cubit_sid_to_ng_fd_[sid];
      RefFace* rf = rf_cache[sid];

      // Triangles on this surface
      std::vector<int> tri_ids = CubitInterface::parse_cubit_list("tri",
          "in surface " + std::to_string(sid));
      for (int tid : tri_ids) {
        std::vector<int> conn = CubitInterface::get_connectivity("tri", tid);
        if ((int)conn.size() < 3) continue;

        ng::Element2d el(ng::TRIG);
        el.SetIndex(ng_fd);
        for (int k = 0; k < 3; k++) {
          auto it = cubit_nid_to_ng_pi_.find(conn[k]);
          if (it == cubit_nid_to_ng_pi_.end()) goto skip_tri;
          el[k] = it->second;

          auto coords = CubitInterface::get_nodal_coordinates(conn[k]);
          CubitVector loc(coords[0], coords[1], coords[2]);
          double u = 0, v = 0;
          if (rf) rf->u_v_from_position(loc, u, v);
          ng::PointGeomInfo gi;
          gi.trignum = ng_fd;
          gi.u = u;
          gi.v = v;
          el.GeomInfo()[k] = gi;
        }
        ng_mesh_->AddSurfaceElement(el);
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
        for (int k = 0; k < 4; k++) {
          auto it = cubit_nid_to_ng_pi_.find(conn[k]);
          if (it == cubit_nid_to_ng_pi_.end()) goto skip_quad;
          el[k] = it->second;

          auto coords = CubitInterface::get_nodal_coordinates(conn[k]);
          CubitVector loc(coords[0], coords[1], coords[2]);
          double u = 0, v = 0;
          if (rf) rf->u_v_from_position(loc, u, v);
          ng::PointGeomInfo gi;
          gi.trignum = ng_fd;
          gi.u = u;
          gi.v = v;
          el.GeomInfo()[k] = gi;
        }
        ng_mesh_->AddSurfaceElement(el);
        nse_added++;
        skip_quad:;
      }
    }
    PRINT_INFO("  Phase1e: %d surface elements added\n", nse_added);
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
      RefEdge* re = GeometryQueryTool::instance()->get_ref_edge(cid);
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
      for (int k = 0; k < 8; k++)
        el[k] = cubit_nid_to_ng_pi_[elem.conn[k]];
      el.SetIndex(1);
      ng_mesh_->AddVolumeElement(el);
    }
    else if (nn == 6 && (et == ::WEDGE6 || et == ::WEDGE)) {
      ng::Element el(ng::PRISM);
      for (int k = 0; k < 6; k++)
        el[k] = cubit_nid_to_ng_pi_[elem.conn[k]];
      el.SetIndex(1);
      ng_mesh_->AddVolumeElement(el);
    }
    else if (nn == 5 && (et == ::PYRAMID5 || et == ::PYRAMID)) {
      ng::Element el(ng::PYRAMID);
      for (int k = 0; k < 5; k++)
        el[k] = cubit_nid_to_ng_pi_[elem.conn[k]];
      el.SetIndex(1);
      ng_mesh_->AddVolumeElement(el);
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
  int num_surfaces = (int)cubit_sid_to_ng_fd_.size();

  // Project callback: (surfnr, x,y,z, u_hint, v_hint) -> (xp,yp,zp, u,v)
  auto project = [this](int surfnr, double x, double y, double z,
                         double u_hint, double v_hint)
    -> std::tuple<double,double,double,double,double>
  {
    auto it = ng_fd_to_cubit_sid_.find(surfnr);
    if (it == ng_fd_to_cubit_sid_.end())
      return {x, y, z, u_hint, v_hint};

    int cubit_sid = it->second;
    RefFace* rf = GeometryQueryTool::instance()->get_ref_face(cubit_sid);
    if (!rf) return {x, y, z, u_hint, v_hint};

    ::Surface* surf = rf->get_surface_ptr();
    CubitVector loc(x, y, z);
    CubitVector closest;

    // Use trimmed closest point (global search, matches bridge.py)
    surf->closest_point_trimmed(loc, closest);
    double u, v;
    rf->u_v_from_position(closest, u, v);
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
    RefFace* rf = GeometryQueryTool::instance()->get_ref_face(cubit_sid);
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

  auto edge_project = [surfpair_to_curve](int surfnr1, int surfnr2,
                                           double x, double y, double z)
    -> std::tuple<double,double,double>
  {
    int64_t key = ((int64_t)surfnr1 << 32) | surfnr2;
    auto it = surfpair_to_curve.find(key);
    if (it != surfpair_to_curve.end()) {
      RefEdge* re = GeometryQueryTool::instance()->get_ref_edge(it->second);
      if (re) {
        CubitVector loc(x, y, z);
        CubitVector closest;
        re->get_curve_ptr()->closest_point_trimmed(loc, closest);
        return {closest.x(), closest.y(), closest.z()};
      }
    }
    return {x, y, z};  // no projection if curve not found
  };

  auto geom = std::make_shared<ng::CallbackGeometry>(
      project, normal, num_surfaces, nullptr, edge_project);
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

  // Extract curved node positions by evaluating element transformations
  // at reference-space coordinates of high-order nodes.
  //
  // For now, we extract all surface element curved data.
  // The mapping from netgen curved elements back to Cubit node IDs
  // is done via the element connectivity.

  int nse = ng_mesh_->GetNSE();
  for (int sei = 0; sei < nse; sei++) {
    if (!curved.IsSurfaceElementCurved(sei)) continue;

    const ng::Element2d & el = ng_mesh_->SurfaceElement(sei+1);
    int np = el.GetNP();

    // Evaluate at equidistant points for the given order
    // For a triangle of order p, there are (p+1)(p+2)/2 nodes
    // For a quad of order p, there are (p+1)^2 nodes
    // Edge nodes: evaluated at equidistant points along edges
    // Face nodes: evaluated at interior equidistant points

    // For each edge of the element, add (order-1) interior points
    // The first `np` points are the vertices (already in the mesh)

    if (np == 3) {  // Triangle
      // Edge points — register in edge_ho_nodes_ for connectivity lookup
      for (int edge = 0; edge < 3; edge++) {
        int v0_ng = el[edge];
        int v1_ng = el[(edge+1) % 3];
        int v0_cubit = ng_pi_to_cubit_nid_[v0_ng];
        int v1_cubit = ng_pi_to_cubit_nid_[v1_ng];

        HoEdgeKey key(v0_cubit, v1_cubit);
        if (edge_ho_nodes_.count(key)) continue;  // already processed

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
        // Store in min→max order
        if (v0_cubit <= v1_cubit)
          edge_ho_nodes_[key] = edge_nodes;
        else
          edge_ho_nodes_[key] = std::vector<int>(edge_nodes.rbegin(), edge_nodes.rend());
      }

      // Face interior points (for order >= 3)
      if (order >= 3) {
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
          }
        }
      }
    }
    else if (np == 4) {  // Quad
      // Edge points with edge_ho_nodes_ registration
      for (int edge = 0; edge < 4; edge++) {
        int v0_ng = el[edge];
        int v1_ng = el[(edge+1) % 4];
        int v0_cubit = ng_pi_to_cubit_nid_[v0_ng];
        int v1_cubit = ng_pi_to_cubit_nid_[v1_ng];

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

      // Face interior points
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

  // Also extract curved volume element data
  int nve = ng_mesh_->GetNE();
  for (int ei = 0; ei < nve; ei++) {
    if (!curved.IsElementCurved(ei)) continue;

    const ng::Element & el = ng_mesh_->VolumeElement(ei+1);
    int np = el.GetNP();

    if (np == 4) {  // Tet
      // Edge points (6 edges) — register in edge_ho_nodes_ for connectivity
      static const int tet_edges[][2] = {{0,1},{0,2},{0,3},{1,2},{1,3},{2,3}};
      for (auto &e : tet_edges) {
        int v0_cubit = ng_pi_to_cubit_nid_[el[e[0]]];
        int v1_cubit = ng_pi_to_cubit_nid_[el[e[1]]];
        HoEdgeKey key(v0_cubit, v1_cubit);
        if (edge_ho_nodes_.count(key)) continue;

        std::vector<int> edge_nodes;
        for (int k = 1; k < order; k++) {
          double t = (double)k / order;
          double lam[4] = {0,0,0,0};
          lam[e[0]] = 1.0 - t;
          lam[e[1]] = t;
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
          for (int i = 1; i < order; i++) {
            for (int j = 1; j < order - i; j++) {
              double lam[4] = {0,0,0,0};
              double a = (double)i / order;
              double b = (double)j / order;
              double c = 1.0 - a - b;
              lam[f[0]] = c;
              lam[f[1]] = a;
              lam[f[2]] = b;
              ng::Point<3> ref(lam[1], lam[2], lam[3]);
              ng::Point<3> phys;
              curved.CalcElementTransformation(ref, ei, phys);

              int new_id = next_node_id_++;
              total_nodes_++;
              ho_node_coords_[new_id] = {phys[0], phys[1], phys[2]};
            }
          }
        }
      }

      // Volume interior points (for order >= 4)
      if (order >= 4) {
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
            }
          }
        }
      }
    }
    else if (np == 8) {  // Hex
      // Hex reference coords: (xi, eta, zeta) in [0,1]^3
      // 12 edges
      static const int hex_edges[][2] = {
        {0,1},{1,2},{2,3},{3,0},  // bottom
        {4,5},{5,6},{6,7},{7,4},  // top
        {0,4},{1,5},{2,6},{3,7}   // vertical
      };
      // Hex vertex reference coords
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
      if (order >= 3) {
        // Face defined by 4 vertices in ref space
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
      // Prism reference: triangle (lam1, lam2) x [0,1] in zeta
      // lam0 = 1 - lam1 - lam2
      // 9 edges
      static const double prism_ref[][3] = {
        {0,0,0},{1,0,0},{0,1,0},  // bottom triangle
        {0,0,1},{1,0,1},{0,1,1}   // top triangle
      };
      static const int prism_edges[][2] = {
        {0,1},{1,2},{2,0},  // bottom
        {3,4},{4,5},{5,3},  // top
        {0,3},{1,4},{2,5}   // vertical
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
      if (order >= 3) {
        // 2 triangular faces (bottom, top)
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
      // Pyramid reference: quad base (xi, eta) at zeta=0, apex at zeta=1
      static const double pyr_ref[][3] = {
        {0,0,0},{1,0,0},{1,1,0},{0,1,0},  // base
        {0.5,0.5,1}                         // apex
      };
      static const int pyr_edges[][2] = {
        {0,1},{1,2},{2,3},{3,0},  // base
        {0,4},{1,4},{2,4},{3,4}   // lateral
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
      if (order >= 3) {
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
void NetgenCurver::write_ho_nodes_gmsh(std::ofstream &) const {}
void NetgenCurver::write_ho_nodes_nastran(std::ofstream &, bool) const {}
int NetgenCurver::gmsh_element_type(int, int, int) { return 0; }

#endif  // HAVE_NETGEN
