#include "MeshData.hpp"
#include "NetgenCurver.hpp"
#include "MeshExportInterface.hpp"
#include "CubitInterface.hpp"
#include "CubitMessage.hpp"

#include <algorithm>
#include <map>
#include <set>
#include <sstream>
#include <memory>
#include <exception>

// ================================================================
// extract() — main entry point
// ================================================================
bool MeshData::extract(int req_order) try
{
  CubitBaseInterface* raw = CubitInterface::get_interface("MeshExport");
  const auto release = [](CubitBaseInterface* p) { CubitInterface::release_interface(p); };
  std::unique_ptr<CubitBaseInterface, decltype(release)> interface_owner(raw, release);
  MeshExportInterface *iface = dynamic_cast<MeshExportInterface*>(raw);
  if (!iface) {
    PRINT_ERROR("Unable to get MeshExport interface.\n");
    return false;
  }

  iface->set_use_sequential_ids(false);
  iface->initialize_export();

  if (iface->get_num_nodes() == 0) {
    PRINT_WARNING("No mesh nodes found. Nothing to export.\n");
    return false;
  }

  order = req_order;

  // Extract linear mesh data first (single pass through iface)
  extract_nodes(iface);
  extract_elements(iface);
  if (!check_block_coverage()) {
    return false;
  }
  extract_sidesets(iface);
  extract_nodesets(iface);

  // Build high-order from extracted MeshData — no more iface iteration
  if (order >= 2) {
    auto nc = std::make_shared<NetgenCurver>();
    netgen_curver_ = nc; // retain failure diagnostics without retaining the interface
    PRINT_INFO("build_high_order: trying NetgenCurver (order %d)...\n", order);
    if (nc->build(*this, order)) {
      PRINT_INFO("build_high_order: NetgenCurver succeeded.\n");
      netgen_curver_ = nc;  // keep alive for the .vol writer
    } else {
      PRINT_ERROR("NetgenCurver failed for order %d. "
                  "No fallback — check mesh and geometry.\n", order);
      return false;
    }
  }

  return true;
}
catch (const std::exception & e) {
  PRINT_ERROR("Mesh extraction failed: %s\n", e.what());
  return false;
}
catch (...) {
  PRINT_ERROR("Mesh extraction failed: unknown exception\n");
  return false;
}

// ================================================================
// extract_nodes
// ================================================================
void MeshData::extract_nodes(MeshExportInterface *iface)
{
  const int buf = 1000;
  std::vector<double> xc(buf), yc(buf), zc(buf);
  std::vector<int> ids(buf);

  int start = 0, found = 0;
  int idx = 0;
  while ((found = iface->get_coords(start, buf, xc, yc, zc, ids)) > 0) {
    for (int i = 0; i < found; i++) {
      nodes.push_back({ids[i], xc[i], yc[i], zc[i]});
      node_id_to_index[ids[i]] = idx++;
      if (ids[i] > max_original_node_id)
        max_original_node_id = ids[i];
    }
    start += found;
  }
  num_original_nodes = (int)nodes.size();
}

// ================================================================
// extract_elements — from blocks
// ================================================================
void MeshData::extract_elements(MeshExportInterface *iface)
{
  // SNAPSHOT user-defined blocks BEFORE calling iface->get_block_list.
  //
  // BUG (2026-05-30, six-turn coil loft model):
  //   `export vol` errors with "ERROR: No block with ID 2 was
  //   found" AND leaves a phantom block 2 in the Cubit session.
  //
  // Cause: MeshExportInterface::get_block_list returns Cubit's internal
  // block iterator handles, INCLUDING a "default block" that holds the
  // elements not assigned to any user block.  Just calling
  // iface->id_from_handle(block) on that handle MATERIALISES the
  // phantom block in Cubit's database as a real block (id ~ N+1 where N
  // is the user's max block id).  So:
  //   1. ExportVolCommand calls `parse_cubit_list("block", "all")`
  //      AFTER extract_elements -- the phantom is already in the list.
  //   2. The skip check at ExportVolCommand:866
  //      (`if (cubit_block_set.find(bid) == end)`) does NOT fire.
  //   3. parse_cubit_list("volume", "in block N+1") errors because the
  //      phantom block has no volume membership.
  //   4. The phantom block survives the export and pollutes the model.
  //
  // Fix: take the user-block snapshot HERE, before the side-effecting
  // get_block_list/id_from_handle calls run, and only register block_ids
  // that were user-defined.  Phantom blocks never enter md.block_ids,
  // so downstream callers (the JSON writer, the FaceDescriptor mapping)
  // never try to look them up.
  std::vector<int> user_block_ids = CubitInterface::parse_cubit_list("block", "all");
  std::set<int> user_block_set(user_block_ids.begin(), user_block_ids.end());

  std::vector<BlockHandle> blocks;
  iface->get_block_list(blocks);

  std::set<int> seen_bids;
  // Deduplicate by sorted connectivity (ElementHandle may differ across blocks)
  std::set<std::vector<int>> seen_conn;

  const int ebuf = 100;
  std::vector<ElementType> etypes(ebuf);
  std::vector<ElementHandle> handles(ebuf);

  int dup_count = 0;
  for (size_t bi = 0; bi < blocks.size(); bi++) {
    BlockHandle block = blocks[bi];
    int bid = iface->id_from_handle(block);
    // Skip the MeshExportInterface "default block" that holds elements
    // not assigned to any user-defined Cubit block.  Without this guard,
    // (a) the phantom id leaks into block_ids, (b) downstream
    // parse_cubit_list("volume", "in block N") errors, (c) the phantom
    // block survives the export and pollutes the Cubit model.  See the
    // comment at the top of this function for the full diagnosis.
    if (user_block_set.find(bid) == user_block_set.end())
      continue;
    if (seen_bids.insert(bid).second)
      block_ids.push_back(bid);

    int estart = 0, count = 0;
    while ((count = iface->get_block_elements(estart, ebuf, block, etypes, handles)) > 0) {
      for (int i = 0; i < count; i++) {
        std::vector<int> conn(27);
        int nn = iface->get_connectivity(handles[i], conn);
        conn.resize(nn);

        // Skip duplicate elements (same connectivity in multiple blocks)
        std::vector<int> key(conn.begin(), conn.begin() + nn);
        std::sort(key.begin(), key.end());
        if (!seen_conn.insert(key).second) {
          dup_count++;
          continue;
        }

        MeshElement elem;
        elem.group_id = bid;
        elem.type = etypes[i];
        elem.nv = nn;
        elem.conn = std::move(conn);
        elements.push_back(std::move(elem));
      }
      estart += count;
    }
  }
  if (dup_count > 0)
    PRINT_WARNING("Skipped %d duplicate elements across blocks.\n", dup_count);
}

// ================================================================
// check_block_coverage — report volume elements outside every block
//
// Only block elements are exported (see extract_elements), so a hex or
// tet that was meshed but never added to a block is dropped without a
// trace.  Warn with per-type counts, and fail when the caller needs a
// volume mesh and no block holds a single volume element (a meshed
// `brick` without `block 1 add volume all` exported 0 elements).
// ================================================================
static const char *volume_kind(ElementType t)
{
  switch (t) {
    case TETRA: case TETRA4: case TETRA8: case TETRA10: case TETRA14: case TETRA15:
      return "tet";
    case HEX: case HEX8: case HEX9: case HEX20: case HEX21: case HEX26: case HEX27:
      return "hex";
    case WEDGE: case WEDGE6: case WEDGE15: case WEDGE16: case WEDGE20: case WEDGE21:
      return "wedge";
    case PYRAMID: case PYRAMID5: case PYRAMID8: case PYRAMID13: case PYRAMID18:
      return "pyramid";
    default:
      return nullptr;
  }
}

bool MeshData::check_block_coverage() const
{
  std::map<std::string, int> in_blocks;
  for (const auto &elem : elements)
    if (const char *kind = volume_kind(elem.type))
      in_blocks[kind]++;

  const std::pair<const char *, int> in_model[] = {
    {"tet", CubitInterface::get_tet_count()},
    {"hex", CubitInterface::get_hex_count()},
    {"wedge", CubitInterface::get_wedge_count()},
    {"pyramid", CubitInterface::get_pyramid_count()},
  };
  int n_volume = 0, n_model = 0;
  std::ostringstream missing;
  for (const auto &[kind, total] : in_model) {
    int exported = in_blocks[kind];
    n_volume += exported;
    n_model += total;
    if (total > exported)
      missing << (missing.tellp() > 0 ? ", " : "") << (total - exported)
              << " of " << total << " " << kind;
  }

  if (require_volume_elements && n_volume == 0) {
    PRINT_ERROR("No block contains volume elements; the model has %d meshed "
                "volume elements outside every block.  Add them to a block "
                "before exporting, e.g. 'block 1 add volume all'.\n", n_model);
    return false;
  }
  if (!missing.str().empty())
    PRINT_WARNING("Volume elements outside every block are not exported: %s.  "
                  "Add their volumes to a block to include them.\n",
                  missing.str().c_str());
  return true;
}

// ================================================================
// extract_sidesets — face elements from sideset surfaces
//
// ExodusII side table (get_sideset_sides + get_nodes_per_side) fails
// for wedge elements (returns nv=0). Instead, get the surfaces in
// each sideset and collect tri/face elements via parse_cubit_list,
// which works for all element types (tet, hex, wedge, pyramid).
// ================================================================
void MeshData::extract_sidesets(MeshExportInterface *iface)
{
  // SNAPSHOT user-defined sidesets BEFORE calling iface->get_sideset_list.
  // Same side-effect pattern as the block phantom-bug fix above:
  // iface->get_sideset_list returns a "default sideset" handle for faces
  // not assigned to any user-defined sideset; calling
  // iface->id_from_handle() on that handle MATERIALISES a phantom
  // sideset in Cubit's database (id ~ N+1).  Without this skip, the
  // phantom sideset leaks into downstream consumers (GMSH
  // $PhysicalNames as a junk entry, Netgen .vol as an extra bcname
  // with the fallback "Surface_K" name).  Pre-snapshot of
  // parse_cubit_list("sideset", "all") gives us the set of
  // user-defined sideset IDs; the phantom is skipped per the same
  // user_block_set pattern.  Bugs of this class were fixed for
  // blocks 2026-05-30 (six-turn coil loft); sidesets +
  // nodesets fixed 2026-05-31.
  std::vector<int> user_sideset_ids = CubitInterface::parse_cubit_list("sideset", "all");
  std::set<int> user_sideset_set(user_sideset_ids.begin(), user_sideset_ids.end());

  std::vector<SidesetHandle> ss_handles;
  iface->get_sideset_list(ss_handles);

  for (auto &ssh : ss_handles) {
    SidesetGroup sg;
    sg.id = iface->id_from_handle(ssh);
    // Skip phantom default-sideset materialised by iface->get_sideset_list.
    if (user_sideset_set.find(sg.id) == user_sideset_set.end())
      continue;
    sg.name = iface->get_sideset_name(ssh);

    // A sideset can own geometry surfaces, individual mesh faces, or both.
    // Sculpt imports its generated hex mesh as free elements, so its
    // gen_sidesets faces have no parent Cubit surface.  Keep one canonical
    // connectivity set and collect both representations; otherwise the
    // exported .vol silently loses the complete exterior skin.
    std::set<std::vector<int>> seen_faces;
    auto append_face = [&](ElementType type, std::vector<int> conn, int nv) {
      if ((int)conn.size() < nv) return;
      conn.resize(nv);
      std::vector<int> key(conn);
      std::sort(key.begin(), key.end());
      if (!seen_faces.insert(key).second) return;

      MeshElement face;
      face.group_id = sg.id;
      face.type = type;
      face.nv = nv;
      face.conn = std::move(conn);
      sg.faces.push_back(std::move(face));
    };

    // Get surfaces in this sideset
    std::vector<int> ss_surfs = CubitInterface::get_sideset_surfaces(sg.id);

    for (int sid : ss_surfs) {
      // Triangles on this surface (tet/wedge boundary faces)
      std::vector<int> tri_ids = CubitInterface::parse_cubit_list(
          "tri", "in surface " + std::to_string(sid));
      for (int tid : tri_ids) {
        append_face(TRI3, CubitInterface::get_connectivity("tri", tid), 3);
      }

      // Quads on this surface (hex/wedge boundary faces)
      std::vector<int> face_ids = CubitInterface::parse_cubit_list(
          "face", "in surface " + std::to_string(sid));
      for (int fid : face_ids) {
        append_face(QUAD4, CubitInterface::get_connectivity("face", fid), 4);
      }
    }

    // Direct mesh entities are essential for Sculpt/free-mesh output.  The
    // Cubit API documents these separately from get_sideset_surfaces().
    for (int tid : CubitInterface::get_sideset_tris(sg.id))
      append_face(TRI3, CubitInterface::get_connectivity("tri", tid), 3);
    for (int qid : CubitInterface::get_sideset_quads(sg.id))
      append_face(QUAD4, CubitInterface::get_connectivity("face", qid), 4);

    sidesets.push_back(std::move(sg));
  }
}

// ================================================================
// extract_nodesets
// ================================================================
void MeshData::extract_nodesets(MeshExportInterface *iface)
{
  // SNAPSHOT user-defined nodesets BEFORE calling iface->get_nodeset_list.
  // Same side-effect pattern as block + sideset above -- phantom
  // default-nodeset gets materialised by iface->id_from_handle() and
  // would pollute the Cubit session if not skipped here.  See full
  // diagnosis at extract_elements.
  std::vector<int> user_nodeset_ids = CubitInterface::parse_cubit_list("nodeset", "all");
  std::set<int> user_nodeset_set(user_nodeset_ids.begin(), user_nodeset_ids.end());

  std::vector<NodesetHandle> ns_handles;
  iface->get_nodeset_list(ns_handles);

  const int nbuf = 1000;
  std::vector<int> nids(nbuf);

  for (auto &nsh : ns_handles) {
    NodesetGroup ng;
    ng.id = iface->id_from_handle(nsh);
    // Skip phantom default-nodeset materialised by iface->get_nodeset_list.
    if (user_nodeset_set.find(ng.id) == user_nodeset_set.end())
      continue;
    ng.name = iface->get_nodeset_name(nsh);

    int nstart = 0, ncount = 0;
    while ((ncount = iface->get_nodeset_node_ids(nsh, nstart, nbuf, nids)) > 0) {
      for (int i = 0; i < ncount; i++)
        ng.node_ids.push_back(nids[i]);
      nstart += ncount;
    }
    nodesets.push_back(std::move(ng));
  }
}

// ================================================================
// total_element_count
// ================================================================
int MeshData::total_element_count() const
{
  int total = (int)elements.size();
  for (auto &sg : sidesets)
    total += (int)sg.faces.size();
  return total;
}

