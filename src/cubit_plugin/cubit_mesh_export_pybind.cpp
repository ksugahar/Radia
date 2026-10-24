// ============================================================
// cubit_mesh_export_pybind — pybind11 module for Cubit→Netgen mesh curving
//
// Cubit-free: receives mesh data and projection/normal callbacks
// from Python (cubit Python API), uses only Netgen C++ for
// high-order curving via BuildCurvedElements.
//
// Usage (external Python 3.12):
//   from cubit_mesh_export import cubit_mesh_curver
//   # Supply mesh arrays and geometry callbacks to build_curved_mesh.
// ============================================================

#include <meshing.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/numpy.h>
#include <pybind11/functional.h>

namespace py = pybind11;
namespace ng = netgen;

#include "NetgenCurverPure.hpp"

// ============================================================
// Module definition
// ============================================================
PYBIND11_MODULE(cubit_mesh_curver, m)
{
    m.doc() = R"doc(
        Cubit mesh curving (C++ accelerated, Cubit-free).

        Receives linear mesh data and geometry callbacks from Python,
        builds a curved Netgen mesh using BuildCurvedElements.
        Does NOT link to any Cubit C++ libraries.

        Use cubit_mesh_export.cubit_mesh_curver.build_curved_mesh with
        node arrays, element connectivity and geometry callbacks.
    )doc";

    // --------------------------------------------------------
    // build_curved_mesh — core function
    //
    // Accepts mesh data as numpy arrays + Python callables for
    // geometry projection/normal.  Returns netgen.meshing.Mesh.
    // --------------------------------------------------------
    m.def("build_curved_mesh",
        [](py::array_t<double, py::array::c_style> node_coords,  // (N, 3)
           py::array_t<int, py::array::c_style>    node_ids,      // (N,)
           py::list vol_elements,        // list of (type_code, [conn...])
           py::list surface_data,        // list of dicts: {id, tris, quads, uvs}
           py::function project_func,    // (surfnr, x,y,z, u,v) -> (xp,yp,zp, u,v)
           py::function normal_func,     // (surfnr, x,y,z) -> (nx,ny,nz)
           int order,
           bool surface_only,
           bool split_quads,
           py::object edge_project_func,
           py::list edge_elements) -> py::object
        {
            // Import netgen.meshing so pybind11 knows the Mesh type
            py::module_::import("netgen.meshing");

            // --- Convert node data ---
            if (node_ids.ndim() != 1 || node_coords.ndim() != 2 ||
                node_coords.shape(1) != 3 || node_coords.shape(0) != node_ids.shape(0))
                throw py::value_error("NetgenCurverPure: node_coords must have shape (len(node_ids), 3) and node_ids must be one-dimensional");
            auto nc_buf = node_coords.unchecked<2>();
            auto ni_buf = node_ids.unchecked<1>();
            int N = (int)ni_buf.shape(0);

            std::vector<double> coords(N * 3);
            std::vector<int> nids(N);
            for (int i = 0; i < N; i++) {
                nids[i] = ni_buf(i);
                coords[3*i]   = nc_buf(i, 0);
                coords[3*i+1] = nc_buf(i, 1);
                coords[3*i+2] = nc_buf(i, 2);
            }

            // --- Convert volume elements ---
            // Each element is (type_code, conn_list) or (type_code, conn_list, domain)
            std::vector<VolumeElement> velems;
            for (auto item : vol_elements) {
                auto tup = item.cast<py::tuple>();
                if (py::len(tup) != 2 && py::len(tup) != 3)
                    throw py::value_error("NetgenCurverPure: volume element must contain type, connectivity and optional domain");
                VolumeElement ve;
                const int type_code = tup[0].cast<int>();
                if (type_code < PURE_TET || type_code > PURE_PYRAMID)
                    throw py::value_error("NetgenCurverPure: unsupported volume element type " +
                                          std::to_string(type_code));
                ve.type = static_cast<PureElemType>(type_code);
                ve.conn = tup[1].cast<std::vector<int>>();
                ve.domain = (py::len(tup) > 2) ? tup[2].cast<int>() : 1;
                velems.push_back(ve);
            }

            // --- Convert surface data ---
            std::vector<SurfaceInfo> surfs;
            for (auto item : surface_data) {
                auto d = item.cast<py::dict>();
                SurfaceInfo si;
                si.id = d["id"].cast<int>();

                // tris: flat list [n0,n1,n2, n0,n1,n2, ...]
                auto tri_list = d["tris"].cast<std::vector<int>>();
                if (tri_list.size() % 3 != 0)
                    throw py::value_error("NetgenCurverPure: surface " + std::to_string(si.id) + " tris length must be divisible by 3");
                for (size_t t = 0; t < tri_list.size(); t += 3) {
                    si.tris.push_back({tri_list[t], tri_list[t+1], tri_list[t+2]});
                }

                // quads: flat list [n0,n1,n2,n3, ...]
                auto quad_list = d["quads"].cast<std::vector<int>>();
                if (quad_list.size() % 4 != 0)
                    throw py::value_error("NetgenCurverPure: surface " + std::to_string(si.id) + " quads length must be divisible by 4");
                for (size_t q = 0; q < quad_list.size(); q += 4) {
                    si.quads.push_back({quad_list[q], quad_list[q+1],
                                        quad_list[q+2], quad_list[q+3]});
                }

                // uvs: dict {node_id: (u, v)}
                auto uv_dict = d["uvs"].cast<py::dict>();
                for (auto [key, val] : uv_dict) {
                    int nid = key.cast<int>();
                    auto uv = val.cast<std::tuple<double, double>>();
                    si.uvs.push_back({nid, {std::get<0>(uv), std::get<1>(uv)}});
                }

                surfs.push_back(std::move(si));
            }

            // --- Wrap Python callables as std::function ---
            // surfnr in callbacks is 1-based FaceDescriptor index.
            // The Python side maps this to Cubit surface IDs.
            ProjectFunc proj = [project_func](int surfnr, double x, double y, double z,
                                               double u_hint, double v_hint, bool has_hint)
                -> std::tuple<double,double,double,double,double>
            {
                py::object result = project_func(surfnr, x, y, z, u_hint, v_hint, has_hint);
                return result.cast<std::tuple<double,double,double,double,double>>();
            };

            NormalFunc norm = [normal_func](int surfnr, double x, double y, double z)
                -> std::tuple<double,double,double>
            {
                py::object result = normal_func(surfnr, x, y, z);
                return result.cast<std::tuple<double,double,double>>();
            };

            // --- Edge projection callback (optional) ---
            EdgeProjectFunc edge_proj = nullptr;
            if (!edge_project_func.is_none()) {
                py::function epf = edge_project_func.cast<py::function>();
                edge_proj = [epf](int edge, int s1, int s2, double x, double y, double z)
                    -> std::tuple<double,double,double>
                {
                    py::object result = epf(edge, s1, s2, x, y, z);
                    return result.cast<std::tuple<double,double,double>>();
                };
            }

            // --- Convert edge elements ---
            std::vector<EdgeInfo> edges;
            for (auto item : edge_elements) {
                auto d = item.cast<py::dict>();
                EdgeInfo ei;
                ei.surfnr1 = d["surfnr1"].cast<int>();
                ei.surfnr2 = d["surfnr2"].cast<int>();
                auto seg_list = d["segments"].cast<py::list>();
                for (auto s : seg_list) {
                    auto pair = s.cast<std::vector<int>>();
                    if (pair.size() != 2)
                        throw py::value_error("NetgenCurverPure: edge " + std::to_string(edges.size() + 1) +
                                              " segment " + std::to_string(ei.segments.size() + 1) + " requires exactly 2 node IDs");
                    ei.segments.push_back({pair[0], pair[1]});
                }
                // dist (normalized curve parameter) per segment endpoint
                if (d.contains("dists")) {
                    auto dist_list = d["dists"].cast<py::list>();
                    for (auto dd : dist_list) {
                        auto dv = dd.cast<std::vector<double>>();
                        if (dv.size() != 2)
                            throw py::value_error("NetgenCurverPure: edge " + std::to_string(edges.size() + 1) +
                                                  " dists entry " + std::to_string(ei.dists.size() + 1) + " requires exactly 2 values");
                        ei.dists.push_back({dv[0], dv[1]});
                    }
                }
                edges.push_back(std::move(ei));
            }

            // --- Build curved mesh ---
            NetgenCurverPure curver;
            auto ng_mesh = curver.build(coords, nids, velems, surfs,
                                         proj, norm, order, edge_proj, edges);

            // --- surface_only: remove volume elements ---
            if (surface_only) {
                ng_mesh->ClearVolumeElements();
            }

            // --- split_quads: split quad surface elements into triangles ---
            if (split_quads) {
                std::vector<ng::Element2d> new_tris;

                for (auto & el : ng_mesh->SurfaceElements()) {
                    if (el.GetType() == ng::QUAD || el.GetType() == ng::QUAD6
                        || el.GetType() == ng::QUAD8) {
                        ng::Element2d tri1(ng::TRIG);
                        tri1[0] = el[0]; tri1[1] = el[1]; tri1[2] = el[2];
                        tri1.SetIndex(el.GetIndex());
                        if (el.GetType() == ng::QUAD) {
                            tri1.GeomInfo()[0] = el.GeomInfo()[0];
                            tri1.GeomInfo()[1] = el.GeomInfo()[1];
                            tri1.GeomInfo()[2] = el.GeomInfo()[2];
                        }
                        new_tris.push_back(tri1);

                        ng::Element2d tri2(ng::TRIG);
                        tri2[0] = el[0]; tri2[1] = el[2]; tri2[2] = el[3];
                        tri2.SetIndex(el.GetIndex());
                        if (el.GetType() == ng::QUAD) {
                            tri2.GeomInfo()[0] = el.GeomInfo()[0];
                            tri2.GeomInfo()[1] = el.GeomInfo()[2];
                            tri2.GeomInfo()[2] = el.GeomInfo()[3];
                        }
                        new_tris.push_back(tri2);

                        el.Delete();
                    }
                }

                for (auto & tri : new_tris)
                    ng_mesh->AddSurfaceElement(tri);

                ng_mesh->Compress();
                ng_mesh->UpdateTopology();
            }

            return py::cast(ng_mesh);
        },
        py::arg("node_coords"),
        py::arg("node_ids"),
        py::arg("vol_elements"),
        py::arg("surface_data"),
        py::arg("project_func"),
        py::arg("normal_func"),
        py::arg("order") = 2,
        py::arg("surface_only") = false,
        py::arg("split_quads") = false,
        py::arg("edge_project_func") = py::none(),
        py::arg("edge_elements") = py::list(),
        R"pbdoc(
            Build a curved Netgen mesh from external mesh data and geometry callbacks.

            This function does NOT call any Cubit C++ APIs. All mesh data and
            geometry queries are provided by the caller (typically via cubit Python API).

            Args:
                node_coords: Node coordinates, shape (N, 3), float64.
                node_ids: Node IDs, shape (N,), int32.
                vol_elements: Volume elements as list of (type_code, conn_list[, domain]).
                    type_code: 0=TET, 1=HEX, 2=PRISM, 3=PYRAMID.
                    domain: 1-based material domain index (default 1).
                surface_data: Per-surface data as list of dicts with keys:
                    id (int), tris (flat list), quads (flat list),
                    uvs (dict {node_id: (u, v)}).
                project_func: Callable (surfnr, x,y,z, u_hint, v_hint, has_hint) -> (xp,yp,zp, u,v).
                    surfnr is 1-based FaceDescriptor index.
                    has_hint=True signals (u_hint, v_hint) is a valid UV
                    (use closest_point_uv_guess); has_hint=False means no
                    hint, do a global projection. Distinguishes legitimate
                    UV=(0,0) from "no hint".
                normal_func: Callable (surfnr, x,y,z) -> (nx,ny,nz).
                edge_project_func: Callable (edgenr, surfnr1, surfnr2, x,y,z)
                    -> (xp,yp,zp). edgenr is the geometric EdgeDescriptor.EdgeNr(),
                    not a surface-pair key or the descriptor table index.
                order: Polynomial order for curving (>= 2).
                surface_only: Remove volume elements (BEM workflow).
                split_quads: Split quad faces into triangles.

            Returns:
                netgen.meshing.Mesh: Curved mesh.
        )pbdoc");

    // --------------------------------------------------------
    // Version and capability info
    // --------------------------------------------------------
    // Must equal packages/cubit-mesh-export/pyproject.toml (contract test).
    m.attr("__version__") = "2.2.6";
    m.attr("has_netgen") = true;
}
