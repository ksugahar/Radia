// _netgen_curving: Python access to netgen_curving/ for the regression tests.
// Takes netgen.meshing.Mesh objects; built against the same Netgen SDK (and its
// bundled pybind11) as the netgen package loaded at run time.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "curved_access.hpp"
#include "geometric_refit.hpp"
#include "periodic_copy.hpp"

namespace py = pybind11;
namespace nc = radia_netgen_curving;

// Test helpers only: range- and storage-check an edge index before touching
// coefficient storage (polynomial orders 1..3, as ValidatePeriodicStorage requires).
static void checked_edge(netgen::Mesh & mesh, int edge)
{
  const auto & top = mesh.GetTopology();
  if (!mesh.GetCurvedElements().IsHighOrder())
    throw std::runtime_error("edge coefficients: the mesh is not curved (call Curve first)");
  if (edge < 0 || size_t(edge) >= top.GetNEdges())
    throw std::out_of_range("edge coefficients: edge index out of range");
  nc::curved_access(mesh).validate(top.GetNEdges(), top.GetNFaces());
}

PYBIND11_MODULE(_netgen_curving, m)
{
  py::module_::import("netgen.meshing");   // registers netgen::Mesh
  m.doc() = "GeometricRefit and periodic coefficient copy, moved out of the Netgen fork";

  m.def("geometric_refit", [](netgen::Mesh & mesh)
        {
          const nc::RefitStats st = nc::geometric_refit(mesh);   // GIL held, as in the fork binding
          py::dict d;
          d["edges_tried"] = st.edges_tried; d["edges_accepted"] = st.edges_accepted;
          d["edges_failed"] = st.edges_failed; d["edges_rejected_distortion"] = st.edges_rejected_distortion;
          d["faces_tried"] = st.faces_tried; d["faces_accepted"] = st.faces_accepted;
          d["faces_failed"] = st.faces_failed; d["faces_rejected_distortion"] = st.faces_rejected_distortion;
          d["edge_dist_before"] = st.edge_dist_before; d["edge_dist_after"] = st.edge_dist_after;
          d["face_dist_before"] = st.face_dist_before; d["face_dist_after"] = st.face_dist_after;
          return d;
        }, py::arg("mesh"),
        "Refit curved-element coefficients of order >= 3 to the mesh geometry; "
        "call after Curve(order). Returns statistics.");

  m.def("copy_periodic_coefficients",
        [](netgen::Mesh & mesh, std::vector<std::array<int, 3>> master,
           std::vector<std::array<int, 3>> slave)
        {
          py::list rows;
          for (const auto & tr : nc::copy_periodic_coefficients(mesh, master, slave))
            {
              py::dict d;
              d["entity"] = tr.face ? "face" : "edge";
              d["master"] = tr.master; d["slave"] = tr.slave;
              d["mode"] = tr.mode; d["sign"] = tr.sign;
              d["original"] = py::make_tuple(tr.original[0], tr.original[1], tr.original[2]);
              d["expected"] = py::make_tuple(tr.expected[0], tr.expected[1], tr.expected[2]);
              rows.append(d);
            }
          return rows;
        }, py::arg("mesh"), py::arg("master"), py::arg("slave"),
        "Copy curved coefficients from master triangles (1-based point numbers) "
        "to their periodic images; returns the applied transfers.");

  // Test-only helpers: locate a topology edge by its 1-based end points, and
  // read or overwrite its coefficients so a test can plant known values.
  // Every edge index is range-checked against the topology and the storage.
  m.def("edge_number", [](netgen::Mesh & mesh, int a, int b)
        {
          const auto & top = mesh.GetTopology();
          if (a > b) std::swap(a, b);
          for (size_t e = 0; e < top.GetNEdges(); e++)
            {
              auto v = top.GetEdgeVertices(netgen::IndexBASE<netgen::EdgeIndex>() + e);
              int p = v[0].Nr1(), q = v[1].Nr1();
              if (p > q) std::swap(p, q);
              if (p == a && q == b) return int(e);
            }
          throw std::runtime_error("edge_number: no topology edge between the given points");
        }, py::arg("mesh"), py::arg("a"), py::arg("b"));
  m.def("edge_coefficients", [](netgen::Mesh & mesh, int edge)
        {
          checked_edge(mesh, edge);
          auto access = nc::curved_access(mesh);
          py::list out;
          const int n = access.edge_count(edge);
          for (int i = 0; i < n; i++)
            {
              const auto & c = access.edge_coeffs(edge)[i];
              out.append(py::make_tuple(c(0), c(1), c(2)));
            }
          return out;
        }, py::arg("mesh"), py::arg("edge"));
  m.def("set_edge_coefficients", [](netgen::Mesh & mesh, int edge,
                                    std::vector<std::array<double, 3>> values)
        {
          checked_edge(mesh, edge);
          auto access = nc::curved_access(mesh);
          if (int(values.size()) != access.edge_count(edge))
            throw std::runtime_error("set_edge_coefficients: wrong number of coefficients");
          for (size_t i = 0; i < values.size(); i++)
            for (int k = 0; k < 3; k++) access.edge_coeffs(edge)[i](k) = values[i][k];
        }, py::arg("mesh"), py::arg("edge"), py::arg("values"));
}
