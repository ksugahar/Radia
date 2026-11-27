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

  m.def("edge_coefficient_count", [](netgen::Mesh & mesh, int edge)
        { return nc::curved_access(mesh).edge_count(edge); }, py::arg("mesh"), py::arg("edge"));
}
