// _cubit_to_vol_mesh: the `export vol` mesh in memory, still attached to Cubit.
// Built against the same Netgen SDK (and its bundled pybind11) as the netgen
// package loaded at run time, so the result is a netgen.meshing.Mesh.  Load it
// in an external Python after `import cubit` and `cubit.init(...)`; Cubit's
// embedded Python never loads it.
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include "ExportVolCommand.hpp"
#include "LiveVolMesh.hpp"
#include "netgen_meshing.hpp"

#include <map>
#include <optional>
#include <stdexcept>
#include <string>

namespace py = pybind11;

static double periodic_number(const std::map<std::string, std::string> &fields, const char *name)
{
  const std::string &text = fields.at(name);
  size_t used = 0;
  const double value = std::stod(text, &used);
  if (used != text.size()) throw std::invalid_argument(std::string(name) + " is not a number: " + text);
  return value;
}

PYBIND11_MODULE(_cubit_to_vol_mesh, m)
{
  py::module_::import("netgen.meshing");   // registers netgen::Mesh
  m.doc() = "The export vol mesh in memory, attached to the live Cubit geometry";

  m.def("build",
        [](int order, std::optional<std::array<double, 3>> kelvin_translation,
           std::optional<double> kelvin_match,
           std::optional<std::map<std::string, std::string>> periodic)
        {
          VolExportRequest request;
          request.order = order;
          if (kelvin_translation.has_value() != kelvin_match.has_value())
            throw std::invalid_argument("Kelvin needs both kelvin_translation and kelvin_match");
          if (kelvin_translation) {
            request.has_kelvin = true;
            request.kelvin_translation = *kelvin_translation;
            request.kelvin_match = *kelvin_match;
          }
          if (periodic) {
            const auto &f = *periodic;
            for (const char *name : {"periodic_name", "periodic_master", "periodic_slave",
                                     "periodic_unit", "periodic_patches", "periodic_tx",
                                     "periodic_ty", "periodic_tz", "periodic_match"})
              if (!f.count(name)) throw std::invalid_argument(std::string("periodic needs ") + name);
            if (f.size() != 9) throw std::invalid_argument("periodic takes exactly the nine fields");
            request.has_periodic = true;
            request.periodic.name = f.at("periodic_name");
            request.periodic.master_label = f.at("periodic_master");
            request.periodic.slave_label = f.at("periodic_slave");
            request.periodic.coordinate_unit = f.at("periodic_unit");
            request.periodic_patch_text = f.at("periodic_patches");
            request.periodic.translation = {periodic_number(f, "periodic_tx"),
                                            periodic_number(f, "periodic_ty"),
                                            periodic_number(f, "periodic_tz")};
            request.periodic.matching_distance = periodic_number(f, "periodic_match");
          }
          std::shared_ptr<netgen::Mesh> mesh;
          if (!ExportVolCommand::run(request, &mesh) || !mesh)
            throw std::runtime_error("the export vol mesh could not be built; the Cubit "
                                     "messages above give the reason");
          return mesh;
        },
        py::arg("order"), py::arg("kelvin_translation") = py::none(),
        py::arg("kelvin_match") = py::none(), py::arg("periodic") = py::none(),
        "Build the labelled, curved export vol mesh from the current Cubit model "
        "without writing a file; it stays attached to the Cubit geometry.");

  m.def("curve",
        [](std::shared_ptr<netgen::Mesh> mesh, int order)
        {
          if (!mesh) throw std::invalid_argument("mesh is None");
          cubit_to_vol_export::recurve(*mesh, order);
        },
        py::arg("mesh"), py::arg("order"),
        "Curve a mesh from build() again at order 1-5 on the live Cubit geometry, "
        "with the same refit and periodic coefficient copy as export vol.");

  m.def("is_live",
        [](std::shared_ptr<netgen::Mesh> mesh)
        {
          return bool(mesh) && bool(std::dynamic_pointer_cast<cubit_to_vol_export::LiveCubitGeometry>(
                                   mesh->GetGeometry()));
        },
        py::arg("mesh"), "True for a mesh from build() that still carries its Cubit geometry.");
}
