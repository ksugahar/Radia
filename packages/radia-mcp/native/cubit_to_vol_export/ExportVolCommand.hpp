#ifndef EXPORT_VOL_COMMAND_HPP
#define EXPORT_VOL_COMMAND_HPP

#include "CubitCommandInterface.hpp"
#include "NamedPeriodicIdentification.hpp"
#include <array>
#include <memory>
#include <string>
#include <vector>

namespace netgen { class Mesh; }

/// One `export vol` request, parsed from the Cubit command or given by the
/// in-memory mesh API.  The named periodic fields are all present or all absent,
/// as are the Kelvin fields; `has_*` says which.
struct VolExportRequest
{
  std::string filename;   // unused for an in-memory mesh
  int order = 2;
  bool has_periodic = false;
  cubit_to_vol_export::NamedPeriodicContract periodic;
  std::string periodic_patch_text;   // "master:slave,..."
  bool has_kelvin = false;
  std::array<double, 3> kelvin_translation{};
  double kelvin_match = 0;
};

/// APREPRO command: export vol "file.vol" [order N]
/// Pure C++ export via MeshData + NetgenCurver (no subprocess).
/// Produces .vol + companion .json with CAD reference values.
class ExportVolCommand : public CubitCommand
{
public:
  ExportVolCommand();
  ~ExportVolCommand();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute(CubitCommandData &data);

  /// The whole export.  With `live` null it writes the .vol and its companion
  /// JSON; otherwise it writes nothing and stores the labelled, curved mesh in
  /// `*live`, still attached to its live Cubit geometry.  False on failure,
  /// reported through Cubit's messages.
  static bool run(const VolExportRequest &request, std::shared_ptr<netgen::Mesh> *live);
};

#endif
