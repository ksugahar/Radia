#ifndef EXPORT_FEMEEM_COMMAND_HPP
#define EXPORT_FEMEEM_COMMAND_HPP

#include "CubitCommandInterface.hpp"
#include <fstream>
#include <string>
#include <vector>

class MeshData;

class ExportFemeemCommand : public CubitCommand
{
public:
  ExportFemeemCommand();
  ~ExportFemeemCommand();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute(CubitCommandData &data);

  bool write_femeem(const std::string &dirname, double scale);

private:
  bool write_in_dat(const std::string &path, const MeshData &mesh, double scale);
  bool write_sin_dat_B(const std::string &path, const MeshData &mesh);
  bool write_sina_dat(const std::string &path);
  bool write_d3(const std::string &path);
};

#endif
