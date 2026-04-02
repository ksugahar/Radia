#ifndef EXPORT_MEG_COMMAND_HPP
#define EXPORT_MEG_COMMAND_HPP

#include "CubitCommandInterface.hpp"
#include <string>
#include <vector>

class MeshData;

class ExportMegCommand : public CubitCommand
{
public:
  ExportMegCommand();
  ~ExportMegCommand();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute(CubitCommandData &data);

  // dim: 'T' (3D), 'K' (2D, z=0), 'R' (axisymmetric, y=0, x>0)
  bool write_meg(const std::string &filename, char dim);

private:
  // ELF element type string: first 4 chars of block name + DIM char
  static std::string elf_type(const std::string &block_name, char dim);
};

#endif
