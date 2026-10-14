#ifndef EXPORT_MEG_COMMAND_HPP
#define EXPORT_MEG_COMMAND_HPP

#include "CubitCommandInterface.hpp"
#include "ElementType.h"
#include <string>
#include <vector>
#include <unordered_map>

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
  // label_map: block_id -> 3-char ELF prefix (from GUI or command line)
  bool write_meg(const std::string &filename, char dim,
                 const std::unordered_map<int, std::string> &label_map = {});

private:
  static std::string elf_type(const std::string &prefix, int nv,
                              ElementType type, char dim);
};

#endif
