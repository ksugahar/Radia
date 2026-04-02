#ifndef EXPORT_NASTRAN_COMMAND_HPP
#define EXPORT_NASTRAN_COMMAND_HPP

#include "CubitCommandInterface.hpp"
#include <fstream>
#include <string>
#include <vector>

class MeshData;
struct MeshElement;

class ExportNastranCommand : public CubitCommand
{
public:
  ExportNastranCommand();
  ~ExportNastranCommand();

  std::vector<std::string> get_syntax();
  std::vector<std::string> get_syntax_help();
  std::vector<std::string> get_help();
  bool execute(CubitCommandData &data);

  bool write_nastran(const std::string &filename, const std::string &dim,
                     bool pyram, int order);

private:

  // Shared helpers operating on MeshData
  void write_header(std::ofstream &fid, const std::string &filename);
  void write_nodes(std::ofstream &fid, const MeshData &mesh, bool is_3d);
  int  write_elements(std::ofstream &fid, const MeshData &mesh,
                      bool is_3d, bool pyram, int start_eid);
  int  write_sidesets(std::ofstream &fid, const MeshData &mesh,
                      bool is_3d, int start_eid);
  void write_nodesets(std::ofstream &fid, const MeshData &mesh);
  void write_properties(std::ofstream &fid, const MeshData &mesh);

  // Write one element (block or sideset face) — returns next eid
  int write_element_card(std::ofstream &fid, const MeshElement &elem,
                         int eid, int pid, bool pyram, int order);
};

#endif
