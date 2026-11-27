#include "ExportNodeCountCommand.hpp"
#include "utf8_path.hpp"

#include "CubitInterface.hpp"
#include "CubitMessage.hpp"
#include "MeshExportInterface.hpp"

#include <cerrno>
#include <cstdio>
#include <cstring>
#include <filesystem>
#include <memory>
#include <set>
#include <string>

namespace {

// Open for writing.  Without `overwrite` the create is exclusive ("x"), so
// a file that appears between any check and the open is never truncated.
std::FILE* open_output(const std::filesystem::path& path, bool overwrite)
{
#ifdef _WIN32
  std::FILE* f = nullptr;
  if (_wfopen_s(&f, path.c_str(), overwrite ? L"wb" : L"wbx") != 0) return nullptr;
  return f;
#else
  return std::fopen(path.c_str(), overwrite ? "wb" : "wbx");
#endif
}

}  // namespace

std::vector<std::string> ExportNodeCountCommand::get_syntax()
{
  // `export <new keyword>` extends the built-in export verb.  Pick a keyword
  // no built-in exporter uses.
  return {
    "export nodecount <string:label='filename',help='<filename>'> [overwrite]"
  };
}

std::vector<std::string> ExportNodeCountCommand::get_syntax_help()
{
  return {"export nodecount \"filename\" [overwrite]"};
}

std::vector<std::string> ExportNodeCountCommand::get_help()
{
  return {"Write node, user-block and element counts of the current mesh.\n"};
}

bool ExportNodeCountCommand::execute(CubitCommandData& data)
{
  std::string filename;
  data.get_string("filename", filename);
  const bool overwrite = data.find_keyword("overwrite");

  // Cubit hands over UTF-8 bytes; build the path from UTF-16 on Windows or
  // a Japanese directory name fails to open.
  std::filesystem::path path;
  if (!u8_string_to_path(filename, path)) {
    PRINT_ERROR("File name is empty or not valid UTF-8: %s\n", filename.c_str());
    return false;
  }

  // Interfaces from CubitInterface must be released; an RAII holder keeps
  // every early return correct.
  std::unique_ptr<CubitBaseInterface, void (*)(CubitBaseInterface*)> raw(
      CubitInterface::get_interface("MeshExport"),
      [](CubitBaseInterface* p) { if (p) CubitInterface::release_interface(p); });
  auto* iface = dynamic_cast<MeshExportInterface*>(raw.get());
  if (!iface) {
    PRINT_ERROR("MeshExport interface is not available\n");
    return false;
  }
  iface->set_use_sequential_ids(false);  // keep Cubit's ids (returns void)
  if (!iface->initialize_export()) {
    PRINT_ERROR("MeshExport initialize_export() failed\n");
    return false;
  }

  const int nodes = iface->get_num_nodes();
  if (nodes <= 0) {
    PRINT_ERROR("No mesh nodes to export\n");
    return false;
  }

  // BlockHandle is a pointer (ElementOutputGroup*), not an id.  Map it with
  // id_from_handle() and keep only the user's blocks: with unassigned mesh
  // the list can also hold a group that is not a user block.
  const std::vector<int> user_ids = CubitInterface::parse_cubit_list("block", "all");
  const std::set<int> wanted(user_ids.begin(), user_ids.end());
  std::vector<BlockHandle> blocks;
  iface->get_block_list(blocks);
  long long elements = 0;
  int user_blocks = 0;
  for (BlockHandle block : blocks) {
    if (!wanted.count(iface->id_from_handle(block))) continue;
    int n = 0;
    if (!iface->get_block_size(block, n) || n < 0) {
      PRINT_ERROR("get_block_size failed for block %d\n", iface->id_from_handle(block));
      return false;
    }
    elements += n;
    ++user_blocks;
  }

  std::FILE* out = open_output(path, overwrite);
  if (!out) {
    if (!overwrite && errno == EEXIST)
      PRINT_ERROR("File exists: %s (add 'overwrite')\n", filename.c_str());
    else
      PRINT_ERROR("Cannot open %s for writing: %s\n", filename.c_str(), std::strerror(errno));
    return false;
  }
  const int written = std::fprintf(out, "nodes %d\nblocks %d\nblock_elements %lld\n",
                                   nodes, user_blocks, elements);
  const bool closed = std::fclose(out) == 0;
  if (written < 0 || !closed) {
    PRINT_ERROR("Write failed: %s\n", filename.c_str());
    return false;
  }
  PRINT_INFO("MYPLUGIN_EXPORT nodes=%d blocks=%d elements=%lld\n",
             nodes, user_blocks, elements);
  return true;
}
