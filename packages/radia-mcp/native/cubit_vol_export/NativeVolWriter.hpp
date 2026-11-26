// Single writer for Netgen .vol text written by the exporter.
//
// Mesh::Save(path) opens a plain .vol in text mode, so on Windows every line
// ends in CRLF. Netgen strips that CR when it reads a plain .vol, but its
// .vol.gz reader keeps it, which turns region names into "body\r": a solver
// that selects "body" then integrates over nothing. The format is LF text, so
// the exporter writes the serialized bytes itself in binary mode.
#pragma once

#include <filesystem>
#include <fstream>
#include <sstream>
#include <stdexcept>
#include <string>

#include <meshing.hpp>

namespace cubit_mesh_export {

inline void save_native_mesh(const netgen::Mesh &mesh, const std::filesystem::path &path)
{
  if (path.extension() != ".vol") {
    // .vol.gz, .vol.bin and Netgen's own suffix rules stay with Netgen.
    mesh.Save(path);
    return;
  }
  std::ostringstream text;
  mesh.Save(text);
  if (!text) throw std::runtime_error("Netgen .vol serialization failed");
  const std::string bytes = text.str();
  if (bytes.find('\r') != std::string::npos)
    throw std::runtime_error("Netgen .vol serialization contains a carriage return");
  std::ofstream out(path, std::ios::binary | std::ios::trunc);
  if (!out) throw std::runtime_error("Cannot open .vol output: " + path.u8string());
  out.write(bytes.data(), static_cast<std::streamsize>(bytes.size()));
  out.close();
  if (!out) throw std::runtime_error("Incomplete .vol output: " + path.u8string());
}

}  // namespace cubit_mesh_export
