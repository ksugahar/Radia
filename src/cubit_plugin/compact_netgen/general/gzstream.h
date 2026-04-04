// Compact Netgen: stub gzstream.h (no zlib dependency)
// .vol.gz not supported, .vol text only
#ifndef GZSTREAM_H_STUB
#define GZSTREAM_H_STUB

#include <fstream>

// Stub: ogzstream/igzstream are just regular file streams
// (Mesh::Save will use the text path, not .vol.gz)
class ogzstream : public std::ofstream {
public:
  ogzstream() {}
  explicit ogzstream(const char* name, std::ios_base::openmode mode = std::ios_base::out)
    : std::ofstream(name, mode) {}
  explicit ogzstream(const std::string& name, std::ios_base::openmode mode = std::ios_base::out)
    : std::ofstream(name, mode) {}
  template<typename T>
  explicit ogzstream(const T& name, std::ios_base::openmode mode = std::ios_base::out)
    : std::ofstream(name, mode) {}
};

class igzstream : public std::ifstream {
public:
  igzstream() {}
  explicit igzstream(const char* name, std::ios_base::openmode mode = std::ios_base::in)
    : std::ifstream(name, mode) {}
  explicit igzstream(const std::string& name, std::ios_base::openmode mode = std::ios_base::in)
    : std::ifstream(name, mode) {}
  template<typename T>
  explicit igzstream(const T& name, std::ios_base::openmode mode = std::ios_base::in)
    : std::ifstream(name, mode) {}
};

#endif // GZSTREAM_H_STUB
