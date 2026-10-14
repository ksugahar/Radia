#ifndef UTF8_PATH_HPP
#define UTF8_PATH_HPP

// u8_string_to_path: UTF-8 std::string -> std::filesystem::path that
// survives Unicode / Japanese directory names on Windows.
//
// Cubit passes filename arguments as std::string carrying UTF-8 bytes.
// On Japanese Windows (system codepage cp932), the implicit
// std::filesystem::path(std::string) constructor decodes the bytes as
// cp932 and mojibakes them — downstream std::ofstream / netgen Save
// calls then fail with:
//
//   "No mapping for the Unicode character exists in the target
//    multi-byte code page."
//
// This helper converts UTF-8 -> UTF-16 via MultiByteToWideChar(CP_UTF8)
// and constructs std::filesystem::path from the wide string.  On
// non-Windows builds the std::string path is passed through unchanged.
//
// Observed 2026-04-21 with path `C:\temp\日本語テスト\3turncoil.stp`
// failing every `export netgen / gmsh / nastran / vtk / meg` call.

#include <filesystem>
#include <string>

#ifdef _WIN32
#include <windows.h>
#endif

inline std::filesystem::path u8_string_to_path(const std::string &s)
{
#ifdef _WIN32
  if (s.empty()) return std::filesystem::path();
  int wlen = MultiByteToWideChar(CP_UTF8, 0, s.c_str(),
                                 static_cast<int>(s.size()), nullptr, 0);
  if (wlen <= 0) return std::filesystem::path(s);  // decode failure fallback
  std::wstring w(static_cast<size_t>(wlen), 0);
  MultiByteToWideChar(CP_UTF8, 0, s.c_str(),
                      static_cast<int>(s.size()), &w[0], wlen);
  return std::filesystem::path(w);
#else
  return std::filesystem::path(s);
#endif
}

#endif
