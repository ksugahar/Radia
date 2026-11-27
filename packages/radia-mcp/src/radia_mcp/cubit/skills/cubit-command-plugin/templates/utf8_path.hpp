#ifndef MY_PLUGIN_UTF8_PATH_HPP
#define MY_PLUGIN_UTF8_PATH_HPP

// Cubit passes file names as std::string holding UTF-8.  On Japanese
// Windows (cp932) std::filesystem::path(std::string) decodes them as cp932,
// so a path like C:\temp\日本語\a.vol fails to open.  Convert via UTF-16.
// UTF-8 is the input contract: an empty or invalid name is rejected rather
// than reinterpreted in the active code page.

#include <filesystem>
#include <string>

#ifdef _WIN32
#ifndef NOMINMAX
#define NOMINMAX
#endif
#include <windows.h>
#endif

inline bool u8_string_to_path(const std::string& s, std::filesystem::path& out)
{
  if (s.empty()) return false;
#ifdef _WIN32
  const int len = static_cast<int>(s.size());
  const int n = MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, s.data(), len, nullptr, 0);
  if (n <= 0) return false;
  std::wstring w(static_cast<size_t>(n), L'\0');
  if (MultiByteToWideChar(CP_UTF8, MB_ERR_INVALID_CHARS, s.data(), len, w.data(), n) != n)
    return false;
  out = std::filesystem::path(w);
#else
  out = std::filesystem::path(s);
#endif
  return true;
}

#endif
