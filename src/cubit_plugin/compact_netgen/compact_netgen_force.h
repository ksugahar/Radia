// Force include: override DLL_HEADER and NGCORE_API to empty for static linking.
// Included via MSVC /FI before any netgen header, so these #defines
// take effect before mydefs.hpp / ngcore_api.hpp are parsed.
// Both netgen_fork headers have #ifndef guards, so these pre-definitions
// prevent them from setting __declspec(dllimport).
#ifndef COMPACT_NETGEN_FORCE_H
#define COMPACT_NETGEN_FORCE_H

#ifndef COMPACT_NETGEN_STATIC
#define COMPACT_NETGEN_STATIC
#endif

#define NGCORE_API_EXPORT
#define NGCORE_API_IMPORT
#define NGCORE_API
#define DLL_HEADER

#endif // COMPACT_NETGEN_FORCE_H
