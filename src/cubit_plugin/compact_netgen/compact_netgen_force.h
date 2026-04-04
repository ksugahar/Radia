// Force include: override DLL_HEADER and NGCORE_API to empty for static linking
#ifndef COMPACT_NETGEN_FORCE_H
#define COMPACT_NETGEN_FORCE_H

#ifndef COMPACT_NETGEN_STATIC
#define COMPACT_NETGEN_STATIC
#endif

// Pre-define these to prevent dllexport/dllimport
#define NGCORE_API_EXPORT
#define NGCORE_API_IMPORT
#define NGCORE_API
#define DLL_HEADER

#endif // COMPACT_NETGEN_FORCE_H
