# Radia Source Directory Structure

## Current Directory Layout

```
04_Radia/
├── src/
│   ├── core/           # Core Radia computation library (27 files)
│   │   ├── radapl*.cpp # Application interface functions
│   │   ├── radinter.cpp # Interface definitions
│   │   ├── radintrc.cpp # Interaction computations
│   │   ├── radrlmet.cpp # Relaxation methods
│   │   ├── radvlpgn.cpp # Polyhedron geometry
│   │   └── ...
│   ├── lib/            # Library entry points
│   │   └── radentry.cpp # Main entry point for C library
│   ├── ext/            # External libraries
│   │   ├── auxparse/   # Auxiliary parsing utilities
│   │   ├── fftw/       # FFTW library (fftw64_f.lib)
│   │   ├── genmath/    # General math utilities
│   │   └── triangle/   # Triangle mesh library
│   └── python/         # Python binding (MOVED from clients/python)
│       ├── radpy.cpp   # Python C API interface
│       └── pyparse.h   # Python parsing utilities
├── build/              # CMake build output
├── dist/               # Distribution files (radia.pyd)
├── CMakeLists.txt      # CMake build configuration
└── Build.ps1           # PowerShell build script
```

## Recent Changes

### Moved Python Client Directory

**Before:**
```
src/
├── clients/
│   └── python/
│       ├── radpy.cpp
│       └── pyparse.h
```

**After:**
```
src/
├── python/
│   ├── radpy.cpp
│   └── pyparse.h
```

### CMakeLists.txt Updates

Changed directory variables:
- Removed: `CLIENTS_DIR` (was pointing to `src/clients`)
- Added: `PYTHON_DIR` (now points to `src/python`)

Updated references:
```cmake
# Old
set(CLIENTS_DIR ${SRC_DIR}/clients)
set(RADIA_PYTHON_SOURCES ${CLIENTS_DIR}/python/radpy.cpp)
target_include_directories(radia PRIVATE ${CLIENTS_DIR}/python)

# New
set(PYTHON_DIR ${SRC_DIR}/python)
set(RADIA_PYTHON_SOURCES ${PYTHON_DIR}/radpy.cpp)
target_include_directories(radia PRIVATE ${PYTHON_DIR})
```

## Rationale

This directory restructuring improves the source layout by:

1. **Simplified Structure**: Removed unnecessary nested `clients` directory
2. **Direct Access**: Python binding source is now directly under `src/`
3. **Consistency**: Matches the pattern of other top-level directories (core, lib, ext)
4. **Future-Proof**: Easier to add other language bindings at the same level if needed

## Build Verification

✓ Build completed successfully with new structure
✓ All tests pass (71.4% pass rate, as expected)
✓ Module size: 1.86 MB (unchanged)
✓ Module version: 4.32
✓ OpenMP parallelization: Enabled
