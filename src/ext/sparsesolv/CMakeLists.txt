cmake_minimum_required(VERSION 3.16)
project(SparseSolv VERSION 2.6.0 LANGUAGES CXX)

# Header-only library
add_library(sparsesolv INTERFACE)
target_include_directories(sparsesolv INTERFACE
    $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
    $<INSTALL_INTERFACE:include>
)
target_compile_features(sparsesolv INTERFACE cxx_std_17)

# OpenMP (optional)
option(SPARSESOLV_USE_OPENMP "Use OpenMP for parallelization" ON)
if(SPARSESOLV_USE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        target_link_libraries(sparsesolv INTERFACE OpenMP::OpenMP_CXX)
    endif()
endif()

# Alias for consistent naming
add_library(SparseSolv::sparsesolv ALIAS sparsesolv)

# ==============================================================================
# NGSolve Python extension module (optional)
# ==============================================================================
option(SPARSESOLV_BUILD_NGSOLVE "Build NGSolve Python extension module" OFF)

if(SPARSESOLV_BUILD_NGSOLVE)
    # Auto-detect NGSolve and Netgen cmake directories from Python.
    # pip-installed NGSolve places NGSolveConfig.cmake in ngsolve/cmake/
    # and NetgenConfig.cmake in netgen/cmake/ (separate packages).
    execute_process(
        COMMAND python -c "
import pathlib
paths = []
try:
    import ngsolve
    d = pathlib.Path(ngsolve.__file__).parent
    for p in [d/'cmake', d.parent, d.parent.parent/'cmake']:
        if p.exists(): paths.append(p.as_posix())
except ImportError: pass
try:
    import netgen
    p = pathlib.Path(netgen.__file__).parent / 'cmake'
    if p.exists(): paths.append(p.as_posix())
except ImportError: pass
print(';'.join(paths))
"
        OUTPUT_VARIABLE _PYTHON_HINTS
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    list(APPEND CMAKE_PREFIX_PATH ${_PYTHON_HINTS} "${NGSOLVE_INSTALL_DIR}")
    find_package(NGSolve CONFIG REQUIRED
        HINTS ${_PYTHON_HINTS} "${NGSOLVE_INSTALL_DIR}"
    )
    message(STATUS "Found NGSolve: ${NGSOLVE_INCLUDE_DIRS}")

    # MKL (BLAS/LAPACK for dense matrix operations)
    # Use LP64 interface (MKL_INT = int, 32-bit) to match sparsesolv::index_t
    set(MKL_INTERFACE_FULL intel_lp64)
    # Auto-detect MKL cmake config: oneAPI toolkit, pip mkl-devel, or pip mkl
    execute_process(
        COMMAND python -c "
import pathlib, sys
paths = []
# pip-installed mkl-devel: sys.prefix/Library/lib/cmake/mkl (Windows)
prefix = pathlib.Path(sys.prefix)
for p in [prefix/'Library'/'lib'/'cmake'/'mkl',
          prefix/'lib'/'cmake'/'mkl']:
    if p.exists(): paths.append(p.as_posix())
# import mkl fallback
try:
    import mkl
    d = pathlib.Path(mkl.__file__).parent
    for p in [d/'lib'/'cmake'/'mkl']:
        if p.exists(): paths.append(p.as_posix())
except ImportError: pass
if paths: print(paths[0])
"
        OUTPUT_VARIABLE _MKL_HINTS
        OUTPUT_STRIP_TRAILING_WHITESPACE
    )
    find_package(MKL CONFIG REQUIRED
        HINTS ${_MKL_HINTS}
              "C:/Program Files (x86)/Intel/oneAPI/mkl/latest/lib/cmake/mkl")
    message(STATUS "Found MKL: ${MKL_ROOT}")

    add_ngsolve_python_module(sparsesolv_ngsolve ngsolve/python_module.cpp)

    target_include_directories(sparsesolv_ngsolve PRIVATE
        ${CMAKE_CURRENT_SOURCE_DIR}/include
    )
    target_compile_features(sparsesolv_ngsolve PRIVATE cxx_std_17)
    if(MSVC)
        target_compile_options(sparsesolv_ngsolve PRIVATE /arch:AVX2)
    else()
        target_compile_options(sparsesolv_ngsolve PRIVATE -mavx2 -mfma)
    endif()
    target_link_libraries(sparsesolv_ngsolve PRIVATE MKL::MKL)

    # Use NGSolve TaskManager for parallel_for/parallel_reduce_sum
    target_compile_definitions(sparsesolv_ngsolve PRIVATE SPARSESOLV_USE_NGSOLVE_TASKMANAGER)

    # MKL BLAS/LAPACK for dense matrix operations
    target_compile_definitions(sparsesolv_ngsolve PRIVATE SPARSESOLV_USE_MKL)

    # Install the Python module
    if(SKBUILD)
        # scikit-build-core: wheel-relative install path
        install(TARGETS sparsesolv_ngsolve
            LIBRARY DESTINATION sparsesolv_ngsolve
            RUNTIME DESTINATION sparsesolv_ngsolve)
        file(WRITE "${CMAKE_CURRENT_BINARY_DIR}/__init__.py"
             "import ngsolve  # ensure DLL paths are set up\nfrom .sparsesolv_ngsolve import *\n")
        install(FILES "${CMAKE_CURRENT_BINARY_DIR}/__init__.py"
                DESTINATION sparsesolv_ngsolve)
        install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/sparsesolv_ngsolve.pyi"
                DESTINATION sparsesolv_ngsolve RENAME __init__.pyi)
        install(FILES "${CMAKE_CURRENT_SOURCE_DIR}/py.typed"
                DESTINATION sparsesolv_ngsolve)
    else()
        # Manual build: install next to NGSolve
        install(TARGETS sparsesolv_ngsolve
            LIBRARY DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/sparsesolv_ngsolve
            RUNTIME DESTINATION ${NGSOLVE_INSTALL_DIR_PYTHON}/sparsesolv_ngsolve)
    endif()
endif()

# ==============================================================================
# Header-only library installation
# ==============================================================================
include(GNUInstallDirs)
install(TARGETS sparsesolv
    EXPORT SparseSolvTargets
    INCLUDES DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
)

install(DIRECTORY include/sparsesolv
    DESTINATION ${CMAKE_INSTALL_INCLUDEDIR}
    FILES_MATCHING PATTERN "*.hpp"
    PATTERN "ngsolve" EXCLUDE
)

install(EXPORT SparseSolvTargets
    FILE SparseSolvTargets.cmake
    NAMESPACE SparseSolv::
    DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/SparseSolv
)
