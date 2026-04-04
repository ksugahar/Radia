cmake_minimum_required(VERSION 3.12)
project(RadiaCubitPlugin)

# --------------------------------------------------------
# Cubit SDK
# --------------------------------------------------------
find_package(Cubit REQUIRED CONFIG)

# --------------------------------------------------------
# Netgen (optional, for high-order curving order >= 2)
# --------------------------------------------------------
set(NETGEN_DIR "" CACHE PATH "Netgen installation directory (full Netgen, optional)")
set(NETGEN_SRC_DIR "" CACHE PATH "Netgen source directory (for compact build)")
set(COMPACT_NETGEN_OVERRIDES "" CACHE PATH "Compact Netgen overrides directory")

# Option 1: Compact Netgen (source compiled directly into plugin)
# No separate build, no DLL, no .lib — just source files.
if(NETGEN_SRC_DIR AND EXISTS "${NETGEN_SRC_DIR}/libsrc/meshing/meshing.hpp")
  if(NOT COMPACT_NETGEN_OVERRIDES)
    set(COMPACT_NETGEN_OVERRIDES "C:/compact_netgen_build/overrides")
  endif()
  if(EXISTS "${COMPACT_NETGEN_OVERRIDES}/meshclass_patched.cpp")
    message(STATUS "Compact Netgen: sources from ${NETGEN_SRC_DIR}")
    set(HAVE_NETGEN ON)
    set(USE_COMPACT_NETGEN ON)

    # Minimal Netgen sources (mesh + curving + save only)
    # NOTE: To bisect static init failures, comment out groups below
    # Bisect static init: add sources incrementally.
    # Group A: core only (stubs + safe sources)
    set(COMPACT_NETGEN_SOURCES
      ${COMPACT_NETGEN_OVERRIDES}/core_stubs.cpp
    )
    # Group B: core originals
    list(APPEND COMPACT_NETGEN_SOURCES
      ${NETGEN_SRC_DIR}/libsrc/core/archive.cpp
      ${NETGEN_SRC_DIR}/libsrc/core/bitarray.cpp
      ${NETGEN_SRC_DIR}/libsrc/core/exception.cpp
      ${COMPACT_NETGEN_OVERRIDES}/flags_patched.cpp
      ${NETGEN_SRC_DIR}/libsrc/core/localheap.cpp
      ${NETGEN_SRC_DIR}/libsrc/core/table.cpp
      ${COMPACT_NETGEN_OVERRIDES}/utils_patched.cpp
      ${COMPACT_NETGEN_OVERRIDES}/version_patched.cpp
      ${COMPACT_NETGEN_OVERRIDES}/statushandler_patched.cpp
    )
    # Group C: general + gprim
    list(APPEND COMPACT_NETGEN_SOURCES
      ${NETGEN_SRC_DIR}/libsrc/general/hashtabl.cpp
      ${NETGEN_SRC_DIR}/libsrc/general/mystring.cpp
      ${NETGEN_SRC_DIR}/libsrc/general/table.cpp
      ${NETGEN_SRC_DIR}/libsrc/general/optmem.cpp
      ${NETGEN_SRC_DIR}/libsrc/gprim/geom3d.cpp
      ${NETGEN_SRC_DIR}/libsrc/gprim/transform3d.cpp
      ${NETGEN_SRC_DIR}/libsrc/linalg/densemat.cpp
    )
    # Group D: meshing
    list(APPEND COMPACT_NETGEN_SOURCES
      ${COMPACT_NETGEN_OVERRIDES}/meshtype_patched.cpp
      ${COMPACT_NETGEN_OVERRIDES}/meshclass_patched.cpp
      ${NETGEN_SRC_DIR}/libsrc/meshing/curvedelems.cpp
      ${NETGEN_SRC_DIR}/libsrc/meshing/topology.cpp
      ${NETGEN_SRC_DIR}/libsrc/meshing/localh.cpp
      ${COMPACT_NETGEN_OVERRIDES}/basegeom_patched.cpp
      ${NETGEN_SRC_DIR}/libsrc/meshing/msghandler.cpp
      ${COMPACT_NETGEN_OVERRIDES}/compact_stubs.cpp
    )

    set(COMPACT_NETGEN_INCLUDES
      ${COMPACT_NETGEN_OVERRIDES}/..    # netgen_config.hpp, netgen_version.hpp
      ${COMPACT_NETGEN_OVERRIDES}       # gzstream.h stub, patched headers
      ${NETGEN_SRC_DIR}/libsrc
      ${NETGEN_SRC_DIR}/libsrc/core
      ${NETGEN_SRC_DIR}/libsrc/include
      ${NETGEN_SRC_DIR}/libsrc/general
      ${NETGEN_SRC_DIR}/libsrc/gprim
      ${NETGEN_SRC_DIR}/libsrc/linalg
      ${NETGEN_SRC_DIR}/libsrc/meshing
    )

    set(COMPACT_NETGEN_DEFS
      COMPACT_NETGEN              # Skip mesher-dependent code paths
      COMPACT_NETGEN_STATIC       # DLL_HEADER / NGCORE_API = empty
      NO_PARALLEL_THREADS         # NgMutex = empty class (TaskManager disabled)
      # OpenMP is still used for BuildCurvedElements parallelization
    )
    # Enable OpenMP for edge/face curving parallelization
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
      message(STATUS "Compact Netgen: OpenMP enabled for BuildCurvedElements")
    endif()
  else()
    message(WARNING "NETGEN_SRC_DIR set but overrides not found at ${COMPACT_NETGEN_OVERRIDES}")
    set(HAVE_NETGEN OFF)
    set(USE_COMPACT_NETGEN OFF)
  endif()

# Option 2: Full Netgen (shared libs, has Python dependency)
elseif(NETGEN_DIR AND EXISTS "${NETGEN_DIR}/lib/nglib.lib")
  message(STATUS "Netgen found: ${NETGEN_DIR}")
  set(HAVE_NETGEN ON)
  set(USE_COMPACT_NETGEN OFF)

  add_library(netgen::nglib SHARED IMPORTED)
  set_target_properties(netgen::nglib PROPERTIES
    IMPORTED_IMPLIB "${NETGEN_DIR}/lib/nglib.lib"
    IMPORTED_LOCATION "${NETGEN_DIR}/nglib.dll"
  )

  add_library(netgen::ngcore SHARED IMPORTED)
  set_target_properties(netgen::ngcore PROPERTIES
    IMPORTED_IMPLIB "${NETGEN_DIR}/lib/ngcore.lib"
    IMPORTED_LOCATION "${NETGEN_DIR}/ngcore.dll"
  )
else()
  message(STATUS "Netgen not found. High-order disabled.")
  set(HAVE_NETGEN OFF)
  set(USE_COMPACT_NETGEN OFF)
endif()

# --------------------------------------------------------
# Qt5 (for GUI Component in .ccl)
# --------------------------------------------------------
set(CMAKE_PREFIX_PATH "${CMAKE_PREFIX_PATH};C:/Qt/5.15.2/msvc2019_64")
find_package(Qt5Core REQUIRED)
find_package(Qt5Widgets REQUIRED)

# --------------------------------------------------------
# Target 1: radia_cubit.ccm (command plugin, goes to plugins/)
#
# Registers APREPRO commands: export gmsh/nastran/vtk/meg/netgen
# Loaded by Python API (CUBIT_PLUGIN_DIR) and by cubitcomp
# (-commandplugindir).
# --------------------------------------------------------
set(CCM_SOURCES
    RadiaPlugin.cpp
    RadiaPlugin.hpp
    MeshData.cpp
    MeshData.hpp
    # HighOrderMesh removed — NetgenCurver is the sole curving path (no fallback)
    NetgenCurver.cpp
    NetgenCurver.hpp
    ExportNastranCommand.cpp
    ExportNastranCommand.hpp
    ExportGmshCommand.cpp
    ExportGmshCommand.hpp
    ExportMegCommand.cpp
    ExportMegCommand.hpp
    ExportVtkCommand.cpp
    ExportVtkCommand.hpp
)

# Netgen-dependent sources: callbackgeom + ExportNetgenCommand
if(HAVE_NETGEN)
  list(APPEND CCM_SOURCES
    callbackgeom.cpp callbackgeom.hpp
    ExportNetgenCommand.cpp ExportNetgenCommand.hpp
  )
  if(USE_COMPACT_NETGEN)
    list(APPEND CCM_SOURCES ${COMPACT_NETGEN_SOURCES})
  endif()
endif()

add_library(radia_cubit_ccm MODULE ${CCM_SOURCES})
set_target_properties(radia_cubit_ccm PROPERTIES
    OUTPUT_NAME "radia_cubit"
    SUFFIX ".ccm"
    CXX_STANDARD 17
)
target_link_libraries(radia_cubit_ccm cubiti cubit_util cubit_geom)

# DELAYLOAD cubit_geom so the .ccm can still load even if
# cubit_geom.dll init is deferred until first geometry call
if(MSVC)
  set_target_properties(radia_cubit_ccm PROPERTIES
    LINK_FLAGS "/DELAYLOAD:cubit_geom.dll"
  )
  target_link_libraries(radia_cubit_ccm delayimp)
endif()

if(HAVE_NETGEN)
  target_compile_definitions(radia_cubit_ccm PRIVATE HAVE_NETGEN)
  if(USE_COMPACT_NETGEN)
    target_include_directories(radia_cubit_ccm PRIVATE ${COMPACT_NETGEN_INCLUDES})
    target_compile_definitions(radia_cubit_ccm PRIVATE ${COMPACT_NETGEN_DEFS})
    if(MSVC)
      target_compile_options(radia_cubit_ccm PRIVATE /bigobj)
    endif()
    if(OpenMP_CXX_FOUND)
      target_link_libraries(radia_cubit_ccm OpenMP::OpenMP_CXX)
    endif()
  else()
    target_include_directories(radia_cubit_ccm PRIVATE
      "${NETGEN_DIR}/include"
      "${NETGEN_DIR}/include/include"
    )
    target_link_libraries(radia_cubit_ccm netgen::nglib netgen::ngcore)
  endif()
endif()

# --------------------------------------------------------
# Target 2: radia_cubit.ccl (GUI component + export, goes to bin/)
#
# Claro Component with Export Mesh menu + all export logic.
# Loaded by clarofw GUI scanner on startup. No -commandplugindir needed.
# --------------------------------------------------------
set(CCL_SOURCES
    RadiaComp.cpp
    RadiaComp.hpp
    MeshData.cpp
    MeshData.hpp
    HighOrderMesh.cpp
    HighOrderMesh.hpp
    ExportNetgenCommand.cpp
    ExportNetgenCommand.hpp
    ExportGmshCommand.cpp
    ExportGmshCommand.hpp
    ExportNastranCommand.cpp
    ExportNastranCommand.hpp
    ExportVtkCommand.cpp
    ExportVtkCommand.hpp
    ExportMegCommand.cpp
    ExportMegCommand.hpp
)

if(HAVE_NETGEN)
  list(APPEND CCL_SOURCES
    NetgenCurver.cpp
    NetgenCurver.hpp
    callbackgeom.cpp
    callbackgeom.hpp
  )
  if(USE_COMPACT_NETGEN)
    list(APPEND CCL_SOURCES ${COMPACT_NETGEN_SOURCES})
  endif()
endif()

add_library(radia_cubit_ccl MODULE ${CCL_SOURCES})
set_target_properties(radia_cubit_ccl PROPERTIES
    OUTPUT_NAME "radia_cubit"
    SUFFIX ".ccl"
    CXX_STANDARD 17
    AUTOMOC ON
)
target_link_libraries(radia_cubit_ccl cubiti cubit_util cubit_geom clarofw clarogui Qt5::Core Qt5::Widgets)

if(MSVC)
  set_target_properties(radia_cubit_ccl PROPERTIES
    LINK_FLAGS "/DELAYLOAD:cubit_geom.dll"
  )
  target_link_libraries(radia_cubit_ccl delayimp)
endif()

if(HAVE_NETGEN)
  target_compile_definitions(radia_cubit_ccl PRIVATE HAVE_NETGEN)
  if(USE_COMPACT_NETGEN)
    target_include_directories(radia_cubit_ccl PRIVATE ${COMPACT_NETGEN_INCLUDES})
    target_compile_definitions(radia_cubit_ccl PRIVATE ${COMPACT_NETGEN_DEFS})
    if(MSVC)
      target_compile_options(radia_cubit_ccl PRIVATE /bigobj)
    endif()
  else()
    target_include_directories(radia_cubit_ccl PRIVATE
      "${NETGEN_DIR}/include"
      "${NETGEN_DIR}/include/include"
    )
    target_link_libraries(radia_cubit_ccl netgen::nglib netgen::ngcore)
    # DELAYLOAD Netgen DLLs — they live in plugins/, not bin/ where .ccl is
    if(MSVC)
      set_property(TARGET radia_cubit_ccl APPEND_STRING PROPERTY
        LINK_FLAGS " /DELAYLOAD:nglib.dll /DELAYLOAD:ngcore.dll")
    endif()
  endif()
endif()

# --------------------------------------------------------
# pybind11 module: radia_cubit_mesh (.pyd) — Cubit-free
# --------------------------------------------------------
find_package(Python3 COMPONENTS Interpreter Development)
execute_process(
  COMMAND ${Python3_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_path('purelib'))"
  OUTPUT_VARIABLE PYTHON_SITE_PACKAGES OUTPUT_STRIP_TRAILING_WHITESPACE
)
set(pybind11_DIR "${PYTHON_SITE_PACKAGES}/pybind11/share/cmake/pybind11")
find_package(pybind11 QUIET)

if(pybind11_FOUND AND HAVE_NETGEN AND NOT USE_COMPACT_NETGEN)
  message(STATUS "pybind11 + Netgen found: building radia_cubit_mesh (Cubit-free)")

  pybind11_add_module(radia_cubit_mesh
      radia_cubit_pybind.cpp
      NetgenCurverPure.cpp
      NetgenCurverPure.hpp
      callbackgeom.cpp
      callbackgeom.hpp
  )
  set_target_properties(radia_cubit_mesh PROPERTIES CXX_STANDARD 17)

  target_include_directories(radia_cubit_mesh PRIVATE
    "${NETGEN_DIR}/include"
    "${NETGEN_DIR}/include/include"
    "${CMAKE_CURRENT_SOURCE_DIR}"
  )
  target_link_libraries(radia_cubit_mesh PRIVATE netgen::nglib netgen::ngcore)

  install(TARGETS radia_cubit_mesh
      LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/plugins"
  )
elseif(pybind11_FOUND AND NOT HAVE_NETGEN)
  message(STATUS "pybind11 found but Netgen missing: skipping radia_cubit_mesh")
else()
  message(STATUS "pybind11 not found: skipping radia_cubit_mesh Python module")
endif()

# --------------------------------------------------------
# Install
# --------------------------------------------------------
install(TARGETS radia_cubit_ccm
    LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/plugins"
)
install(TARGETS radia_cubit_ccl
    LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/bin"
)

if(HAVE_NETGEN AND NOT USE_COMPACT_NETGEN)
  install(FILES
    "${NETGEN_DIR}/nglib.dll"
    "${NETGEN_DIR}/ngcore.dll"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/plugins"
  )
endif()
