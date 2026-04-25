cmake_minimum_required(VERSION 3.12)
project(RadiaCubitPlugin)

# --------------------------------------------------------
# Emit PDBs for Release builds so crash dumps can be symbolicated.
# /Zi: compiler-side debug info (PDB)
# /DEBUG:FULL: linker emits full PDB (not stripped)
# Optimization (/O2) is unchanged — this is like RelWithDebInfo.
# --------------------------------------------------------
if(MSVC)
  add_compile_options($<$<CONFIG:Release>:/Zi>)
  add_link_options($<$<CONFIG:Release>:/DEBUG:FULL>)
endif()

# --------------------------------------------------------
# Cubit SDK
# --------------------------------------------------------
find_package(Cubit REQUIRED CONFIG)

# --------------------------------------------------------
# Netgen (optional, for high-order curving order >= 2)
# --------------------------------------------------------
set(NETGEN_DIR "" CACHE PATH "Netgen installation directory (full Netgen, optional)")
set(COMPACT_NETGEN_OVERRIDES "" CACHE PATH "Compact Netgen overrides directory")

# Option 1: Compact Netgen (source compiled directly into plugin)
# Netgen sources are snapshotted in compact_netgen/netgen_src/ (no external dependency).
# To update: copy libsrc/ files from upstream netgen, apply COMPACT_NETGEN_STATIC patch
# to ngcore_api.hpp and mydefs.hpp.
if(NOT COMPACT_NETGEN_OVERRIDES)
  set(COMPACT_NETGEN_OVERRIDES "${CMAKE_CURRENT_SOURCE_DIR}/compact_netgen")
endif()
set(NETGEN_SNAP "${COMPACT_NETGEN_OVERRIDES}/netgen_src")
if(EXISTS "${NETGEN_SNAP}/meshing/meshing.hpp" AND EXISTS "${COMPACT_NETGEN_OVERRIDES}/meshclass_patched.cpp")
    message(STATUS "Compact Netgen: sources from ${NETGEN_SNAP}")
    set(HAVE_NETGEN ON)
    set(USE_COMPACT_NETGEN ON)

    set(COMPACT_NETGEN_SOURCES
      ${COMPACT_NETGEN_OVERRIDES}/core_stubs.cpp
    )
    list(APPEND COMPACT_NETGEN_SOURCES
      ${NETGEN_SNAP}/core/archive.cpp
      ${NETGEN_SNAP}/core/bitarray.cpp
      ${NETGEN_SNAP}/core/exception.cpp
      ${COMPACT_NETGEN_OVERRIDES}/flags_patched.cpp
      ${NETGEN_SNAP}/core/localheap.cpp
      ${NETGEN_SNAP}/core/table.cpp
      ${COMPACT_NETGEN_OVERRIDES}/utils_patched.cpp
      ${COMPACT_NETGEN_OVERRIDES}/version_patched.cpp
      ${COMPACT_NETGEN_OVERRIDES}/statushandler_patched.cpp
    )
    list(APPEND COMPACT_NETGEN_SOURCES
      ${NETGEN_SNAP}/general/hashtabl.cpp
      ${NETGEN_SNAP}/general/mystring.cpp
      ${NETGEN_SNAP}/general/table.cpp
      ${NETGEN_SNAP}/general/optmem.cpp
      ${NETGEN_SNAP}/gprim/geom3d.cpp
      ${NETGEN_SNAP}/gprim/transform3d.cpp
      ${NETGEN_SNAP}/linalg/densemat.cpp
    )
    list(APPEND COMPACT_NETGEN_SOURCES
      ${COMPACT_NETGEN_OVERRIDES}/meshtype_patched.cpp
      ${COMPACT_NETGEN_OVERRIDES}/meshclass_patched.cpp
      ${NETGEN_SNAP}/meshing/curvedelems.cpp
      ${NETGEN_SNAP}/meshing/topology.cpp
      ${NETGEN_SNAP}/meshing/localh.cpp
      ${COMPACT_NETGEN_OVERRIDES}/basegeom_patched.cpp
      ${NETGEN_SNAP}/meshing/msghandler.cpp
      ${COMPACT_NETGEN_OVERRIDES}/compact_stubs.cpp
    )

    set(COMPACT_NETGEN_INCLUDES
      ${COMPACT_NETGEN_OVERRIDES}/..    # netgen_config.hpp, netgen_version.hpp
      ${COMPACT_NETGEN_OVERRIDES}       # gzstream.h stub, patched headers
      ${NETGEN_SNAP}
      ${NETGEN_SNAP}/core
      ${NETGEN_SNAP}/include
      ${NETGEN_SNAP}/general
      ${NETGEN_SNAP}/gprim
      ${NETGEN_SNAP}/linalg
      ${NETGEN_SNAP}/meshing
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
# Use minimal Qt5 SDK bundled in repo (headers + .lib only).
# Runtime DLLs come from the Cubit installation.
# This avoids requiring a full Qt5 installation or C:\Qt.
set(QT5_SDK_DIR "${CMAKE_CURRENT_SOURCE_DIR}/qt5_sdk")

# Qt5::Core
add_library(Qt5::Core SHARED IMPORTED)
set_target_properties(Qt5::Core PROPERTIES
    IMPORTED_LOCATION "${QT5_SDK_DIR}/bin/Qt5Core.dll"
    IMPORTED_IMPLIB   "${QT5_SDK_DIR}/lib/Qt5Core.lib"
    INTERFACE_INCLUDE_DIRECTORIES "${QT5_SDK_DIR}/include;${QT5_SDK_DIR}/include/QtCore"
)
# Qt5::Gui (transitive dependency of Widgets)
add_library(Qt5::Gui SHARED IMPORTED)
set_target_properties(Qt5::Gui PROPERTIES
    IMPORTED_LOCATION "${QT5_SDK_DIR}/bin/Qt5Gui.dll"
    IMPORTED_IMPLIB   "${QT5_SDK_DIR}/lib/Qt5Gui.lib"
    INTERFACE_INCLUDE_DIRECTORIES "${QT5_SDK_DIR}/include;${QT5_SDK_DIR}/include/QtGui"
)
# Qt5::Widgets
add_library(Qt5::Widgets SHARED IMPORTED)
set_target_properties(Qt5::Widgets PROPERTIES
    IMPORTED_LOCATION "${QT5_SDK_DIR}/bin/Qt5Widgets.dll"
    IMPORTED_IMPLIB   "${QT5_SDK_DIR}/lib/Qt5Widgets.lib"
    INTERFACE_INCLUDE_DIRECTORIES "${QT5_SDK_DIR}/include;${QT5_SDK_DIR}/include/QtWidgets"
    INTERFACE_LINK_LIBRARIES "Qt5::Core;Qt5::Gui"
)
# moc for AUTOMOC
set(CMAKE_AUTOMOC ON)
set(QT_VERSION_MAJOR 5)
set(Qt5Core_VERSION_MAJOR 5)
set(Qt5Core_VERSION_MINOR 15)
set(Qt5Core_VERSION "5.15.2")
set(QT_MOC_EXECUTABLE "${QT5_SDK_DIR}/bin/moc.exe")
if(NOT TARGET Qt5::moc)
    add_executable(Qt5::moc IMPORTED)
    set_target_properties(Qt5::moc PROPERTIES IMPORTED_LOCATION "${QT5_SDK_DIR}/bin/moc.exe")
endif()

# --------------------------------------------------------
# Target 1: radia_cubit.ccm (command plugin, goes to plugins/)
#
# Registers APREPRO commands: export gmsh/nastran/vtk/netgen
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
    ExportVtkCommand.cpp
    ExportVtkCommand.hpp
    ExportFemeemCommand.cpp
    ExportFemeemCommand.hpp
    ExportMegCommand.cpp
    ExportMegCommand.hpp
    CoilCommand.cpp
    CoilCommand.hpp
    VerifyLauncherCommand.cpp
    VerifyLauncherCommand.hpp
    LauncherLogic.hpp
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
    # Separate PDB per target so .ccm and .ccl don't clobber each other.
    PDB_NAME "radia_cubit_ccm"
    COMPILE_PDB_NAME "radia_cubit_ccm"
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

option(DEBUG_HO_NODES "Print debug info for high-order node registration" OFF)

if(HAVE_NETGEN)
  target_compile_definitions(radia_cubit_ccm PRIVATE HAVE_NETGEN)
  if(DEBUG_HO_NODES)
    target_compile_definitions(radia_cubit_ccm PRIVATE DEBUG_HO_NODES)
  endif()
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
    ExportFemeemCommand.cpp
    ExportFemeemCommand.hpp
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
    # Separate PDB per target so .ccm and .ccl don't clobber each other.
    PDB_NAME "radia_cubit_ccl"
    COMPILE_PDB_NAME "radia_cubit_ccl"
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
    # No DELAYLOAD for nglib/ngcore: data symbols (netgen::Mesh::default_bc)
    # cannot be delay-loaded. ensure_netgen_dll_path() sets SetDllDirectoryA
    # to plugins/ at startup, so the DLLs are found via normal load.
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
