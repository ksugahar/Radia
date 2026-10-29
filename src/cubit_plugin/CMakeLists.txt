cmake_minimum_required(VERSION 3.12)
project(CubitMeshExportPlugin)

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
if(NOT TARGET cubit_geom)
  # Cubit 2026.8's CGMPlugin SDK example uses the split CGM core target.
  find_package(CubitCGM REQUIRED CONFIG HINTS "${Cubit_DIR}")
  add_library(cubit_geom ALIAS cubit_cgmcore)
endif()

# --------------------------------------------------------
# Netgen (optional, for high-order curving order >= 2)
# --------------------------------------------------------
set(NETGEN_DIR "" CACHE PATH "Netgen installation directory (full Netgen, optional)")
set(COMPACT_NETGEN_OVERRIDES "" CACHE PATH "Compact Netgen overrides directory")

# Option 1: Compact Netgen (source compiled directly into plugin)
# Netgen sources are snapshotted in compact_netgen/netgen_src/ (no external dependency).
# Snapshot: upstream v6.2.2607 plus local patches (static-library macros in
# ngcore_api.hpp/mydefs.hpp, BlockAllocator mutex lifetime, CurvedElements
# coefficient accessors, pyramid base-edge blending in EvaluateMapping).
# To update, run compact_netgen/update_snapshot.py (three-way merge).
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
      CUBIT_MESH_EXPORT_NETGEN_EDGE_DESCRIPTOR_API  # snapshot >= 6.2.2605
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
# Retired Qt5 target removed.
# --------------------------------------------------------
# The old Qt5 .ccl GUI component is gone; all GUI work is now in
# packages/cubit-mesh-export/src/cubit_mesh_export/cubit_gui/cubit_export_menu.py
# (PySide6).  The .ccm (APREPRO commands) and .pyd (cubit_mesh_curver)
# below have no Qt dependency.

# --------------------------------------------------------
# Target 1: cubit_mesh_export.ccm (command plugin, goes to plugins/)
#
# Registers APREPRO commands: export gmsh/nastran/vtk/netgen
# Loaded by Python API (CUBIT_PLUGIN_DIR) and by cubitcomp
# (-commandplugindir).
# --------------------------------------------------------
set(CCM_SOURCES
    CubitMeshExportPlugin.cpp
    CubitMeshExportPlugin.hpp
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
)

# Netgen-dependent sources: callbackgeom + ExportNetgenCommand
if(HAVE_NETGEN)
  list(APPEND CCM_SOURCES
    callbackgeom.cpp callbackgeom.hpp
    NetgenGeometricRefit.cpp NetgenGeometricRefit.hpp
    ExportNetgenCommand.cpp ExportNetgenCommand.hpp NamedPeriodicIdentification.hpp NamedPeriodicMeshAdapter.hpp
  )
  if(USE_COMPACT_NETGEN)
    list(APPEND CCM_SOURCES ${COMPACT_NETGEN_SOURCES})
  endif()
endif()

add_library(cubit_mesh_export_ccm MODULE ${CCM_SOURCES})
set_target_properties(cubit_mesh_export_ccm PROPERTIES
    OUTPUT_NAME "cubit_mesh_export"
    SUFFIX ".ccm"
    CXX_STANDARD 17
    # Separate PDB per target so .ccm and .ccl don't clobber each other.
    PDB_NAME "cubit_mesh_export_ccm"
    COMPILE_PDB_NAME "cubit_mesh_export_ccm"
)
target_link_libraries(cubit_mesh_export_ccm cubiti cubit_util cubit_geom)
# Cubit 2025.12 loads its bundled MSVCP140.dll 14.31.  Since MSVC 14.40 the
# std::mutex constructor is constexpr, and a mutex built that way crashes in
# lock() against an older runtime.  Keep the pre-14.40 constructor.
target_compile_definitions(cubit_mesh_export_ccm PRIVATE _DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR)

# DELAYLOAD cubit_geom so the .ccm can still load even if
# cubit_geom.dll init is deferred until first geometry call
if(MSVC)
  # CodeView records contain a filename; full local PDBs remain available.
  target_link_options(cubit_mesh_export_ccm PRIVATE
    "$<$<CONFIG:Release>:/PDBALTPATH:$<TARGET_PDB_FILE_NAME:cubit_mesh_export_ccm>>")
  target_link_options(cubit_mesh_export_ccm PRIVATE
    "/DELAYLOAD:$<TARGET_FILE_NAME:cubit_geom>")
  target_link_libraries(cubit_mesh_export_ccm delayimp)
endif()

option(DEBUG_HO_NODES "Print debug info for high-order node registration" OFF)

if(HAVE_NETGEN)
  target_compile_definitions(cubit_mesh_export_ccm PRIVATE HAVE_NETGEN)
  if(DEBUG_HO_NODES)
    target_compile_definitions(cubit_mesh_export_ccm PRIVATE DEBUG_HO_NODES)
  endif()
  if(USE_COMPACT_NETGEN)
    target_include_directories(cubit_mesh_export_ccm PRIVATE ${COMPACT_NETGEN_INCLUDES})
    target_compile_definitions(cubit_mesh_export_ccm PRIVATE ${COMPACT_NETGEN_DEFS})
    if(MSVC)
      target_compile_options(cubit_mesh_export_ccm PRIVATE /bigobj)
    endif()
    if(OpenMP_CXX_FOUND)
      target_link_libraries(cubit_mesh_export_ccm OpenMP::OpenMP_CXX)
    endif()
  else()
    target_include_directories(cubit_mesh_export_ccm PRIVATE
      "${NETGEN_DIR}/include"
      "${NETGEN_DIR}/include/include"
    )
    target_link_libraries(cubit_mesh_export_ccm netgen::nglib netgen::ngcore)
  endif()
endif()

# --------------------------------------------------------
# The old cubit_mesh_export.ccl target is intentionally absent.
# --------------------------------------------------------
# The Qt5 GUI Claro Component is replaced by the PySide6 toolbar at
# packages/cubit-mesh-export/src/cubit_mesh_export/cubit_gui/cubit_export_menu.py.

# --------------------------------------------------------
# pybind11 module: cubit_mesh_curver (.pyd) — Cubit-free
# --------------------------------------------------------
find_package(Python3 COMPONENTS Interpreter Development)
# Both compact-writer and full-Netgen contracts load the same installed reader.
# Discover package paths without executing netgen.__init__ or native CAD code.
function(cubit_contract_runtime_path output)
  execute_process(COMMAND "${Python3_EXECUTABLE}" -c
    "import importlib.metadata as m, importlib.util as u, pathlib, sys; spec=u.find_spec('netgen'); assert spec and spec.origin; ng=pathlib.Path(spec.origin).resolve().parent; roots={str(pathlib.Path(f.locate()).resolve().parent) for f in m.files('netgen-occt') if f.name.lower().endswith('.dll')}; roots.add(sys.base_prefix); print(';'.join([str(ng),str(ng.parent/'netgen_mesher.libs')]+sorted(roots)))"
    OUTPUT_VARIABLE runtime_dirs OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE runtime_status)
  if(NOT runtime_status EQUAL 0 OR NOT runtime_dirs)
    message(FATAL_ERROR "Cannot discover Netgen/OCC contract runtime DLL directories")
  endif()
  if(NETGEN_DIR)
    list(GET runtime_dirs 0 discovered_netgen)
    get_filename_component(discovered_netgen "${discovered_netgen}" REALPATH)
    get_filename_component(configured_netgen "${NETGEN_DIR}" REALPATH)
    if(NOT discovered_netgen STREQUAL configured_netgen)
      message(FATAL_ERROR "Configured Netgen and Python reader package paths differ")
    endif()
  endif()
  set(runtime_path "${runtime_dirs};$ENV{PATH}")
  string(REPLACE ";" "\\;" runtime_path "${runtime_path}")
  set(${output} "${runtime_path}" PARENT_SCOPE)
endfunction()
option(CUBIT_PERIODIC_CURVED_CONTRACT_TESTS "Build compact periodic curved geometry contracts" OFF)
if(CUBIT_PERIODIC_CURVED_CONTRACT_TESTS)
  if(NOT USE_COMPACT_NETGEN)
    message(FATAL_ERROR "Curved coefficient contract requires the compact build variant")
  endif()
  enable_testing()
  add_executable(cubit_periodic_curved_contract
    tests/periodic_curved_contract.cpp ${COMPACT_NETGEN_SOURCES})
  target_compile_features(cubit_periodic_curved_contract PRIVATE cxx_std_17)
  target_include_directories(cubit_periodic_curved_contract PRIVATE
    "${CMAKE_CURRENT_SOURCE_DIR}" ${COMPACT_NETGEN_INCLUDES})
  target_compile_definitions(cubit_periodic_curved_contract PRIVATE ${COMPACT_NETGEN_DEFS})
  if(MSVC)
    target_compile_options(cubit_periodic_curved_contract PRIVATE /bigobj)
  endif()
  if(OpenMP_CXX_FOUND)
    target_link_libraries(cubit_periodic_curved_contract PRIVATE OpenMP::OpenMP_CXX)
  endif()
  add_test(NAME cubit_periodic_curved_coefficients COMMAND cubit_periodic_curved_contract)
  set_tests_properties(cubit_periodic_curved_coefficients PROPERTIES
    FIXTURES_SETUP periodic_curved_fixture TIMEOUT 120 FAIL_REGULAR_EXPRESSION "[Ii]llegal face")
  add_test(NAME cubit_periodic_curved_independent_load COMMAND "${Python3_EXECUTABLE}"
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/verify_periodic_curved.py" "${CMAKE_CURRENT_BINARY_DIR}")
  set_tests_properties(cubit_periodic_curved_independent_load PROPERTIES
    FIXTURES_REQUIRED periodic_curved_fixture TIMEOUT 120)
  cubit_contract_runtime_path(COMPACT_READER_RUNTIME_PATH)
  set_tests_properties(cubit_periodic_curved_independent_load PROPERTIES
    ENVIRONMENT "PATH=${COMPACT_READER_RUNTIME_PATH};PYTHONUNBUFFERED=1;PYTHONDONTWRITEBYTECODE=1")
endif()
# Honor the caller's explicit config directory; Python purelib may differ from
# the importable dependency location in an isolated environment.
find_package(pybind11 CONFIG REQUIRED)

if(pybind11_FOUND AND HAVE_NETGEN AND NOT USE_COMPACT_NETGEN)
  message(STATUS "pybind11 + Netgen found: building cubit_mesh_curver (Cubit-free)")

  pybind11_add_module(cubit_mesh_curver
      cubit_mesh_export_pybind.cpp
      NetgenCurverPure.cpp
      NetgenCurverPure.hpp
      callbackgeom.cpp
      callbackgeom.hpp
  )
  set_target_properties(cubit_mesh_curver PROPERTIES CXX_STANDARD 17)
  if(MSVC)
    target_link_options(cubit_mesh_curver PRIVATE
      "$<$<CONFIG:Release>:/PDBALTPATH:$<TARGET_PDB_FILE_NAME:cubit_mesh_curver>>")
  endif()

  target_include_directories(cubit_mesh_curver PRIVATE
    "${NETGEN_DIR}/include"
    "${NETGEN_DIR}/include/include"
    "${CMAKE_CURRENT_SOURCE_DIR}"
  )
  target_compile_definitions(cubit_mesh_curver PRIVATE CUBIT_MESH_EXPORT_NETGEN_EDGE_DESCRIPTOR_API)
  target_link_libraries(cubit_mesh_curver PRIVATE netgen::nglib netgen::ngcore)

  install(TARGETS cubit_mesh_curver
      LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/plugins"
  )
  option(CUBIT_CALLBACK_CONTRACT_TESTS "Build native geometry callback contracts" OFF)
  if(CUBIT_CALLBACK_CONTRACT_TESTS)
    enable_testing()
    add_library(cubit_named_periodic_header_compile OBJECT tests/named_periodic_header_compile.cpp)
    target_compile_features(cubit_named_periodic_header_compile PRIVATE cxx_std_17)
    target_include_directories(cubit_named_periodic_header_compile PRIVATE
      "${CMAKE_CURRENT_SOURCE_DIR}" "${NETGEN_DIR}/include" "${NETGEN_DIR}/include/include")
    add_executable(cubit_callback_contract tests/callback_curve_identity.cpp callbackgeom.cpp NetgenCurver.cpp)
    add_dependencies(cubit_callback_contract cubit_named_periodic_header_compile)
    target_compile_features(cubit_callback_contract PRIVATE cxx_std_17)
    target_include_directories(cubit_callback_contract PRIVATE
      "${CMAKE_CURRENT_SOURCE_DIR}" "${NETGEN_DIR}/include" "${NETGEN_DIR}/include/include")
    target_compile_definitions(cubit_callback_contract PRIVATE HAVE_NETGEN CUBIT_MESH_EXPORT_MAPPING_CONTRACT_ONLY)
    target_link_libraries(cubit_callback_contract PRIVATE netgen::nglib netgen::ngcore)
    add_test(NAME cubit_callback_curve_identity COMMAND cubit_callback_contract)
    set_tests_properties(cubit_callback_curve_identity PROPERTIES FIXTURES_SETUP named_periodic_fixture)
    add_test(NAME cubit_named_periodic_vol_load COMMAND "${Python3_EXECUTABLE}"
      "${CMAKE_CURRENT_SOURCE_DIR}/tests/verify_named_periodic_vol.py"
      "${CMAKE_CURRENT_BINARY_DIR}/named-periodic-contract.vol")
    set_tests_properties(cubit_named_periodic_vol_load PROPERTIES FIXTURES_REQUIRED named_periodic_fixture)
    # Native executables do not run netgen.__init__'s OCC/delvewheel preload.
    # Discover the exact installed package DLL locations, without importing CAD.
    cubit_contract_runtime_path(CALLBACK_RUNTIME_PATH)
    set_tests_properties(cubit_callback_curve_identity cubit_named_periodic_vol_load PROPERTIES
      ENVIRONMENT "PATH=${CALLBACK_RUNTIME_PATH}" TIMEOUT 120)
  endif()
elseif(pybind11_FOUND AND NOT HAVE_NETGEN)
  message(STATUS "pybind11 found but Netgen missing: skipping cubit_mesh_curver")
else()
  message(STATUS "pybind11 not found: skipping cubit_mesh_curver Python module")
endif()

# --------------------------------------------------------
# Install
# --------------------------------------------------------
install(TARGETS cubit_mesh_export_ccm
    LIBRARY DESTINATION "${CMAKE_INSTALL_PREFIX}/plugins"
)
# No .ccl is installed; the packaged PySide6 toolbar is the sole GUI.

if(HAVE_NETGEN AND NOT USE_COMPACT_NETGEN)
  install(FILES
    "${NETGEN_DIR}/nglib.dll"
    "${NETGEN_DIR}/ngcore.dll"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/plugins"
  )
endif()
