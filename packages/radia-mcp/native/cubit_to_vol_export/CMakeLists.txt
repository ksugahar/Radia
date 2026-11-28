cmake_minimum_required(VERSION 3.20)
project(cubit_to_vol_export CXX)
set(CMAKE_CXX_STANDARD 17)
# CUBIT_TO_VOL_EXPORT_PLUGIN builds the Cubit command plugin (needs the Cubit SDK).
# NETGEN_CURVING_PYMODULE builds _netgen_curving, the Python test module for
# netgen_curving/ (needs only the Netgen SDK, no Cubit).
option(CUBIT_TO_VOL_EXPORT_PLUGIN "Build cubit_to_vol_export.dll" ON)
option(NETGEN_CURVING_PYMODULE "Build the _netgen_curving test module" OFF)
# CUBIT_TO_VOL_EXPORT_CONTRACT_TESTS builds the license-free native contracts in
# tests/ (Netgen SDK and Python with netgen/ngsolve only, no Cubit) and
# registers them with CTest.
option(CUBIT_TO_VOL_EXPORT_CONTRACT_TESTS "Build the native contract tests" OFF)
find_package(Netgen REQUIRED CONFIG)

add_library(netgen_curving STATIC netgen_curving/geometric_refit.cpp netgen_curving/periodic_copy.cpp)
set_target_properties(netgen_curving PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(netgen_curving PUBLIC nglib ngcore)
if(MSVC)
  target_compile_options(netgen_curving PRIVATE /bigobj /utf-8)
  target_compile_definitions(netgen_curving PUBLIC _DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR)
endif()

if(CUBIT_TO_VOL_EXPORT_PLUGIN)
  find_package(Cubit REQUIRED CONFIG)
  if(NOT TARGET cubit_geom)
    find_package(CubitCGM REQUIRED CONFIG HINTS "${Cubit_DIR}")
    add_library(cubit_geom ALIAS cubit_cgmcore)
  endif()
  # The export itself, shared by the command plugin and the in-memory mesh module.
  add_library(cubit_to_vol_export_core OBJECT
    MeshData.cpp NetgenCurver.cpp NetgenCurverIdentity.cpp callbackgeom.cpp
    ExportVolCommand.cpp LiveVolMesh.cpp)
  target_link_libraries(cubit_to_vol_export_core PUBLIC cubiti cubit_util cubit_geom netgen_curving nglib ngcore)
  target_compile_definitions(cubit_to_vol_export_core PUBLIC _DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR)
  target_compile_options(cubit_to_vol_export_core PUBLIC /bigobj /utf-8)
  add_library(cubit_to_vol_export_plugin MODULE CubitToVolExportPlugin.cpp)
  set_target_properties(cubit_to_vol_export_plugin PROPERTIES OUTPUT_NAME "cubit_to_vol_export" PREFIX "" SUFFIX ".dll")
  target_link_libraries(cubit_to_vol_export_plugin PRIVATE cubit_to_vol_export_core)
  # _cubit_to_vol_mesh: the same export kept in memory as a netgen.meshing.Mesh,
  # for an external Python that imported cubit.  Uses the pybind11 headers the
  # Netgen SDK ships, as netgen.libngpy does, so the Mesh casts across modules.
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
  Python3_add_library(_cubit_to_vol_mesh MODULE WITH_SOABI python_live_mesh.cpp)
  target_link_libraries(_cubit_to_vol_mesh PRIVATE cubit_to_vol_export_core)
endif()

if(NETGEN_CURVING_PYMODULE)
  # netgen.libngpy is built with the pybind11 headers the SDK ships; use the
  # same ones so a netgen.meshing.Mesh argument casts across modules.
  find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
  Python3_add_library(_netgen_curving MODULE WITH_SOABI netgen_curving/python_module.cpp)
  target_link_libraries(_netgen_curving PRIVATE netgen_curving nglib ngcore)
  if(MSVC)
    target_compile_options(_netgen_curving PRIVATE /bigobj /utf-8)
  endif()
endif()

if(CUBIT_TO_VOL_EXPORT_CONTRACT_TESTS)
  # The verifiers read the written .vol with the Netgen/NGSolve that this build
  # links; a missing interpreter or package is a configure error, not a skip.
  find_package(Python3 REQUIRED COMPONENTS Interpreter)
  execute_process(COMMAND "${Python3_EXECUTABLE}" -I -c
    "import importlib.metadata as m, importlib.util as u, pathlib, sys; spec=u.find_spec('netgen'); assert spec and spec.origin, 'netgen is not importable'; assert u.find_spec('ngsolve'), 'ngsolve is not importable'; ng=pathlib.Path(spec.origin).resolve().parent; roots={str(pathlib.Path(f.locate()).resolve().parent) for f in (m.files('netgen-occt') or []) if f.name.lower().endswith('.dll')}; roots.add(sys.base_prefix); print(';'.join([str(ng), str(ng.parent/'netgen_mesher.libs')]+sorted(roots)))"
    OUTPUT_VARIABLE contract_runtime_dirs OUTPUT_STRIP_TRAILING_WHITESPACE
    RESULT_VARIABLE contract_runtime_status)
  if(NOT contract_runtime_status EQUAL 0 OR contract_runtime_dirs STREQUAL "")
    message(FATAL_ERROR "Cannot find the Netgen/NGSolve runtime of ${Python3_EXECUTABLE}")
  endif()
  set(contract_path "${contract_runtime_dirs};$ENV{PATH}")
  string(REPLACE ";" "\\;" contract_path "${contract_path}")
  enable_testing()
  # The contracts refuse to overwrite their outputs, so every CTest run starts
  # from an empty output directory.
  set(contract_out "${CMAKE_CURRENT_BINARY_DIR}/contract-output")
  add_test(NAME cubit_to_vol_export_contract_output_clean
    COMMAND "${CMAKE_COMMAND}" -E rm -rf "${contract_out}")
  add_test(NAME cubit_to_vol_export_contract_output_create
    COMMAND "${CMAKE_COMMAND}" -E make_directory "${contract_out}")
  set_tests_properties(cubit_to_vol_export_contract_output_clean PROPERTIES
    FIXTURES_SETUP contract_output_clean)
  set_tests_properties(cubit_to_vol_export_contract_output_create PROPERTIES
    FIXTURES_REQUIRED contract_output_clean FIXTURES_SETUP contract_output)

  function(cubit_to_vol_export_contract_target target)
    target_compile_features(${target} PRIVATE cxx_std_17)
    target_include_directories(${target} PRIVATE "${CMAKE_CURRENT_SOURCE_DIR}")
    target_link_libraries(${target} PRIVATE netgen_curving nglib ngcore)
    if(MSVC)
      target_compile_options(${target} PRIVATE /bigobj /utf-8)
    endif()
  endfunction()

  # The production matcher must compile on its own before any contract links.
  add_library(cubit_to_vol_export_named_periodic_header_compile OBJECT
    tests/named_periodic_header_compile.cpp)
  cubit_to_vol_export_contract_target(cubit_to_vol_export_named_periodic_header_compile)

  # Same production mapping and callback geometry as cubit_to_vol_export.dll.
  add_executable(cubit_to_vol_export_callback_contract
    tests/callback_curve_identity.cpp callbackgeom.cpp NetgenCurverIdentity.cpp)
  add_dependencies(cubit_to_vol_export_callback_contract cubit_to_vol_export_named_periodic_header_compile)
  cubit_to_vol_export_contract_target(cubit_to_vol_export_callback_contract)
  add_test(NAME cubit_to_vol_export_callback_curve_identity COMMAND cubit_to_vol_export_callback_contract)
  set_tests_properties(cubit_to_vol_export_callback_curve_identity PROPERTIES
    FIXTURES_SETUP named_periodic_vol FIXTURES_REQUIRED contract_output)
  add_test(NAME cubit_to_vol_export_named_periodic_vol_load COMMAND "${Python3_EXECUTABLE}" -I
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/verify_named_periodic_vol.py"
    "${contract_out}/named-periodic-contract.vol")
  set_tests_properties(cubit_to_vol_export_named_periodic_vol_load PROPERTIES
    FIXTURES_REQUIRED "named_periodic_vol;contract_output")

  add_executable(cubit_to_vol_export_periodic_curved_contract tests/periodic_curved_contract.cpp)
  cubit_to_vol_export_contract_target(cubit_to_vol_export_periodic_curved_contract)
  add_test(NAME cubit_to_vol_export_periodic_curved_coefficients COMMAND cubit_to_vol_export_periodic_curved_contract)
  set_tests_properties(cubit_to_vol_export_periodic_curved_coefficients PROPERTIES
    FIXTURES_SETUP periodic_curved_vol FIXTURES_REQUIRED contract_output
    FAIL_REGULAR_EXPRESSION "[Ii]llegal face")
  add_test(NAME cubit_to_vol_export_periodic_curved_vol_load COMMAND "${Python3_EXECUTABLE}" -I
    "${CMAKE_CURRENT_SOURCE_DIR}/tests/verify_periodic_curved.py" "${contract_out}")
  set_tests_properties(cubit_to_vol_export_periodic_curved_vol_load PROPERTIES
    FIXTURES_REQUIRED "periodic_curved_vol;contract_output")

  set_tests_properties(
    cubit_to_vol_export_callback_curve_identity cubit_to_vol_export_named_periodic_vol_load
    cubit_to_vol_export_periodic_curved_coefficients cubit_to_vol_export_periodic_curved_vol_load
    PROPERTIES TIMEOUT 120 WORKING_DIRECTORY "${contract_out}"
    ENVIRONMENT "PATH=${contract_path};PYTHONUNBUFFERED=1;PYTHONDONTWRITEBYTECODE=1")
endif()
