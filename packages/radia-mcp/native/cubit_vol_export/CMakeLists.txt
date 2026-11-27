cmake_minimum_required(VERSION 3.20)
project(radia_cubit_export CXX)
set(CMAKE_CXX_STANDARD 17)
# CUBIT_VOL_EXPORT_PLUGIN builds the Cubit command plugin (needs the Cubit SDK).
# NETGEN_CURVING_PYMODULE builds _netgen_curving, the Python test module for
# netgen_curving/ (needs only the Netgen SDK, no Cubit).
option(CUBIT_VOL_EXPORT_PLUGIN "Build radia_cubit_export.dll" ON)
option(NETGEN_CURVING_PYMODULE "Build the _netgen_curving test module" OFF)
find_package(Netgen REQUIRED CONFIG)

add_library(netgen_curving STATIC netgen_curving/geometric_refit.cpp netgen_curving/periodic_copy.cpp)
set_target_properties(netgen_curving PROPERTIES POSITION_INDEPENDENT_CODE ON)
target_link_libraries(netgen_curving PUBLIC nglib ngcore)
if(MSVC)
  target_compile_options(netgen_curving PRIVATE /bigobj /utf-8)
  target_compile_definitions(netgen_curving PUBLIC _DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR)
endif()

if(CUBIT_VOL_EXPORT_PLUGIN)
  find_package(Cubit REQUIRED CONFIG)
  if(NOT TARGET cubit_geom)
    find_package(CubitCGM REQUIRED CONFIG HINTS "${Cubit_DIR}")
    add_library(cubit_geom ALIAS cubit_cgmcore)
  endif()
  add_library(radia_cubit_export_plugin MODULE
    CubitMeshExportPlugin.cpp MeshData.cpp NetgenCurver.cpp callbackgeom.cpp
    ExportNetgenCommand.cpp)
  set_target_properties(radia_cubit_export_plugin PROPERTIES OUTPUT_NAME "radia_cubit_export" PREFIX "" SUFFIX ".dll")
  target_link_libraries(radia_cubit_export_plugin cubiti cubit_util cubit_geom netgen_curving nglib ngcore)
  target_compile_definitions(radia_cubit_export_plugin PRIVATE _DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR)
  target_compile_options(radia_cubit_export_plugin PRIVATE /bigobj /utf-8)
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
