cmake_minimum_required(VERSION 3.26)
project(MyCubitPlugin CXX)

# Cubit_DIR = "<Cubit install>/cmake" (pass with -DCubit_DIR=...).
find_package(Cubit REQUIRED CONFIG)

# Release only: Cubit ships Release binaries and a Debug plugin will not
# load (the verification marker says debug=yes) or will crash on the CRT.
if(NOT CMAKE_BUILD_TYPE)
  set(CMAKE_BUILD_TYPE Release CACHE STRING "" FORCE)
endif()

if(MSVC)
  # RelWithDebInfo-style PDBs so a crash dump in Cubit can be symbolicated.
  add_compile_options($<$<CONFIG:Release>:/Zi> /utf-8)
  add_link_options($<$<CONFIG:Release>:/DEBUG:FULL>)
endif()

add_library(my_cubit_plugin MODULE
  MyPlugin.cpp MyPlugin.hpp
  HelloCommand.cpp HelloCommand.hpp
  ExportNodeCountCommand.cpp ExportNodeCountCommand.hpp
  utf8_path.hpp)
set_target_properties(my_cubit_plugin PROPERTIES CXX_STANDARD 17 CXX_STANDARD_REQUIRED ON)
target_link_libraries(my_cubit_plugin PRIVATE cubiti cubit_util)

if(MSVC)
  # Cubit loads its own bundled MSVCP140.dll, which can be older than the
  # compiler.  Since MSVC 14.40 std::mutex has a constexpr constructor; a
  # mutex built that way crashes in lock() against an older runtime.
  target_compile_definitions(my_cubit_plugin PRIVATE _DISABLE_CONSTEXPR_MUTEX_CONSTRUCTOR)
  # Embed only the PDB file name, not the build machine's path.
  target_link_options(my_cubit_plugin PRIVATE
    "$<$<CONFIG:Release>:/PDBALTPATH:$<TARGET_PDB_FILE_NAME:my_cubit_plugin>>")
endif()

# If the plugin needs geometry (ACIS/CGM) calls, link cubit_geom and
# delay-load it so the plugin still loads before the geometry engine is up:
#   if(NOT TARGET cubit_geom)
#     find_package(CubitCGM REQUIRED CONFIG HINTS "${Cubit_DIR}")
#     add_library(cubit_geom ALIAS cubit_cgmcore)   # 2026.8 split target
#   endif()
#   target_link_libraries(my_cubit_plugin PRIVATE cubit_geom)
#   target_link_options(my_cubit_plugin PRIVATE "/DELAYLOAD:$<TARGET_FILE_NAME:cubit_geom>")
#   target_link_libraries(my_cubit_plugin PRIVATE delayimp)
