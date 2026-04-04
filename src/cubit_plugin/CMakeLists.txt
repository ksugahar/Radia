cmake_minimum_required(VERSION 3.12)
project(RadiaCubitPlugin)

# --------------------------------------------------------
# Cubit SDK
# --------------------------------------------------------
find_package(Cubit REQUIRED CONFIG)

# --------------------------------------------------------
# Netgen (optional, for high-order curving order >= 2)
# --------------------------------------------------------
set(NETGEN_DIR "" CACHE PATH "Netgen installation directory (optional)")

if(NETGEN_DIR AND EXISTS "${NETGEN_DIR}/lib/nglib.lib")
  message(STATUS "Netgen found: ${NETGEN_DIR}")
  set(HAVE_NETGEN ON)

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
  message(STATUS "Netgen not found (NETGEN_DIR=${NETGEN_DIR}). High-order disabled.")
  set(HAVE_NETGEN OFF)
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
# Registers APREPRO commands: export gmsh/nastran/vtk/meg
# Loaded by Python API (CUBIT_PLUGIN_DIR) and by cubitcomp
# (-commandplugindir). May link Netgen for high-order curving.
# --------------------------------------------------------
set(CCM_SOURCES
    RadiaPlugin.cpp
    RadiaPlugin.hpp
    MeshData.cpp
    MeshData.hpp
    HighOrderMesh.cpp
    HighOrderMesh.hpp
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
# ExportNetgenCommand uses MeshData+NetgenCurver, requires Netgen DLLs
if(HAVE_NETGEN)
  list(APPEND CCM_SOURCES
    callbackgeom.cpp callbackgeom.hpp
    ExportNetgenCommand.cpp ExportNetgenCommand.hpp
  )
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
  target_include_directories(radia_cubit_ccm PRIVATE
    "${NETGEN_DIR}/include"
    "${NETGEN_DIR}/include/include"
  )
  target_link_libraries(radia_cubit_ccm netgen::nglib netgen::ngcore)
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

if(pybind11_FOUND AND HAVE_NETGEN)
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

if(HAVE_NETGEN)
  install(FILES
    "${NETGEN_DIR}/nglib.dll"
    "${NETGEN_DIR}/ngcore.dll"
    DESTINATION "${CMAKE_INSTALL_PREFIX}/plugins"
  )
endif()
