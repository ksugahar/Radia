cmake_minimum_required(VERSION 3.16)

# Addon name MUST match the argument of PYBIND11_MODULE() below.
set(addon_name radia_axifemm)

project(${addon_name})
include(ngsolve_addon.cmake)

# Phase 2-B/C are merged: AxiHenrotteFE/FESpace + DiffOps build by default.
# Toggle off only to rebuild the bootstrap smoke-test variant.
option(RADIA_AXIFEMM_PHASE_2B "Build the AxiHenrotteFE/FESpace classes" ON)

set(SOURCES src/radia_axifemm.cpp)
if(RADIA_AXIFEMM_PHASE_2B)
    list(APPEND SOURCES
         src/axi_henrotte_fe.cpp
         src/axi_henrotte_fespace.cpp
         src/axi_henrotte_integrators.cpp)
    add_compile_definitions(RADIA_AXIFEMM_PHASE_2B)
endif()

add_ngsolve_addon(${addon_name} ${SOURCES})

# Install compiled module + Python init/demos
set(install_dir ${ADDON_INSTALL_DIR_PYTHON}/${addon_name})
install(TARGETS ${addon_name} DESTINATION ${install_dir})
install(FILES src/__init__.py DESTINATION ${install_dir})

# Generate IDE stubs (must be last so it picks up installed modules).
ngsolve_generate_stub_files(${addon_name})
