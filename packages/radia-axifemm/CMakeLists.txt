cmake_minimum_required(VERSION 3.16)

# Addon name MUST match the argument of PYBIND11_MODULE() below.
set(addon_name radia_axifemm)

project(${addon_name})
include(ngsolve_addon.cmake)

# Phase 2-A: build only the smoke-test entry; FE/FESpace land in 2-B.
# Set RADIA_AXIFEMM_PHASE_2B=ON once the abstract methods are implemented.
option(RADIA_AXIFEMM_PHASE_2B "Build the AxiHenrotteFE/FESpace classes" OFF)

set(SOURCES src/radia_axifemm.cpp)
if(RADIA_AXIFEMM_PHASE_2B)
    list(APPEND SOURCES src/axi_henrotte_fe.cpp src/axi_henrotte_fespace.cpp)
    add_compile_definitions(RADIA_AXIFEMM_PHASE_2B)
endif()

add_ngsolve_addon(${addon_name} ${SOURCES})

# Install compiled module + Python init/demos
set(install_dir ${ADDON_INSTALL_DIR_PYTHON}/${addon_name})
install(TARGETS ${addon_name} DESTINATION ${install_dir})
install(FILES src/__init__.py DESTINATION ${install_dir})

# Generate IDE stubs (must be last so it picks up installed modules).
ngsolve_generate_stub_files(${addon_name})
