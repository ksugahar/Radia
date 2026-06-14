cmake_minimum_required(VERSION 3.16)

# Addon name MUST match the argument of PYBIND11_MODULE() below.
set(addon_name radia_vim)

project(${addon_name})
include(ngsolve_addon.cmake)

# Phase F-1 bootstrap: skeleton with version()/hello()/precision_info() only.
# Subsequent phases gated by RADIA_VIM_PHASE_F2 and friends:
#   F-2: HDiv div-free hex basis (any order p)
#   F-3: Self-coincident Sauter-Schwab quadrature
#   F-4: K matrix assembler (OpenMP parallel)
#   F-5: Foster solver + Cauer pipeline (Python wrapper)
option(RADIA_VIM_PHASE_F2 "Build the HDiv div-free hex basis" ON)
option(RADIA_VIM_PHASE_F3 "Build the Sauter-Schwab quadrature routines" ON)
option(RADIA_VIM_PHASE_F4 "Build the K matrix assembler" ON)

# Optional high-precision support via Boost.Multiprecision. Quad and MPFR are
# only needed for orders >= 5 (see precision_sensitivity_study.py). Defaulting
# off keeps F-1 dependency-light.
option(RADIA_VIM_WITH_QUAD "Enable __float128 / cpp_bin_float_quad type"  OFF)
option(RADIA_VIM_WITH_MPFR "Enable mpfr_float arbitrary-precision type"   OFF)

set(SOURCES src/radia_vim.cpp)

if(RADIA_VIM_PHASE_F2)
    list(APPEND SOURCES src/hdiv_divfree_hex_basis.cpp)
    add_compile_definitions(RADIA_VIM_PHASE_F2)
endif()
if(RADIA_VIM_PHASE_F3)
    list(APPEND SOURCES src/sauter_schwab_hex.cpp)
    add_compile_definitions(RADIA_VIM_PHASE_F3)
endif()
if(RADIA_VIM_PHASE_F4)
    list(APPEND SOURCES src/volume_galerkin.cpp)
    add_compile_definitions(RADIA_VIM_PHASE_F4)
endif()

add_ngsolve_addon(${addon_name} ${SOURCES})

if(RADIA_VIM_WITH_QUAD)
    target_compile_definitions(${addon_name} PRIVATE RADIA_VIM_WITH_QUAD)
    if(NOT MSVC)
        # GCC/Clang: link quadmath for __float128 support
        target_link_libraries(${addon_name} PRIVATE quadmath)
    endif()
endif()
if(RADIA_VIM_WITH_MPFR)
    target_compile_definitions(${addon_name} PRIVATE RADIA_VIM_WITH_MPFR)
    target_link_libraries(${addon_name} PRIVATE mpfr gmp)
endif()

# OpenMP for K matrix parallel assembly (F-4 and onwards)
find_package(OpenMP)
if(OpenMP_CXX_FOUND AND RADIA_VIM_PHASE_F4)
    target_link_libraries(${addon_name} PRIVATE OpenMP::OpenMP_CXX)
endif()

# Install compiled module + Python init
set(install_dir ${ADDON_INSTALL_DIR_PYTHON}/${addon_name})
install(TARGETS ${addon_name} DESTINATION ${install_dir})
install(FILES src/__init__.py DESTINATION ${install_dir})

# Generate IDE stubs (must be last so it picks up installed modules).
ngsolve_generate_stub_files(${addon_name})
