cmake_minimum_required(VERSION 3.21)
project(radia C CXX)

option(RADIA_BUILD_MATLAB_MEX "Build the MATLAB MEX gateway" OFF)

# C++17 required for structured bindings and inline variables
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# FFTW is replaced by Intel MKL DFTI (direct integration in fftw.h)
# No FFTW library linking required - MKL DFTI provides FFT functionality

# Parallelization: NGSolve TaskManager (replaces OpenMP)
# NGSolve is a hard prerequisite for this fork.
# TaskManager provides work-stealing parallel loops via ngstd library.
message(STATUS "Parallelization: NGSolve TaskManager (OpenMP removed)")

# Find Python3 first to get the installation path
find_package(Python3 REQUIRED COMPONENTS Interpreter Development)
message("Python3_EXECUTABLE: ${Python3_EXECUTABLE}")

# Get Python's site-packages directory
execute_process(
	COMMAND ${Python3_EXECUTABLE} -c "import sysconfig; print(sysconfig.get_path('purelib'))"
	OUTPUT_VARIABLE PYTHON_SITE_PACKAGES
	OUTPUT_STRIP_TRAILING_WHITESPACE
)
message("Python site-packages: ${PYTHON_SITE_PACKAGES}")

# Set pybind11 directory
set(pybind11_DIR ${PYTHON_SITE_PACKAGES}/pybind11/share/cmake/pybind11)

message("Finding pybind11...")
find_package(pybind11 REQUIRED)

# NGSolve is always required (provides TaskManager for parallelization)
# Auto-detect via Python import, override with -DNGSolve_DIR=... -DNetgen_DIR=...
if(NOT NGSolve_DIR)
	# Try to find NGSolve cmake config via Python import
	execute_process(
		COMMAND ${Python3_EXECUTABLE} -c "import ngsolve, os; print(os.path.dirname(ngsolve.__file__))"
		OUTPUT_VARIABLE _NGSOLVE_PY_PATH
		OUTPUT_STRIP_TRAILING_WHITESPACE
		RESULT_VARIABLE _NGSOLVE_PY_RESULT
		ERROR_QUIET
	)
	if(_NGSOLVE_PY_RESULT EQUAL 0)
		# pip-installed: site-packages/ngsolve/cmake/
		if(EXISTS "${_NGSOLVE_PY_PATH}/cmake/NGSolveConfig.cmake")
			set(NGSolve_DIR "${_NGSOLVE_PY_PATH}/cmake")
		else()
			# Source-installed: go up from .../Lib/site-packages/ngsolve to install root
			cmake_path(GET _NGSOLVE_PY_PATH PARENT_PATH _SP_DIR)
			cmake_path(GET _SP_DIR PARENT_PATH _LIB_DIR)
			cmake_path(GET _LIB_DIR PARENT_PATH _INSTALL_ROOT)
			if(EXISTS "${_INSTALL_ROOT}/cmake/NGSolveConfig.cmake")
				set(NGSolve_DIR "${_INSTALL_ROOT}/cmake")
			else()
				set(NGSolve_DIR ${PYTHON_SITE_PACKAGES}/ngsolve/cmake)
			endif()
		endif()
	else()
		set(NGSolve_DIR ${PYTHON_SITE_PACKAGES}/ngsolve/cmake)
	endif()
endif()
if(NOT Netgen_DIR)
	# Netgen config is typically in the same directory as NGSolve config
	if(EXISTS "${NGSolve_DIR}/NetgenConfig.cmake")
		set(Netgen_DIR "${NGSolve_DIR}")
	else()
		execute_process(
			COMMAND ${Python3_EXECUTABLE} -c "import netgen, os; print(os.path.dirname(netgen.__file__))"
			OUTPUT_VARIABLE _NETGEN_PY_PATH
			OUTPUT_STRIP_TRAILING_WHITESPACE
			RESULT_VARIABLE _NETGEN_PY_RESULT
			ERROR_QUIET
		)
		if(_NETGEN_PY_RESULT EQUAL 0)
			if(EXISTS "${_NETGEN_PY_PATH}/cmake/NetgenConfig.cmake")
				set(Netgen_DIR "${_NETGEN_PY_PATH}/cmake")
			else()
				cmake_path(GET _NETGEN_PY_PATH PARENT_PATH _SP_DIR2)
				cmake_path(GET _SP_DIR2 PARENT_PATH _LIB_DIR2)
				cmake_path(GET _LIB_DIR2 PARENT_PATH _INSTALL_ROOT2)
				if(EXISTS "${_INSTALL_ROOT2}/cmake/NetgenConfig.cmake")
					set(Netgen_DIR "${_INSTALL_ROOT2}/cmake")
				else()
					set(Netgen_DIR ${PYTHON_SITE_PACKAGES}/netgen/cmake)
				endif()
			endif()
		else()
			set(Netgen_DIR ${PYTHON_SITE_PACKAGES}/netgen/cmake)
		endif()
	endif()
endif()
message(STATUS "NGSolve_DIR: ${NGSolve_DIR}")
message(STATUS "Netgen_DIR: ${Netgen_DIR}")

# Find NGSolve (required for TaskManager parallelization in all targets)
find_package(NGSolve CONFIG REQUIRED)

# Compute the netgen include dir explicitly.  ngsolve 2604 ships a
# bundled pybind11 (e.g. 3.1.0a0 alpha) under
# `<site-packages>/netgen/include/pybind11/`; the system pybind11 from
# PyPI is 3.0.2 stable — their internals layouts differ.  Building
# `_radia_pybind.pyd` against system pybind11 then loading under
# ngsolve 2604 produces `unknown base type "ngfem::CoefficientFunction"`
# because the type registries are not shared across pybind11 versions.
# We need to put the netgen include dir AHEAD of pybind11's system
# include so `#include <pybind11/pybind11.h>` resolves to ngsolve's
# bundled headers.
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import netgen, os; print(os.path.dirname(netgen.__file__))"
    OUTPUT_VARIABLE _NETGEN_PKG_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
)
set(NGSOLVE_BUNDLED_INC "${_NETGEN_PKG_DIR}/include")
if(EXISTS "${NGSOLVE_BUNDLED_INC}/pybind11/pybind11.h")
    message(STATUS "ngsolve bundled pybind11: ${NGSOLVE_BUNDLED_INC}/pybind11")
else()
    message(WARNING "ngsolve bundled pybind11 NOT FOUND at ${NGSOLVE_BUNDLED_INC}/pybind11 -- ABI mismatch with system pybind11 likely")
endif()

# Radia source directories
set(CORE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/core")
set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/lib")
set(EXT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/ext")
set(ACOUSTIC_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/acoustic")
set(PYTHON_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/radia")
# FFTW directory removed - using MKL DFTI via fftw.h compatibility layer


# BLAS/LAPACK: Intel MKL 2026 only
# Use either a full oneAPI installation or the pip mkl-devel package.

# Detect Intel compiler
if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM" OR CMAKE_C_COMPILER_ID MATCHES "IntelLLVM")
	set(USING_INTEL_COMPILER TRUE)
	message(STATUS "Intel oneAPI compiler detected (icx/icx-cl)")
else()
	set(USING_INTEL_COMPILER FALSE)
endif()

if(WIN32)
	# Intel oneAPI MKL (required).  Support both full oneAPI installs and
	# pip mkl-devel installs under Python's Library directory.
	set(_MKL_CANDIDATES)
	if(DEFINED ENV{MKLROOT})
		list(APPEND _MKL_CANDIDATES "$ENV{MKLROOT}")
	endif()
	list(APPEND _MKL_CANDIDATES "C:/Program Files (x86)/Intel/oneAPI/mkl/latest")
	execute_process(
		COMMAND ${Python3_EXECUTABLE} -c "import pathlib, sys; print(pathlib.Path(sys.prefix) / 'Library')"
		OUTPUT_VARIABLE _PYTHON_MKL_ROOT
		OUTPUT_STRIP_TRAILING_WHITESPACE
		ERROR_QUIET
	)
	if(_PYTHON_MKL_ROOT)
		list(APPEND _MKL_CANDIDATES "${_PYTHON_MKL_ROOT}")
	endif()
	set(MKL_DIR "")
	foreach(_MKL_CANDIDATE IN LISTS _MKL_CANDIDATES)
		file(TO_CMAKE_PATH "${_MKL_CANDIDATE}" _MKL_CANDIDATE_NORM)
		if(EXISTS "${_MKL_CANDIDATE_NORM}/lib/mkl_rt.lib")
			set(MKL_DIR "${_MKL_CANDIDATE_NORM}")
			break()
		endif()
	endforeach()
	set(MKL_INCLUDE_DIR "${MKL_DIR}/include")
	set(MKL_LIB_DIR "${MKL_DIR}/lib")
	set(MKL_BIN_DIR "${MKL_DIR}/bin")
	set(MKL_RUNTIME_ABI "3")
	set(MKL_RT_DLL "${MKL_BIN_DIR}/mkl_rt.${MKL_RUNTIME_ABI}.dll")
	set(MKL_RUNTIME_DLLS)
	foreach(_MKL_RUNTIME_NAME IN ITEMS
		mkl_rt
		mkl_core
		mkl_intel_thread
		mkl_avx2
		mkl_def
	)
		set(_MKL_RUNTIME_DLL
			"${MKL_BIN_DIR}/${_MKL_RUNTIME_NAME}.${MKL_RUNTIME_ABI}.dll")
		if(NOT EXISTS "${_MKL_RUNTIME_DLL}")
			message(FATAL_ERROR
				"Intel MKL 2026 runtime not found: ${_MKL_RUNTIME_DLL}. "
				"Install pip package mkl-devel>=2026,<2027.")
		endif()
		list(APPEND MKL_RUNTIME_DLLS "${_MKL_RUNTIME_DLL}")
	endforeach()

	# Intel compiler runtime directory
	set(INTEL_COMPILER_DIR "C:/Program Files (x86)/Intel/oneAPI/compiler/latest")
	set(INTEL_COMPILER_BIN "${INTEL_COMPILER_DIR}/bin")

	# Use mkl_rt.lib (dynamic linking) - works with both Intel and MSVC compilers
	set(MKL_LIB "${MKL_LIB_DIR}/mkl_rt.lib")
	if(EXISTS ${MKL_LIB})
		message(STATUS "Intel MKL 2026 found: ${MKL_LIB}")
		message(STATUS "Intel MKL runtime: ${MKL_RT_DLL}")
		set(BLAS_LIB ${MKL_LIB})
		set(BLAS_INCLUDE_DIR ${MKL_INCLUDE_DIR})
		set(USE_MKL_FOUND TRUE)
	else()
		message(FATAL_ERROR "Intel MKL not found. Install pip package mkl-devel>=2026,<2027.")
	endif()
endif()

# Radia core sources (organized by API functionality)
set(RADIA_CORE_SOURCES
	# API-organized files (new naming convention)
	${CORE_DIR}/rad_obj_creation.cpp          # Object creation (ObjRecMag, ObjCylMag, etc.)
	${CORE_DIR}/rad_material_impl.cpp         # Materials (MatLin, MatSatIsoFrm, MatApl, etc.)
	${CORE_DIR}/rad_transform_impl.cpp        # Transformations (TrfTrsl, TrfRot, etc.)
	${CORE_DIR}/rad_field_compute.cpp         # Field computation (Fld, FldLst, FldEnr, etc.)
	${CORE_DIR}/rad_obj_subdivision.cpp       # Legacy Triangle subdivision stubs
	${CORE_DIR}/rad_point_classify.cpp        # Point classification (inside/outside)
	${CORE_DIR}/rad_field_unified.cpp         # Unified field computation (inside/outside + batch)

	# Geometry files (descriptive names)
	${CORE_DIR}/rad_arc_current.cpp           # Arc with azimuthal current
	${CORE_DIR}/rad_elliptic_integral.cpp     # Elliptic integrals for coil field
	${CORE_DIR}/rad_coefficient_field.cpp     # CoefficientFunction-based field source
	${CORE_DIR}/rad_type_cast.cpp             # Dynamic type casting
	${CORE_DIR}/rad_highorder_nodes.cpp       # High-order BND nodes for GMSH export
	${CORE_DIR}/rad_extruded_polygon.cpp      # Extruded polygon (prism)
	${CORE_DIR}/rad_filament.cpp              # Filament conductor
	${CORE_DIR}/rad_geometry_3d.cpp           # 3D geometry base class
	${CORE_DIR}/rad_group.cpp                 # Container/group of field sources
	${CORE_DIR}/rad_c_interface.cpp           # C-style interface functions
	${CORE_DIR}/rad_interaction.cpp           # Magnetic interaction between objects
	${CORE_DIR}/rad_planar_charges.cpp        # Shared 2D planar exterior field + Maxwell torque
	${CORE_DIR}/rad_io_buffer.cpp             # I/O buffer (errors and warnings)
	${CORE_DIR}/rad_math_methods.cpp          # Mathematical/numerical methods
	${CORE_DIR}/rad_material_def.cpp          # Material relaxation auxiliary
	${CORE_DIR}/rad_planar_2d_part1.cpp       # 2D planar objects (part 1)
	${CORE_DIR}/rad_planar_2d_part2.cpp       # 2D planar objects (part 2)
	${CORE_DIR}/rad_poly_analytical.cpp       # Analytical polygon field formulas
	${CORE_DIR}/rad_rectangular_block.cpp     # Rectangular parallelepiped
	${CORE_DIR}/rad_relaxation_methods.cpp    # Relaxation methods
	# rad_subdivided_*.cpp REMOVED (Phase C, 2026-04-16) - Subdivide API gone, use Cubit/Netgen meshing
	${CORE_DIR}/rad_serialization.cpp         # Interface functions
	${CORE_DIR}/rad_polyhedron.cpp            # Polyhedron with magnetization
	${CORE_DIR}/rad_hacapk.cpp                # HACApK H-matrix manager (Method 2)

	# PEEC Loop-Star matrix construction (Laplace kernel only)
	${CORE_DIR}/rad_peec_matrices.cpp         # PEEC Loop-Star matrix construction
	${CORE_DIR}/rad_hacapk_peec.cpp           # HACApK PEEC adapter (Ruehli filament kernel)
	${CORE_DIR}/rad_hacapk_bem.cpp            # HACApK scalar BEM adapter (Laplace SL/DL Galerkin)
	${CORE_DIR}/rad_hdiv_vim.cpp              # Symmetric HDiv-type VIM demag operator (N = B^T G B)
	${CORE_DIR}/rad_hdiv_hysteresis.cpp       # Energy-based B-input vector Stop material
	${CORE_DIR}/rad_hdiv_field_evaluator.cpp  # Persistent direct/tree BDM1 exterior-field evaluator
	${CORE_DIR}/rad_lie_map_batch.cpp         # Batched Dragt-Finn Lie-map ensemble tracking
	${CORE_DIR}/rad_lie_map_kernel.cpp        # Fourth-order Lie-map construction (jets + stage RK4 + composition)
	${CORE_DIR}/rad_orbit_tracker.cpp         # Native 3D reference-orbit RK4 with measured planarity gate
	${CORE_DIR}/rad_hacapk_hdiv.cpp           # HACApK H-matrix for the HDiv-type VIM demag operator
	${CORE_DIR}/rad_bem_galerkin.cpp          # Fast Galerkin SL/DL assembler (Sauter-Schwab + OpenMP)
	${CORE_DIR}/rad_biot_savart_filaments.cpp # Fast Biot-Savart H/A from finite-segment filaments (complex I)
	${CORE_DIR}/rad_biot_savart_surface.cpp   # Fast Biot-Savart B/A from triangulated surface (complex J_s)
	${CORE_DIR}/rad_equivalence_source.cpp    # Equivalence-theorem near-field source (Schelkunoff/Love static H, Phase A)
	${CORE_DIR}/rad_average_field.cpp         # Closed-form average B over a target box (Wakao Part 6 §7)
	${CORE_DIR}/rad_stream_function.cpp       # (ACA+)+TSVD stream-function solver (HACApK ACA+ + TSVD recompression)
	${CORE_DIR}/rad_evrs_tmethod.cpp          # EVRS/T-method algebra kernel (T --curl--> J, reduced VIM/CLN projection)
	${CORE_DIR}/rad_hybrid_vim_schur.cpp      # Complex Schur/SIBC helper kernels for hybrid VIM / mixed Galerkin reduction
	${CORE_DIR}/rad_motor_rom_c.cpp           # Angle-periodic motor ROM C ABI
	${CORE_DIR}/rad_beam_transfer.cpp         # Distributed canonical R/T/U transfer attribution
	${CORE_DIR}/rad_beam_dynamics.cpp         # Native EarlyTimes migration: fields, Lorentz RHS, RK4/Boris tracking
)

# Fast analytic regression for the dependency-free beam-transfer kernel.
# Solver-heavy accelerator validation remains under validation_test/.
include(CTest)
if(BUILD_TESTING)
	add_executable(test_rad_beam_transfer
		${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp/test_rad_beam_transfer.cpp
		${CORE_DIR}/rad_beam_transfer.cpp
	)
	target_include_directories(test_rad_beam_transfer PRIVATE ${CORE_DIR})
	if(MSVC)
		target_compile_options(test_rad_beam_transfer PRIVATE /W4 /permissive-)
	else()
		target_compile_options(test_rad_beam_transfer PRIVATE -Wall -Wextra -Wpedantic)
	endif()
	add_test(NAME radia.beam_transfer.cpp COMMAND test_rad_beam_transfer)

	add_executable(test_rad_beam_dynamics
		${CMAKE_CURRENT_SOURCE_DIR}/tests/cpp/test_rad_beam_dynamics.cpp
		${CORE_DIR}/rad_beam_dynamics.cpp
	)
	target_include_directories(test_rad_beam_dynamics PRIVATE ${CORE_DIR})
	if(MSVC)
		target_compile_options(test_rad_beam_dynamics PRIVATE /W4 /permissive-)
	else()
		target_compile_options(test_rad_beam_dynamics PRIVATE -Wall -Wextra -Wpedantic)
	endif()
	add_test(NAME radia.beam_dynamics.cpp COMMAND test_rad_beam_dynamics)
endif()

# Acoustic application kernels are deliberately separate from Radia's
# Laplace-only MQS core. The same source is linked into pybind11 and MEX.
set(RADIA_ACOUSTIC_SOURCES
	${ACOUSTIC_DIR}/rad_acoustic_analytic.cpp
)

# Platform-neutral C ABI for Simulink S-functions, Modelica, and FMI wrappers.
# The same source is also compiled into _radia_pybind above so the ABI cannot
# silently drift away from the Python package build.
add_library(radia_motor_rom SHARED
	${CORE_DIR}/rad_motor_rom_c.cpp
)
if(MSVC)
	# This source is also compiled into _radia_pybind and owns the C ABI's
	# exception-to-status translation boundary.
	set_source_files_properties(${CORE_DIR}/rad_motor_rom_c.cpp PROPERTIES
		COMPILE_OPTIONS /EHsc
	)
endif()
target_include_directories(radia_motor_rom PUBLIC ${CORE_DIR})
set_target_properties(radia_motor_rom PROPERTIES
	PREFIX ""
	OUTPUT_NAME "radia_motor_rom"
)
if(MSVC)
	target_compile_options(radia_motor_rom PRIVATE /O2 /W4)
else()
	target_compile_options(radia_motor_rom PRIVATE -O3 -Wall -Wextra)
endif()
add_custom_command(TARGET radia_motor_rom POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
		"$<TARGET_FILE:radia_motor_rom>"
		"${CMAKE_SOURCE_DIR}/src/radia/$<TARGET_FILE_NAME:radia_motor_rom>"
	COMMENT "Copying motor ROM C ABI library to src/radia"
)

# Radia ext sources
set(RADIA_EXT_SOURCES
	${EXT_DIR}/auxparse/auxparse.cpp
	${EXT_DIR}/genmath/gmfft.cpp
	${EXT_DIR}/genmath/gmtrans.cpp
)

# HACApK sources (H-matrix with ACA+ compression)
# Note: C files compiled as C (not C++) to avoid typedef struct * compatibility issues
# Using standard HACApK 1.0.0 (same ACA+ implementation as ELF)
# MPI replaced with stub (single-process), TaskManager parallelization, unified logging
set(HACAPK_SOURCES
	${EXT_DIR}/HACApK/cHACApK_base.c
	${EXT_DIR}/HACApK/cHACApK_lib.c
	${EXT_DIR}/HACApK/cHACApK_cpp_impl.c
	${EXT_DIR}/HACApK/cHACApK_pca_cluster.c  # PCA cluster split for flat meshes (Phase 0.5 of H-LU project)
	${EXT_DIR}/HACApK/cHACApK_block_tree.c   # Block-tree view on flat leaves (Phase 0.7 of H-LU)
	${EXT_DIR}/HACApK/cHACApK_harith.c        # H-arithmetic Phase 1/2 minimal (dense-leaf H-LU)
	${EXT_DIR}/HACApK/cHACApK_harith_par.cpp  # ngcore TaskManager bridge for block-parallel H-LU
	${CORE_DIR}/rad_hacapk_parallel.cpp  # C++ TaskManager wrappers for HACApK C code
)

# Radia lib sources (main radia module)
set(RADIA_LIB_SOURCES
	${LIB_DIR}/radentry.cpp
)

# PEEC matrices pybind11 module sources (separate module)
set(PEEC_MATRICES_API_SOURCE
	${LIB_DIR}/rad_peec_matrices_api.cpp  # pybind11 PEEC matrix construction API
)

# ========================================
# Legacy modules (REMOVED)
# ========================================
# radia_ngsolve: Merged into _radia_pybind (RadiaFieldCF integrated, 2026-03)
# _radia (radpy_pyapi.cpp): Removed 2026-01, replaced by pybind11
#
# For the new pybind11 module, see: _radia_pybind target
# The Python package (__init__.py) imports from _radia_pybind.

# ========================================
# Test Programs
# ========================================

# Test for radTHMatrixFieldSource (Phase 2 validation)
# NOTE: Temporarily disabled - test_radhmat.cpp renamed to test_geometry_extraction.cpp
# add_executable(test_radhmat test_radhmat.cpp)
#
# # Link with Radia core libraries
# target_link_libraries(test_radhmat PRIVATE
# 	hacapk
# )
#
# # Link with FFTW and system libraries
# if(WIN32)
# 	if(EXISTS ${FFTW_LIB})
# 		target_link_libraries(test_radhmat PRIVATE ${FFTW_LIB})
# 	endif()
# 	target_link_libraries(test_radhmat PRIVATE version.lib)
# endif()
#
# # Add Radia core source files directly to test (since radia is a Python module)
# target_sources(test_radhmat PRIVATE
# 	${RADIA_CORE_SOURCES}
# 	${EXT_DIR}/auxparse/auxparse.cpp
# 	${EXT_DIR}/genmath/gmfft.cpp
# 	${EXT_DIR}/genmath/gmtrans.cpp
# 	${RADIA_LIB_SOURCES}
# )
#
# # Add include directories
# target_include_directories(test_radhmat PRIVATE
# 	${CORE_DIR}
# 	${LIB_DIR}
# 	${EXT_DIR}/auxparse
# 	${EXT_DIR}/genmath
# 	${EXT_DIR}/HACApK
# 	${FFTW_DIR}
# )
#
# # Add compile definitions
# target_compile_definitions(test_radhmat PRIVATE
# 	ALPHA__LIB__
# 	FFTW_ENABLE_FLOAT
# 	NO_TIMER
# 	ANSI_DECLARATORS
# 	_GM_WITHOUT_BASE
# 	_CRT_SECURE_NO_WARNINGS
# 	$<$<CONFIG:Debug>:_DEBUG>
# 	$<$<CONFIG:Release>:NDEBUG>
# 	$<$<PLATFORM_ID:Windows>:WIN32>
# 	$<$<PLATFORM_ID:Windows>:_WINDOWS>
# )
#
# # Compiler-specific flags
# if(MSVC)
# 	target_compile_options(test_radhmat PRIVATE
# 		/wd4244  # Disable conversion warnings
# 		/wd4267  # Disable size_t conversion warnings
# 	)
# elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
# 	target_compile_options(test_radhmat PRIVATE
# 		-Wno-conversion
# 		-Wno-sign-conversion
# 		-Wno-unused-variable
# 	)
# elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
# 	target_compile_options(test_radhmat PRIVATE
# 		-Wno-conversion
# 		-Wno-sign-conversion
# 		-Wno-unused-variable
# 	)
# endif()
#
# message(STATUS "test_radhmat executable configured")

# test_polygon_field removed (no longer needed)

# ========================================
# PEEC Matrices Python Module (peec_matrices.pyd)
# ========================================
# Independent pybind11 module for PEEC matrix construction
# Style follows radia_ngsolve.cpp for consistency

pybind11_add_module(peec_matrices
	${PEEC_MATRICES_API_SOURCE}
	${CORE_DIR}/rad_peec_matrices.cpp
)

# Add include directories
target_include_directories(peec_matrices PRIVATE
	${CORE_DIR}
	${EXT_DIR}/genmath
)

# Add compile definitions
target_compile_definitions(peec_matrices PRIVATE
	_CRT_SECURE_NO_WARNINGS
	$<$<CONFIG:Debug>:_DEBUG>
	$<$<CONFIG:Release>:NDEBUG>
	$<$<PLATFORM_ID:Windows>:WIN32>
)

# Intel MKL for BLAS operations (if available)
if(WIN32 AND USE_MKL_FOUND)
	target_link_libraries(peec_matrices PRIVATE ${BLAS_LIB})
	target_include_directories(peec_matrices PRIVATE ${BLAS_INCLUDE_DIR})
	target_compile_definitions(peec_matrices PRIVATE HAVE_LAPACK)
	message(STATUS "Using Intel MKL for peec_matrices: ${BLAS_LIB}")
endif()

# Link Windows system libraries
if(WIN32)
	target_link_libraries(peec_matrices PRIVATE version.lib)
	if(MSVC)
		target_link_options(peec_matrices PRIVATE /NODEFAULTLIB:LIBCMT)
	endif()
endif()

# Parallelization: NGSolve TaskManager (via ngstd)
target_link_libraries(peec_matrices PRIVATE ngstd)
target_include_directories(peec_matrices BEFORE PRIVATE ${NGSOLVE_BUNDLED_INC} ${Netgen_INCLUDE_DIRS} ${NGSolve_INCLUDE_DIRS})

# Compiler-specific flags
if(MSVC)
	target_compile_options(peec_matrices PRIVATE
		/O2
		/fp:fast
		/W0
		/wd4244
		/wd4267
	)
	message(STATUS "Using MSVC compiler for peec_matrices")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
	target_compile_options(peec_matrices PRIVATE /O3 /fp:fast /W0)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	target_compile_options(peec_matrices PRIVATE
		-Wno-conversion
		-Wno-sign-conversion
	)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	target_compile_options(peec_matrices PRIVATE
		-Wno-conversion
		-Wno-sign-conversion
	)
endif()

# Post-build: Copy required DLLs
if(WIN32 AND USE_MKL_FOUND)
	add_custom_command(TARGET peec_matrices POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			"${MKL_RT_DLL}"
			"$<TARGET_FILE_DIR:peec_matrices>"
		COMMENT "Copying mkl_rt.3.dll for peec_matrices"
	)
endif()

message(STATUS "peec_matrices module configured")

# ========================================
# CLN Core Python Module (cln_core.pyd)
# ========================================
# Lanczos algorithm for PEEC model order reduction
# CLN I form: R_diag (diagonal) + L_tridiag (tridiagonal)

set(CLN_CORE_SOURCES
	${LIB_DIR}/rad_cln_api.cpp
	${CORE_DIR}/rad_cln.cpp
)

pybind11_add_module(cln_core
	${CLN_CORE_SOURCES}
)

# Add include directories
target_include_directories(cln_core PRIVATE
	${CORE_DIR}
)

# Add compile definitions
target_compile_definitions(cln_core PRIVATE
	_CRT_SECURE_NO_WARNINGS
	$<$<CONFIG:Debug>:_DEBUG>
	$<$<CONFIG:Release>:NDEBUG>
	$<$<PLATFORM_ID:Windows>:WIN32>
)

# Intel MKL for LAPACK operations (required for Lanczos)
if(WIN32 AND USE_MKL_FOUND)
	target_link_libraries(cln_core PRIVATE ${BLAS_LIB})
	target_include_directories(cln_core PRIVATE ${BLAS_INCLUDE_DIR})
	target_compile_definitions(cln_core PRIVATE HAVE_LAPACK)
	message(STATUS "Using Intel MKL for cln_core: ${BLAS_LIB}")
endif()

# Link Windows system libraries
if(WIN32)
	target_link_libraries(cln_core PRIVATE version.lib)
	if(MSVC)
		target_link_options(cln_core PRIVATE /NODEFAULTLIB:LIBCMT)
	endif()
endif()

# Parallelization: NGSolve TaskManager (via ngstd)
target_link_libraries(cln_core PRIVATE ngstd)
target_include_directories(cln_core BEFORE PRIVATE ${NGSOLVE_BUNDLED_INC} ${Netgen_INCLUDE_DIRS} ${NGSolve_INCLUDE_DIRS})

# Compiler-specific flags
if(MSVC)
	target_compile_options(cln_core PRIVATE
		/O2
		/fp:fast
		/W0
		/wd4244
		/wd4267
	)
	message(STATUS "Using MSVC compiler for cln_core")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
	target_compile_options(cln_core PRIVATE /O3 /fp:fast /W0)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	target_compile_options(cln_core PRIVATE
		-Wno-conversion
		-Wno-sign-conversion
	)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	target_compile_options(cln_core PRIVATE
		-Wno-conversion
		-Wno-sign-conversion
	)
endif()

# Post-build: Copy required DLLs
if(WIN32 AND USE_MKL_FOUND)
	add_custom_command(TARGET cln_core POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			"${MKL_RT_DLL}"
			"$<TARGET_FILE_DIR:cln_core>"
		COMMENT "Copying mkl_rt.3.dll for cln_core"
	)
endif()

message(STATUS "cln_core module configured")

# ========================================
# Radia pybind11 Module (_radia_pybind.pyd)
# ========================================
# pybind11 bindings following NGSolve patterns
# This is the sole Python binding (radpy_pyapi.cpp removed 2026-01)

set(RADIA_PYBIND_SOURCES
	${CMAKE_CURRENT_SOURCE_DIR}/src/lib/radia_pybind.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/lib/rad_beam_transfer_pybind.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/lib/rad_beam_dynamics_pybind.cpp
	${CMAKE_CURRENT_SOURCE_DIR}/src/lib/rad_beam_ngsolve.cpp
)

pybind11_add_module(_radia_pybind
	${RADIA_CORE_SOURCES}
	${RADIA_ACOUSTIC_SOURCES}

	${EXT_DIR}/auxparse/auxparse.cpp
	${EXT_DIR}/genmath/gmfft.cpp
	${EXT_DIR}/genmath/gmtrans.cpp
	${HACAPK_SOURCES}
	${RADIA_LIB_SOURCES}
	${RADIA_PYBIND_SOURCES}
)

# Add include directories
target_include_directories(_radia_pybind PRIVATE
	${CORE_DIR}
	${LIB_DIR}
	${ACOUSTIC_DIR}
	${EXT_DIR}/auxparse
	${EXT_DIR}/genmath
	${EXT_DIR}/HACApK
)

# Add preprocessor definitions
target_compile_definitions(_radia_pybind PRIVATE
	_CRT_SECURE_NO_WARNINGS
	_GM_WITHOUT_BASE
	ALPHA__DLL__
	ANSI_DECLARATORS
	RADIA_USE_HACAPK
	HACAPK_RADIA_LOGGING
)

# Intel MKL for _radia_pybind
if(WIN32 AND USE_MKL_FOUND)
	target_link_libraries(_radia_pybind PRIVATE ${BLAS_LIB})
	target_include_directories(_radia_pybind PRIVATE ${BLAS_INCLUDE_DIR})
	target_compile_definitions(_radia_pybind PRIVATE HAVE_LAPACK)
	message(STATUS "Using Intel MKL for _radia_pybind: ${BLAS_LIB}")
endif()

# Link Windows system libraries
if(WIN32)
	target_link_libraries(_radia_pybind PRIVATE version.lib)
	if(MSVC)
		target_link_options(_radia_pybind PRIVATE /NODEFAULTLIB:LIBCMT)
		# Increase stack size to 16MB for deep HACApK cluster tree recursion
		# Default is 1MB which causes stack overflow at ~30,000+ elements
		target_link_options(_radia_pybind PRIVATE /STACK:16777216)
	endif()
endif()

# NGSolve integration: TaskManager (ngstd) + CoefficientFunction (ngfem, ngcomp)
target_link_libraries(_radia_pybind PRIVATE ngstd ngfem ngcomp ngbla)
# Use BEFORE so the bundled pybind11 (in netgen/include/pybind11/) is
# found ahead of the system pybind11.  ngsolve 2604 ships pybind11
# 3.1.0a0 (not on PyPI; only 3.0.x is) and pybind11's `internals` type
# registry has a different layout between 3.0.x and 3.1.0a0 — building
# `_radia_pybind` against system pybind11 3.0.2 then loading under
# ngsolve 2604 fails with `generic_type: type "RadiaField" referenced
# unknown base type "ngfem::CoefficientFunction"`.  Putting the netgen
# include dir BEFORE the system pybind11 path makes
# `#include <pybind11/pybind11.h>` resolve to ngsolve's bundled
# 3.1.0a0 headers, which matches ngsolve's compiled-in pybind11.
target_include_directories(_radia_pybind BEFORE PRIVATE
    ${NGSOLVE_BUNDLED_INC}
    ${Netgen_INCLUDE_DIRS}
    ${NGSolve_INCLUDE_DIRS}
)

# Compiler-specific flags
if(MSVC)
	target_compile_options(_radia_pybind PRIVATE
		/O2
		/fp:fast
		/W0
		/wd4244
		/wd4267
	)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
	target_compile_options(_radia_pybind PRIVATE /O3 /fp:fast /W0)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	target_compile_options(_radia_pybind PRIVATE
		-Wno-conversion
		-Wno-sign-conversion
	)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	target_compile_options(_radia_pybind PRIVATE
		-Wno-conversion
		-Wno-sign-conversion
	)
endif()

# Post-build: Copy to src/radia
add_custom_command(TARGET _radia_pybind POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
		"$<TARGET_FILE:_radia_pybind>"
		"${CMAKE_SOURCE_DIR}/src/radia/_radia_pybind.pyd"
	COMMENT "Copying _radia_pybind.pyd to src/radia"
)

message(STATUS "_radia_pybind module configured (pybind11 bindings)")

# ========================================
# MATLAB MEX gateway (optional)
# ========================================
# Stateful C++ objects use checked uint64 handles.  NGSolve remains the owner
# of TaskManager, mesh loading, finite-element orientation, and space updates.
if(RADIA_BUILD_MATLAB_MEX)
	find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY MEX_LIBRARY)
	matlab_add_mex(
		NAME radia_mex
		R2018a
		SRC
			${CMAKE_CURRENT_SOURCE_DIR}/src/matlab/radia_mex.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/src/matlab/radia_ih_mex_commands.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/src/matlab/radia_reactor_mex_commands.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/src/matlab/radia_beam_mex_commands.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/src/lib/rad_beam_ngsolve.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/src/radia/simulink/radia_ih_runtime.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/src/radia/simulink/radia_ih_transport.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/src/radia/simulink/radia_ih_thermal.cpp
			${CMAKE_CURRENT_SOURCE_DIR}/src/radia/simulink/radia_nonlinear_reactor_runtime.cpp
			${CORE_DIR}/rad_cln.cpp
			${RADIA_CORE_SOURCES}
			${RADIA_ACOUSTIC_SOURCES}
			${EXT_DIR}/axifem/axi_henrotte_numeric.cpp
			${EXT_DIR}/auxparse/auxparse.cpp
			${EXT_DIR}/genmath/gmfft.cpp
			${EXT_DIR}/genmath/gmtrans.cpp
			${HACAPK_SOURCES}
			# The legacy C ABI is still the native Radia compatibility layer.  Its
			# MEX build is Python-free because radentry.cpp is compiled with the
			# RADIA_MEX_NO_PYTHON feature below.
			${RADIA_LIB_SOURCES}
	)
	target_include_directories(radia_mex PRIVATE
		${CMAKE_CURRENT_SOURCE_DIR}/src/matlab
		${CMAKE_CURRENT_SOURCE_DIR}/src/radia/simulink
		${CORE_DIR}
		${LIB_DIR}
		${ACOUSTIC_DIR}
		${EXT_DIR}/axifem
		${EXT_DIR}/auxparse
		${EXT_DIR}/genmath
		${EXT_DIR}/HACApK
		${NGSOLVE_BUNDLED_INC}
		${Netgen_INCLUDE_DIRS}
		${NGSolve_INCLUDE_DIRS}
		${Python3_INCLUDE_DIRS}
	)
	target_compile_definitions(radia_mex PRIVATE
		_CRT_SECURE_NO_WARNINGS
		_GM_WITHOUT_BASE
		ALPHA__DLL__
		ANSI_DECLARATORS
		RADIA_USE_HACAPK
		HACAPK_RADIA_LOGGING
		NOMINMAX
		RADIA_MEX_NO_PYTHON
	)
	target_compile_features(radia_mex PRIVATE cxx_std_20)
	# FindMatlab uses the plain signature internally, so this call must too.
	# The MEX gateway uses NGSolve's native C++ API directly.  Do not link the
	# Python runtime here: MATLAB/Simulink must be able to load radia_mex
	# without starting or embedding Python.
	target_link_libraries(radia_mex ngstd ngfem ngcomp ngbla)
	if(WIN32 AND USE_MKL_FOUND)
		target_link_libraries(radia_mex ${BLAS_LIB})
		target_include_directories(radia_mex PRIVATE ${BLAS_INCLUDE_DIR})
		target_compile_definitions(radia_mex PRIVATE HAVE_LAPACK)
	endif()
	if(MSVC)
		# The bundled Windows NGSolve package is built with Python support and
		# propagates these feature macros through its imported targets.  The MEX
		# gateway uses only the native C++ API, so suppress Python-only header
		# instantiations here; otherwise radia_mex.obj requests python312.lib.
		target_compile_options(radia_mex PRIVATE
			/O2 /fp:fast /W0
			/UNGS_PYTHON /UNETGEN_PYTHON /UNG_PYTHON
			/UPYBIND11_SIMPLE_GIL_MANAGEMENT
		)
		target_link_options(radia_mex PRIVATE /NODEFAULTLIB:LIBCMT)
		target_link_libraries(radia_mex version.lib)
	endif()
	add_custom_command(TARGET radia_mex POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E make_directory "${CMAKE_SOURCE_DIR}/matlab"
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			"$<TARGET_FILE:radia_mex>"
			"${CMAKE_SOURCE_DIR}/matlab/$<TARGET_FILE_NAME:radia_mex>"
		COMMENT "Copying radia_mex to matlab"
	)
	if(WIN32 AND USE_MKL_FOUND)
		# MATLAB resolves MEX dependencies from the module directory.  Keep the
		# oneMKL dispatch DLLs beside radia_mex so it does not depend on the
		# host Python installation.  Do not copy libiomp5md.dll: MATLAB already
		# owns an OpenMP runtime and loading a second copy aborts the process.
		foreach(_MKL_RUNTIME_DLL IN LISTS MKL_RUNTIME_DLLS)
			get_filename_component(_MKL_RUNTIME_FILENAME "${_MKL_RUNTIME_DLL}" NAME)
			add_custom_command(TARGET radia_mex POST_BUILD
				COMMAND ${CMAKE_COMMAND} -E copy_if_different
					"${_MKL_RUNTIME_DLL}"
					"${CMAKE_SOURCE_DIR}/matlab/${_MKL_RUNTIME_FILENAME}"
				COMMENT "Copying ${_MKL_RUNTIME_FILENAME} beside radia_mex"
			)
		endforeach()
	endif()
	message(STATUS "radia_mex configured (MATLAB R2018a interleaved complex API)")
endif()

# ========================================
# sparsesolv_ngsolve module (NGSolve preconditioners + Krylov solvers)
# ========================================
# Compact AMS + COCR / ICCG / SGS for HCurl eddy-current systems.
# Source: src/ext/sparsesolv/ (formerly the standalone ngsolve-sparsesolv
# PyPI package; absorbed into the Radia wheel 2026-05-08 to fix the
# 100号機 ImportError that was blocking calc_fem_kelvin and
# vector_potential_solver).
#
# Imported as `radia.sparsesolv_ngsolve` (the .pyd lives at
# src/radia/sparsesolv_ngsolve.pyd alongside _radia_pybind.pyd).

set(SPARSESOLV_DIR "${EXT_DIR}/sparsesolv")

add_ngsolve_python_module(sparsesolv_ngsolve
	${SPARSESOLV_DIR}/ngsolve/python_module.cpp
)

target_include_directories(sparsesolv_ngsolve PRIVATE
	${SPARSESOLV_DIR}/include
)

# pybind11 ABI fix (same rationale as _radia_pybind line ~575): force
# ngsolve's bundled pybind11 headers ahead of system pybind11.
target_include_directories(sparsesolv_ngsolve BEFORE PRIVATE ${NGSOLVE_BUNDLED_INC})

target_compile_definitions(sparsesolv_ngsolve PRIVATE
	SPARSESOLV_USE_NGSOLVE_TASKMANAGER
	SPARSESOLV_USE_MKL
	_CRT_SECURE_NO_WARNINGS
	$<$<CONFIG:Debug>:_DEBUG>
	$<$<CONFIG:Release>:NDEBUG>
)

if(MSVC)
	target_compile_options(sparsesolv_ngsolve PRIVATE
		/arch:AVX2
		/O2
		/fp:fast
		/W0
		/wd4244
		/wd4267
	)
else()
	target_compile_options(sparsesolv_ngsolve PRIVATE -mavx2 -mfma)
endif()

# Intel MKL for BLAS/LAPACK (sparsesolv uses dense ops in AMG coarse solve)
if(WIN32 AND USE_MKL_FOUND)
	target_link_libraries(sparsesolv_ngsolve PRIVATE ${BLAS_LIB})
	target_include_directories(sparsesolv_ngsolve PRIVATE ${BLAS_INCLUDE_DIR})
endif()

# Post-build: Copy to src/radia
add_custom_command(TARGET sparsesolv_ngsolve POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
		"$<TARGET_FILE:sparsesolv_ngsolve>"
		"${CMAKE_SOURCE_DIR}/src/radia/sparsesolv_ngsolve.pyd"
	COMMENT "Copying sparsesolv_ngsolve.pyd to src/radia"
)

message(STATUS "sparsesolv_ngsolve module configured (NGSolve preconditioners)")

# ========================================
# axifem module (Henrotte/Meeker axisymmetric FE for NGSolve)
# ========================================
# FEMM-style axisymmetric finite elements: shape functions linear in
# {1, r^2, z}, giving B_z constant per element and B_r ~ 1/r per element
# without the standard P1 axis singularity.  P1 triangles + Q1/Q2 axis-
# aligned quads + bilinear form integrators (Stiffness, SigmaMass).
#
# Source: src/ext/axifem/ (formerly a standalone axifem add-on; absorbed
# into the radia-core wheel 2026-05-10 to fix the
# 100号機 / mdx ImportError that was blocking IH axisymmetric work).
#
# Imported as `radia.axifem` (the .pyd lives at
# src/radia/axifem.pyd alongside _radia_pybind.pyd).
#
# Reference: F. Henrotte et al., IEEE TMag 9(2):1352-1355, 1993.
#            D. Meeker, FEMM 4.2 axisymmetric formulation notes.

set(AXIFEM_DIR "${EXT_DIR}/axifem")

add_ngsolve_python_module(axifem
	${AXIFEM_DIR}/axifem.cpp
	${AXIFEM_DIR}/axi_henrotte_fe.cpp
	${AXIFEM_DIR}/axi_henrotte_fespace.cpp
	${AXIFEM_DIR}/axi_henrotte_integrators.cpp
	${AXIFEM_DIR}/axi_henrotte_numeric.cpp
)

target_include_directories(axifem PRIVATE
	${AXIFEM_DIR}
)

# pybind11 ABI fix (same rationale as _radia_pybind line ~575): force
# ngsolve's bundled pybind11 headers ahead of system pybind11.  Without
# this, axifem.pyd is built against system pybind11 3.0.2 and
# fails at import with `ImportError: generic_type: ... referenced unknown
# base type "ngfem::FiniteElement"` when ngsolve 2604 uses bundled
# 3.1.0a0 at runtime.  (Caught by v4.80.2 CI 2026-05-28 in
# test_heat_henrotte_smoke.)
target_include_directories(axifem BEFORE PRIVATE ${NGSOLVE_BUNDLED_INC})

target_compile_definitions(axifem PRIVATE
	RADIA_AXIFEM_PHASE_2B
	_CRT_SECURE_NO_WARNINGS
	$<$<CONFIG:Debug>:_DEBUG>
	$<$<CONFIG:Release>:NDEBUG>
)

if(MSVC)
	target_compile_options(axifem PRIVATE
		/O2
		/fp:fast
		/W0
		/wd4244
		/wd4267
	)
endif()

# Post-build: Copy to src/radia
add_custom_command(TARGET axifem POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
		"$<TARGET_FILE:axifem>"
		"${CMAKE_SOURCE_DIR}/src/radia/axifem.pyd"
	COMMENT "Copying axifem.pyd to src/radia"
)

message(STATUS "axifem module configured (Henrotte axisymmetric FE)")

# ========================================
# _equation module (LaTeX equations: edit, render, put into Office)
# ========================================
# Source: src/ext/equation/.  Imported as `radia.equation`, which re-exports
# this .pyd plus the Office writers in src/radia/equation/office.py.
#
# An equation is stored as LaTeX, usually inside a Markdown file, and every
# output is produced from that: OMML for Word and PowerPoint (a native,
# editable equation, so nothing has to be installed on the reader's machine),
# SVG for a picture, and a layout the editor draws with.  MTEF, the Equation
# Editor 3.x binary format, is here only to read equations out of documents
# that already contain them.
#
# No NGSolve and no MKL: this is text and font metrics, not numerics.  GDI
# supplies the glyph boxes the layout needs.

set(EQUATION_DIR "${EXT_DIR}/equation")

pybind11_add_module(_equation
	${EQUATION_DIR}/equation_pybind.cpp
	${EQUATION_DIR}/mtef_node.cpp
	${EQUATION_DIR}/mtef_parser.cpp
	${EQUATION_DIR}/line_pass.cpp
	${EQUATION_DIR}/latex_emitter.cpp
	${EQUATION_DIR}/mtef2tex.cpp
	${EQUATION_DIR}/tex2mtef.cpp
	${EQUATION_DIR}/tex_parser.cpp
	${EQUATION_DIR}/math_writer.cpp
	${EQUATION_DIR}/mtef_omml.cpp
	${EQUATION_DIR}/mtef_rtf.cpp
	${EQUATION_DIR}/mtef_mathml.cpp
	${EQUATION_DIR}/mtef_svg.cpp
	${EQUATION_DIR}/mtef_gdi.cpp
	${EQUATION_DIR}/mtef_dump.cpp
	${EQUATION_DIR}/eq_edit.cpp
	${EQUATION_DIR}/md_doc.cpp
	${EQUATION_DIR}/md_blocks.cpp
	${EQUATION_DIR}/md_layout.cpp
)

target_include_directories(_equation PRIVATE
	${EQUATION_DIR}
)

target_compile_definitions(_equation PRIVATE
	_CRT_SECURE_NO_WARNINGS
	$<$<CONFIG:Debug>:_DEBUG>
	$<$<CONFIG:Release>:NDEBUG>
	$<$<PLATFORM_ID:Windows>:WIN32>
)

if(WIN32)
	# Glyph metrics come from GDI: the ink box of each glyph, not the font's
	# global ascent and descent, is what a script has to clear.
	target_link_libraries(_equation PRIVATE gdi32 user32 gdiplus ole32)
endif()

if(MSVC)
	target_compile_options(_equation PRIVATE /O2 /W0 /wd4244 /wd4267)
endif()

# Post-build: Copy to src/radia
add_custom_command(TARGET _equation POST_BUILD
	COMMAND ${CMAKE_COMMAND} -E copy_if_different
		"$<TARGET_FILE:_equation>"
		"${CMAKE_SOURCE_DIR}/src/radia/_equation.pyd"
	COMMENT "Copying _equation.pyd to src/radia"
)

message(STATUS "_equation module configured (LaTeX equations -> OMML / SVG)")
