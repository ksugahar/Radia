cmake_minimum_required(VERSION 3.21)
project(radia_ngsolve C CXX)

# C++17 required for inline variables in header-only ExaFMM library
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# Option to disable OpenMP (for debugging)
option(RADIA_ENABLE_OPENMP "Enable OpenMP parallelization" ON)

# Option to build radia_ngsolve module (requires NGSolve to be installed)
option(RADIA_BUILD_NGSOLVE "Build radia_ngsolve integration module" ON)

# FFTW is replaced by Intel MKL DFTI (direct integration in fftw.h)
# No FFTW library linking required - MKL DFTI provides FFT functionality

# Find OpenMP for parallel matrix construction
if(RADIA_ENABLE_OPENMP)
    find_package(OpenMP)
    if(OpenMP_CXX_FOUND)
        message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
    else()
        message(WARNING "OpenMP not found - parallel matrix construction disabled")
    endif()
else()
    message(STATUS "OpenMP disabled by RADIA_ENABLE_OPENMP=OFF")
endif()

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

# Set NGSolve and Netgen directories (only if building radia_ngsolve)
if(RADIA_BUILD_NGSOLVE)
	# Use pip-installed NGSolve from site-packages (version-agnostic)
	# Override with -DNGSolve_DIR=... -DNetgen_DIR=... if needed
	if(NOT NGSolve_DIR)
		set(NGSolve_DIR ${PYTHON_SITE_PACKAGES}/ngsolve/cmake)
	endif()
	if(NOT Netgen_DIR)
		set(Netgen_DIR ${PYTHON_SITE_PACKAGES}/netgen/cmake)
	endif()
	message(STATUS "NGSolve_DIR: ${NGSolve_DIR}")
	message(STATUS "Netgen_DIR: ${Netgen_DIR}")

	# Find NGSolve
	find_package(NGSolve CONFIG REQUIRED)
endif()

# Radia source directories
set(CORE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/core")
set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/lib")
set(EXT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/ext")
set(PYTHON_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/radia")
# FFTW directory removed - using MKL DFTI via fftw.h compatibility layer


# BLAS/LAPACK: Intel MKL only
# Requires Intel oneAPI to be installed

# Detect Intel compiler
if(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM" OR CMAKE_C_COMPILER_ID MATCHES "IntelLLVM")
	set(USING_INTEL_COMPILER TRUE)
	message(STATUS "Intel oneAPI compiler detected (icx/icx-cl)")
else()
	set(USING_INTEL_COMPILER FALSE)
endif()

if(WIN32)
	# Intel oneAPI MKL (required)
	set(MKL_DIR "C:/Program Files (x86)/Intel/oneAPI/mkl/latest")
	set(MKL_INCLUDE_DIR "${MKL_DIR}/include")
	set(MKL_LIB_DIR "${MKL_DIR}/lib")
	set(MKL_BIN_DIR "${MKL_DIR}/bin")

	# Intel compiler runtime directory
	set(INTEL_COMPILER_DIR "C:/Program Files (x86)/Intel/oneAPI/compiler/latest")
	set(INTEL_COMPILER_BIN "${INTEL_COMPILER_DIR}/bin")

	# Use mkl_rt.lib (dynamic linking) - works with both Intel and MSVC compilers
	set(MKL_LIB "${MKL_LIB_DIR}/mkl_rt.lib")
	if(EXISTS ${MKL_LIB})
		message(STATUS "Intel MKL found: ${MKL_LIB}")
		set(BLAS_LIB ${MKL_LIB})
		set(BLAS_INCLUDE_DIR ${MKL_INCLUDE_DIR})
		# MKL FFTW3 wrapper for ExaFMM-t compatibility
		set(MKL_FFTW3_INCLUDE_DIR "${MKL_INCLUDE_DIR}/fftw")
		set(USE_MKL_FOUND TRUE)
	else()
		message(FATAL_ERROR "Intel MKL not found. Please install Intel oneAPI Base Toolkit.")
	endif()
endif()

# Option to enable ExaFMM acceleration for fast dipole field computation
# ExaFMM-t is a pure C++ FMM library (no Fortran dependency)
option(RADIA_ENABLE_EXAFMM "Enable ExaFMM library for fast field computation" ON)

# Radia core sources (organized by API functionality)
set(RADIA_CORE_SOURCES
	# API-organized files (new naming convention)
	${CORE_DIR}/rad_obj_creation.cpp          # Object creation (ObjRecMag, ObjCylMag, etc.)
	${CORE_DIR}/rad_material_impl.cpp         # Materials (MatLin, MatSatIsoFrm, MatApl, etc.)
	${CORE_DIR}/rad_transform_impl.cpp        # Transformations (TrfTrsl, TrfRot, etc.)
	${CORE_DIR}/rad_field_compute.cpp         # Field computation (Fld, FldLst, FldEnr, etc.)
	${CORE_DIR}/rad_obj_subdivision.cpp       # Subdivision (ObjDivMag, multi-generation extrusion)
	${CORE_DIR}/rad_exafmm.cpp                # ExaFMM wrapper for fast multipole acceleration
	${CORE_DIR}/rad_point_classify.cpp        # Point classification for FMM near/far
	${CORE_DIR}/rad_dipole_collect.cpp        # Collect dipoles from Radia for FMM
	${CORE_DIR}/rad_field_unified.cpp         # Unified field computation (inside/outside + FMM)

	# Geometry files (descriptive names)
	${CORE_DIR}/rad_arc_current.cpp           # Arc with azimuthal current
	${CORE_DIR}/rad_elliptic_integral.cpp     # Elliptic integrals for coil field
	${CORE_DIR}/rad_coefficient_field.cpp     # CoefficientFunction-based field source
	${CORE_DIR}/rad_type_cast.cpp             # Dynamic type casting
	${CORE_DIR}/rad_extruded_polygon.cpp      # Extruded polygon (prism)
	${CORE_DIR}/rad_filament.cpp              # Filament conductor
	${CORE_DIR}/rad_geometry_3d.cpp           # 3D geometry base class
	${CORE_DIR}/rad_graphics_3d.cpp           # 3D graphics/visualization
	${CORE_DIR}/rad_group.cpp                 # Container/group of field sources
	${CORE_DIR}/rad_c_interface.cpp           # C-style interface functions
	${CORE_DIR}/rad_interaction.cpp           # Magnetic interaction between objects
	${CORE_DIR}/rad_io_buffer.cpp             # I/O buffer (errors and warnings)
	${CORE_DIR}/rad_math_methods.cpp          # Mathematical/numerical methods
	${CORE_DIR}/rad_material_def.cpp          # Material relaxation auxiliary
	${CORE_DIR}/rad_planar_2d_part1.cpp       # 2D planar objects (part 1)
	${CORE_DIR}/rad_planar_2d_part2.cpp       # 2D planar objects (part 2)
	${CORE_DIR}/rad_poly_analytical.cpp       # Analytical polygon field formulas
	${CORE_DIR}/rad_particle_trajectory.cpp   # Particle trajectory/dynamics
	${CORE_DIR}/rad_rectangular_block.cpp     # Rectangular parallelepiped
	${CORE_DIR}/rad_relaxation_methods.cpp    # Relaxation methods
	${CORE_DIR}/rad_subdivided_arc_current.cpp        # Subdivided arc with current
	${CORE_DIR}/rad_subdivided_extruded_polygon.cpp   # Subdivided extruded polygon
	${CORE_DIR}/rad_subdivided_rectangle.cpp          # Subdivided rectangle
	${CORE_DIR}/rad_subdivided_polyhedron.cpp         # Subdivided polyhedron
	${CORE_DIR}/rad_serialization.cpp         # Interface functions
	${CORE_DIR}/rad_polyhedron.cpp            # Polyhedron with magnetization
	${CORE_DIR}/rad_hacapk.cpp                # HACApK H-matrix manager (Method 2)

	# PEEC Loop-Star matrix construction (Laplace kernel only)
	${CORE_DIR}/rad_peec_matrices.cpp         # PEEC Loop-Star matrix construction
)

# Radia ext sources
set(RADIA_EXT_SOURCES
	${EXT_DIR}/auxparse/auxparse.cpp
	${EXT_DIR}/genmath/gmfft.cpp
	${EXT_DIR}/genmath/gmtrans.cpp
	${EXT_DIR}/triangle/triangle.c
)

# HACApK sources (H-matrix with ACA+ compression)
# Note: C files compiled as C (not C++) to avoid typedef struct * compatibility issues
# Using standard HACApK 1.0.0 (same ACA+ implementation as ELF)
# MPI replaced with stub (single-process), OpenMP only, unified logging
set(HACAPK_SOURCES
	${EXT_DIR}/HACApK/cHACApK_base.c
	${EXT_DIR}/HACApK/cHACApK_lib.c
	${EXT_DIR}/HACApK/cHACApK_cpp_impl.c
)

# Radia lib sources (main radia module)
set(RADIA_LIB_SOURCES
	${LIB_DIR}/radentry.cpp
)

# PEEC matrices pybind11 module sources (separate module)
set(PEEC_MATRICES_API_SOURCE
	${LIB_DIR}/rad_peec_matrices_api.cpp  # pybind11 PEEC matrix construction API
)

# Triangle library (C source needs special handling)
# Create as static library first
add_library(triangle_lib STATIC ${EXT_DIR}/triangle/triangle.c)
target_compile_definitions(triangle_lib PRIVATE TRILIBRARY ANSI_DECLARATORS NO_TIMER)
target_include_directories(triangle_lib PUBLIC ${EXT_DIR}/triangle)

# Create NGSolve Python module using add_ngsolve_python_module
# (following EMPY_Field pattern exactly)
# Only build if RADIA_BUILD_NGSOLVE is ON (requires NGSolve to be installed)
if(RADIA_BUILD_NGSOLVE)
add_ngsolve_python_module(radia_ngsolve
	${RADIA_CORE_SOURCES}

	${EXT_DIR}/auxparse/auxparse.cpp
	${EXT_DIR}/genmath/gmfft.cpp
	${EXT_DIR}/genmath/gmtrans.cpp
	${HACAPK_SOURCES}
	${RADIA_LIB_SOURCES}
	${PYTHON_DIR}/radia_ngsolve.cpp
)

# Disable IPO/LTO for Intel compiler (produces LLVM bitcode, not COFF objects)
# This is required when linking with MSVC-compiled NGSolve libraries
if(CMAKE_CXX_COMPILER_ID STREQUAL "IntelLLVM")
	set_target_properties(radia_ngsolve PROPERTIES
		INTERPROCEDURAL_OPTIMIZATION OFF
	)
	# Force object file format to COFF (disable LLVM bitcode)
	# $<$<COMPILE_LANGUAGE:CXX>:-fno-lto> for C++ only
	# $<$<COMPILE_LANGUAGE:C>:-fno-lto> for C only
	target_compile_options(radia_ngsolve PRIVATE
		$<$<COMPILE_LANGUAGE:CXX>:-fno-lto>
		$<$<COMPILE_LANGUAGE:C>:-fno-lto>
	)
	message(STATUS "IPO/LTO disabled for radia_ngsolve (Intel compiler)")
endif()

# Link triangle library
target_link_libraries(radia_ngsolve PRIVATE triangle_lib)

# Add include directories
target_include_directories(radia_ngsolve PRIVATE
	${CORE_DIR}
	${LIB_DIR}
	${EXT_DIR}/auxparse
	${EXT_DIR}/genmath
	${EXT_DIR}/triangle
	${EXT_DIR}/HACApK
	${CMAKE_CURRENT_SOURCE_DIR}/external/exafmm-t/include
	${PYTHON_DIR}
	${MKL_FFTW3_INCLUDE_DIR}
)

# Add compile definitions (from Radia)
target_compile_definitions(radia_ngsolve PRIVATE
	ALPHA__LIB__
	FFTW_ENABLE_FLOAT
	NO_TIMER
	ANSI_DECLARATORS
	TRILIBRARY
	_GM_WITHOUT_BASE
	RADIA_USE_HACAPK
	HACAPK_RADIA_LOGGING
	$<$<BOOL:${RADIA_ENABLE_EXAFMM}>:RADIA_USE_EXAFMM>
	$<$<CONFIG:Debug>:_DEBUG>
	$<$<CONFIG:Release>:NDEBUG>
	$<$<PLATFORM_ID:Windows>:WIN32>
	$<$<PLATFORM_ID:Windows>:_WINDOWS>
)

# FFT functionality provided by Intel MKL DFTI (via fftw.h compatibility layer)
# No external FFTW library required - MKL is already linked for BLAS/LAPACK
message(STATUS "Using Intel MKL DFTI for FFT operations (radia_ngsolve)")

# Intel MKL for radia_ngsolve (required)
if(WIN32 AND USE_MKL_FOUND)
	target_link_libraries(radia_ngsolve PRIVATE ${BLAS_LIB})
	target_include_directories(radia_ngsolve PRIVATE ${BLAS_INCLUDE_DIR})
	target_compile_definitions(radia_ngsolve PRIVATE HAVE_LAPACK)
	message(STATUS "Using Intel MKL for radia_ngsolve: ${BLAS_LIB}")
endif()

# ExaFMM is header-only, no library linking needed
if(RADIA_ENABLE_EXAFMM)
	message(STATUS "ExaFMM enabled for radia_ngsolve (header-only library)")
endif()

# Link Windows system libraries
if(WIN32)
	target_link_libraries(radia_ngsolve PRIVATE version.lib)
	# Fix LIBCMT conflict warning
	if(MSVC)
		target_link_options(radia_ngsolve PRIVATE /NODEFAULTLIB:LIBCMT)
		# Increase stack size to 16MB for deep HACApK cluster tree recursion
		# Default is 1MB which causes stack overflow at ~30,000+ elements
		target_link_options(radia_ngsolve PRIVATE /STACK:16777216)
	endif()
endif()

# Link OpenMP for radia_ngsolve
# Use Intel OpenMP (libiomp5md.dll) for better performance with Intel MKL
if(RADIA_ENABLE_OPENMP)
	if(WIN32 AND MSVC AND NOT CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
		# MSVC: Link Intel OpenMP directly instead of using OpenMP::OpenMP_CXX
		# OpenMP::OpenMP_CXX would link to vcomp140.dll (MSVC OpenMP)
		set(INTEL_OPENMP_LIB "${INTEL_COMPILER_DIR}/lib/libiomp5md.lib")
		if(EXISTS ${INTEL_OPENMP_LIB})
			target_link_libraries(radia_ngsolve PRIVATE ${INTEL_OPENMP_LIB})
			target_compile_definitions(radia_ngsolve PRIVATE _OPENMP)
			message(STATUS "OpenMP enabled for radia_ngsolve (Intel OpenMP: ${INTEL_OPENMP_LIB})")
		elseif(OpenMP_CXX_FOUND)
			target_link_libraries(radia_ngsolve PRIVATE OpenMP::OpenMP_CXX)
			message(STATUS "OpenMP enabled for radia_ngsolve (MSVC OpenMP - fallback)")
		endif()
	elseif(OpenMP_CXX_FOUND)
		target_link_libraries(radia_ngsolve PRIVATE OpenMP::OpenMP_CXX)
		message(STATUS "OpenMP enabled for radia_ngsolve")
	endif()
endif()

# Compiler-specific flags for radia_ngsolve
# MSVC is preferred for radia_ngsolve due to NGSolve library compatibility
if(MSVC)
	# MSVC with Intel MKL (required for NGSolve compatibility)
	target_compile_options(radia_ngsolve PRIVATE
		/O2                # Optimization
		/fp:fast           # Fast floating-point
		/W0                # Disable warnings
		/wd4244            # Disable conversion warnings
		/wd4267            # Disable size_t conversion warnings
	)
	if(RADIA_ENABLE_OPENMP)
		# Use standard /openmp flag - Intel OpenMP runtime is linked separately
		# The libiomp5md.lib is linked via target_link_libraries above
		target_compile_options(radia_ngsolve PRIVATE /openmp)
	endif()
	message(STATUS "Using MSVC compiler with Intel MKL for radia_ngsolve")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
	# Intel oneAPI compiler - may have ABI issues with MSVC-compiled NGSolve
	set(RADIA_NGSOLVE_INTEL_FLAGS
		/O3                # Maximum optimization
		/fp:fast           # Fast floating-point model
		/W0                # Disable warnings
	)
	if(RADIA_ENABLE_OPENMP)
		list(APPEND RADIA_NGSOLVE_INTEL_FLAGS /Qopenmp)
	endif()
	target_compile_options(radia_ngsolve PRIVATE ${RADIA_NGSOLVE_INTEL_FLAGS})
	if(RADIA_ENABLE_OPENMP)
		target_link_options(radia_ngsolve PRIVATE /Qopenmp)
	endif()
	message(STATUS "Using Intel compiler for radia_ngsolve (may have ABI issues)")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	# GCC
	target_compile_options(radia_ngsolve PRIVATE
		-Wno-conversion
		-Wno-sign-conversion
		-Wno-unused-variable
	)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	# Clang (including Apple Clang)
	target_compile_options(radia_ngsolve PRIVATE
		-Wno-conversion
		-Wno-sign-conversion
		-Wno-unused-variable
	)
endif()

message(STATUS "radia_ngsolve module configured")
endif() # RADIA_BUILD_NGSOLVE

# ========================================
# Legacy _radia Python Module (REMOVED)
# ========================================
# The legacy Python C API module (_radia.pyd using radpy_pyapi.cpp) has been
# removed as of 2026-01. All Python bindings now use pybind11 exclusively
# via the _radia_pybind module (see below at line ~956).
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
# 	triangle_lib
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
# 	${EXT_DIR}/triangle
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
# 	TRILIBRARY
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

# Link OpenMP for parallel matrix construction
if(RADIA_ENABLE_OPENMP)
	if(WIN32 AND MSVC AND NOT CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
		set(INTEL_OPENMP_LIB "${INTEL_COMPILER_DIR}/lib/libiomp5md.lib")
		if(EXISTS ${INTEL_OPENMP_LIB})
			target_link_libraries(peec_matrices PRIVATE ${INTEL_OPENMP_LIB})
			target_compile_definitions(peec_matrices PRIVATE _OPENMP)
			message(STATUS "OpenMP enabled for peec_matrices (Intel OpenMP)")
		elseif(OpenMP_CXX_FOUND)
			target_link_libraries(peec_matrices PRIVATE OpenMP::OpenMP_CXX)
			message(STATUS "OpenMP enabled for peec_matrices (MSVC OpenMP)")
		endif()
	elseif(OpenMP_CXX_FOUND)
		target_link_libraries(peec_matrices PRIVATE OpenMP::OpenMP_CXX)
		message(STATUS "OpenMP enabled for peec_matrices")
	endif()
endif()

# Compiler-specific flags
if(MSVC)
	target_compile_options(peec_matrices PRIVATE
		/O2
		/fp:fast
		/W0
		/wd4244
		/wd4267
	)
	if(RADIA_ENABLE_OPENMP)
		target_compile_options(peec_matrices PRIVATE /openmp)
	endif()
	message(STATUS "Using MSVC compiler for peec_matrices")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
	target_compile_options(peec_matrices PRIVATE /O3 /fp:fast /W0)
	if(RADIA_ENABLE_OPENMP)
		target_compile_options(peec_matrices PRIVATE /Qopenmp)
		target_link_options(peec_matrices PRIVATE /Qopenmp)
	endif()
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
			"${MKL_BIN_DIR}/mkl_rt.2.dll"
			"$<TARGET_FILE_DIR:peec_matrices>"
		COMMENT "Copying mkl_rt.2.dll for peec_matrices"
	)
	if(EXISTS "${INTEL_COMPILER_BIN}/libiomp5md.dll")
		add_custom_command(TARGET peec_matrices POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy_if_different
				"${INTEL_COMPILER_BIN}/libiomp5md.dll"
				"$<TARGET_FILE_DIR:peec_matrices>"
			COMMENT "Copying libiomp5md.dll for peec_matrices"
		)
	endif()
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

# Link OpenMP for parallel frequency sweep
if(RADIA_ENABLE_OPENMP)
	if(WIN32 AND MSVC AND NOT CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
		set(INTEL_OPENMP_LIB "${INTEL_COMPILER_DIR}/lib/libiomp5md.lib")
		if(EXISTS ${INTEL_OPENMP_LIB})
			target_link_libraries(cln_core PRIVATE ${INTEL_OPENMP_LIB})
			target_compile_definitions(cln_core PRIVATE _OPENMP)
			message(STATUS "OpenMP enabled for cln_core (Intel OpenMP)")
		elseif(OpenMP_CXX_FOUND)
			target_link_libraries(cln_core PRIVATE OpenMP::OpenMP_CXX)
			message(STATUS "OpenMP enabled for cln_core (MSVC OpenMP)")
		endif()
	elseif(OpenMP_CXX_FOUND)
		target_link_libraries(cln_core PRIVATE OpenMP::OpenMP_CXX)
		message(STATUS "OpenMP enabled for cln_core")
	endif()
endif()

# Compiler-specific flags
if(MSVC)
	target_compile_options(cln_core PRIVATE
		/O2
		/fp:fast
		/W0
		/wd4244
		/wd4267
	)
	if(RADIA_ENABLE_OPENMP)
		target_compile_options(cln_core PRIVATE /openmp)
	endif()
	message(STATUS "Using MSVC compiler for cln_core")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
	target_compile_options(cln_core PRIVATE /O3 /fp:fast /W0)
	if(RADIA_ENABLE_OPENMP)
		target_compile_options(cln_core PRIVATE /Qopenmp)
		target_link_options(cln_core PRIVATE /Qopenmp)
	endif()
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
			"${MKL_BIN_DIR}/mkl_rt.2.dll"
			"$<TARGET_FILE_DIR:cln_core>"
		COMMENT "Copying mkl_rt.2.dll for cln_core"
	)
	if(EXISTS "${INTEL_COMPILER_BIN}/libiomp5md.dll")
		add_custom_command(TARGET cln_core POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy_if_different
				"${INTEL_COMPILER_BIN}/libiomp5md.dll"
				"$<TARGET_FILE_DIR:cln_core>"
			COMMENT "Copying libiomp5md.dll for cln_core"
		)
	endif()
endif()

message(STATUS "cln_core module configured")

# ========================================
# MMM Core Python Module (mmm_core.pyd)
# ========================================
# Standalone MMM (Magnetic Moment Method) matrix construction
# for NGSolve mesh integration and external solvers

set(MMM_CORE_SOURCES
	${LIB_DIR}/rad_mmm_matrices_api.cpp
	${CORE_DIR}/rad_mmm_matrices.cpp
)

pybind11_add_module(mmm_core
	${MMM_CORE_SOURCES}
)

# Add include directories
target_include_directories(mmm_core PRIVATE
	${CORE_DIR}
)

# Add compile definitions
target_compile_definitions(mmm_core PRIVATE
	_CRT_SECURE_NO_WARNINGS
	$<$<CONFIG:Debug>:_DEBUG>
	$<$<CONFIG:Release>:NDEBUG>
	$<$<PLATFORM_ID:Windows>:WIN32>
)

# Intel MKL for LAPACK operations (required for LU solver)
if(WIN32 AND USE_MKL_FOUND)
	target_link_libraries(mmm_core PRIVATE ${BLAS_LIB})
	target_include_directories(mmm_core PRIVATE ${BLAS_INCLUDE_DIR})
	target_compile_definitions(mmm_core PRIVATE HAVE_LAPACK)
	message(STATUS "Using Intel MKL for mmm_core: ${BLAS_LIB}")
endif()

# Link Windows system libraries
if(WIN32)
	target_link_libraries(mmm_core PRIVATE version.lib)
	if(MSVC)
		target_link_options(mmm_core PRIVATE /NODEFAULTLIB:LIBCMT)
	endif()
endif()

# Link OpenMP for parallel matrix construction
if(RADIA_ENABLE_OPENMP)
	if(WIN32 AND MSVC AND NOT CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
		set(INTEL_OPENMP_LIB "${INTEL_COMPILER_DIR}/lib/libiomp5md.lib")
		if(EXISTS ${INTEL_OPENMP_LIB})
			target_link_libraries(mmm_core PRIVATE ${INTEL_OPENMP_LIB})
			target_compile_definitions(mmm_core PRIVATE _OPENMP)
			message(STATUS "OpenMP enabled for mmm_core (Intel OpenMP)")
		elseif(OpenMP_CXX_FOUND)
			target_link_libraries(mmm_core PRIVATE OpenMP::OpenMP_CXX)
			message(STATUS "OpenMP enabled for mmm_core (MSVC OpenMP)")
		endif()
	elseif(OpenMP_CXX_FOUND)
		target_link_libraries(mmm_core PRIVATE OpenMP::OpenMP_CXX)
		message(STATUS "OpenMP enabled for mmm_core")
	endif()
endif()

# Compiler-specific flags
if(MSVC)
	target_compile_options(mmm_core PRIVATE
		/O2
		/fp:fast
		/W0
		/wd4244
		/wd4267
	)
	if(RADIA_ENABLE_OPENMP)
		target_compile_options(mmm_core PRIVATE /openmp)
	endif()
	message(STATUS "Using MSVC compiler for mmm_core")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
	target_compile_options(mmm_core PRIVATE /O3 /fp:fast /W0)
	if(RADIA_ENABLE_OPENMP)
		target_compile_options(mmm_core PRIVATE /Qopenmp)
		target_link_options(mmm_core PRIVATE /Qopenmp)
	endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	target_compile_options(mmm_core PRIVATE
		-Wno-conversion
		-Wno-sign-conversion
	)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	target_compile_options(mmm_core PRIVATE
		-Wno-conversion
		-Wno-sign-conversion
	)
endif()

# Post-build: Copy required DLLs
if(WIN32 AND USE_MKL_FOUND)
	add_custom_command(TARGET mmm_core POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			"${MKL_BIN_DIR}/mkl_rt.2.dll"
			"$<TARGET_FILE_DIR:mmm_core>"
		COMMENT "Copying mkl_rt.2.dll for mmm_core"
	)
	if(EXISTS "${INTEL_COMPILER_BIN}/libiomp5md.dll")
		add_custom_command(TARGET mmm_core POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy_if_different
				"${INTEL_COMPILER_BIN}/libiomp5md.dll"
				"$<TARGET_FILE_DIR:mmm_core>"
			COMMENT "Copying libiomp5md.dll for mmm_core"
		)
	endif()
endif()

message(STATUS "mmm_core module configured")

# ========================================
# Radia pybind11 Module (_radia_pybind.pyd)
# ========================================
# pybind11 bindings following NGSolve patterns
# This is the sole Python binding (radpy_pyapi.cpp removed 2026-01)

set(RADIA_PYBIND_SOURCES
	${PYTHON_DIR}/radia_pybind.cpp
)

pybind11_add_module(_radia_pybind
	${RADIA_CORE_SOURCES}

	${EXT_DIR}/auxparse/auxparse.cpp
	${EXT_DIR}/genmath/gmfft.cpp
	${EXT_DIR}/genmath/gmtrans.cpp
	${HACAPK_SOURCES}
	${RADIA_LIB_SOURCES}
	${RADIA_PYBIND_SOURCES}
)

# Link triangle library
target_link_libraries(_radia_pybind PRIVATE triangle_lib)

# Add include directories
target_include_directories(_radia_pybind PRIVATE
	${CORE_DIR}
	${LIB_DIR}
	${EXT_DIR}/auxparse
	${EXT_DIR}/genmath
	${EXT_DIR}/HACApK
	${CMAKE_CURRENT_SOURCE_DIR}/external/exafmm-t/include
	${MKL_FFTW3_INCLUDE_DIR}
)

# Add preprocessor definitions
target_compile_definitions(_radia_pybind PRIVATE
	_CRT_SECURE_NO_WARNINGS
	_GM_WITHOUT_BASE
	ALPHA__DLL__
	ANSI_DECLARATORS
	RADIA_USE_HACAPK
	HACAPK_RADIA_LOGGING
	$<$<BOOL:${RADIA_ENABLE_EXAFMM}>:RADIA_USE_EXAFMM>
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
	endif()
endif()

# Link OpenMP for parallel operations
if(RADIA_ENABLE_OPENMP)
	if(WIN32 AND MSVC AND NOT CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
		set(INTEL_OPENMP_LIB "${INTEL_COMPILER_DIR}/lib/libiomp5md.lib")
		if(EXISTS ${INTEL_OPENMP_LIB})
			target_link_libraries(_radia_pybind PRIVATE ${INTEL_OPENMP_LIB})
			target_compile_definitions(_radia_pybind PRIVATE _OPENMP)
			message(STATUS "OpenMP enabled for _radia_pybind (Intel OpenMP)")
		elseif(OpenMP_CXX_FOUND)
			target_link_libraries(_radia_pybind PRIVATE OpenMP::OpenMP_CXX)
			message(STATUS "OpenMP enabled for _radia_pybind (MSVC OpenMP)")
		endif()
	elseif(OpenMP_CXX_FOUND)
		target_link_libraries(_radia_pybind PRIVATE OpenMP::OpenMP_CXX)
		message(STATUS "OpenMP enabled for _radia_pybind")
	endif()
endif()

# Compiler-specific flags
if(MSVC)
	target_compile_options(_radia_pybind PRIVATE
		/O2
		/fp:fast
		/W0
		/wd4244
		/wd4267
	)
	if(RADIA_ENABLE_OPENMP)
		target_compile_options(_radia_pybind PRIVATE /openmp)
	endif()
elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
	target_compile_options(_radia_pybind PRIVATE /O3 /fp:fast /W0)
	if(RADIA_ENABLE_OPENMP)
		target_compile_options(_radia_pybind PRIVATE /Qopenmp)
		target_link_options(_radia_pybind PRIVATE /Qopenmp)
	endif()
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
