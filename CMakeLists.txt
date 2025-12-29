cmake_minimum_required(VERSION 3.21)
project(radia_ngsolve C CXX)

# Option to disable OpenMP (for debugging)
option(RADIA_ENABLE_OPENMP "Enable OpenMP parallelization" ON)

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

# Set NGSolve and Netgen directories
set(NGSolve_DIR ${PYTHON_SITE_PACKAGES}/ngsolve/cmake)
set(Netgen_DIR ${PYTHON_SITE_PACKAGES}/netgen/cmake)

# Find NGSolve
find_package(NGSolve CONFIG REQUIRED)

# Radia source directories
set(CORE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/core")
set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/lib")
set(EXT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/ext")
set(PYTHON_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/radia")
set(FFTW_DIR "${EXT_DIR}/fftw")


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
		set(USE_MKL_FOUND TRUE)
	else()
		message(FATAL_ERROR "Intel MKL not found. Please install Intel oneAPI Base Toolkit.")
	endif()
endif()

# Option to enable FMM3D acceleration
# Note: FMM3D library requires Fortran compilation - temporarily OFF until library is built
option(RADIA_ENABLE_FMM3D "Enable FMM3D library for fast field computation" OFF)

# Radia core sources (organized by API functionality)
set(RADIA_CORE_SOURCES
	# API-organized files (new naming convention)
	${CORE_DIR}/rad_obj_creation.cpp          # Object creation (ObjRecMag, ObjCylMag, etc.)
	${CORE_DIR}/rad_material_impl.cpp         # Materials (MatLin, MatSatIsoFrm, MatApl, etc.)
	${CORE_DIR}/rad_transform_impl.cpp        # Transformations (TrfTrsl, TrfRot, etc.)
	${CORE_DIR}/rad_field_compute.cpp         # Field computation (Fld, FldLst, FldEnr, etc.)
	${CORE_DIR}/rad_obj_subdivision.cpp       # Subdivision (ObjDivMag, multi-generation extrusion)
	${CORE_DIR}/rad_fmm3d.cpp                 # FMM3D wrapper for fast multipole acceleration
	${CORE_DIR}/rad_point_classify.cpp        # Point classification for FMM near/far

	# Geometry files (descriptive names)
	${CORE_DIR}/rad_arc_current.cpp           # Arc with azimuthal current
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
	${CORE_DIR}/rad_extruded_polygon_msc.cpp  # MSC extruded polygon (6 DOF hexahedra)
	${CORE_DIR}/rad_relaxation_methods.cpp    # Relaxation methods
	${CORE_DIR}/rad_subdivided_arc_current.cpp        # Subdivided arc with current
	${CORE_DIR}/rad_subdivided_extruded_polygon.cpp   # Subdivided extruded polygon
	${CORE_DIR}/rad_subdivided_rectangle.cpp          # Subdivided rectangle
	${CORE_DIR}/rad_subdivided_polyhedron.cpp         # Subdivided polyhedron
	${CORE_DIR}/rad_serialization.cpp         # Interface functions
	${CORE_DIR}/rad_polyhedron.cpp            # Polyhedron with magnetization
	${CORE_DIR}/rad_hacapk.cpp                # HACApK H-matrix manager (Method 2)
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

# Radia lib sources
set(RADIA_LIB_SOURCES
	${LIB_DIR}/radentry.cpp
)

# Triangle library (C source needs special handling)
# Create as static library first
add_library(triangle_lib STATIC ${EXT_DIR}/triangle/triangle.c)
target_compile_definitions(triangle_lib PRIVATE TRILIBRARY ANSI_DECLARATORS NO_TIMER)
target_include_directories(triangle_lib PUBLIC ${EXT_DIR}/triangle)

# Create NGSolve Python module using add_ngsolve_python_module
# (following EMPY_Field pattern exactly)
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
	${EXT_DIR}/fmm3d/c
	${PYTHON_DIR}
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
	RADIA_MSC_SUPPORT
	$<$<BOOL:${RADIA_ENABLE_FMM3D}>:RADIA_FMM3D_ENABLED>
	$<$<CONFIG:Debug>:_DEBUG>
	$<$<CONFIG:Release>:NDEBUG>
	$<$<PLATFORM_ID:Windows>:WIN32>
	$<$<PLATFORM_ID:Windows>:_WINDOWS>
)

# FFTW library - platform-specific detection
if(WIN32)
	# Windows: Look for precompiled FFTW library
	set(FFTW_LIB ${FFTW_DIR}/fftw64_f.lib)
	if(EXISTS ${FFTW_LIB})
		target_link_libraries(radia_ngsolve PRIVATE ${FFTW_LIB})
		message(STATUS "Found FFTW library: ${FFTW_LIB}")
	else()
		message(WARNING "FFTW library not found: ${FFTW_LIB}")
	endif()
elseif(APPLE)
	# macOS: Try to find FFTW via pkg-config or Homebrew
	find_library(FFTW_LIBRARY NAMES fftw3f fftw3
		PATHS /usr/local/lib /opt/homebrew/lib
		NO_DEFAULT_PATH
	)
	if(FFTW_LIBRARY)
		target_link_libraries(radia_ngsolve PRIVATE ${FFTW_LIBRARY})
		find_path(FFTW_INCLUDE_DIR fftw3.h
			PATHS /usr/local/include /opt/homebrew/include
			NO_DEFAULT_PATH
		)
		if(FFTW_INCLUDE_DIR)
			target_include_directories(radia_ngsolve PRIVATE ${FFTW_INCLUDE_DIR})
		endif()
		message(STATUS "Found FFTW library: ${FFTW_LIBRARY}")
	else()
		message(WARNING "FFTW library not found. Install with: brew install fftw")
	endif()
else()
	# Linux: Try to find FFTW via pkg-config
	find_library(FFTW_LIBRARY NAMES fftw3f fftw3
		PATHS /usr/lib /usr/local/lib /usr/lib/x86_64-linux-gnu
	)
	if(FFTW_LIBRARY)
		target_link_libraries(radia_ngsolve PRIVATE ${FFTW_LIBRARY})
		find_path(FFTW_INCLUDE_DIR fftw3.h
			PATHS /usr/include /usr/local/include
		)
		if(FFTW_INCLUDE_DIR)
			target_include_directories(radia_ngsolve PRIVATE ${FFTW_INCLUDE_DIR})
		endif()
		message(STATUS "Found FFTW library: ${FFTW_LIBRARY}")
	else()
		message(WARNING "FFTW library not found. Install with: sudo apt-get install libfftw3-dev")
	endif()
endif()

# Intel MKL for radia_ngsolve (required)
if(WIN32 AND USE_MKL_FOUND)
	target_link_libraries(radia_ngsolve PRIVATE ${BLAS_LIB})
	target_include_directories(radia_ngsolve PRIVATE ${BLAS_INCLUDE_DIR})
	target_compile_definitions(radia_ngsolve PRIVATE HAVE_LAPACK)
	message(STATUS "Using Intel MKL for radia_ngsolve: ${BLAS_LIB}")
endif()

# FMM3D library for fast multipole method (optional)
if(RADIA_ENABLE_FMM3D)
	set(FMM3D_LIB "${EXT_DIR}/fmm3d/lib/fmm3d.lib")
	if(WIN32 AND EXISTS ${FMM3D_LIB})
		target_link_libraries(radia_ngsolve PRIVATE ${FMM3D_LIB})
		message(STATUS "FMM3D enabled for radia_ngsolve: ${FMM3D_LIB}")
	else()
		message(STATUS "FMM3D library not found for radia_ngsolve, using direct computation fallback")
	endif()
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

# ========================================
# Radia Python Module (radia.pyd)
# ========================================

# Create radia Python module using pybind11
pybind11_add_module(radia
	${RADIA_CORE_SOURCES}
	${EXT_DIR}/auxparse/auxparse.cpp
	${EXT_DIR}/genmath/gmfft.cpp
	${EXT_DIR}/genmath/gmtrans.cpp
	${HACAPK_SOURCES}
	${RADIA_LIB_SOURCES}
	${PYTHON_DIR}/radpy_pyapi.cpp
)

# Link triangle library
target_link_libraries(radia PRIVATE triangle_lib)

# Add include directories
target_include_directories(radia PRIVATE
	${CORE_DIR}
	${LIB_DIR}
	${EXT_DIR}/auxparse
	${EXT_DIR}/genmath
	${EXT_DIR}/HACApK
	${EXT_DIR}/fmm3d/c
	${FFTW_DIR}
)

# Add preprocessor definitions
target_compile_definitions(radia PRIVATE
	_CRT_SECURE_NO_WARNINGS
	_GM_WITHOUT_BASE
	ALPHA__DLL__
	ANSI_DECLARATORS
	RADIA_USE_HACAPK
	HACAPK_RADIA_LOGGING
	RADIA_MSC_SUPPORT
	$<$<BOOL:${RADIA_ENABLE_FMM3D}>:RADIA_FMM3D_ENABLED>
)

# FFTW library - platform-specific detection (same as radia_ngsolve)
if(WIN32)
	# Windows: Look for precompiled FFTW library
	if(EXISTS ${FFTW_LIB})
		target_link_libraries(radia PRIVATE ${FFTW_LIB})
		message(STATUS "Found FFTW library for radia: ${FFTW_LIB}")
	else()
		message(WARNING "FFTW library not found for radia: ${FFTW_LIB}")
	endif()
elseif(APPLE)
	# macOS: Try to find FFTW via pkg-config or Homebrew
	if(FFTW_LIBRARY)
		target_link_libraries(radia PRIVATE ${FFTW_LIBRARY})
		if(FFTW_INCLUDE_DIR)
			target_include_directories(radia PRIVATE ${FFTW_INCLUDE_DIR})
		endif()
		message(STATUS "Found FFTW library for radia: ${FFTW_LIBRARY}")
	else()
		message(WARNING "FFTW library not found for radia. Install with: brew install fftw")
	endif()
else()
	# Linux: Try to find FFTW via pkg-config
	if(FFTW_LIBRARY)
		target_link_libraries(radia PRIVATE ${FFTW_LIBRARY})
		if(FFTW_INCLUDE_DIR)
			target_include_directories(radia PRIVATE ${FFTW_INCLUDE_DIR})
		endif()
		message(STATUS "Found FFTW library for radia: ${FFTW_LIBRARY}")
	else()
		message(WARNING "FFTW library not found for radia. Install with: sudo apt-get install libfftw3-dev")
	endif()
endif()

# Intel MKL for radia (required for fast LU decomposition and matrix operations)
if(WIN32 AND USE_MKL_FOUND)
	target_link_libraries(radia PRIVATE ${BLAS_LIB})
	target_include_directories(radia PRIVATE ${BLAS_INCLUDE_DIR})
	target_compile_definitions(radia PRIVATE HAVE_LAPACK)
	message(STATUS "Using Intel MKL for radia: ${BLAS_LIB}")
endif()

# FMM3D library for fast multipole method (optional)
if(RADIA_ENABLE_FMM3D)
	set(FMM3D_LIB "${EXT_DIR}/fmm3d/lib/fmm3d.lib")
	if(WIN32 AND EXISTS ${FMM3D_LIB})
		target_link_libraries(radia PRIVATE ${FMM3D_LIB})
		# FMM3D requires Intel Fortran runtime
		# These are linked via MKL and Intel compiler runtime
		message(STATUS "FMM3D enabled for radia: ${FMM3D_LIB}")
	else()
		message(STATUS "FMM3D library not found, using direct computation fallback")
	endif()
endif()

# Link Windows system libraries
if(WIN32)
	target_link_libraries(radia PRIVATE version.lib)
	# Fix LIBCMT conflict warning
	if(MSVC)
		target_link_options(radia PRIVATE /NODEFAULTLIB:LIBCMT)
		# Increase stack size to 16MB for deep HACApK cluster tree recursion
		# Default is 1MB which causes stack overflow at ~30,000+ elements
		target_link_options(radia PRIVATE /STACK:16777216)
	endif()
endif()

# Link OpenMP for parallel matrix construction
# Use Intel OpenMP (libiomp5md.dll) for better performance with Intel MKL
if(RADIA_ENABLE_OPENMP)
	if(WIN32 AND MSVC AND NOT CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
		# MSVC: Link Intel OpenMP directly instead of using OpenMP::OpenMP_CXX
		# OpenMP::OpenMP_CXX would link to vcomp140.dll (MSVC OpenMP)
		set(INTEL_OPENMP_LIB "${INTEL_COMPILER_DIR}/lib/libiomp5md.lib")
		if(EXISTS ${INTEL_OPENMP_LIB})
			target_link_libraries(radia PRIVATE ${INTEL_OPENMP_LIB})
			target_compile_definitions(radia PRIVATE _OPENMP)
			message(STATUS "OpenMP enabled for radia (Intel OpenMP: ${INTEL_OPENMP_LIB})")
		elseif(OpenMP_CXX_FOUND)
			target_link_libraries(radia PRIVATE OpenMP::OpenMP_CXX)
			message(STATUS "OpenMP enabled for radia (MSVC OpenMP - fallback)")
		endif()
	else()
		if(OpenMP_CXX_FOUND)
			target_link_libraries(radia PRIVATE OpenMP::OpenMP_CXX)
			message(STATUS "OpenMP enabled for radia (C++)")
		endif()
		if(OpenMP_C_FOUND)
			target_link_libraries(radia PRIVATE OpenMP::OpenMP_C)
			message(STATUS "OpenMP enabled for radia (C)")
		endif()
	endif()
endif()

# Compiler-specific flags for radia
# MSVC with Intel MKL or Intel oneAPI compiler
if(MSVC)
	# MSVC with Intel MKL
	target_compile_options(radia PRIVATE
		/O2                # Optimization
		/fp:fast           # Fast floating-point
		/W0                # Disable warnings
		/wd4244            # Disable conversion warnings
		/wd4267            # Disable size_t conversion warnings
	)
	if(RADIA_ENABLE_OPENMP)
		# Use standard /openmp flag - Intel OpenMP runtime is linked separately
		# The libiomp5md.lib is linked via target_link_libraries above
		target_compile_options(radia PRIVATE /openmp)
	endif()
	message(STATUS "Using MSVC compiler with Intel MKL for radia")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "IntelLLVM")
	# Intel oneAPI compiler (icx-cl) - optimal performance
	set(RADIA_INTEL_FLAGS
		/O3                # Maximum optimization
		/Qipo              # Interprocedural optimization
		/fp:fast           # Fast floating-point model
		/Qopt-report:2     # Optimization report
		/W0                # Disable warnings
	)
	if(RADIA_ENABLE_OPENMP)
		list(APPEND RADIA_INTEL_FLAGS /Qopenmp)  # Enable OpenMP (Intel style)
	endif()
	target_compile_options(radia PRIVATE ${RADIA_INTEL_FLAGS})
	# Link with Intel OpenMP
	if(RADIA_ENABLE_OPENMP)
		target_link_options(radia PRIVATE /Qopenmp)
	endif()
	message(STATUS "Using Intel oneAPI compiler for radia (optimal performance)")
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	# GCC
	target_compile_options(radia PRIVATE
		-Wno-conversion
		-Wno-sign-conversion
		-Wno-unused-variable
	)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	# Clang (including Apple Clang)
	target_compile_options(radia PRIVATE
		-Wno-conversion
		-Wno-sign-conversion
		-Wno-unused-variable
	)
endif()

message(STATUS "radia module configured")

# Post-build: Copy required DLLs for distribution
if(WIN32 AND USE_MKL_FOUND)
	# Copy MKL runtime DLL
	add_custom_command(TARGET radia POST_BUILD
		COMMAND ${CMAKE_COMMAND} -E copy_if_different
			"${MKL_BIN_DIR}/mkl_rt.2.dll"
			"$<TARGET_FILE_DIR:radia>"
		COMMENT "Copying mkl_rt.2.dll for distribution"
	)

	# Copy Intel OpenMP runtime if available
	if(EXISTS "${INTEL_COMPILER_BIN}/libiomp5md.dll")
		add_custom_command(TARGET radia POST_BUILD
			COMMAND ${CMAKE_COMMAND} -E copy_if_different
				"${INTEL_COMPILER_BIN}/libiomp5md.dll"
				"$<TARGET_FILE_DIR:radia>"
			COMMENT "Copying libiomp5md.dll for distribution"
		)
	endif()
endif()

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
