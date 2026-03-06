cmake_minimum_required(VERSION 3.21)
project(radia_ngsolve C CXX)

# Find OpenMP for parallel matrix construction
find_package(OpenMP)
if(OpenMP_CXX_FOUND)
    message(STATUS "OpenMP found: ${OpenMP_CXX_VERSION}")
else()
    message(WARNING "OpenMP not found - parallel matrix construction disabled")
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


# OpenBLAS directory and library (for optimized LAPACK LU solver)
set(OPENBLAS_DIR "${EXT_DIR}/openblas")
if(WIN32)
	set(OPENBLAS_LIB ${OPENBLAS_DIR}/lib/libopenblas.lib)
endif()

# Radia core sources (organized by API functionality)
set(RADIA_CORE_SOURCES
	# API-organized files (new naming convention)
	${CORE_DIR}/rad_obj_creation.cpp          # Object creation (ObjRecMag, ObjCylMag, etc.)
	${CORE_DIR}/rad_material_impl.cpp         # Materials (MatLin, MatSatIsoFrm, MatApl, etc.)
	${CORE_DIR}/rad_transform_impl.cpp        # Transformations (TrfTrsl, TrfRot, etc.)
	${CORE_DIR}/rad_field_compute.cpp         # Field computation (Fld, FldLst, FldEnr, etc.)
	${CORE_DIR}/rad_obj_subdivision.cpp       # Subdivision (ObjDivMag, multi-generation extrusion)

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
	${CORE_DIR}/rad_pyramid_msc.cpp           # MMM with MSC (pyramid elements)
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

# OpenBLAS library for radia_ngsolve (LAPACK for HACApK)
if(WIN32)
	if(EXISTS ${OPENBLAS_LIB})
		target_link_libraries(radia_ngsolve PRIVATE ${OPENBLAS_LIB})
		target_include_directories(radia_ngsolve PRIVATE ${OPENBLAS_DIR}/include)
		target_compile_definitions(radia_ngsolve PRIVATE HAVE_LAPACK)
		message(STATUS "Found OpenBLAS library for radia_ngsolve: ${OPENBLAS_LIB}")
	endif()
endif()

# Link Windows system libraries
if(WIN32)
	target_link_libraries(radia_ngsolve PRIVATE version.lib)
	# Fix LIBCMT conflict warning
	if(MSVC)
		target_link_options(radia_ngsolve PRIVATE /NODEFAULTLIB:LIBCMT)
	endif()
endif()

# Link OpenMP for radia_ngsolve
if(OpenMP_CXX_FOUND)
	target_link_libraries(radia_ngsolve PRIVATE OpenMP::OpenMP_CXX)
	message(STATUS "OpenMP enabled for radia_ngsolve")
endif()

# Compiler-specific flags
if(MSVC)
	# MSVC (Visual Studio)
	target_compile_options(radia_ngsolve PRIVATE
		/wd4244  # Disable conversion warnings
		/wd4267  # Disable size_t conversion warnings
	)
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

# OpenBLAS library (LAPACK for fast LU decomposition)
if(WIN32)
	if(EXISTS ${OPENBLAS_LIB})
		target_link_libraries(radia PRIVATE ${OPENBLAS_LIB})
		target_include_directories(radia PRIVATE ${OPENBLAS_DIR}/include)
		target_compile_definitions(radia PRIVATE HAVE_LAPACK)
		message(STATUS "Found OpenBLAS library: ${OPENBLAS_LIB}")
	else()
		message(WARNING "OpenBLAS library not found: ${OPENBLAS_LIB}")
	endif()
endif()

# Link Windows system libraries
if(WIN32)
	target_link_libraries(radia PRIVATE version.lib)
	# Fix LIBCMT conflict warning
	if(MSVC)
		target_link_options(radia PRIVATE /NODEFAULTLIB:LIBCMT)
	endif()
endif()

# Link OpenMP for parallel matrix construction
if(OpenMP_CXX_FOUND)
	target_link_libraries(radia PRIVATE OpenMP::OpenMP_CXX)
	message(STATUS "OpenMP enabled for radia (C++)")
endif()
if(OpenMP_C_FOUND)
	target_link_libraries(radia PRIVATE OpenMP::OpenMP_C)
	message(STATUS "OpenMP enabled for radia (C)")
endif()

# Compiler-specific flags
if(MSVC)
	# MSVC (Visual Studio)
	target_compile_options(radia PRIVATE
		/wd4244  # Disable conversion warnings
		/wd4267  # Disable size_t conversion warnings
	)
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
