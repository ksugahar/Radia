cmake_minimum_required(VERSION 3.21)
project(HACApK_OpenMP_Tests C CXX)

# C++17 standard (same as Radia)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# C11 standard for H-matrix implementation
set(CMAKE_C_STANDARD 11)
set(CMAKE_C_STANDARD_REQUIRED ON)

# ========================================
# OpenMP-based H-Matrix Library
# ========================================

# H-matrix source files (our OpenMP implementation)
set(HMATRIX_SOURCES
	hacapk_openmp.c
	hacapk_openmp.h
)

# Create H-matrix static library
add_library(hmatrix_openmp STATIC ${HMATRIX_SOURCES})

# Include directories
target_include_directories(hmatrix_openmp PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}
)

# Find and link OpenMP
find_package(OpenMP REQUIRED)

# Compiler-specific flags for H-matrix library
if(MSVC)
	# MSVC (Visual Studio)
	target_compile_options(hmatrix_openmp PRIVATE
		/W3          # Warning level 3
		/wd4244      # Disable conversion warnings
		/wd4267      # Disable size_t conversion warnings
	)
	target_link_libraries(hmatrix_openmp PUBLIC OpenMP::OpenMP_C)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	# GCC
	target_compile_options(hmatrix_openmp PRIVATE
		-Wall
		-Wextra
		-Wno-unused-parameter
		-Wno-unused-variable
	)
	target_link_libraries(hmatrix_openmp PUBLIC OpenMP::OpenMP_C)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	# Clang
	target_compile_options(hmatrix_openmp PRIVATE
		-Wall
		-Wextra
		-Wno-unused-parameter
		-Wno-unused-variable
	)
	target_link_libraries(hmatrix_openmp PUBLIC OpenMP::OpenMP_C)
endif()

# Link math library (required for sqrt, etc.)
if(UNIX)
	target_link_libraries(hmatrix_openmp PUBLIC m)
endif()

# ========================================
# Test Executables
# ========================================

# Test 1: Basic H-matrix functionality
add_executable(test_hmatrix_basic test_hmatrix_basic.cpp)
target_link_libraries(test_hmatrix_basic PRIVATE hmatrix_openmp)
target_include_directories(test_hmatrix_basic PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}
)

# Test 2: Radia integration with OpenMP H-matrix
add_executable(test_hmatrix_radia test_hmatrix_radia.cpp)
target_link_libraries(test_hmatrix_radia PRIVATE hmatrix_openmp)
target_include_directories(test_hmatrix_radia PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}
)

# Compiler-specific flags for test executables
if(MSVC)
	# MSVC
	target_compile_options(test_hmatrix_basic PRIVATE
		/W4          # Warning level 4
		/wd4244      # Disable conversion warnings
		/wd4267      # Disable size_t conversion warnings
		/D_USE_MATH_DEFINES  # Define M_PI
	)
	target_compile_options(test_hmatrix_radia PRIVATE
		/W4
		/wd4244
		/wd4267
		/D_USE_MATH_DEFINES  # Define M_PI
	)
	target_link_libraries(test_hmatrix_basic PRIVATE OpenMP::OpenMP_CXX)
	target_link_libraries(test_hmatrix_radia PRIVATE OpenMP::OpenMP_CXX)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	# GCC
	target_compile_options(test_hmatrix_basic PRIVATE
		-Wall
		-Wextra
		-Wpedantic
		-Wno-unused-parameter
	)
	target_compile_options(test_hmatrix_radia PRIVATE
		-Wall
		-Wextra
		-Wpedantic
		-Wno-unused-parameter
	)
	target_link_libraries(test_hmatrix_basic PRIVATE OpenMP::OpenMP_CXX)
	target_link_libraries(test_hmatrix_radia PRIVATE OpenMP::OpenMP_CXX)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	# Clang
	target_compile_options(test_hmatrix_basic PRIVATE
		-Wall
		-Wextra
		-Wpedantic
		-Wno-unused-parameter
	)
	target_compile_options(test_hmatrix_radia PRIVATE
		-Wall
		-Wextra
		-Wpedantic
		-Wno-unused-parameter
	)
	target_link_libraries(test_hmatrix_basic PRIVATE OpenMP::OpenMP_CXX)
	target_link_libraries(test_hmatrix_radia PRIVATE OpenMP::OpenMP_CXX)
endif()

# ========================================
# Installation
# ========================================

install(TARGETS test_hmatrix_basic test_hmatrix_radia
	RUNTIME DESTINATION bin
)

# ========================================
# Testing with CTest
# ========================================

enable_testing()
add_test(NAME HMatrix_Basic COMMAND test_hmatrix_basic)
add_test(NAME HMatrix_Radia COMMAND test_hmatrix_radia)

# ========================================
# Configuration Summary
# ========================================

message(STATUS "==============================================")
message(STATUS "H-Matrix OpenMP Tests Configuration Summary")
message(STATUS "==============================================")
message(STATUS "CMAKE_BUILD_TYPE:     ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_COMPILER:   ${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_C_COMPILER:     ${CMAKE_C_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "OpenMP_C_FLAGS:       ${OpenMP_C_FLAGS}")
message(STATUS "OpenMP_CXX_FLAGS:     ${OpenMP_CXX_FLAGS}")
message(STATUS "Build targets:")
message(STATUS "  - hmatrix_openmp (static library)")
message(STATUS "  - test_hmatrix_basic")
message(STATUS "  - test_hmatrix_radia")
message(STATUS "==============================================")
