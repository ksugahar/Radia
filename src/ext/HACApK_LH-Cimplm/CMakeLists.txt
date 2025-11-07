cmake_minimum_required(VERSION 3.21)
project(HACApK_CPP CXX)

# C++17 standard (same as Radia)
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# ========================================
# HACApK C++ Library
# ========================================

# HACApK source files (C++ only, no Fortran)
set(HACAPK_SOURCES
	hacapk.cpp
	hacapk.hpp
)

# Create HACApK static library
add_library(hacapk STATIC ${HACAPK_SOURCES})

# Include directories
target_include_directories(hacapk PUBLIC
	${CMAKE_CURRENT_SOURCE_DIR}
)

# Find and link OpenMP
find_package(OpenMP REQUIRED)

# Compiler-specific flags
if(MSVC)
	# MSVC (Visual Studio)
	target_compile_options(hacapk PRIVATE
		/W3          # Warning level 3
		/wd4244      # Disable conversion warnings
		/wd4267      # Disable size_t conversion warnings
	)
	target_link_libraries(hacapk PUBLIC OpenMP::OpenMP_CXX)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	# GCC
	target_compile_options(hacapk PRIVATE
		-Wall
		-Wextra
		-Wno-unused-parameter
		-Wno-unused-variable
	)
	target_link_libraries(hacapk PUBLIC OpenMP::OpenMP_CXX)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	# Clang
	target_compile_options(hacapk PRIVATE
		-Wall
		-Wextra
		-Wno-unused-parameter
		-Wno-unused-variable
	)
	target_link_libraries(hacapk PUBLIC OpenMP::OpenMP_CXX)
endif()

# Link math library (required for sqrt, etc.)
if(UNIX)
	target_link_libraries(hacapk PUBLIC m)
endif()

# ========================================
# Test Executable
# ========================================

add_executable(test_hacapk test_hacapk.cpp)
target_link_libraries(test_hacapk PRIVATE hacapk)
target_include_directories(test_hacapk PRIVATE
	${CMAKE_CURRENT_SOURCE_DIR}
)

# Compiler-specific flags for test
if(MSVC)
	target_compile_options(test_hacapk PRIVATE
		/W4          # Warning level 4
		/wd4244      # Disable conversion warnings
		/wd4267      # Disable size_t conversion warnings
		/D_USE_MATH_DEFINES  # Define M_PI
	)
	target_link_libraries(test_hacapk PRIVATE OpenMP::OpenMP_CXX)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "GNU")
	target_compile_options(test_hacapk PRIVATE
		-Wall
		-Wextra
		-Wpedantic
		-Wno-unused-parameter
	)
	target_link_libraries(test_hacapk PRIVATE OpenMP::OpenMP_CXX)
elseif(CMAKE_CXX_COMPILER_ID MATCHES "Clang")
	target_compile_options(test_hacapk PRIVATE
		-Wall
		-Wextra
		-Wpedantic
		-Wno-unused-parameter
	)
	target_link_libraries(test_hacapk PRIVATE OpenMP::OpenMP_CXX)
endif()

# ========================================
# Testing with CTest
# ========================================

enable_testing()
add_test(NAME HACApK_CPP COMMAND test_hacapk)

# ========================================
# Installation
# ========================================

install(TARGETS hacapk test_hacapk
	ARCHIVE DESTINATION lib
	LIBRARY DESTINATION lib
	RUNTIME DESTINATION bin
)

install(FILES hacapk.hpp
	DESTINATION include
)

# ========================================
# Configuration Summary
# ========================================

message(STATUS "=============================================")
message(STATUS "HACApK C++ Configuration Summary")
message(STATUS "=============================================")
message(STATUS "CMAKE_BUILD_TYPE:     ${CMAKE_BUILD_TYPE}")
message(STATUS "CMAKE_CXX_COMPILER:   ${CMAKE_CXX_COMPILER}")
message(STATUS "CMAKE_CXX_COMPILER_ID: ${CMAKE_CXX_COMPILER_ID}")
message(STATUS "OpenMP_CXX_FLAGS:     ${OpenMP_CXX_FLAGS}")
message(STATUS "Build targets:")
message(STATUS "  - hacapk (static library)")
message(STATUS "  - test_hacapk (test executable)")
message(STATUS "=============================================")
