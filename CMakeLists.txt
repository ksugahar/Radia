cmake_minimum_required(VERSION 3.20)

project(RadiaPython VERSION 1.0 LANGUAGES C CXX)

# Force 64-bit build
if(NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
    message(FATAL_ERROR "This project requires 64-bit build. Please use -A x64 for Visual Studio generator.")
endif()

# Set C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# Force static runtime library for MSVC
if(MSVC)
    set(CMAKE_MSVC_RUNTIME_LIBRARY "MultiThreaded$<$<CONFIG:Debug>:Debug>")
    add_compile_options(/W3 /MP)
    add_compile_definitions(_CRT_SECURE_NO_WARNINGS _CRT_SECURE_NO_DEPRECATE)
endif()

# Set output directories
set(CMAKE_LIBRARY_OUTPUT_DIRECTORY ${CMAKE_BINARY_DIR}/lib)

# Project directories
set(PROJECT_ROOT ${CMAKE_CURRENT_SOURCE_DIR})
set(SRC_DIR ${PROJECT_ROOT}/src)
set(CORE_DIR ${SRC_DIR}/core)
set(LIB_DIR ${SRC_DIR}/lib)
set(EXT_DIR ${SRC_DIR}/ext)
set(PYTHON_DIR ${SRC_DIR}/python)
set(FFTW_DIR ${EXT_DIR}/fftw)

# Find Python 3.12 (64-bit only)
find_package(Python3 3.12 EXACT COMPONENTS Interpreter Development REQUIRED)

# Find OpenMP for CPU core parallelization
find_package(OpenMP REQUIRED)

# Verify Python is 64-bit
execute_process(
    COMMAND ${Python3_EXECUTABLE} -c "import struct; print(struct.calcsize('P'))"
    OUTPUT_VARIABLE PYTHON_POINTER_SIZE
    OUTPUT_STRIP_TRAILING_WHITESPACE
)

if(NOT PYTHON_POINTER_SIZE EQUAL 8)
    message(FATAL_ERROR "Python must be 64-bit. Found: ${PYTHON_POINTER_SIZE} bytes")
endif()

# Verify Python version is exactly 3.12
if(NOT Python3_VERSION_MAJOR EQUAL 3 OR NOT Python3_VERSION_MINOR EQUAL 12)
    message(FATAL_ERROR "Python 3.12 is required. Found: ${Python3_VERSION}")
endif()

# Collect Radia core source files
set(RADIA_CORE_SOURCES
    ${CORE_DIR}/radapl1.cpp
    ${CORE_DIR}/radapl2.cpp
    ${CORE_DIR}/radapl3.cpp
    ${CORE_DIR}/radapl4.cpp
    ${CORE_DIR}/radapl5.cpp
    ${CORE_DIR}/radapl6.cpp
    ${CORE_DIR}/radapl7.cpp
    ${CORE_DIR}/radarccu.cpp
    ${CORE_DIR}/radcast.cpp
    ${CORE_DIR}/radexpgn.cpp
    ${CORE_DIR}/radflm.cpp
    ${CORE_DIR}/radg3d.cpp
    ${CORE_DIR}/radg3dgr.cpp
    ${CORE_DIR}/radgroup.cpp
    ${CORE_DIR}/radinter.cpp
    ${CORE_DIR}/radintrc.cpp
    ${CORE_DIR}/radiobuf.cpp
    ${CORE_DIR}/radmamet.cpp
    ${CORE_DIR}/radmater.cpp
    ${CORE_DIR}/radplnr1.cpp
    ${CORE_DIR}/radplnr2.cpp
    ${CORE_DIR}/radptrj.cpp
    ${CORE_DIR}/radrec.cpp
    ${CORE_DIR}/radrlmet.cpp
    ${CORE_DIR}/radsbdac.cpp
    ${CORE_DIR}/radsbdep.cpp
    ${CORE_DIR}/radsbdrc.cpp
    ${CORE_DIR}/radsbdvp.cpp
    ${CORE_DIR}/radsend.cpp
    ${CORE_DIR}/radvlpgn.cpp
)

# Extension library sources (no GLUT viewer)
set(RADIA_EXT_SOURCES
    ${EXT_DIR}/auxparse/auxparse.cpp
    ${EXT_DIR}/genmath/gmfft.cpp
    ${EXT_DIR}/genmath/gmtrans.cpp
    ${EXT_DIR}/triangle/triangle.c
)

# Library entry point
set(RADIA_LIB_SOURCES
    ${LIB_DIR}/radentry.cpp
)

# Python client source
set(RADIA_PYTHON_SOURCES
    ${PYTHON_DIR}/radpy.cpp
)

# Create Python extension module (.pyd)
add_library(radia MODULE
    ${RADIA_CORE_SOURCES}
    ${RADIA_EXT_SOURCES}
    ${RADIA_LIB_SOURCES}
    ${RADIA_PYTHON_SOURCES}
)

# Set output name and properties
set_target_properties(radia PROPERTIES
    PREFIX ""
    OUTPUT_NAME "radia"
    SUFFIX ".pyd"
)

# Include directories
target_include_directories(radia PRIVATE
    ${CORE_DIR}
    ${LIB_DIR}
    ${EXT_DIR}/auxparse
    ${EXT_DIR}/genmath
    ${EXT_DIR}/triangle
    ${PYTHON_DIR}
    ${Python3_INCLUDE_DIRS}
)

# Preprocessor definitions (no GLUT, no viewer)
target_compile_definitions(radia PRIVATE
    WIN32
    _WINDOWS
    ALPHA__LIB__
    FFTW_ENABLE_FLOAT
    NO_TIMER
    ANSI_DECLARATORS
    TRILIBRARY
    _GM_WITHOUT_BASE
    $<$<CONFIG:Debug>:_DEBUG>
    $<$<CONFIG:Release>:NDEBUG>
)

# Link Python libraries
target_link_libraries(radia PRIVATE
    Python3::Python
    OpenMP::OpenMP_CXX
)

# FFTW library (64-bit MSVC)
set(FFTW_LIB ${FFTW_DIR}/fftw64_f.lib)
if(EXISTS ${FFTW_LIB})
    target_link_libraries(radia PRIVATE ${FFTW_LIB})
else()
    message(FATAL_ERROR "FFTW library not found: ${FFTW_LIB}")
endif()

# Windows system libraries
if(WIN32)
    target_link_libraries(radia PRIVATE version.lib)
endif()

# Compiler-specific flags
if(MSVC)
    target_compile_options(radia PRIVATE
        /wd4244  # Disable conversion warnings
        /wd4267  # Disable size_t conversion warnings
        $<$<CONFIG:Release>:/O2 /GL>
        $<$<CONFIG:Debug>:/Od /Zi>
    )

    # Link-time code generation for Release
    if(CMAKE_BUILD_TYPE STREQUAL "Release")
        set_target_properties(radia PROPERTIES
            LINK_FLAGS "/LTCG"
        )
    endif()
endif()

# Print configuration summary
message(STATUS "")
message(STATUS "========================================")
message(STATUS "  Radia Python 3.12 Module")
message(STATUS "========================================")
message(STATUS "Build type:        ${CMAKE_BUILD_TYPE}")
message(STATUS "C++ Compiler:      ${CMAKE_CXX_COMPILER_ID} ${CMAKE_CXX_COMPILER_VERSION}")
message(STATUS "Platform:          ${CMAKE_SYSTEM_NAME} ${CMAKE_SYSTEM_PROCESSOR}")
message(STATUS "Architecture:      64-bit only")
message(STATUS "")
message(STATUS "Python:")
message(STATUS "  Version:         ${Python3_VERSION} (3.12 REQUIRED)")
message(STATUS "  Executable:      ${Python3_EXECUTABLE}")
message(STATUS "")
message(STATUS "Features:")
message(STATUS "  GLUT Viewer:     DISABLED (removed)")
message(STATUS "  MPI:             DISABLED")
message(STATUS "  OpenMP:          ENABLED (CPU parallelization)")
message(STATUS "  Igor XOP:        DISABLED (removed)")
message(STATUS "  Mathematica:     DISABLED (removed)")
message(STATUS "  Multi-Python:    DISABLED (3.12 only)")
message(STATUS "")
message(STATUS "Output:")
message(STATUS "  Module:          radia.pyd")
message(STATUS "  Location:        ${CMAKE_LIBRARY_OUTPUT_DIRECTORY}")
message(STATUS "========================================")
message(STATUS "")
