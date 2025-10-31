cmake_minimum_required(VERSION 3.21)
project(rad_ngsolve C CXX)

# Set Python environment path (following EMPY_Field pattern)
set(ENVIRONMENT_PATH "C:/Program Files/Python312/")
message("ENVIRONMENT_PATH:" ${ENVIRONMENT_PATH})

# Set pybind11 directory
set(pybind11_DIR ${ENVIRONMENT_PATH}Lib/site-packages/pybind11/share/cmake/pybind11)

message("Finding pybind11...")
find_package(pybind11 REQUIRED HINTS ${ENVIRONMENT_PATH})

# Set NGSolve and Netgen directories
set(NGSolve_DIR ${ENVIRONMENT_PATH}Lib/site-packages/ngsolve/cmake)
set(Netgen_DIR ${ENVIRONMENT_PATH}Lib/site-packages/netgen/cmake)

# Find NGSolve
find_package(NGSolve CONFIG REQUIRED HINTS ${ENVIRONMENT_PATH})

# Radia source directories
set(CORE_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/core")
set(LIB_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/lib")
set(EXT_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/ext")
set(PYTHON_DIR "${CMAKE_CURRENT_SOURCE_DIR}/src/python")
set(FFTW_DIR "${EXT_DIR}/fftw")

# Radia core sources (explicit list from original CMakeLists.txt)
set(RADIA_CORE_SOURCES
    ${CORE_DIR}/radapl1.cpp
    ${CORE_DIR}/radapl2.cpp
    ${CORE_DIR}/radapl3.cpp
    ${CORE_DIR}/radapl4.cpp
    ${CORE_DIR}/radapl5.cpp
    ${CORE_DIR}/radapl6.cpp
    ${CORE_DIR}/radapl7.cpp
    ${CORE_DIR}/radarccu.cpp
    ${CORE_DIR}/radcffld.cpp
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

# Radia ext sources
set(RADIA_EXT_SOURCES
    ${EXT_DIR}/auxparse/auxparse.cpp
    ${EXT_DIR}/genmath/gmfft.cpp
    ${EXT_DIR}/genmath/gmtrans.cpp
    ${EXT_DIR}/triangle/triangle.c
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
add_ngsolve_python_module(rad_ngsolve
    ${RADIA_CORE_SOURCES}
    ${EXT_DIR}/auxparse/auxparse.cpp
    ${EXT_DIR}/genmath/gmfft.cpp
    ${EXT_DIR}/genmath/gmtrans.cpp
    ${RADIA_LIB_SOURCES}
    ${PYTHON_DIR}/rad_ngsolve.cpp
)

# Link triangle library
target_link_libraries(rad_ngsolve PRIVATE triangle_lib)

# Add include directories
target_include_directories(rad_ngsolve PRIVATE
    ${CORE_DIR}
    ${LIB_DIR}
    ${EXT_DIR}/auxparse
    ${EXT_DIR}/genmath
    ${EXT_DIR}/triangle
    ${PYTHON_DIR}
)

# Add compile definitions (from Radia)
target_compile_definitions(rad_ngsolve PRIVATE
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

# FFTW library (64-bit MSVC)
set(FFTW_LIB ${FFTW_DIR}/fftw64_f.lib)
if(EXISTS ${FFTW_LIB})
    target_link_libraries(rad_ngsolve PRIVATE ${FFTW_LIB})
else()
    message(WARNING "FFTW library not found: ${FFTW_LIB}")
endif()

# Link Windows system libraries
if(WIN32)
    target_link_libraries(rad_ngsolve PRIVATE version.lib)
endif()

# Compiler-specific flags (following EMPY_Field pattern - NO OpenMP)
if(MSVC)
    target_compile_options(rad_ngsolve PRIVATE
        /wd4244  # Disable conversion warnings
        /wd4267  # Disable size_t conversion warnings
    )
endif()

message(STATUS "rad_ngsolve module configured")

# ========================================
# Radia Python Module (radia.pyd)
# ========================================

# Create radia Python module using pybind11
pybind11_add_module(radia
    ${RADIA_CORE_SOURCES}
    ${EXT_DIR}/auxparse/auxparse.cpp
    ${EXT_DIR}/genmath/gmfft.cpp
    ${EXT_DIR}/genmath/gmtrans.cpp
    ${RADIA_LIB_SOURCES}
    ${PYTHON_DIR}/radpy.cpp
)

# Link triangle library
target_link_libraries(radia PRIVATE triangle_lib)

# Add include directories
target_include_directories(radia PRIVATE
    ${CORE_DIR}
    ${LIB_DIR}
    ${EXT_DIR}/auxparse
    ${EXT_DIR}/genmath
    ${FFTW_DIR}
)

# Add preprocessor definitions
target_compile_definitions(radia PRIVATE
    _CRT_SECURE_NO_WARNINGS
    _GM_WITHOUT_BASE
    ALPHA__DLL__
    ANSI_DECLARATORS
)

# FFTW library (64-bit MSVC)
if(EXISTS ${FFTW_LIB})
    target_link_libraries(radia PRIVATE ${FFTW_LIB})
else()
    message(WARNING "FFTW library not found: ${FFTW_LIB}")
endif()

# Link Windows system libraries
if(WIN32)
    target_link_libraries(radia PRIVATE version.lib)
endif()

# Compiler-specific flags
if(MSVC)
    target_compile_options(radia PRIVATE
        /wd4244  # Disable conversion warnings
        /wd4267  # Disable size_t conversion warnings
    )
endif()

message(STATUS "radia module configured")
