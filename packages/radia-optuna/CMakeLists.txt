cmake_minimum_required(VERSION 3.20)
project(radia_optuna_kernels LANGUAGES CXX)
if(NOT CMAKE_SIZEOF_VOID_P EQUAL 8)
  message(FATAL_ERROR "Windows x64 pointer ABI required")
endif()
# nvcc object is built explicitly by the package entry; no SDK/compiler fallback.
if(NOT WIN32 OR NOT MSVC OR NOT MSVC_VERSION EQUAL 1944)
  message(FATAL_ERROR "Windows MSVC 14.44 baseline required")
endif()
if(NOT EXISTS "${OPTUNA_CUDA_OBJECT}" OR NOT EXISTS "${OPTUNA_CUDA_ROOT}/include/cuda_runtime.h")
  message(FATAL_ERROR "Explicit nvcc object and bundled CUDA headers required")
endif()
set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_MSVC_RUNTIME_LIBRARY MultiThreadedDLL)
set(CMAKE_EXPORT_COMPILE_COMMANDS ON)
find_package(Python3 REQUIRED COMPONENTS Interpreter Development.Module)
find_package(pybind11 3.0.2 EXACT CONFIG REQUIRED)
find_package(OpenMP REQUIRED COMPONENTS CXX)
find_library(OPTUNA_CUDART NAMES cudart PATHS "${OPTUNA_CUDA_ROOT}/lib/x64" NO_DEFAULT_PATH REQUIRED)

set_source_files_properties("${OPTUNA_CUDA_OBJECT}" PROPERTIES EXTERNAL_OBJECT TRUE GENERATED TRUE)
add_library(optuna_kernels STATIC "${OPTUNA_CUDA_OBJECT}")
set_target_properties(optuna_kernels PROPERTIES LINKER_LANGUAGE CXX)
target_link_libraries(optuna_kernels PUBLIC "${OPTUNA_CUDART}" OpenMP::OpenMP_CXX)
target_include_directories(optuna_kernels PUBLIC native)
file(GENERATE OUTPUT "${CMAKE_BINARY_DIR}/object-list.txt" CONTENT "${OPTUNA_CUDA_OBJECT}\n")

pybind11_add_module(_optuna_kernels NO_EXTRAS native/optuna_python.cpp native/optuna_python_runtime.cpp)
target_link_libraries(_optuna_kernels PRIVATE optuna_kernels)
target_compile_definitions(_optuna_kernels PRIVATE NOMINMAX)
target_include_directories(_optuna_kernels PRIVATE "${OPTUNA_CUDA_ROOT}/include")
set_target_properties(_optuna_kernels PROPERTIES
  LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/radia_optuna"
  RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/src/radia_optuna")

option(OPTUNA_BUILD_MATLAB "Link plain MEX only after Python kernel verification" OFF)
if(OPTUNA_BUILD_MATLAB)
  find_package(Matlab REQUIRED COMPONENTS MX_LIBRARY)
  find_library(OPTUNA_MATLAB_GPU_LIBRARY NAMES gpumexbinder
    PATHS "${Matlab_ROOT_DIR}/extern/lib/win64/microsoft" NO_DEFAULT_PATH REQUIRED)
  matlab_add_mex(NAME optuna_mex SRC native/optuna_mex.cpp native/optuna_matlab_runtime.cpp
    LINK_TO optuna_kernels "${OPTUNA_MATLAB_GPU_LIBRARY}" R2018a)
  target_compile_definitions(optuna_mex PRIVATE NOMINMAX
    _CRT_SECURE_NO_DEPRECATE _SCL_SECURE_NO_DEPRECATE _SECURE_SCL=0)
  target_compile_options(optuna_mex PRIVATE /Zp8 /GR /W3 /EHs /MD /O2 /Oy- /openmp)
  target_include_directories(optuna_mex PRIVATE "${Matlab_ROOT_DIR}/toolbox/parallel/gpu/extern/include")
  set_target_properties(optuna_mex PROPERTIES
    RUNTIME_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/matlab"
    LIBRARY_OUTPUT_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}/matlab")
endif()
