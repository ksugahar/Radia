cmake_minimum_required(VERSION 3.21)
project(radia_simulink_native_acceptance LANGUAGES CXX)
enable_testing()
get_filename_component(RADIA_ROOT "${CMAKE_CURRENT_LIST_DIR}/../../../.." ABSOLUTE)
add_executable(test_radia_mex_sparse_symmetry
    "${RADIA_ROOT}/tests/cpp/test_radia_mex_sparse_symmetry.cpp")
target_include_directories(test_radia_mex_sparse_symmetry PRIVATE "${RADIA_ROOT}/src/core")
target_compile_features(test_radia_mex_sparse_symmetry PRIVATE cxx_std_17)
if(MSVC)
    target_compile_options(test_radia_mex_sparse_symmetry PRIVATE /W4 /permissive-)
else()
    target_compile_options(test_radia_mex_sparse_symmetry PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME radia.mex_sparse_symmetry.cpp COMMAND test_radia_mex_sparse_symmetry)
add_executable(test_radia_ih_thermal_acceptance
    "${RADIA_ROOT}/tests/cpp/test_radia_ih_thermal_acceptance.cpp"
    "${RADIA_ROOT}/src/runtime/ih/radia_ih_runtime.cpp"
    "${RADIA_ROOT}/src/runtime/ih/radia_ih_thermal.cpp"
    "${RADIA_ROOT}/src/runtime/ih/radia_ih_transport.cpp")
target_include_directories(test_radia_ih_thermal_acceptance PRIVATE "${RADIA_ROOT}/src/runtime/ih")
target_compile_features(test_radia_ih_thermal_acceptance PRIVATE cxx_std_17)
add_test(NAME radia.ih_thermal_acceptance.cpp COMMAND test_radia_ih_thermal_acceptance)

add_executable(test_ih_coupled_config test_ih_coupled_config.cpp)
target_compile_features(test_ih_coupled_config PRIVATE cxx_std_17)
if(MSVC)
    target_compile_options(test_ih_coupled_config PRIVATE /W4 /permissive-)
else()
    target_compile_options(test_ih_coupled_config PRIVATE -Wall -Wextra -Wpedantic)
endif()
add_test(NAME radia.ih_coupled_config.cpp COMMAND test_ih_coupled_config)
