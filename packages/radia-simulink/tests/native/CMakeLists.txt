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

# Host link checks compile production sources in this test tree. Synthetic
# transaction physics is linked only into its separate scenario executable.
add_library(radia_ih_coupled_domain_host OBJECT "../../native/ih_coupled_domain.cpp")
target_compile_features(radia_ih_coupled_domain_host PRIVATE cxx_std_17)
add_library(radia_ih_coupled_domain STATIC "../../native/ih_coupled_domain.cpp"
    "../../native/ih_coupled_backend.cpp")
target_compile_features(radia_ih_coupled_domain PRIVATE cxx_std_17)
add_executable(test_ih_coupled_production_unavailable test_ih_coupled_production_unavailable.cpp)
target_link_libraries(test_ih_coupled_production_unavailable PRIVATE radia_ih_coupled_domain)
target_compile_features(test_ih_coupled_production_unavailable PRIVATE cxx_std_17)
add_test(NAME radia.ih_domain.production_unavailable COMMAND test_ih_coupled_production_unavailable)
add_executable(test_ih_coupled_domain_scenarios test_ih_coupled_domain_scenarios.cpp
    $<TARGET_OBJECTS:radia_ih_coupled_domain_host>)
target_compile_features(test_ih_coupled_domain_scenarios PRIVATE cxx_std_17)
foreach(target radia_ih_coupled_domain_host radia_ih_coupled_domain
        test_ih_coupled_production_unavailable test_ih_coupled_domain_scenarios)
    if(MSVC)
        target_compile_options(${target} PRIVATE /W4 /permissive-)
    else()
        target_compile_options(${target} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
endforeach()
foreach(scenario eddy_solve heat_projection thermal_solve port_gate thermal_gate
        host_snapshot_allocation reset_copy_or_sync
        device_failure trial_sync malformed_snapshot oversized_snapshot
        outputs_repeated reset_success closed_snapshot
        owner_idempotency owner_failed_advance owner_reset_generation owner_poison)
    add_test(NAME radia.ih_domain.${scenario} COMMAND test_ih_coupled_domain_scenarios ${scenario})
endforeach()

# Peer-agreed registry/runtime contracts use no production physics factory.
find_package(Threads REQUIRED)
foreach(name IN ITEMS ih_coupled_registry ih_coupled_runtime)
    add_executable(test_${name} test_${name}.cpp)
    target_compile_features(test_${name} PRIVATE cxx_std_17)
    if(MSVC)
        target_compile_options(test_${name} PRIVATE /W4 /permissive-)
    else()
        target_compile_options(test_${name} PRIVATE -Wall -Wextra -Wpedantic)
    endif()
    target_link_libraries(test_${name} PRIVATE Threads::Threads)
    add_test(NAME radia.${name}.cpp COMMAND test_${name})
endforeach()
