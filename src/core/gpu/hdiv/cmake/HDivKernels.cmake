# Reusable numerical archive: no Python, MATLAB or NGSolve dependency.
# RADIA_CUDA_OBJECT is compiled once by the reviewed nvcc entrypoint.
foreach(required RADIA_CUDA_OBJECT RADIA_OUTPUT RADIA_CUDART RADIA_CUBLAS RADIA_CUSOLVER)
  if(NOT DEFINED ${required} OR "${${required}}" STREQUAL "")
    message(FATAL_ERROR "Missing explicit ${required}")
  endif()
endforeach()
if(NOT EXISTS "${RADIA_CUDA_OBJECT}")
  message(FATAL_ERROR "Missing reviewed CUDA object")
endif()
if(TARGET radia_hdiv_kernels)
  message(FATAL_ERROR "HDiv numerical archive target already defined")
endif()
set_source_files_properties("${RADIA_CUDA_OBJECT}" PROPERTIES EXTERNAL_OBJECT TRUE GENERATED TRUE)
add_library(radia_hdiv_kernels STATIC "${RADIA_CUDA_OBJECT}")
set_target_properties(radia_hdiv_kernels PROPERTIES LINKER_LANGUAGE CXX
  POSITION_INDEPENDENT_CODE ON ARCHIVE_OUTPUT_DIRECTORY "${RADIA_OUTPUT}")
target_link_libraries(radia_hdiv_kernels PUBLIC
  "${RADIA_CUDART}" "${RADIA_CUBLAS}" "${RADIA_CUSOLVER}")
