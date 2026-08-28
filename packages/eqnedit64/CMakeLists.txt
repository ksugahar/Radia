cmake_minimum_required(VERSION 3.20)
project(eqnedit64_python LANGUAGES CXX RC)

if(NOT WIN32)
    message(FATAL_ERROR "eqnedit64 wheels currently support Windows only")
endif()

find_package(Python COMPONENTS Interpreter Development.Module REQUIRED)
execute_process(
    COMMAND "${Python_EXECUTABLE}" -m pybind11 --cmakedir
    OUTPUT_VARIABLE PYBIND11_CMAKE_DIR
    OUTPUT_STRIP_TRAILING_WHITESPACE
    COMMAND_ERROR_IS_FATAL ANY
)
find_package(pybind11 CONFIG REQUIRED PATHS "${PYBIND11_CMAKE_DIR}" NO_DEFAULT_PATH)

get_filename_component(RADIA_ROOT "${CMAKE_CURRENT_SOURCE_DIR}/../.." ABSOLUTE)
set(EQNEDIT64_ROOT "${RADIA_ROOT}/tools/eqnedit64")
set(EQNEDIT64_SRC "${EQNEDIT64_ROOT}/src")
set(EQNEDIT64_ASSETS "${EQNEDIT64_ROOT}/assets")

foreach(required
        "${EQNEDIT64_SRC}/eqnedit_pybind.cpp"
        "${EQNEDIT64_ASSETS}/latinmodern-math.otf"
        "${EQNEDIT64_ROOT}/web/equation-editor.js")
    if(NOT EXISTS "${required}")
        message(FATAL_ERROR "canonical Eqnedit64 source is missing: ${required}")
    endif()
endforeach()

# The extension uses the same process-private embedded math font as the GUI.
# Generate the resource in the build tree so the repository keeps one asset
# copy under tools/eqnedit64.
file(MAKE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/assets")
configure_file(
    "${EQNEDIT64_ASSETS}/latinmodern-math.otf"
    "${CMAKE_CURRENT_BINARY_DIR}/assets/latinmodern-math.otf"
    COPYONLY
)
set(EQNEDIT64_PY_RESOURCE "${CMAKE_CURRENT_BINARY_DIR}/eqnedit64_python.rc")
file(WRITE "${EQNEDIT64_PY_RESOURCE}"
    "#include <windows.h>\n200 RCDATA \"assets\\\\latinmodern-math.otf\"\n")

pybind11_add_module(eqnedit_core MODULE
    "${EQNEDIT64_SRC}/eqnedit_pybind.cpp"
    "${EQNEDIT64_SRC}/equation_render.cpp"
    "${EQNEDIT64_SRC}/tex_parser.cpp"
    "${EQNEDIT64_SRC}/equation_edit.cpp"
    "${EQNEDIT64_SRC}/equation_node.cpp"
    "${EQNEDIT64_SRC}/latex_emitter.cpp"
    "${EQNEDIT64_SRC}/mathml_emitter.cpp"
    "${EQNEDIT64_SRC}/math_symbols.cpp"
    "${EQNEDIT64_SRC}/palettes.cpp"
    "${EQNEDIT64_SRC}/tex_document.cpp"
    "${EQNEDIT64_PY_RESOURCE}"
)
target_include_directories(eqnedit_core PRIVATE "${EQNEDIT64_SRC}")
target_compile_features(eqnedit_core PRIVATE cxx_std_17)
target_compile_definitions(eqnedit_core PRIVATE
    _CRT_SECURE_NO_WARNINGS UNICODE _UNICODE)
target_compile_options(eqnedit_core PRIVATE /utf-8 /W4)
target_link_libraries(eqnedit_core PRIVATE gdi32 user32)
install(TARGETS eqnedit_core
    LIBRARY DESTINATION eqnedit64
    RUNTIME DESTINATION eqnedit64)

if(NOT DEFINED ENV{EQNEDIT64_RELEASE_EXE})
    message(FATAL_ERROR
        "EQNEDIT64_RELEASE_EXE must name the Eqnedit64.exe to package")
endif()
file(TO_CMAKE_PATH "$ENV{EQNEDIT64_RELEASE_EXE}" EQNEDIT64_RELEASE_EXE)
if(NOT EXISTS "${EQNEDIT64_RELEASE_EXE}")
    message(FATAL_ERROR "Eqnedit64 release executable is missing: ${EQNEDIT64_RELEASE_EXE}")
endif()
install(FILES "${EQNEDIT64_RELEASE_EXE}" DESTINATION eqnedit64)
install(FILES
    "${EQNEDIT64_ROOT}/web/equation-editor.js"
    "${EQNEDIT64_ROOT}/web/equation-editor.fragment.html"
    DESTINATION eqnedit64/web)
install(FILES
    "${EQNEDIT64_ASSETS}/GUST-FONT-LICENSE.txt"
    "${EQNEDIT64_ASSETS}/THIRD_PARTY_NOTICES.md"
    DESTINATION eqnedit64/licenses)
