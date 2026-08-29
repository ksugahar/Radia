# Convenience build for IDE navigation and the C++ unit test.  It is NOT the
# release path: `build\build_eqnedt64.bat` produces dist\Eqnedit64.exe and
# `build\build_pymodule.bat` the headless module the Python tests import, and
# `build\build_tests.bat` runs the required checks.  Nothing here writes to
# dist, so a build from this file can never become the shipped executable.
cmake_minimum_required(VERSION 3.20)
project(Eqnedit64 VERSION 3.0.8 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(EQNEDIT_CORE_SOURCES
    src/equation_edit.cpp
    src/equation_node.cpp
    src/equation_render.cpp
    src/latex_emitter.cpp
    src/mathml_emitter.cpp
    src/math_symbols.cpp
    src/palettes.cpp
    src/tex_document.cpp
    src/tex_parser.cpp
)

add_library(eqnedit_core STATIC ${EQNEDIT_CORE_SOURCES})
target_include_directories(eqnedit_core PUBLIC src)
target_compile_definitions(eqnedit_core PRIVATE _CRT_SECURE_NO_WARNINGS)
if(MSVC)
    target_compile_options(eqnedit_core PRIVATE /W4 /utf-8)
endif()

if(WIN32)
    add_executable(Eqnedit64 WIN32
        src/eqnedt64_app.cpp
        src/eqnedit64.rc
    )
    target_include_directories(Eqnedit64 PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_definitions(Eqnedit64 PRIVATE
        _CRT_SECURE_NO_WARNINGS UNICODE _UNICODE
    )
    target_compile_options(Eqnedit64 PRIVATE /W4 /utf-8)
    target_link_libraries(Eqnedit64 PRIVATE
        eqnedit_core user32 gdi32 comctl32 comdlg32 shell32 ole32 imm32
        windowscodecs
    )
endif()

enable_testing()
add_executable(test_tex_document
    tests/test_tex_document.cpp
    src/tex_document.cpp
)
target_include_directories(test_tex_document PRIVATE src)
if(MSVC)
    target_compile_options(test_tex_document PRIVATE /W4 /utf-8)
endif()
add_test(NAME tex_document COMMAND test_tex_document)

add_executable(test_math_alphabets tests/test_math_alphabets.cpp)
target_link_libraries(test_math_alphabets PRIVATE eqnedit_core)
target_include_directories(test_math_alphabets PRIVATE src)
if(MSVC)
    target_compile_options(test_math_alphabets PRIVATE /W4 /utf-8)
endif()
add_test(NAME math_alphabets COMMAND test_math_alphabets)
