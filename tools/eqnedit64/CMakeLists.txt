# Convenience build for IDE navigation and the C++ unit test.  It is NOT the
# release path: `build\build_eqnedt64.bat` produces dist\Eqnedit64.exe and
# `build\build_pymodule.bat` the headless module the Python tests import, and
# `build\build_tests.bat` runs the required checks.  Nothing here writes to
# dist, so a build from this file can never become the shipped executable.
cmake_minimum_required(VERSION 3.20)
project(Eqnedit64 VERSION 3.1.2 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

set(EQNEDIT_CORE_SOURCES
    src/equation_edit.cpp
    src/equation_node.cpp
    src/equation_render.cpp
    src/latex_emitter.cpp
    src/mathml_emitter.cpp
    src/office_html.cpp
    src/math_symbols.cpp
    src/palettes.cpp
    src/tex_document.cpp
    src/tex_parser.cpp
)

add_library(eqnedit_core STATIC ${EQNEDIT_CORE_SOURCES})
target_include_directories(eqnedit_core PUBLIC src)
target_compile_definitions(eqnedit_core PRIVATE _CRT_SECURE_NO_WARNINGS)
if(MSVC)
    target_compile_options(eqnedit_core PRIVATE /W4 /WX /w14062 /utf-8)
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
    # Proof images call the real owner-draw implementation, not a duplicate
    # UI renderer. Run only in an isolated font session; never on LAB desktop.
    add_executable(test_palette_visual tests/test_palette_visual.cpp src/eqnedit64.rc)
    target_include_directories(test_palette_visual PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_compile_definitions(test_palette_visual PRIVATE _CRT_SECURE_NO_WARNINGS UNICODE _UNICODE)
    target_compile_options(test_palette_visual PRIVATE /W4 /WX /w14062 /utf-8)
    target_link_libraries(test_palette_visual PRIVATE eqnedit_core user32 gdi32
        comctl32 comdlg32 shell32 ole32 imm32 windowscodecs)
endif()

enable_testing()
add_executable(test_color tests/test_color.cpp)
target_link_libraries(test_color PRIVATE eqnedit_core)
if(MSVC)
    target_compile_options(test_color PRIVATE /W4 /WX /utf-8)
endif()
add_test(NAME color COMMAND test_color)
if(WIN32)
    # Manual disposable-CI diagnostic, intentionally not a CTest/release target.
    add_executable(test_color_render tests/test_color_render.cpp src/eqnedit64.rc)
    target_include_directories(test_color_render PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
    target_link_libraries(test_color_render PRIVATE eqnedit_core gdi32 user32)
    target_compile_definitions(test_color_render PRIVATE NOMINMAX WIN32_LEAN_AND_MEAN _CRT_SECURE_NO_WARNINGS)
    target_compile_options(test_color_render PRIVATE /W4 /WX /utf-8)
    add_test(NAME color_render COMMAND test_color_render)
    add_executable(font_lifecycle_probe EXCLUDE_FROM_ALL tests/font_lifecycle_probe.cpp)
    target_compile_options(font_lifecycle_probe PRIVATE /W4 /WX /utf-8)
    target_link_libraries(font_lifecycle_probe PRIVATE gdi32 user32)
    add_executable(test_font_trace tests/test_font_trace.cpp)
    target_include_directories(test_font_trace PRIVATE src)
    target_compile_options(test_font_trace PRIVATE /W4 /WX /utf-8)
    add_test(NAME font_trace COMMAND test_font_trace "${CMAKE_CURRENT_BINARY_DIR}/font-trace-test")
endif()
add_executable(test_prime_base tests/test_prime_base.cpp)
target_link_libraries(test_prime_base PRIVATE eqnedit_core)
if(MSVC)
    target_compile_options(test_prime_base PRIVATE /W4 /WX /utf-8)
endif()
add_test(NAME prime_base COMMAND test_prime_base)
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

add_executable(test_node_dispatch tests/test_node_dispatch.cpp)
target_link_libraries(test_node_dispatch PRIVATE eqnedit_core)
if(WIN32)
    target_link_libraries(test_node_dispatch PRIVATE gdi32 user32)
endif()
if(MSVC)
    target_compile_options(test_node_dispatch PRIVATE /W4 /WX /utf-8)
endif()
add_test(NAME node_dispatch COMMAND test_node_dispatch)

# Pure parser/serializer probe: no font registration, rendering or Office.
add_executable(office_payload_probe tests/office_payload_probe.cpp)
target_link_libraries(office_payload_probe PRIVATE eqnedit_core)
if(MSVC)
    target_compile_options(office_payload_probe PRIVATE /W4 /WX /utf-8)
endif()
