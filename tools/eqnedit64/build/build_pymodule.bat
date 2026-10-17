@echo off
setlocal enabledelayedexpansion
REM Build the TeX-only Eqnedit64 core module for background tests.

call "%~dp0select_vs2022.bat"
if errorlevel 1 exit /b 1

cd /d "%~dp0.."
for /f "usebackq delims=" %%I in (`python -c "import sysconfig;print(sysconfig.get_paths()['include'])"`) do set "PYINC=%%I"
for /f "usebackq delims=" %%I in (`python -c "import pybind11;print(pybind11.get_include())"`) do set "PYBIND=%%I"
for /f "usebackq delims=" %%I in (`python -c "import sysconfig,os;print(os.path.join(sysconfig.get_config_var('installed_base'),'libs'))"`) do set "PYLIBS=%%I"
for /f "usebackq delims=" %%I in (`python -c "import sysconfig;print(sysconfig.get_config_var('EXT_SUFFIX'))"`) do set "EXTSUF=%%I"

REM Embed the same resources the executable carries -- above all the font.
REM Loading it from the file on the share instead was unreliable: private
REM font registrations are torn down when a process exits, and a process
REM starting during that teardown intermittently got nothing, so every glyph
REM measured zero wide and three unrelated layout checks failed at once.
rc /nologo /fo build\Eqnedit64.res src\eqnedit64.rc
if errorlevel 1 ( echo [ERROR] resources failed & exit /b 1 )

echo [INFO] building build\eqnedit_core!EXTSUF!
cl /nologo /O2 /W4 /WX /EHsc /MD /std:c++17 /utf-8 ^
   /D_CRT_SECURE_NO_WARNINGS /I src /I "!PYINC!" /I "!PYBIND!" ^
   /TP src\eqnedit_pybind.cpp src\equation_render.cpp src\tex_parser.cpp ^
       src\equation_edit.cpp src\equation_node.cpp src\latex_emitter.cpp ^
       src\mathml_emitter.cpp src\office_html.cpp src\math_symbols.cpp src\palettes.cpp src\tex_document.cpp ^
   /LD /Fo:build\ /Fe:build\eqnedit_core!EXTSUF! ^
   /link /LIBPATH:"!PYLIBS!" build\Eqnedit64.res gdi32.lib user32.lib
if errorlevel 1 ( echo [ERROR] eqnedit_core module build failed & exit /b 1 )

if not exist "build\eqnedit_core!EXTSUF!" (
    echo [ERROR] eqnedit_core module missing
    exit /b 1
)
echo [OK] build\eqnedit_core!EXTSUF!
