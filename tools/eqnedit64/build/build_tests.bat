@echo off
setlocal
cd /d "%~dp0.."

if /I not "%EQNEDIT64_ISOLATED_TEST_SESSION%"=="1" (
    echo [ERROR] Full Eqnedit64 tests register and release a private math font.
    echo [ERROR] Run them in a disposable CI/VM/user session and set
    echo [ERROR] EQNEDIT64_ISOLATED_TEST_SESSION=1 there.
    exit /b 90
)

call build\build_eqnedt64.bat || exit /b 1
call build\build_pymodule.bat || exit /b 1
build\test_tex_document.exe || exit /b 1
python tests\run_model_tests.py || exit /b 1
pwsh -NoProfile -ExecutionPolicy Bypass -File build\test_background.ps1 || exit /b 1
pwsh -NoProfile -ExecutionPolicy Bypass -File build\test_ui_fuzz.ps1 || exit /b 1
pwsh -Sta -NoProfile -ExecutionPolicy Bypass -File build\test_external_paste.ps1 || exit /b 1

echo [OK] Eqnedit64 portable executable, background GUI fuzz, and external paste tests passed.
