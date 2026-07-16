@echo off
setlocal
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0build_launcher_exe.ps1"
if errorlevel 1 (
    echo.
    echo Failed to build Ophiuchus.exe.
    pause
    exit /b 1
)
echo.
echo Ophiuchus.exe is ready.
pause
