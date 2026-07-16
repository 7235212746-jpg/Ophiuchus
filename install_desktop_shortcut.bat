@echo off
setlocal
powershell.exe -NoProfile -ExecutionPolicy Bypass -File "%~dp0install_desktop_shortcut.ps1"
if errorlevel 1 (
    echo.
    echo Failed to create the Ophiuchus desktop shortcut.
    pause
    exit /b 1
)
echo.
echo The Ophiuchus desktop shortcut is ready.
pause
