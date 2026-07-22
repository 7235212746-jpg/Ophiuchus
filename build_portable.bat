@echo off
setlocal
cd /d "%~dp0"
powershell -NoProfile -ExecutionPolicy Bypass -File "%~dp0build_portable.ps1" -Python "%USERPROFILE%\anaconda3\envs\ophi\python.exe"
set "OPHI_EXIT=%ERRORLEVEL%"
if not "%OPHI_EXIT%"=="0" pause
exit /b %OPHI_EXIT%
