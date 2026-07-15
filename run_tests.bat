@echo off
setlocal
cd /d "%~dp0"
set "OPHI_PY_CMD="

if defined OPHI_PYTHON if exist "%OPHI_PYTHON%" set OPHI_PY_CMD="%OPHI_PYTHON%"

if not defined OPHI_PY_CMD if exist "%USERPROFILE%\anaconda3\envs\ophi\python.exe" (
    set OPHI_PY_CMD="%USERPROFILE%\anaconda3\envs\ophi\python.exe"
)

if not defined OPHI_PY_CMD if exist "%USERPROFILE%\miniconda3\envs\ophi\python.exe" (
    set OPHI_PY_CMD="%USERPROFILE%\miniconda3\envs\ophi\python.exe"
)

if not defined OPHI_PY_CMD (
    where conda >nul 2>nul
    if not errorlevel 1 (
        conda run -n ophi python -c "import sys" >nul 2>nul
        if not errorlevel 1 set "OPHI_PY_CMD=conda run -n ophi python"
    )
)

if not defined OPHI_PY_CMD (
    where py >nul 2>nul
    if not errorlevel 1 set "OPHI_PY_CMD=py -3"
)

if not defined OPHI_PY_CMD (
    where python >nul 2>nul
    if not errorlevel 1 set "OPHI_PY_CMD=python"
)

if not defined OPHI_PY_CMD (
    echo Python was not found. Install the Ophiuchus environment described in docs\Ophiuchus_操作手册.md.
    pause
    exit /b 1
)

%OPHI_PY_CMD% -m unittest discover -s tests -v
set "OPHI_EXIT=%ERRORLEVEL%"
if not "%OPHI_EXIT%"=="0" pause
exit /b %OPHI_EXIT%
