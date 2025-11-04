@echo off

REM Prompt user to browse to python.exe (defaults to per-user Python39 if it exists)
for /f "usebackq delims=" %%I in (`
  powershell -NoProfile -STA -Command "Add-Type -AssemblyName System.Windows.Forms; $p = Join-Path $env:LOCALAPPDATA 'Programs\Python'; $dlg = New-Object System.Windows.Forms.OpenFileDialog; $dlg.Title = 'Select python.exe'; $dlg.Filter = 'Python executable|python.exe'; if (Test-Path $p) { $dlg.InitialDirectory = $p }; if ($dlg.ShowDialog() -eq [System.Windows.Forms.DialogResult]::OK) { [Console]::WriteLine($dlg.FileName) }"
`) do set "PYTHON_PATH=%%I"

REM Abort if dialog was closed or path is invalid
if not defined PYTHON_PATH (
    echo No file selected. Exiting.
    pause
    exit /b 1
)
if not exist "%PYTHON_PATH%" (
    echo Selected file does not exist: "%PYTHON_PATH%"
    pause
    exit /b 1
)

REM Verify Python version is >= 3.12
set "PY_VER="
"%PYTHON_PATH%" -c "import sys;print('.'.join(map(str, sys.version_info[:3])))" > "%TEMP%\_pyver.txt" 2>NUL
set /p PY_VER=<"%TEMP%\_pyver.txt"
del "%TEMP%\_pyver.txt" >NUL 2>&1

if not defined PY_VER (
    echo Failed to read Python version from "%PYTHON_PATH%".
    pause
    exit /b 1
)

for /f "tokens=1,2 delims=." %%a in ("%PY_VER%") do (
    set "MAJOR=%%a"
    set "MINOR=%%b"
)

if %MAJOR% LSS 3 (
    echo Python version %PY_VER% is too old. Need 3.12 or higher.
    pause
    exit /b 1
)

if %MAJOR% EQU 3 if %MINOR% LSS 12 (
    echo Python version %PY_VER% is too old. Need 3.12 or higher.
    pause
    exit /b 1
)

echo Detected Python %PY_VER% OK.



set ENV_NAME=pypbee
REM Store the virtual environment in "venv\<ENV_NAME>"
set ENV_FOLDER=venv\%ENV_NAME%
set ENV_PYTHON=%ENV_FOLDER%\Scripts\python.exe

REM Check if the specified Python executable exists
if not exist "%PYTHON_PATH%" (
    echo Python executable not found at %PYTHON_PATH%.
    pause
    exit /b 1
)

REM Remove any existing virtual environment and create a fresh one
if exist "%ENV_FOLDER%" (
    echo Removing existing virtual environment at "%ENV_FOLDER%"...
    rmdir /s /q "%ENV_FOLDER%"
)

echo Creating virtual environment...
"%PYTHON_PATH%" -m venv "%ENV_FOLDER%"

REM Activate the virtual environment
echo Activating virtual environment...
call "%ENV_FOLDER%\Scripts\activate"

REM Check if requirements.txt exists
if not exist requirements.txt (
    echo requirements.txt not found in the current directory.
    pause
    exit /b 1
)

REM Install requirements
echo Installing requirements from requirements.txt...
python -m pip install --upgrade --no-cache-dir pip
pip install --no-cache-dir -r requirements.txt

echo Setup complete.
pause