@echo off
REM Set the path to the desired Python executable
set PYTHON_PATH=C:\Users\joel-students\AppData\Local\Programs\Python\Python39\python.exe
set ENV_NAME=pypbee

REM Store the virtual environment in "venv\<ENV_NAME>"
set ENV_FOLDER=venv\%ENV_NAME%

REM Check if the specified Python executable exists
if not exist "%PYTHON_PATH%" (
    echo Python executable not found at %PYTHON_PATH%.
    pause
    exit /b 1
)

REM Create the "venv" folder if it doesn't exist
if not exist "venv" (
    mkdir venv
)

REM Create a virtual environment
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
pip install -r requirements.txt

echo Setup complete.
pause
