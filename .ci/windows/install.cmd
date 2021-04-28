REM Assume a working conda environment
conda install -c conda-forge cmake ninja pybind11 xsimd
if %errorlevel% neq 0 exit /b %errorlevel%

%PYTHON% -m pip install --upgrade pip
if %errorlevel% neq 0 exit /b %errorlevel%

REM Assume Python later than 3.5
%PYTHON% -m pip install -U numpy
if %errorlevel% neq 0 exit /b %errorlevel%

cd %WORKSPACE%
