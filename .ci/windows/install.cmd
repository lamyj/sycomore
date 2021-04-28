REM Assume a working conda environment
conda install --yes -c conda-forge cmake ninja numpy pybind11 xsimd
if %errorlevel% neq 0 exit /b %errorlevel%
