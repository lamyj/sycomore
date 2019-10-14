REM WARNING: assume boost, cmake, and python are installed (this is the case
REM with Appveyor). Otherwise, "choco install" works for CMake and Python.

REM NOTE: OpenMP support is native in Visual C++

cd C:\projects
git clone https://github.com/pybind/pybind11.git
cd pybind11
git checkout v2.4.2
cmake -DPYBIND11_TEST=OFF -DCMAKE_INSTALL_PREFIX=C:\Libraries\pybind11 -DCMAKE_GENERATOR_PLATFORM=%PLATFORM% .
cmake --build . --config Release --target install

%PYTHON% -m pip install --upgrade pip

REM Assume Python later than 3.5
%PYTHON% -m pip install -U numpy

cd %WORKSPACE%
