set BUILD_DIR=%WORKSPACE%\build
set INSTALL_DIR=%WORKSPACE%\install

mkdir %BUILD_DIR%
mkdir %INSTALL_DIR%

cd %BUILD_DIR%

cmake ^
  -DCMAKE_INSTALL_PREFIX=%INSTALL_DIR% ^
  -DCMAKE_BUILD_TYPE=Release ^
  -Dpybind11_DIR="C:/Libraries/pybind11/share/cmake/pybind11" ^
  -DBOOST_ROOT="C:/Libraries/boost_1_67_0" ^
  -DCMAKE_GENERATOR_PLATFORM=%PLATFORM% ^
  ${CMAKE_OPTIONS} ^
  ..
if %errorlevel% neq 0 exit /b %errorlevel%

cmake --build . --config Release --target install --parallel 3
if %errorlevel% neq 0 exit /b %errorlevel%

cd %WORKSPACE%
