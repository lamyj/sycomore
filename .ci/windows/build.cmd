set BUILD_DIR=%WORKSPACE%\build
set INSTALL_DIR=%WORKSPACE%\install

mkdir %BUILD_DIR%
mkdir %INSTALL_DIR%

cd %BUILD_DIR%

cmake ^
  -G Ninja ^
  -DCMAKE_INSTALL_PREFIX=%INSTALL_DIR% ^
  -DCMAKE_BUILD_TYPE=Release ^
  -DBOOST_ROOT="C:/Libraries/boost_1_67_0" ^
  %CMAKE_OPTIONS% ^
  ..
if %errorlevel% neq 0 exit /b %errorlevel%

cmake --build . --config Release --target install
if %errorlevel% neq 0 exit /b %errorlevel%

cd %WORKSPACE%
