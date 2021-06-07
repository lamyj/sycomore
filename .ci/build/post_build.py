import os
import subprocess
import sys

import distutils.sysconfig

workspace = os.environ["WORKSPACE"]
build_dir = os.environ.get("BUILD_DIR", os.path.join(workspace, "build"))
install_dir = os.environ.get("BUILD_DIR", os.path.join(workspace, "install"))

bin_dir = os.path.join(install_dir, "bin")
lib_dir = os.path.join(install_dir, "lib")
python_lib_dir = os.path.join(
    install_dir, distutils.sysconfig.get_python_lib(True, prefix=""))
python_tests_dir = os.path.join(workspace, "tests", "python")

print(lib_dir, python_lib_dir, python_tests_dir)

# Set-up environment: C++ library, Python module and test data location.
for name in ["DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"]:
    os.environ[name] = os.pathsep.join([
        *os.environ.get(name, "").split(os.pathsep), lib_dir])
os.environ["PATH"] = os.pathsep.join([
    *os.environ.get("PATH", "").split(os.pathsep), bin_dir])
os.environ["PYTHONPATH"] = os.pathsep.join([
    *os.environ.get("PYTHONPATH", "").split(os.pathsep), python_lib_dir])
os.environ["SYCOMORE_TEST_DATA"] = os.path.join(workspace, "tests", "data")

# Run C++ and Python tests even if the former fails, return non-zero if any
# failed.
cpp_tests_return_code = subprocess.call(["ctest"], cwd=build_dir)
python_tests_return_code = subprocess.call(
    [sys.executable, "-m", "unittest", "discover", "-s", python_tests_dir], 
    cwd=build_dir)
sys.exit(max(cpp_tests_return_code, python_tests_return_code))
