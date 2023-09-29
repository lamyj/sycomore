import os
import subprocess
import sys
import sysconfig

workspace = os.environ["WORKSPACE"]
build_dir = os.environ.get("BUILD_DIR", os.path.join(workspace, "build"))
install_dir = os.environ.get("INSTALL_DIR", os.path.join(workspace, "install"))

bin_dir = os.path.join(install_dir, "bin")
lib_dir = os.path.join(install_dir, "lib")
python_lib_dir = os.path.join(
    install_dir,
    sysconfig.get_path("platlib", os.name+"_user", {"userbase": "."}))
python_tests_dir = os.path.join(workspace, "tests", "python")

# Set-up common environment: test data location
os.environ["SYCOMORE_TEST_DATA"] = os.path.join(workspace, "tests", "data")

# Run C++ and Python tests even if the former fails, return non-zero if any
# failed.

# C++ tests: only library path is needed
environment = os.environ.copy()
for name in ["DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"]:
    environment[name] = os.pathsep.join([
        lib_dir, *os.environ.get(name, "").split(os.pathsep)])
environment["PATH"] = os.pathsep.join([
    bin_dir, *os.environ.get("PATH", "").split(os.pathsep)])
cpp_tests_return_code = subprocess.call(
    ["ctest"], cwd=build_dir, env=environment)

# Python tests: library path is needed since the Python extension links to it
# Add PYTHONPATH for the pure-python module
environment["PYTHONPATH"] = os.pathsep.join([
    python_lib_dir, *environment.get("PYTHONPATH", "").split(os.pathsep)])
python_tests_return_code = subprocess.call(
    [sys.executable, "-m", "unittest", "discover", "-s", python_tests_dir], 
    cwd=build_dir, env=environment)

sys.exit(max(cpp_tests_return_code, python_tests_return_code))
