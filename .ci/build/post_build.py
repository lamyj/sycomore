import os
import subprocess
import sys
import sysconfig

workspace = os.environ["WORKSPACE"]
build_dir = os.environ.get("BUILD_DIR", os.path.join(workspace, "build"))
install_dir = os.environ.get("BUILD_DIR", os.path.join(workspace, "install"))

bin_dir = os.path.join(install_dir, "bin")
lib_dir = os.path.join(install_dir, "lib")
python_lib_dir = os.path.join(
    install_dir,
    sysconfig.get_path("platlib", os.name+"_user", {"userbase": "."}))
python_tests_dir = os.path.join(workspace, "tests", "python")

# Set-up environment: C++ library, Python module and test data location.
for name in ["DYLD_LIBRARY_PATH", "LD_LIBRARY_PATH"]:
    os.environ[name] = os.pathsep.join([
        lib_dir, *os.environ.get(name, "").split(os.pathsep)])
os.environ["PATH"] = os.pathsep.join([
    bin_dir, *os.environ.get("PATH", "").split(os.pathsep)])
os.environ["PYTHONPATH"] = os.pathsep.join([
    python_lib_dir, *os.environ.get("PYTHONPATH", "").split(os.pathsep)])
os.environ["SYCOMORE_TEST_DATA"] = os.path.join(workspace, "tests", "data")

# Run C++ and Python tests even if the former fails, return non-zero if any
# failed.
cpp_tests_return_code = subprocess.call(["ctest"], cwd=build_dir)
python_tests_return_code = subprocess.call(
    [sys.executable, "-m", "unittest", "discover", "-s", python_tests_dir], 
    cwd=build_dir)
sys.exit(max(cpp_tests_return_code, python_tests_return_code))
