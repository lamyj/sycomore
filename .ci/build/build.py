import multiprocessing
import os
import re
import subprocess
import sys

workspace = os.environ["WORKSPACE"]
build_dir = os.environ.get("BUILD_DIR", os.path.join(workspace, "build"))
install_dir = os.environ.get("INSTALL_DIR", os.path.join(workspace, "install"))

for dir in [build_dir, install_dir]:
    if not os.path.isdir(dir):
        os.makedirs(dir)

os.environ["CMAKE_PREFIX_PATH"] = os.pathsep.join([
    *os.environ.get("CMAKE_PREFIX_PATH", "").split(os.pathsep),
    os.path.join(os.path.expanduser("~"), "local")])

subprocess.check_call(
    [
        "cmake", f"-DPYTHON_EXECUTABLE={sys.executable}",
        "-DCMAKE_BUILD_TYPE=Release", f"-DCMAKE_INSTALL_PREFIX={install_dir}", 
        workspace],
    cwd=build_dir)

cmake_version = {}
with open(os.path.join(build_dir, "CMakeCache.txt")) as fd:
    cmake_version = re.findall(
        r"CMAKE_CACHE_(.+)_VERSION:[^=]+=(\d+)", fd.read())
cmake_version = dict(cmake_version)
cmake_version = tuple(
    [int(cmake_version[x]) for x in ["MAJOR", "MINOR", "PATCH"]])

if cmake_version >= (3, 12, 0):
    parallel = ["--parallel", str(multiprocessing.cpu_count())]
else:
    parallel = []
subprocess.check_call(
    [
        "cmake", "--build", ".", "--target", "install", "--config", "Release",
        *parallel],
    cwd=build_dir)
