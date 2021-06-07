import multiprocessing
import os
import subprocess
import sys

workspace = os.environ["WORKSPACE"]
build_dir = os.environ.get("BUILD_DIR", os.path.join(workspace, "build"))
install_dir = os.environ.get("BUILD_DIR", os.path.join(workspace, "install"))

for dir in [build_dir, install_dir]:
    if not os.path.isdir(dir):
        os.makedirs(dir)

os.environ["CMAKE_PREFIX_PATH"] = os.pathsep.join([
    *os.environ.get("CMAKE_PREFIX_PATH", "").split(os.pathsep),
    os.path.join(os.path.expanduser("~"), "local")
])
subprocess.check_call(
    [
        "cmake", 
        "-DPYTHON_EXECUTABLE={}".format(sys.executable),
        "-DCMAKE_BUILD_TYPE=Release",
        "-DCMAKE_INSTALL_PREFIX={}".format(install_dir), 
        workspace],
    cwd=build_dir)

subprocess.check_call(
    [
        "cmake", "--build", ".", "--target", "install", "--config", "Release",
        "--parallel", str(multiprocessing.cpu_count())],
    cwd=build_dir)
