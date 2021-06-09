import subprocess

subprocess.check_call([
    "conda", "install", "--yes", "-c", "conda-forge",
    "boost", "cmake", "ninja", "numpy", "pybind11", "xsimd"])
