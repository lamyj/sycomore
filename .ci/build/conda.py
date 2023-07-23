import subprocess
import sys

subprocess.check_call([
    "conda", "install", "--yes", "-c", "conda-forge",
    "boost", "cmake", "ninja", "numpy", "pybind11", "scipy", "xsimd",
    "xtensor", "xtensor-python"])
