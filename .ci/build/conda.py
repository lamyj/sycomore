import subprocess
import sys

conda = sys.argv[1] if len(sys.argv) >= 2 else "conda"

subprocess.check_call([
    conda, "install", "--yes", "-c", "conda-forge",
    "boost", "boost-cpp", "cmake", "make", "numpy", "pybind11", "scipy",
    "xsimd", "xtensor", "xtensor-python"])
